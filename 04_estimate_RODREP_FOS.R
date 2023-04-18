library(tidyverse)
library(cmdstanr)

set.seed(1234)

# load retest
df_response <- read_csv("results/df_response.csv")
df_sequence <- read_csv("results/df_sequence.csv")
df_retest <- read_csv("results/retest_data_RODREP.csv")

df_retest %>%
  mutate(tt = factor(tt, levels = 0:35),
         THRESHOLD = factor(THRESHOLD, levels = -1:40)) %>%
  count(tt, THRESHOLD, .drop = FALSE) %>% 
  pivot_wider(names_from = THRESHOLD, 
              values_from = n) %>% 
  as.data.frame() -> Y_retest

stan_data <- list(
  S = 41,
  K = 42,
  X = 36,
  N_response = nrow(df_response),
  seen = df_response$seen,
  db = df_response$db,
  N_sequence = nrow(df_sequence),
  final = df_sequence$final,
  start_row = df_sequence$start_row,
  stop_row = df_sequence$stop_row,
  Y = Y_retest[, -1],
  gamma = 0.03
)


stan_model <- cmdstan_model("biaced_die_model.stan")

stan_fit <- stan_model$variational(
  data = stan_data,
  seed = 1234
)

stan_fit$cmdstan_diagnose()

stan_fit$summary(c("mu", "sigma", "lamda", "asym")) %>% 
  select(variable, mean, q5, q95) %>%
  extract(variable, 
          c("name", "tt"), "([a-z]+)\\[([0-9]+)", 
          convert = TRUE) %>% 
  mutate(tt = tt - 1) -> res_RODREP

res_RODREP %>% 
  filter(name == "sigma") %>% 
  print(n=36)


# for ARVO poster

psi <- function(x, mu, sigma, lamda, gamma) {
  y = gamma + (1 - lamda - gamma) * (1 - pnorm(x, mu, sigma))
}

res_RODREP %>% 
  select(name, tt, mean) %>% 
  pivot_wider(names_from = name, values_from = mean) %>% 
  mutate(gamma = 0.03) %>% 
  crossing(x = seq(0, 40, 0.1)) %>% 
  mutate(y = psi(x, mu, sigma, lamda, gamma)) %>% 
  ggplot(aes(x = x, y = y)) +
  geom_line(aes(group = tt), size = 0.2) +
  scale_x_continuous(limits = c(0, 40), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0), 
                     labels = scales::percent) +
  labs(x = "Stimulus luminance (dB)",
       y = "Probability of seen") +
  theme_linedraw() +
  theme(legend.position = "none")


res_RODREP %>% 
  filter(name == "asym") %>% 
  ggplot(aes(x = tt, y = mean)) +
  geom_ribbon(aes(ymin = q5, ymax = q95), fill = "blue", alpha = 0.2) +
  geom_line() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0),
                     labels = scales::percent) +
  labs(x = "",
       y = "Asymptotic maximum \n response probability") +
  theme_linedraw() +
  theme(panel.grid.minor = element_blank(),
        plot.caption = element_text(hjust = 0)) -> p1

res_RODREP %>% 
  filter(name == "mu") %>% 
  ggplot(aes(x = tt, y = mean)) +
  geom_ribbon(aes(ymin = q5, ymax = q95), fill = "blue", alpha = 0.2) +
  geom_line() +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 40), expand = c(0, 0)) +
  labs(x = "Perimetric Sensitivity (dB)",
       y = "Luminance of \n inflection point (dB)") +
  theme_linedraw() +
  theme(panel.grid.minor = element_blank()) -> p2

res_RODREP %>% 
  filter(name == "sigma") %>% 
  mutate(Henson = map_dbl(tt, ~ min(exp(-0.081 * . + 3.27), 6))) %>%
  ggplot(aes(x = tt, y = mean)) +
  geom_ribbon(aes(ymin = q5, ymax = q95), fill = "blue", alpha = 0.2) +
  geom_line() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 8), expand = c(0, 0)) +
  labs(x = "",
       y = "Slope \n (standard deviation)") +
  theme_linedraw() +
  theme(panel.grid.minor = element_blank()) -> p3

library(cowplot)

plot_grid(p1, p2, p3, labels = "AUTO", 
          axis = "b",
          align = "h",
          nrow = 1)


res_RODREP %>% 
  rename(parameter = name,
         sensitivity = tt) %>% 
  write_csv(file = "estimated_parameters.csv")
