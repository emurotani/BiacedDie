library(tidyverse)
library(OPI)
library(cmdstanr)

source("full_threshold_2.R")

set.seed(1234)

# initialize
chooseOpi("SimHenson")
opiInitialise(type = "C", cap = 6)

makeStim <- function(db, n) { 
  s <- list(level = dbTocd(db))
  class(s) <- "opiStaticStimulus"
  return(s)
}

# set parameters
fpr <- 0.03
fnr <- 0.15

# create test-retest data using conventional FOS curves
crossing(tt = 0:35, n_times = 1:100) %>% 
  mutate(FT2_res = map(tt, ~ FT2(tt = .x, 
                                 fpr = fpr,
                                 fnr = fnr,
                                 makeStim = makeStim))) %>% 
  mutate(final = map_dbl(FT2_res, ~ pluck(.x, "final"))) %>% 
  mutate(final = factor(final, levels = -1:40)) %>% 
  count(tt, final, .drop = FALSE) %>% 
  pivot_wider(names_from = final, values_from = n) %>% 
  as.matrix() -> df_Henson

df_Henson
df_Henson %>% dim()
# 36 rows: true perimetric sensitivity 0, ..., 35dB
# 43 columns: 1, true perimetric sensitivity
#             2 ~ 43, count of final estimates <0, 0, ..., 40dB

# setup data as a list
df_response <- read_csv("results/df_response.csv")
df_sequence <- read_csv("results/df_sequence.csv")

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
  Y = df_Henson[, -1], # remove 1st column
  gamma = 0.03
)

stan_model <- cmdstan_model("biaced_die_model.stan")

fit_vb <- stan_model$variational(
  data = stan_data,
  seed = 1234
)

fit_vb$cmdstan_diagnose()

# extract
fit_vb$summary(c("mu", "sigma", "lamda")) %>% 
  select(variable, mean, q5, q95) %>%
  extract(variable, 
          c("name", "tt"), "([a-z]+)\\[([0-9]+)", 
          convert = TRUE) %>% 
  mutate(tt = tt - 1) -> res_Henson

res_Henson

# visualize results
library(cowplot)

res_Henson %>% 
  filter(name == "mu") %>% 
  mutate(ground_truth = tt) %>% 
  ggplot(aes(x = tt, y = mean)) +
  geom_ribbon(aes(ymin = q5, ymax = q95), 
              alpha = 0.2, 
              fill = "blue") +
  geom_line() +
  geom_line(aes(y = ground_truth), linetype = "dashed") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "",
       y = "Luminance of \n inflection point (dB)") +
  theme_linedraw() +
  theme(panel.grid.minor = element_blank()) -> p_mu

res_Henson %>% 
  filter(name == "sigma") %>% 
  mutate(ground_truth = map_dbl(tt, ~ min(exp(-0.081 * . + 3.27), 6))) %>% 
  ggplot(aes(x = tt, y = mean)) +
  geom_ribbon(aes(ymin = q5, ymax = q95), 
              alpha = 0.2, 
              fill = "blue") +
  geom_line() +
  geom_line(aes(y = ground_truth), linetype = "dashed") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 9), 
                     expand = c(0, 0),
                     breaks = c(2, 4, 6, 8)) +
  labs(x = "Perimetric Sensitivity (dB)",
       y = "Slope \n (standard deviation)") +
  theme_linedraw() +
  theme(panel.grid.minor = element_blank()) -> p_sigma

res_Henson %>% 
  filter(name == "lamda") %>% 
  mutate(ground_truth = fnr) %>%
  ggplot(aes(x = tt, y = mean)) +
  geom_ribbon(aes(ymin = q5, ymax = q95), 
              alpha = 0.2, 
              fill = "blue") +
  geom_line() +
  geom_line(aes(y = ground_truth), linetype = "dashed") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0),
                     labels = scales::percent) +
  labs(x = "",
       y = "Lapse rate") +
  theme_linedraw() +
  theme(panel.grid.minor = element_blank()) -> p_lamda

plot_grid(p_mu, p_sigma, p_lamda, 
          labels = "AUTO", axis = "b", align = "h", nrow = 1)

