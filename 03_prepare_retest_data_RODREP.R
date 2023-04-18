library(tidyverse)
library(broom)

# IOVS
# 130 patients x both eyes x 15 series of visual field = 3900 records
bryan <- read_csv(file = "original_data/Bryan2013.csv")
bryan

# calculate follow up year
read_csv(file = "original_data/VisualFields.csv") %>% 
  select(-IOP) %>% 
  semi_join(bryan, by = "FIELD_ID") %>%
  group_by(STUDY_ID, SITE) %>% 
  mutate(FOLLOW_YR = (AGE - min(AGE)) / 365.25) %>% 
  ungroup() %>% 
  select(STUDY_ID, SITE, FOLLOW_YR, FIELD_ID) %>% 
  arrange(STUDY_ID, SITE, FOLLOW_YR) -> df_vf

df_vf

# choose data at primary points (X, Y) = (+/-9, +/-9)
crossing(X = c(9, -9), 
         Y = c(9, -9)) -> primary_points

primary_points

read_csv(file = "original_data/VFPoints.csv") %>% 
  select(-TOTAL_DEVIATION) %>% 
  semi_join(bryan, by = "FIELD_ID") %>%
  semi_join(primary_points, by = c("X", "Y")) %>% 
  left_join(df_vf, by = "FIELD_ID") %>% 
  group_nest(STUDY_ID, SITE, X, Y, .key = "DATA") %>% 
  mutate(POINT_ID = row_number()) %>% 
  select(POINT_ID, STUDY_ID, DATA) -> df_vfp

df_vfp
df_vfp$DATA[[1]] # THRESHOLD = measured sensitivity


# fit loess smoothing
df_vfp %>% 
  mutate(fit = map(DATA, ~ loess(THRESHOLD ~ FOLLOW_YR, 
                                 data = ., span = 1)),
         aug = map(fit, augment)) -> df_loess

df_loess
df_loess$aug[[1]]

df_loess %>% 
  sample_n(100) %>% 
  select(POINT_ID, aug) %>% 
  unnest(aug) %>% 
  ggplot(aes(x = FOLLOW_YR, group = POINT_ID)) +
  geom_line(aes(y = .fitted), alpha = 0.2)


# set true perimetric sensitivity(tt)
df_loess %>%
  select(POINT_ID, STUDY_ID, aug) %>%
  unnest(aug) %>%
  mutate(tt = round(.fitted)) -> df_retest

df_retest

# remove unreliable test points
df_retest %>% 
  select(POINT_ID, THRESHOLD, tt) %>%
  mutate(IS_TRIGGER_HAPPY = if_else(THRESHOLD > 40, TRUE, FALSE),
         OVER_35 = if_else(tt > 35, TRUE, FALSE),
         BELOW_0 = if_else(tt < 0, TRUE, FALSE)) %>% 
  filter(IS_TRIGGER_HAPPY | OVER_35 | BELOW_0) %>% 
  distinct(POINT_ID) -> unreliable_point
  
unreliable_point

df_retest %>% 
  semi_join(unreliable_point, by = "POINT_ID") %>% 
  ggplot(aes(x = FOLLOW_YR, y = .fitted, group = POINT_ID)) +
    geom_line(alpha = 0.2)


df_retest %>%
  anti_join(unreliable_point, by = "POINT_ID") -> df_retest

df_retest

# retest data summary
table(factor(df_retest$tt)) # true perimetric sensitivity
table(factor(df_retest$THRESHOLD, levels = -1:40)) # measured sensitivity

df_retest %>% 
  ggplot(aes(x = tt)) +
    geom_histogram(binwidth = 1)
  
df_retest %>%
  ggplot(aes(x = THRESHOLD)) +
  geom_histogram(aes(y = ..density..), binwidth = 1) +
  facet_wrap(~ tt, nrow = 4) +
  labs(title = "Retest distribution at primary points")

# save
df_retest %>% 
  write_excel_csv(file = "results/retest_data_RODREP.csv")

# clear environment
rm(list = ls())
