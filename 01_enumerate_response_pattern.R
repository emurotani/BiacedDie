library(tidyverse)
library(OPI)

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

# test FT2 function
crossing(tt = 0, n_times = 1:10) %>% 
  mutate(FT_res = map(tt, ~ FT2(tt = .x, makeStim = makeStim))) %>% 
  mutate(final = map_dbl(FT_res, ~ pluck(.x, "final")))

# copy from 00_explore_OPI_package.R
flatten_respSeq <- function(FT_res) {
  FT_res$respSeq %>% 
    bind_rows() %>% 
    mutate(db_seen = paste0(db, "(", seen, ")")) %>% 
    pull(db_seen) %>% 
    reduce(paste, sep = "-")
}

# to detect rare response patterns
fpr <- 0.30
fnr <- 0.30

# it takes a few minutes
crossing(tt = 0:40, n_times = 1:1000) %>% 
  mutate(FT_res = map(tt, ~ FT2(tt = .x, fpr = fpr, fnr = fnr,
                                makeStim = makeStim))) %>% 
  mutate(npres = map_dbl(FT_res, ~ pluck(.x, "npres")),
         first = map_dbl(FT_res, ~ pluck(.x, "first")),
         final = map_dbl(FT_res, ~ pluck(.x, "final"))) %>% 
  mutate(flat_respSeq = map_chr(FT_res, ~ flatten_respSeq(.x))) %>% 
  mutate(raw_respSeq = map(FT_res, ~ pluck(.x, "respSeq") %>% 
                          bind_rows)) -> df_enum

df_enum

# number of distinct patterns
df_enum %>% 
  distinct(flat_respSeq) %>% 
  count()

# extract distinct patterns
df_enum %>% 
  distinct(flat_respSeq, .keep_all = TRUE) %>%
  select(final, npres, raw_respSeq) %>% 
  arrange(final) %>% 
  mutate(seq_id = row_number(), .before = final) %>% 
  mutate(stop_row = cumsum(npres),
         start_row = stop_row - npres + 1) -> df_enum_distinct

df_enum_distinct

# divide into sequence and response
df_enum_distinct %>% 
  select(seq_id, final, start_row, stop_row) -> df_sequence

df_sequence

df_enum_distinct %>% 
  select(raw_respSeq) %>% 
  unnest(cols = c(raw_respSeq)) -> df_response

df_response

df_sequence %>% 
  ggplot(aes(x = final)) +
  geom_histogram(binwidth = 1) +
  labs(x = "Final estimate (dB)",
       y = "Counts of sequences")

# save
write_excel_csv(df_sequence, file = "results/df_sequence.csv")
write_excel_csv(df_response, file = "results/df_response.csv")

# reset environment
rm(list = ls())
