library(tidyverse)
library(OPI)

set.seed(1234)

# initialize
chooseOpi("SimHenson")
opiInitialise(type = "C", cap = 6)

makeStim <- function(db, n) { 
  s <- list(level=dbTocd(db))
  class(s) <- "opiStaticStimulus"
  return(s)
}

# full threshold: a run
# tt = true perimetric sensitivity
FT_res <- FT(est = 20, tt = 20, 
             fpr = 0.03, fnr = 0.01, makeStim = makeStim)

FT_res
str(FT_res) # list type

# full threshold: repeat
crossing(tt = 10:13, n_times = 1:3) %>% 
  mutate(FT_res = map(tt, ~ FT(tt = .x, makeStim = makeStim))) %>% 
  mutate(final = map_dbl(FT_res, ~ pluck(.x, "final")))

# convert test-retest data into a matrix
crossing(tt = 10:13, n_times = 1:20) %>% 
  mutate(FT_res = map(tt, ~ FT(tt = .x, makeStim = makeStim))) %>% 
  mutate(final = map_dbl(FT_res, ~ pluck(.x, "final"))) %>% 
  group_by(tt, final) %>% 
  count() %>% 
  pivot_wider(names_from = final, 
              values_from = n, values_fill = 0) %>% 
  as.matrix()

# response pattern
# seen = patient's response
# db = stimulus luminance
FT_res$respSeq

# response pattern in one row
FT_res$respSeq %>% 
  bind_rows()

FT_res$respSeq %>% 
  bind_rows() %>% 
  mutate(db_seen = paste0(db, "(", seen, ")")) %>% 
  pull(db_seen) %>%
  reduce(paste, sep = "-")

flatten_respSeq <- function(FT_res) {
  FT_res$respSeq %>% 
    bind_rows() %>% 
    mutate(db_seen = paste0(db, "(", seen, ")")) %>% 
    pull(db_seen) %>% 
    reduce(paste, sep = "-")
}

flatten_respSeq(FT_res)

# convert FT_res to a tidy form
crossing(tt = 10:13, n_times = 1:2) %>% 
  mutate(FT_res = map(tt, ~ FT(tt = .x, makeStim = makeStim))) %>% 
  mutate(npres = map_dbl(FT_res, ~ pluck(.x, "npres")),
         first = map_dbl(FT_res, ~ pluck(.x, "first")),
         final = map_dbl(FT_res, ~ pluck(.x, "final"))) %>% 
  mutate(flat_respSeq = map_chr(FT_res, ~ flatten_respSeq(.x))) %>% 
  mutate(raw_respSeq = map(FT_res, ~ pluck(.x, "respSeq") %>% 
                          bind_rows)) -> FT_res_tidy

FT_res_tidy
FT_res_tidy$raw_respSeq[[1]]

# reset environment
rm(list = ls())
