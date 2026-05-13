library(tidyverse)

cohort <- read_csv("cohort.csv") |>
  replace_na(list(
    bc.age = 0,
    bc.bil = FALSE,
    bc.tn = FALSE,
    oc.age = 0
  )) |>
  mutate(bc.m = (bc.age > 0) & sex == "M") |>
  mutate(sp = (bc.age > 0) & (oc.age > 0)) |>
  mutate(brca = (brca1|brca2), .keep = "unused")

train.ratio <- 0.75
test.ratio  <- 1.0 - train.ratio

train <- slice_head(cohort, prop = train.ratio)
test  <- slice_tail(cohort, prop = test.ratio)

write_csv(train, "train.csv")
write_csv(test, "test.csv")
