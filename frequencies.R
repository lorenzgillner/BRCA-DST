library(tidyverse)

cohort <- read_csv("cohort.csv", show_col_types = FALSE) |>
  replace_na(list(
    bc.age = 0,
    bc.bil = FALSE,
    bc.tn = FALSE,
    oc.age = 0
  )) |>
  mutate(bc.m = (bc.age > 0) & sex == "M") |>
  mutate(sp = (bc.age > 0) & (oc.age > 0)) |>
  mutate(brca = (brca1|brca2), .keep = "unused")

age.bins <- c(0, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85)
age.labels <- c("≤20", "21-25", "26-30", "31-35", "36-40", "41-45", "46-50",
                "51-55", "56-60", "61-65", "66-70", "71-75", "76-80", ">80")

group_df <- function (df) {
  df |> group_by(age.group, .drop = FALSE) |> summarize(n=n())
}

summarize_df <- function(df.total, df.positive) {
  df.total.grouped <- group_df(df.total)
  df.positive.grouped <- group_df(df.positive)
  # df.ratio.grouped <- df.positive.grouped$n/df.total.grouped$n |> replace_na(0)
  # df.reversed.cumulative <- rev(cumsum(rev(df.positive.grouped$n))/nrow(df.total))
  
  return(data.frame(age.group = age.labels,
              positive = df.positive.grouped$n,
              total = df.total.grouped$n))
}

count_cases <- function(cohort) {
  # Female breast cancer ----------------------------------------------------
  
  fbc.total <- cohort |>
    filter((bc.age > 0) & (oc.age == 0) & sex == "F") |>
    rename(age = bc.age) |>
    mutate(age.group = cut(age, breaks = age.bins, labels = age.labels))
  fbc.positive <- fbc.total |> filter(brca)
  
  fbc.summary <- summarize_df(fbc.total, fbc.positive)
  
  
  # Bilateral breast cancer -------------------------------------------------
  
  bbc.total <- fbc.total |> filter(bc.bil)
  bbc.positive <- bbc.total |> filter(brca)
  
  bbc.summary <- summarize_df(bbc.total, bbc.positive)
  
  
  # Triple negative breast cancer -------------------------------------------
  
  tnbc.total <- fbc.total |> filter(bc.tn)
  tnbc.positive <- tnbc.total |> filter(brca)
  
  tnbc.summary <- summarize_df(tnbc.total, tnbc.positive)
  
  
  # Male breast cancer ------------------------------------------------------
  
  mbc.total <- cohort |>
    filter((bc.age > 0) & (oc.age == 0) & sex == "M") |>
    rename(age = bc.age) |>
    mutate(age.group = cut(age, breaks = age.bins, labels = age.labels))
  mbc.positive <- mbc.total |> filter(brca)
  
  mbc.summary <- summarize_df(mbc.total, mbc.positive)
  
  
  # Ovarian cancer ----------------------------------------------------------
  
  oc.total <- cohort |>
    filter((oc.age > 0) & (bc.age == 0) & sex == "F") |>
    rename(age = oc.age) |>
    mutate(age.group = cut(age, breaks = age.bins, labels = age.labels))
  oc.positive <- oc.total |> filter(brca)
  
  oc.summary <- summarize_df(oc.total, oc.positive)
  
  
  # Breast and ovarian cancer in the same person ----------------------------
  
  sp.total <- cohort |>
    filter((bc.age > 0) & (oc.age > 0) & sex == "F") |>
    mutate(age.group = cut(pmax(bc.age, oc.age), breaks = age.bins, labels = age.labels))
  sp.positive <- sp.total |> filter(brca)
  
  sp.summary <- summarize_df(sp.total, sp.positive)
  
  
  # Combine all frequencies into one single data frame ----------------------
  
  return(data.frame(
    age.group = age.labels,
    fbc = fbc.summary |> select(!age.group),
    bbc = bbc.summary |> select(!age.group),
    tnbc = tnbc.summary |> select(!age.group),
    mbc = mbc.summary |> select(!age.group),
    oc = oc.summary |> select(!age.group),
    sp = sp.summary |> select(!age.group)))
}

frequencies <- count_cases(cohort)
write_csv(frequencies, "frequencies.csv")
