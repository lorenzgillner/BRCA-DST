library(tidyverse)

seed <- 23012026
set.seed(seed)

N <- 500000

n.sp  <- round(N * 0.05)

# Based on the 10-year average cancer incidence in Germany
# https://www.krebsdaten.de/Krebs/EN/Database/databasequery_step1_node.html
n.fbc <- round((N - n.sp) * 0.8991)
n.mbc <- round((N - n.sp) * 0.0087)
n.oc  <- round((N - n.sp) * 0.0922)

n.sp = n.sp + (N - (n.fbc + n.mbc + n.oc + n.sp))

if ((n.fbc + n.mbc + n.oc + n.sp) != N) {
  throw("Number of samples doesn't match the total number of cases!")
}

# Based on Buys et al., 2017
# https://doi.org/10.1002/cncr.30498
pv.fbc <- 0.05
pv.mbc <- 0.04
pv.oc  <- 0.08
pv.sp  <- 0.11


# Generator function ------------------------------------------------------

generate <- function(n,
                     cancer,
                     crude.prob,
                     sex,
                     age.range,
                     c.age.mean.sd,
                     nc.age.mean.sd,
                     brca.probs = c(0.55, 0.45)) {
  sex <- rep(sex, n)
  
  reset_stat_defaults()
  brca.status <- rbinom(n, 1, prob = crude.prob)
  brca.type <- ifelse(brca.status, sample(c(1, 2), n, replace = TRUE, prob = brca.probs), 0)
  brca1 <- brca.type == 1
  brca2 <- brca.type == 2
  
  c.age.mean <- first(c.age.mean.sd)
  c.age.sd <- last(c.age.mean.sd)
  c.age.dist <- sample(age.range, n, replace = TRUE, prob = dnorm(age.range, mean = c.age.mean, sd = c.age.sd)^0.9)
  
  nc.age.mean <- first(nc.age.mean.sd)
  nc.age.sd <- last(nc.age.mean.sd)
  nc.age.dist <- sample(age.range, n, replace = TRUE, prob = dnorm(age.range, mean = nc.age.mean, sd = nc.age.sd)^0.2)
  
  age <- ifelse(brca.status, c.age.dist, nc.age.dist)
  
  fh <- sample(c(TRUE, FALSE), n, replace = TRUE)
  
  data.frame(
    cancer = cancer,
    sex = sex,
    fh = fh,
    age.at.onset = age,
    brca1 = brca1,
    brca2 = brca2
  )
}


# Female breast cancer ----------------------------------------------------

df.fbc <- generate(n.fbc + n.sp, "BC", pv.fbc, "F", 16:85, c(42, 12), c(62, 12))
df.fbc$bc.tn <- ifelse((df.fbc$brca1|df.fbc$brca2),
                       sample(c(TRUE, FALSE), n.fbc + n.sp, replace = TRUE, prob = c(0.35, 0.65)),
                       sample(c(TRUE, FALSE), n.fbc + n.sp, replace = TRUE, prob = c(0.15, 0.85)))
df.fbc$bc.bil <- sample(c(TRUE, FALSE), n.fbc + n.sp, replace = TRUE, prob = c(0.22, 0.78))

slice.fbc <- slice_tail(df.fbc, n = n.sp)
df.fbc <- slice_head(df.fbc, n = n.fbc)


# Male breast cancer ------------------------------------------------------

df.mbc <- generate(n.mbc, "BC", pv.mbc, "M", 35:85, c(63, 9), c(66, 9))
df.mbc$bc.tn <- rep(NA, n.mbc)
df.mbc$bc.bil <- rep(NA, n.mbc)


# Ovarian cancer ----------------------------------------------------------

df.oc <- generate(n.oc + n.sp, "OC", pv.oc, "F", 16:85, c(44, 11), c(63, 11), c(0.55, 0.45))

slice.oc <- slice_tail(df.oc, n = n.sp)
df.oc <- slice_head(df.oc, n = n.oc)

df.all <- bind_rows(df.fbc, df.mbc) |>
  dplyr::select(sex, fh, brca1, brca2, age.at.onset, bc.bil, bc.tn)|>
  dplyr::rename(bc.age = age.at.onset)


# Merge all data frames ---------------------------------------------------

df.oc <- dplyr::select(df.oc, !cancer) |>  dplyr::rename(oc.age = age.at.onset)

df.all <- bind_rows(df.all, df.oc)

slice.fbc <- slice.fbc |>
  dplyr::rename(bc.age = age.at.onset, bc.brca1 = brca1, bc.brca2 = brca2) |>
  dplyr::select(sex, fh, bc.brca1, bc.brca2, bc.age, bc.bil, bc.tn)

slice.oc <- slice.oc |>
  dplyr::rename(oc.age = age.at.onset, oc.brca1 = brca1, oc.brca2 = brca2) |>
  dplyr::select(oc.brca1, oc.brca2, oc.age)

df.sp <- bind_cols(slice.fbc, slice.oc) |>
  dplyr::mutate(brca1 = bc.brca1 | oc.brca1, .after = sex, .keep = "unused") |>
  dplyr::mutate(brca2 = bc.brca2 | oc.brca2, .after = brca1, .keep = "unused")

df.all <- bind_rows(df.all, df.sp)

df.all <- sample_n(df.all, N) |> mutate(ID = 1:N, .before = sex)


# Write the final CSV file ------------------------------------------------

write_csv(df.all, "cohort.csv")
