### Purpose: Conduct APC supplementary analyses of time spent alone
### Author:  S Bauldry 
### Date:    Aug 14, 2025

setwd("~/desktop")
library(tidyverse)
library(weights)
library(patchwork)


### Read helper functions
source("isoapc-functions.R")


### Read prepared ATUS data
atus <- read_csv("isoapc-data.csv", col_types = list(a = "f", p = "f", c = "f", day = "f", hol = "f")) |>
  mutate(
    a = fct_relevel(a, "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"),
    p = fct_relevel(p, "1", "2", "3", "4"),
    c = fct_relevel(c, "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16"),
    wke = ifelse(day %in% c(1,7), 1, 0))


### Supplementary analysis 1 -- stratify by weekday vs weekend
atus_day <- atus |> filter(wke == 0)
atus_end <- atus |> filter(wke == 1)

cns_day <- apc_model_mmmapg(atus_day, "nonwork_alone_min")
cns_end <- apc_model_mmmapg(atus_end, "nonwork_alone_min")

te_day <- total_effects_range(atus_day, "nonwork_alone_min", 5, cns_day$alpha[1], cns_day$alpha[2], cns_day$alpha[3])
te_end <- total_effects_range(atus_end, "nonwork_alone_min", 5, cns_end$alpha[1], cns_end$alpha[2], cns_end$alpha[3])

sf1_day_a <- te_fig(subset(te_day, apc == "age"), 100, 600, "Weekday", xn = "age", c("15", "75"), c("15", "79"))
sf1_day_p <- te_fig(subset(te_day, apc == "year"), 100, 600, "", xn = "year", c("2003", "2018"), c("2003", "2022"))
sf1_day_c <- te_fig(subset(te_day, apc == "cohort"), 100, 600, "", xn = "cohort", c("1924", "1999"), c("1924", "2007"))

sf1_end_a <- te_fig(subset(te_end, apc == "age"), 100, 600, "Weekend", xn = "age", c("15", "75"), c("15", "79"))
sf1_end_p <- te_fig(subset(te_end, apc == "year"), 100, 600, "", xn = "year", c("2003", "2018"), c("2003", "2022"))
sf1_end_c <- te_fig(subset(te_end, apc == "cohort"), 100, 600, "", xn = "cohort", c("1924", "1999"), c("1924", "2007"))

sf1 <- sf1_day_a + sf1_day_p + sf1_day_c + sf1_end_a + sf1_end_p + sf1_end_c
sf1
ggsave("sf1.jpg", sf1, dpi = 300)



### Supplementary analysis 2 -- leisure activities
atus_w <- atus |> filter(fem == 1)
atus_m <- atus |> filter(fem == 0)

cns_w <- apc_model_mmmapg(atus_w, "leisure_alone_min")
cns_m <- apc_model_mmmapg(atus_m, "leisure_alone_min")

te_w <- total_effects_range(atus_w, "leisure_alone_min", 5, cns_w$alpha[1], cns_w$alpha[2], cns_w$alpha[3])
te_m <- total_effects_range(atus_m, "leisure_alone_min", 5, cns_m$alpha[1], cns_m$alpha[2], cns_m$alpha[3])

sf2_w_a <- te_fig(subset(te_w, apc == "age"), 0, 300, "Women", xn = "age", c("15", "75"), c("15", "79"))
sf2_w_p <- te_fig(subset(te_w, apc == "year"), 0, 300, "", xn = "year", c("2003", "2018"), c("2003", "2022"))
sf2_w_c <- te_fig(subset(te_w, apc == "cohort"), 0, 300, "", xn = "cohort", c("1924", "1999"), c("1924", "2007"))

sf2_m_a <- te_fig(subset(te_m, apc == "age"), 0, 300, "Men", xn = "age", c("15", "75"), c("15", "79"))
sf2_m_p <- te_fig(subset(te_m, apc == "year"), 0, 300, "", xn = "year", c("2003", "2018"), c("2003", "2022"))
sf2_m_c <- te_fig(subset(te_m, apc == "cohort"), 0, 300, "", xn = "cohort", c("1924", "1999"), c("1924", "2007"))

sf2 <- sf2_w_a + sf2_w_p + sf2_w_c + sf2_m_a + sf2_m_p + sf2_m_c
sf2
ggsave("sf2.jpg", sf2, dpi = 300)



### Supplementary analysis 3 -- check exclude post-pandemic period
atus_wp <- atus |> filter(fem == 1 & year < 2020)
atus_mp <- atus |> filter(fem == 0 & year < 2020)

cns_wp <- apc_model_mmmapg(atus_wp, "nonwork_alone_min")
cns_mp <- apc_model_mmmapg(atus_mp, "nonwork_alone_min")

te_wp <- total_effects_range(atus_wp, "nonwork_alone_min", 5, cns_wp$alpha[1], cns_wp$alpha[2], cns_wp$alpha[3])
te_mp <- total_effects_range(atus_mp, "nonwork_alone_min", 5, cns_mp$alpha[1], cns_mp$alpha[2], cns_mp$alpha[3])

sf3_wp_a <- te_fig(subset(te_wp, apc == "age"), 0, 600, "Women", xn = "age", c("15", "75"), c("15", "79"))
sf3_wp_p <- te_fig(subset(te_wp, apc == "year"), 0, 600, "", xn = "year", c("2003", "2018"), c("2003", "2019"))
sf3_wp_c <- te_fig(subset(te_wp, apc == "cohort"), 0, 600, "", xn = "cohort", c("1924", "1999"), c("1924", "2007"))

sf3_mp_a <- te_fig(subset(te_mp, apc == "age"), 0, 600, "Men", xn = "age", c("15", "75"), c("15", "79"))
sf3_mp_p <- te_fig(subset(te_mp, apc == "year"), 0, 600, "", xn = "year", c("2003", "2018"), c("2003", "2019"))
sf3_mp_c <- te_fig(subset(te_mp, apc == "cohort"), 0, 600, "", xn = "cohort", c("1924", "1999"), c("1924", "2007"))

sf3 <- sf3_wp_a + sf3_wp_p + sf3_wp_c + sf3_mp_a + sf3_mp_p + sf3_mp_c
sf3
ggsave("sf3.jpg", sf3, dpi = 300)

