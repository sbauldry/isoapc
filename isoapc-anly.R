### Purpose: Conduct APC analysis of time spent alone
### Author:  S Bauldry 
### Date:    Oct 17, 2024

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
    c = fct_relevel(c, "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16"))

# select subsamples for women and men
atus_w <- atus |> filter(fem == 1)
atus_m <- atus |> filter(fem == 0)


### Figure 1: Univariate APC trends
f1_w_a <- figUniAPC(atus_w, "age", "nonwork_alone_min", 100, 600, "Women", c(15, 79), c("15", "79"))
f1_w_p <- figUniAPC(atus_w, "year", "nonwork_alone_min", 100, 600, "", c(2003, 2022), c("2003", "2022"))
f1_w_c <- figUniAPC(atus_w, "cohort", "nonwork_alone_min", 100, 600, "", c(1924, 2007), c("1924", "2007"))

f1_m_a <- figUniAPC(atus_m, "age", "nonwork_alone_min", 100, 600, "Men", c(15, 79), c("15", "79"))
f1_m_p <- figUniAPC(atus_m, "year", "nonwork_alone_min", 100, 600, "", c(2003, 2022), c("2003", "2022"))
f1_m_c <- figUniAPC(atus_m, "cohort", "nonwork_alone_min", 100, 600, "", c(1924, 2007), c("1924", "2007"))

f1 <- f1_w_a + f1_w_p + f1_w_c + f1_m_a + f1_m_p + f1_m_c + plot_layout(ncol = 3)
f1
ggsave("fig1.tif", f1, height = 5, width = 7.5, units = "in", dpi = 300, compression = "lzw")


### Figure 2 - Canonical solution line with constrained region noted
cns_w <- apc_model_mmmapg(atus_w, "nonwork_alone_min")
cns_m <- apc_model_mmmapg(atus_m, "nonwork_alone_min")

# report min and max for figure notes
c( cns_w$alpha[1], cns_w$alpha[3], cns_w$pi[1], cns_w$pi[3], cns_w$gamma[1], cns_w$gamma[3] )
c( cns_m$alpha[1], cns_m$alpha[3], cns_m$pi[1], cns_m$pi[3], cns_m$gamma[1], cns_m$gamma[3] )

f2_w <- fig2D(cns_w$theta[1], cns_w$theta[2], -3, 6, cns_w$alpha[1], cns_w$alpha[2], cns_w$alpha[3], "Women")
f2_m <- fig2D(cns_m$theta[1], cns_m$theta[2], -3, 6, cns_m$alpha[1], cns_m$alpha[2], cns_m$alpha[3], "Men")

f2 <- f2_w + f2_m 
f2
ggsave("fig2.tif", f2, height = 5, width = 7.5, units = "in", dpi = 300, compression = "lzw")



### Figure 3 - Range of total effects based on constrained region of canonical solution line
te_w <- total_effects_range(atus_w, "nonwork_alone_min", 5, cns_w$alpha[1], cns_w$alpha[2], cns_w$alpha[3])
te_m <- total_effects_range(atus_m, "nonwork_alone_min", 5, cns_m$alpha[1], cns_m$alpha[2], cns_m$alpha[3])

f3_w_a <- te_fig(subset(te_w, apc == "age"), 100, 500, "Women", xn = "age", c("15", "75"), c("15", "79"))
f3_w_p <- te_fig(subset(te_w, apc == "year"), 100, 500, "", xn = "year", c("2003", "2018"), c("2003", "2022"))
f3_w_c <- te_fig(subset(te_w, apc == "cohort"), 100, 500, "", xn = "cohort", c("1924", "1999"), c("1924", "2007"))

f3_m_a <- te_fig(subset(te_m, apc == "age"), 100, 500, "Men", xn = "age", c("15", "75"), c("15", "79"))
f3_m_p <- te_fig(subset(te_m, apc == "year"), 100, 500, "", xn = "year", c("2003", "2018"), c("2003", "2022"))
f3_m_c <- te_fig(subset(te_m, apc == "cohort"), 100, 500, "", xn = "cohort", c("1924", "1999"), c("1924", "2007"))

f3 <- f3_w_a + f3_w_p + f3_w_c + f3_m_a + f3_m_p + f3_m_c
f3
ggsave("fig3.tif", f3, height = 5, width = 7.5, units = "in", dpi = 300, compression = "lzw")
