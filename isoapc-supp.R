### Purpose: Conduct APC supplementary analyses of time spent alone
### Author:  S Bauldry 
### Date:    Sep 14, 2024

### Set working directory and load libraries
setwd("~/desktop")
library(tidyverse)
library(ggpubr)
library(weights)

### Read helper functions
source("isoapc-functions.R")

### Read prepared ATUS data
atus <- read_csv("isoapc-data.csv", 
                 col_types = list(a = "f", p = "f", c = "f", 
                                  fem = "f", rce = "f", edu = "f", 
                                  day = "f", hol = "f")) |>
  mutate(
    a = fct_relevel(a, "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", 
                    "12", "13"),
    p = fct_relevel(p, "1", "2", "3", "4"),
    c = fct_relevel(c, "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", 
                    "12", "13", "14", "15", "16"))

# select stratified samples
atus_nla <- atus |> filter(lal == 0)
atus_lal <- atus |> filter(lal == 1)

atus_nft <- atus |> filter(fte == 0)
atus_fte <- atus |> filter(fte == 1)


### Table S1 - descriptive statistics
wtd.mean(atus$nonwork_alone_min, weights = atus$rwt)
wpct(atus$lal, weight = atus$rwt)
wpct(atus$fte, weight = atus$rwt)

old <- options(pillar.sigfig = 4)
lal_means <- atus |>
  group_by(lal) |>
  summarise(n = n(),
            wm = sum(nonwork_alone_min*rwt)/sum(rwt))
lal_means

emp_means <- atus |>
  group_by(fte) |>
  summarise(n = n(),
            wm = sum(nonwork_alone_min*rwt)/sum(rwt))
emp_means


### Figure S1 - solution line
theta_1 <- apc_model(atus, "nonwork_alone_min")
figs1 <- fig2D(theta_1[1], theta_1[2], "")
ggsave("figS1.jpg", figs1)


### Figure S2 - solution line with additional covariates
theta_2 <- apc_model_cov(atus, "nonwork_alone_min")
figs2 <- fig2D(theta_2[1], theta_2[2], "")
ggsave("figS2.jpg", figs2)


### Figure S3 - solution line
theta_3 <- apc_model(atus_lal, "nonwork_alone_min")
theta_3
figs3 <- fig2D(theta_3[1], theta_3[2], "")
ggsave("figS3.jpg", figs3)


### Figure S4 - solution line
theta_4 <- apc_model(atus_nla, "nonwork_alone_min")
theta_4
figs4 <- fig2D(theta_4[1], theta_4[2], "")
ggsave("figS4.jpg", figs4)


### Figure S5 - solution line
theta_5 <- apc_model(atus_nft, "nonwork_alone_min")
theta_5
figs5 <- fig2D(theta_5[1], theta_5[2], "")
ggsave("figS5.jpg", figs5)


### Figure S6 - solution line
theta_6 <- apc_model(atus_fte, "nonwork_alone_min")
theta_6
figs6 <- fig2D(theta_6[1], theta_6[2], "")
ggsave("figS6.jpg", figs6)


### Figure S7 - Lexis diagrams by living alone
fs7_lal <- figLexis(atus_lal, "nonwork_alone_min", 175, 700, "living alone")
fs7_nla <- figLexis(atus_nla, "nonwork_alone_min", 175, 700, "not living alone")
figs7 <- ggarrange(fs7_lal, fs7_nla, nrow = 1)
ggsave("figS7.jpg", figs7)


### Figure S8 - Lexis diagrams by working full time
fs8_fte <- figLexis(atus_fte, "nonwork_alone_min", 175, 500, "working full time")
fs8_nft <- figLexis(atus_nft, "nonwork_alone_min", 175, 500, "not working full time")
figs8 <- ggarrange(fs8_fte, fs8_nft, nrow = 1)
ggsave("figS8.jpg", figs8)


### Figure S9 - total effects including additional covariates
ub <- 650
te_all   <- total_effects_range_cov(atus, "nonwork_alone_min")
te_all_a <- te_age_fig(subset(te_all, apc == "age"), 0, ub)
te_all_p <- te_per_fig(subset(te_all, apc == "year"), 0, ub)
te_all_c <- te_coh_fig(subset(te_all, apc == "cohort"), 0, ub)
figs9    <- ggarrange(te_all_a, te_all_p, te_all_c, nrow = 1)
figs9
ggsave("figS9.jpg", figs9)


### Figure S10 - Lexis diagram for leisure activities
fs10 <- figLexis(atus_lal, "leisure_alone_min", 50, 450, "")
ggsave("figS10.jpg", fs10)


### Figure S11 - solution line
theta_11 <- apc_model(atus, "leisure_alone_min")
theta_11
figs11 <- fig2D(theta_11[1], theta_11[2], "")
ggsave("figS11.jpg", figs11)


### Figure S12 - total effects leisure activities
ub <- 450
te_all   <- total_effects_range_cov(atus, "leisure_alone_min")
te_all_a <- te_age_fig(subset(te_all, apc == "age"), 0, ub)
te_all_p <- te_per_fig(subset(te_all, apc == "year"), 0, ub)
te_all_c <- te_coh_fig(subset(te_all, apc == "cohort"), 0, ub)
figs12   <- ggarrange(te_all_a, te_all_p, te_all_c, nrow = 1)
figs12
ggsave("figS12.jpg", figs12)




