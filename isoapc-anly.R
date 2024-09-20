### Purpose: Conduct primary APC analysis of time spent alone
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


### Figure 1 - Univariate trend in time spent alone over time
f1_all_p <- figUniAPC(atus, "year", "nonwork_alone_min", 100, 500)
f1_all_p <- f1_all_p + scale_x_continuous(name = "year", 
                                      breaks = c(2003, 2022), 
                                      labels = c(2003, 2022))

f1_all_a <- figUniAPC(atus, "age", "nonwork_alone_min", 100, 500)
f1_all_a <- f1_all_a + scale_x_continuous(name = "age", 
                                        breaks = c(15,79), 
                                        labels = c(15,79))

f1_all_c <- figUniAPC(atus, "cohort", "nonwork_alone_min", 100, 500)
f1_all_c <- f1_all_c + scale_x_continuous(name = "cohort", 
                                          breaks = c(1924, 2007), 
                                          labels = c(1924, 2007))

f1_all <- ggarrange(f1_all_a, f1_all_p, f1_all_c, nrow = 1)
f1_all
ggsave("f1_all.jpg", f1_all)


### Figure 2 - Lexis diagram in overall sample
f2_all <- figLexis(atus, "nonwork_alone_min", 175, 525, "")
f2_all
ggsave("f2_all.jpg", f2_all)


### Figure 3 - Total net APC effects in overall sample
ub <- 600
te_all   <- total_effects_range(atus, "nonwork_alone_min")
te_all_a <- te_age_fig(subset(te_all, apc == "age"), 0, ub)
te_all_p <- te_per_fig(subset(te_all, apc == "year"), 0, ub)
te_all_c <- te_coh_fig(subset(te_all, apc == "cohort"), 0, ub)
f3_all   <- ggarrange(te_all_a, te_all_p, te_all_c, nrow = 1)
f3_all
ggsave("f3_all.jpg", f3_all)


### Figure 4 - Total net APC effects and marginal effects for living alone
ub <- 900
te_nla   <- total_effects_range(atus_nla, "nonwork_alone_min")
te_nla_a <- te_age_fig(subset(te_nla, apc == "age"), 0, ub)
te_nla_p <- te_per_fig(subset(te_nla, apc == "year"), 0, ub)
te_nla_c <- te_coh_fig(subset(te_nla, apc == "cohort"), 0, ub)
f4_nla   <- ggarrange(te_nla_a, te_nla_p, te_nla_c, nrow = 1)
f4_nla   <- annotate_figure(f4_nla, top = text_grob("Not living alone", size = 16, face = "bold"))

te_lal   <- total_effects_range(atus_lal, "nonwork_alone_min")
te_lal_a <- te_age_fig(subset(te_lal, apc == "age"), 0, ub)
te_lal_p <- te_per_fig(subset(te_lal, apc == "year"), 0, ub)
te_lal_c <- te_coh_fig(subset(te_lal, apc == "cohort"), 0, ub)
f4_lal   <- ggarrange(te_lal_a, te_lal_p, te_lal_c, nrow = 1)
f4_lal   <- annotate_figure(f4_lal, top = text_grob("Living alone", size = 16, face = "bold"))

mte_lal <- te_lal |>
  left_join(te_nla, by = c("apc", "index")) |>
  mutate(dmid = mid.x - mid.y,
         dmin = cmin.x - cmin.y,
         dmax = cmax.x - cmax.y) |>
  select(apc, index, dmid, dmin, dmax)
mte_lal_a <- mte_age_fig(subset(mte_lal, apc == "age"), -200, 500)
mte_lal_p <- mte_per_fig(subset(mte_lal, apc == "year"), -200, 500)
mte_lal_c <- mte_coh_fig(subset(mte_lal, apc == "cohort"), -200, 500)

f4_mlal   <- ggarrange(mte_lal_a, mte_lal_p, mte_lal_c, nrow = 1)
f4_mlal   <- annotate_figure(f4_mlal, top = text_grob("Difference between living alone and not living alone", size = 16, face = "bold"))

f4_alo <- ggarrange(f4_nla, f4_lal, f4_mlal, nrow = 3)
ggsave("f4_alo.jpg", f4_alo, width = 18, height = 18)




### Figure 5 - Total net APC effects and marginal effects for working full-time
ub <- 700
te_nft   <- total_effects_range(atus_nft, "nonwork_alone_min")
te_nft_a <- te_age_fig(subset(te_nft, apc == "age"), 0, ub)
te_nft_p <- te_per_fig(subset(te_nft, apc == "year"), 0, ub)
te_nft_c <- te_coh_fig(subset(te_nft, apc == "cohort"), 0, ub)
f5_nft   <- ggarrange(te_nft_a, te_nft_p, te_nft_c, nrow = 1)
f5_nft   <- annotate_figure(f5_nft, top = text_grob("Not working full time", size = 16, face = "bold"))

te_fte   <- total_effects_range(atus_fte, "nonwork_alone_min")
te_fte_a <- te_age_fig(subset(te_fte, apc == "age"), 0, ub)
te_fte_p <- te_per_fig(subset(te_fte, apc == "year"), 0, ub)
te_fte_c <- te_coh_fig(subset(te_fte, apc == "cohort"), 0, ub)
f5_fte   <- ggarrange(te_fte_a, te_fte_p, te_fte_c, nrow = 1)
f5_fte   <- annotate_figure(f5_fte, top = text_grob("Working full time", size = 16, face = "bold"))

mte_nft <- te_nft |>
  left_join(te_fte, by = c("apc", "index")) |>
  mutate(dmid = mid.x - mid.y,
         dmin = cmin.x - cmin.y,
         dmax = cmax.x - cmax.y) |>
  select(apc, index, dmid, dmin, dmax)
mte_nft_a <- mte_age_fig(subset(mte_nft, apc == "age"), -100, 300)
mte_nft_p <- mte_per_fig(subset(mte_nft, apc == "year"), -100, 300)
mte_nft_c <- mte_coh_fig(subset(mte_nft, apc == "cohort"), -100, 300)

f5_mnft   <- ggarrange(mte_nft_a, mte_nft_p, mte_nft_c, nrow = 1)
f5_mnft   <- annotate_figure(f5_mnft, top = text_grob("Difference between not working full time and working full time", size = 16, face = "bold"))

f5_emp <- ggarrange(f5_nft, f5_fte, f5_mnft, nrow = 3)
ggsave("f5_emp.jpg", f5_emp, width = 18, height = 18)


