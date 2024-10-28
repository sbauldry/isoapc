### Purpose: Conduct APC analysis of time spent alone
### Author:  S Bauldry 
### Date:    Oct 17, 2024

### Set working directory and load libraries
setwd("~/desktop")
library(tidyverse)
library(ggpubr)
library(weights)

### Read helper functions
source("isoapc-functions.R")

### Read prepared ATUS data
atus <- read_csv("isoapc-data.csv", col_types = list(a = "f", p = "f", c = "f", day = "f", hol = "f")) |>
  mutate(
    a = fct_relevel(a, "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"),
    p = fct_relevel(p, "1", "2", "3", "4"),
    c = fct_relevel(c, "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16"))



### Figure 1 - Univariate trend in time spent alone over time
f1_a <- figUniAPC(atus, "age", "nonwork_alone_min", 100, 500, "Age")
f1_a <- f1_a + scale_x_continuous(name = "age", breaks = c(15,79), labels = c(15,79))

f1_p <- figUniAPC(atus, "year", "nonwork_alone_min", 100, 500, "Period")
f1_p <- f1_p + scale_x_continuous(name = "year", breaks = c(2003, 2022), labels = c(2003, 2022))

f1_c <- figUniAPC(atus, "cohort", "nonwork_alone_min", 100, 500, "Cohort")
f1_c <- f1_c + scale_x_continuous(name = "cohort", breaks = c(1924, 2007), labels = c(1924, 2007))

f1 <- ggarrange(f1_a, f1_p, f1_c, nrow = 1)
f1
ggsave("f1.jpg", f1, width = 10, height = 5.38)



### Work out constraints on canonical solution line
# canonical solution line and nonlinear deviations
theta <- apc_model_theta(atus, "nonwork_alone_min", 5)
nld   <- apc_model_nld(atus, "nonwork_alone_min", 5)

# min and max alpha based on non-monotonic age effect with minimum at age 35
contrasts(atus$a) <- contr.poly.weighted(atus$a, width = 5)
nld_a    <- nld[nld$apc == "age", ]
nld_a$lc <- contrasts(atus$a)[, 1]

find_alpha <- function(alpha) {
  nld_a$te <- nld_a$nld + nld_a$lc*(alpha)
  ind <- which.min(nld_a$te)
  age <- nld_a$index[ind]
  return(age)
}

am <- round( theta[1], 1 )
alpha_range <- matrix(NA, length(seq(0, am, 0.1)), 2)
r <- 0
for(i in seq(0, am, 0.1)) {
  r <- r + 1
  alpha_range[r,1] <- i
  alpha_range[r,2] <- find_alpha(i)
}
min_alpha <- min(alpha_range[alpha_range[, 2] == 35, 1])
max_alpha <- max(alpha_range[alpha_range[, 2] == 35, 1])
mid_alpha <- (max_alpha - min_alpha)/2 + min_alpha

min_pi <- theta[1] - max_alpha
mid_pi <- theta[1] - mid_alpha
max_pi <- theta[1] - min_alpha

min_gamma <- theta[2] - max_pi
mid_gamma <- theta[2] - mid_pi
max_gamma <- theta[2] - min_pi

c(min_alpha, mid_alpha, max_alpha)
c(min_pi, mid_pi, max_pi)
c(min_gamma, mid_gamma, max_gamma)

### Figure 2 - Canonical solution line with constrained region noted
f2 <- fig2D(theta[1], theta[2], -2, 5, min_alpha, mid_alpha, max_alpha)
f2
ggsave("f2.jpg", f2, width = 8, height = 8)



### Figure 3 - Range of total effects based on constrained region of canonical solution line
te <- total_effects_range(atus, "nonwork_alone_min", 5, min_alpha, mid_alpha, max_alpha)
te

f3a <- te_fig(subset(te, apc == "age"), 100, 500, "Age")
f3a <- f3a + scale_x_discrete(name = "age", breaks = c("15", "75"), labels = c("15", "79"))

f3p <- te_fig(subset(te, apc == "year"), 100, 500, "Period")
f3p <- f3p + scale_x_discrete(name = "year", breaks = c("2003", "2018"), labels = c("2003", "2022"))

f3c <- te_fig(subset(te, apc == "cohort"), 100, 500, "Cohort")
f3c <- f3c + scale_x_discrete(name = "cohort", breaks = c("1924", "1999"), labels = c("1924", "2007"))

f3  <- ggarrange(f3a, f3p, f3c, nrow = 1)
f3
ggsave("f3.jpg", f3, width = 10, height = 5.38)

