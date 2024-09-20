### Purpose: Prepare ATUS data for APC analysis of time spent alone
### Author:  S Bauldry (based on S Peng Stata code)
### Date:    Sep 14, 2024

### Set working directory and load libraries
setwd("~/desktop")
library(tidyverse)
library(haven)

### Load various ATUS data files
act <- read_dta("atusact_0323.dta")
who <- read_dta("atuswho_0323.dta")
rsp <- read_dta("atusresp_0323.dta")
sum <- read_dta("atussum_0323.dta")

### Prepare variables for analysis
d1 <- left_join(act, who, by = c("tucaseid","tuactivity_n")) |>
  filter(trcodep < 50101 | trcodep > 59999) |>  # remove work-related activities
  mutate(
    # identify leisure activities (eating and drinking; socializing, relaxing, 
    # and leisure; sports, exercise, recreation)
    leisure = ifelse( (trcodep >= 110101 & trcodep <= 119999) |
                      (trcodep >= 120101 & trcodep <= 129999) | 
                      (trcodep >= 130101 & trcodep <= 139999), 1, 0),
    
    # identify activities while alone
    alone = case_when (
      tuwho_code < 0 ~ NA,
      tuwho_code == 18 | tuwho_code == 19 ~ 1,
      tuwho_code >= 20 ~ 0
    ),
    
    # minutes spent in all activities and leisure activities while alone
    nonwork_alone_min = ifelse(alone == 1, tuactdur24, NA),
    leisure_alone_min = ifelse(alone == 1 & leisure == 1, tuactdur24, NA)
  ) |>
  select(tucaseid, nonwork_alone_min, leisure_alone_min) |>
  group_by(tucaseid) |>
  summarise(
    nonwork_alone_min = sum(nonwork_alone_min, na.rm = T),
    leisure_alone_min = sum(leisure_alone_min, na.rm = T)
  )

# merge in respondent characteristics
d2 <- left_join(d1, rsp, by = "tucaseid") |>
  left_join(sum, by = "tucaseid") |>
  filter(teage < 80 & tuyear.x < 2023) |>
  mutate(
    # Continuous APC variables
    age    = teage,
    year   = tuyear.x,
    cohort = year - age,
    
    # Categorical APC variables
    a = as.numeric( cut(age, breaks = seq(15, 80, 5), right = F) ),
    p = as.numeric( cut(year, breaks = seq(2003, 2023, 5), right = F) ),
    c = p - a + 13,
    
    # day of the week and holidays
    day = tudiaryday.x,
    hol = trholiday.x,
    
    # identify R living alone
    lal = ifelse(trnumhou == 1, 1, 0),
    
    # sociodemographics
    fem = ifelse(tesex == 2, 1, 0),
    rce  = case_when (
      pehspnon == 1 ~ 1, # Hispanic
      ptdtrace == 1 ~ 2, # non-Hispanic White
      ptdtrace == 2 ~ 3, # non-Hispanic Black
      ptdtrace == 4 ~ 4, # Asian
      ptdtrace == 3 | ptdtrace > 4 ~ 5 # other race/ethnicity
    ),
    edu = case_when (
      peeduca <= 38 ~ 1, # less than HS
      peeduca == 39 ~ 2, # HS or GED
      peeduca <= 42 ~ 3, # some college/technical
      peeduca <= 46 ~ 4  # college
    ),
    fte = ifelse(telfs.x == 1, 1, 0),
    
    # weights
    wgt = ifelse(tufnwgtp.x != -1, tufnwgtp.x, tu20fwgt.x),
    rwt = wgt/sum(wgt)*n()
  ) |>
  select(tucaseid, nonwork_alone_min, leisure_alone_min, age, year, cohort, 
         a, p, c, day, hol, lal, fem, rce, edu, fte, wgt, rwt)
  
### Save data for analysis
write_csv(d2, "isoapc-data.csv")
