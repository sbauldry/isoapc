### Purpose: Various helper functions for APC analyses
### Author:  S Bauldry 
### Date:    Sep 17, 2024

### Function to calculate weighted orthogonal polynomial contrasts
### Adapted from Elbers (https://github.com/elbersb/weightedcontrasts)
contr.poly.weighted <- function (f, width = 1) {
  
  n       <- length(levels(f))
  weights <- as.numeric(table(f)) / length(f)
  
  y     <- 1:n - sum(1:n * weights)
  X     <- sqrt(weights) * outer(y, seq_len(n) - 1, "^")
  QR    <- qr(X)
  z     <- QR$qr
  z     <- z * (row(z) == col(z))
  raw   <- qr.qy(QR, z) / sqrt(weights)
  contr <- sweep(raw, 2L, apply(raw, 2L, function(x) sqrt(sum(x^2))), "/", check.margin = F)
  
  scores     <- seq(1, width * n, by = width)
  scores     <- scores - sum(scores * weights)
  contr[, 2] <- contr[, 2] * sqrt(sum(scores^2))
  
  dn              <- paste0("^", 1L:n - 1L)
  dn[2:min(4, n)] <- c(".L", ".Q", ".C")[1:min(3, n - 1)]
  colnames(contr) <- dn
  contr[, -1, drop = FALSE]
}


### Functions to extract nonlinear effects and compute total effects
### Adapted from Elbers (https://github.com/elbersb/weightedcontrasts)
deviations_w_intercept <- function(model, set, contrasts) {
  nonlinear_coefs <- coef(model)[grepl(set, names(coef(model)))]
  nonlinear_coefs <- nonlinear_coefs[2:length(nonlinear_coefs)]
  nonlinear_coefs[is.na(nonlinear_coefs)] <- 0
  deviations <- contrasts[, 2:(1 + length(nonlinear_coefs))] %*% nonlinear_coefs
  ne <- coef(model)[1] + deviations[, 1]
  return(ne)
}

total_effect <- function(model, set, contrasts, linear_coef) {
  nonlinear_coefs <- coef(model)[grepl(set, names(coef(model)))]
  nonlinear_coefs <- nonlinear_coefs[2:length(nonlinear_coefs)]
  nonlinear_coefs[is.na(nonlinear_coefs)] <- 0
  linear_contrast <- contrasts[, 1]
  nonlinear_contrasts <- contrasts[, 2:(1 + length(nonlinear_coefs))]
  ne <- nonlinear_contrasts %*% nonlinear_coefs      
  te <- coef(model)[1] + linear_contrast*linear_coef + ne[, 1]
  return(te)
}


### Function for univariate APC plots
figUniAPC <- function(df, v, dv, lb, ub) {
  mn <- df |>
    group_by(!!sym(v)) |>
    summarise(wm = sum(!!sym(dv)*rwt)/sum(rwt))
  f <- ggplot(mn, aes(x = !!sym(v), y = wm)) +
    geom_point() +
    geom_smooth(method = "gam", formula = y ~ s(x), se = F) +
    scale_y_continuous(name = "mean minutes", limits = c(lb, ub)) +
    theme_light() +
    theme(
      axis.text  = element_text(size = 10),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 16)
    )
  return(f)
}


### Function to create Lexis plot
figLexis <- function(df, dv, lb, ub, tit) {
  d <- df |>
    group_by(a, p) |>
    summarise(wm = sum(!!sym(dv)*rwt)/sum(rwt))
  f <- ggplot(data = d, mapping = aes(x = p, y = a, fill = wm)) +
    geom_tile() +
    scale_fill_gradient(low="white", high="darkred", name = "mean \nminutes", 
                        limits = c(lb, ub)) +
    scale_y_discrete(labels = c("15-19", "20-24", "25-29", "30-34", "35-39", 
                                "40-44", "45-49", "50-54", "55-59", "60-64", 
                                "65-69", "70-74", "75-79"), name = "age") +
    scale_x_discrete(labels = c("2003-07", "2008-12", "2013-17", "2018-22"),
                     name = "year") +
    ggtitle(tit) +
    theme_light() +
    theme(
      axis.text  = element_text(size = 10),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 16)
    )
  return(f)
}


### Function to define contrasts, fit APC model, and return theta
### Note 1: model adjusts for day of week and holidays
### Note 2: model smooths nonlinearities for a and c beyond 4th order polynomial
apc_model <- function(df, dv) {
  contrasts(df$a)      <- contr.poly.weighted(df$a, width = 5)
  contrasts(df$p)      <- contr.poly.weighted(df$p, width = 5)
  contrasts(df$c)      <- contr.poly.weighted(df$c, width = 10)
  contrasts(df$p)[, 1] <- 0
  
  for(i in seq(5, 15, 1)) {
    contrasts(df$c)[, i] <- 0
    if(i < 13) {
      contrasts(df$a)[, i] <- 0
    }
  }
  
  model  <- as.formula( paste0(dv, " ~ a + p + c + day + hol") )
  m1     <- lm(model, data = df, weights = rwt)
  theta1 <- coef(m1)["a.L"]
  theta2 <- coef(m1)["c.L"]
  theta  <- c(theta1, theta2)
  return(theta)
}


### Function to define contrasts, fit APC model, and return theta
### Note 1: model adjusts for gender, race, education, day of week, and holidays
### Note 2: model smooths nonlinearities for a and c beyond 4th order polynomial
apc_model_cov <- function(df, dv) {
  contrasts(df$a)      <- contr.poly.weighted(df$a, width = 5)
  contrasts(df$p)      <- contr.poly.weighted(df$p, width = 5)
  contrasts(df$c)      <- contr.poly.weighted(df$c, width = 10)
  contrasts(df$p)[, 1] <- 0
  
  for(i in seq(5, 15, 1)) {
    contrasts(df$c)[, i] <- 0
    if(i < 13) {
      contrasts(df$a)[, i] <- 0
    }
  }
  
  model  <- as.formula( paste0(dv, " ~ a + p + c + fem + rce + edu + day + hol") )
  m1     <- lm(model, data = df, weights = rwt)
  theta1 <- coef(m1)["a.L"]
  theta2 <- coef(m1)["c.L"]
  theta  <- c(theta1, theta2)
  return(theta)
}


### Function for 2D APC plot
fig2D <- function(theta1, theta2, tit) {
  sum_thetas = abs(theta2 - theta1)
  limits     = c(-1.1*sum_thetas, 1.1*sum_thetas)
  
  f <- ggplot() +
    geom_abline(aes(intercept = theta1, slope = -1), linewidth = 1, 
                color = "darkblue") +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_hline(yintercept = -(theta2 - theta1), linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0, ymax = -(theta2 - theta1)), 
              fill = "lightgreen", alpha = 0.3) +
    scale_y_continuous(expression(alpha), 
                       limits = limits,
                       sec.axis = sec_axis(~.+(theta2 - theta1), 
                                           name = expression(gamma))) +
    scale_x_continuous(expression(pi), limits = limits, sec.axis = dup_axis()) +
    ggtitle(tit) +
    theme_light() +
    theme(text = element_text(size = 14),
          axis.title.y.left = element_text(angle = 0, vjust = 0.5),
          axis.title.y.right = element_text(angle = 0, vjust = 0.5))
  return(f)
}


### Function to calculate total effects across range of parameter values
### Note 1: model adjusts for day of week and holidays
### Note 2: model smooths nonlinearities for a and c beyond 5th order polynomial
### Note 3: min, mid, and max values based on pi > 0, alpha > 0, & gamma < 0
total_effects_range <- function(df, dv) {
  contrasts(df$a)      <- contr.poly.weighted(df$a, width = 5)
  contrasts(df$p)      <- contr.poly.weighted(df$p, width = 5)
  contrasts(df$c)      <- contr.poly.weighted(df$c, width = 9)
  contrasts(df$p)[, 1] <- 0
  
  for(i in seq(5, 15, 1)) {
    contrasts(df$c)[, i] <- 0
    if(i < 13) {
      contrasts(df$a)[, i] <- 0
    }
  }
  
  model  <- as.formula( paste0(dv, " ~ a + p + c + day + hol") )
  m1     <- lm(model, data = df, weights = rwt)
  theta1 <- coef(m1)["a.L"]
  theta2 <- coef(m1)["c.L"]
  
  min_pi <- theta2
  max_pi <- theta1
  mid_pi <- (max_pi - min_pi)/2
  
  min_alpha <- theta1 - max_pi
  max_alpha <- theta1 - min_pi
  mid_alpha <- (max_alpha - min_alpha)/2
  
  min_gamma <- theta2 - max_pi
  max_gamma <- theta2 - min_pi
  mid_gamma <- (min_gamma - max_gamma)/2
  
  pi_min_gamma    <- theta2 - min_gamma
  alpha_min_gamma <- theta1 - pi_min_gamma
  
  pi_max_gamma    <- theta2 - max_gamma
  alpha_max_gamma <- theta1 - pi_max_gamma
  
  a_min  <- total_effect(m1, "^a", contrasts(df$a), min_alpha)
  a_mid  <- total_effect(m1, "^a", contrasts(df$a), mid_alpha)
  a_max  <- total_effect(m1, "^a", contrasts(df$a), max_alpha)
  ac_min <- total_effect(m1, "^a", contrasts(df$a), alpha_min_gamma)
  ac_max <- total_effect(m1, "^a", contrasts(df$a), alpha_max_gamma)
  
  adf <- as.data.frame(cbind(rep("age", 13), seq(15, 75, 5), a_min, a_mid, 
                             a_max, ac_min, ac_max))
  colnames(adf) <- c("apc", "index", "min", "mid", "max", "cmin", "cmax")
  
  contrasts(df$p) <- contr.poly.weighted(df$p, width = 5)
  p_min  <- total_effect(m1, "^p", contrasts(df$p), min_pi)
  p_mid  <- total_effect(m1, "^p", contrasts(df$p), mid_pi)
  p_max  <- total_effect(m1, "^p", contrasts(df$p), max_pi)
  pc_min <- total_effect(m1, "^p", contrasts(df$p), pi_min_gamma)
  pc_max <- total_effect(m1, "^p", contrasts(df$p), pi_max_gamma)
  
  pdf <- as.data.frame(cbind(rep("year", 4), seq(2003, 2018, 5), p_min, p_mid, 
                             p_max, pc_min, pc_max))
  colnames(pdf) <- c("apc", "index", "min", "mid", "max", "cmin", "cmax")
  
  c_min <- total_effect(m1, "^c", contrasts(df$c), min_gamma)
  c_mid <- total_effect(m1, "^c", contrasts(df$c), mid_gamma)
  c_max <- total_effect(m1, "^c", contrasts(df$c), max_gamma)
  
  cdf <- as.data.frame(cbind(rep("cohort", 16), seq(1924, 1999, 5), c_min, 
                             c_mid, c_max, c_min, c_max))
  colnames(cdf) <- c("apc", "index", "min", "mid", "max", "cmin", "cmax")
  
  df2 <- rbind(adf, pdf, cdf) |>
    mutate(index = as.numeric(index),
           min   = ifelse(as.numeric(min) > 0, as.numeric(min), 0),
           mid   = as.numeric(mid),
           max   = as.numeric(max),
           cmin  = ifelse(as.numeric(cmin) > 0, as.numeric(cmin), 0),
           cmax  = ifelse(as.numeric(cmax) > 0, as.numeric(cmax), 0))
  return(df2)
}


### Functions for total effect plots
te_age_fig <- function(d, ll, ul) {
  f <- ggplot(data = d, mapping = aes(x = index)) +
    geom_point(mapping = aes(y = mid)) +
    geom_line(mapping = aes(y = mid)) +
    geom_line(mapping = aes(y = cmin), color = "red") +
    geom_line(mapping = aes(y = cmax), color = "green") +
    geom_ribbon(mapping = aes(ymin = min, ymax = max), fill = "lightblue", 
                alpha = 0.5) +
    scale_x_continuous(name = "age", breaks = c(15,75)) +
    scale_y_continuous(name = "minutes alone", limits = c(ll, ul)) +
    ggtitle("Age") +
    theme_light() +
    theme(
      axis.text  = element_text(size = 10),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 16)
    )
  return(f)
}

te_per_fig <- function(d, ll, ul) {
  f <- ggplot(data = d, mapping = aes(x = index)) +
    geom_point(mapping = aes(y = mid)) +
    geom_line(mapping = aes(y = mid)) +
    geom_line(mapping = aes(y = cmin), color = "red") +
    geom_line(mapping = aes(y = cmax), color = "green") +
    geom_ribbon(mapping = aes(ymin = min, ymax = max), fill = "lightblue", 
                alpha = 0.5) +
    scale_x_continuous(name = "year", breaks = c(2003, 2018)) +
    scale_y_continuous(name = "minutes alone", limits = c(ll, ul)) +
    ggtitle("Period") +
    theme_light() +
    theme(
      axis.text  = element_text(size = 10),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 16)
    )
  return(f)
}

te_coh_fig <- function(d, ll, ul) {
  f <- ggplot(data = d, mapping = aes(x = index)) +
    geom_point(mapping = aes(y = mid)) +
    geom_line(mapping = aes(y = mid)) +
    geom_line(mapping = aes(y = cmin), color = "red") +
    geom_line(mapping = aes(y = cmax), color = "green") +
    geom_ribbon(mapping = aes(ymin = min, ymax = max), fill = "lightblue", 
                alpha = 0.5) +
    scale_x_continuous(name = "cohort", breaks = c(1924, 1999)) +
    scale_y_continuous(name = "minutes alone", limits = c(ll, ul)) +
    ggtitle("Cohort") +
    theme_light() +
    theme(
      axis.text  = element_text(size = 10),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 16)
    )
  return(f)
}


### Functions for marginal total effect plots
mte_age_fig <- function(d, ll, ul) {
  f <- ggplot(data = d, mapping = aes(x = index)) +
    geom_point(mapping = aes(y = dmid)) +
    geom_line(mapping = aes(y = dmid)) +
    geom_line(mapping = aes(y = dmin), color = "red") +
    geom_line(mapping = aes(y = dmax), color = "green") +
    geom_ribbon(mapping = aes(ymin = dmin, ymax = dmax), fill = "lightblue", alpha = 0.5) +
    scale_x_continuous(name = "age", breaks = seq(15, 80, 10)) +
    scale_y_continuous(name = "difference in minutes alone", limits = c(ll, ul)) +
    ggtitle("Age") +
    theme_light()
  return(f)
}

mte_per_fig <- function(d, ll, ul) {
  f <- ggplot(data = d, mapping = aes(x = index)) +
    geom_point(mapping = aes(y = dmid)) +
    geom_line(mapping = aes(y = dmid)) +
    geom_line(mapping = aes(y = dmin), color = "red") +
    geom_line(mapping = aes(y = dmax), color = "green") +
    geom_ribbon(mapping = aes(ymin = dmin, ymax = dmax), fill = "lightblue", alpha = 0.5) +
    scale_x_continuous(name = "year", breaks = seq(2003, 2023, 5)) +
    scale_y_continuous(name = "difference in minutes alone", limits = c(ll, ul)) +
    ggtitle("Period") +
    theme_light()
  return(f)
}

mte_coh_fig <- function(d, ll, ul) {
  f <- ggplot(data = d, mapping = aes(x = index)) +
    geom_point(mapping = aes(y = dmid)) +
    geom_line(mapping = aes(y = dmid)) +
    geom_line(mapping = aes(y = dmin), color = "red") +
    geom_line(mapping = aes(y = dmax), color = "green") +
    geom_ribbon(mapping = aes(ymin = dmin, ymax = dmax), fill = "lightblue", alpha = 0.5) +
    scale_x_continuous(name = "cohort", breaks = seq(1923, 2008, 20)) +
    scale_y_continuous(name = "difference in minutes alone", limits = c(ll, ul)) +
    ggtitle("Cohort") +
    theme_light()
  return(f)
}



### Function to calculate total effects across range of parameter values
### Note 1: model adjusts for gender, race, education, day of week, and holidays
### Note 2: model smooths nonlinearities for a and c beyond 5th order polynomial
### Note 3: min, mid, and max values based on pi > 0, alpha > 0, & gamma < 0
total_effects_range_cov <- function(df, dv) {
  contrasts(df$a)      <- contr.poly.weighted(df$a, width = 5)
  contrasts(df$p)      <- contr.poly.weighted(df$p, width = 5)
  contrasts(df$c)      <- contr.poly.weighted(df$c, width = 9)
  contrasts(df$p)[, 1] <- 0
  
  for(i in seq(5, 15, 1)) {
    contrasts(df$c)[, i] <- 0
    if(i < 13) {
      contrasts(df$a)[, i] <- 0
    }
  }
  
  model  <- as.formula( paste0(dv, " ~ a + p + c + fem + rce + edu + day + hol") )
  m1     <- lm(model, data = df, weights = rwt)
  theta1 <- coef(m1)["a.L"]
  theta2 <- coef(m1)["c.L"]
  
  min_pi <- theta2
  max_pi <- theta1
  mid_pi <- (max_pi - min_pi)/2
  
  min_alpha <- theta1 - max_pi
  max_alpha <- theta1 - min_pi
  mid_alpha <- (max_alpha - min_alpha)/2
  
  min_gamma <- theta2 - max_pi
  max_gamma <- theta2 - min_pi
  mid_gamma <- (min_gamma - max_gamma)/2
  
  pi_min_gamma    <- theta2 - min_gamma
  alpha_min_gamma <- theta1 - pi_min_gamma
  
  pi_max_gamma    <- theta2 - max_gamma
  alpha_max_gamma <- theta1 - pi_max_gamma
  
  a_min  <- total_effect(m1, "^a", contrasts(df$a), min_alpha)
  a_mid  <- total_effect(m1, "^a", contrasts(df$a), mid_alpha)
  a_max  <- total_effect(m1, "^a", contrasts(df$a), max_alpha)
  ac_min <- total_effect(m1, "^a", contrasts(df$a), alpha_min_gamma)
  ac_max <- total_effect(m1, "^a", contrasts(df$a), alpha_max_gamma)
  
  adf <- as.data.frame(cbind(rep("age", 13), seq(15, 75, 5), a_min, a_mid, 
                             a_max, ac_min, ac_max))
  colnames(adf) <- c("apc", "index", "min", "mid", "max", "cmin", "cmax")
  
  contrasts(df$p) <- contr.poly.weighted(df$p, width = 5)
  p_min  <- total_effect(m1, "^p", contrasts(df$p), min_pi)
  p_mid  <- total_effect(m1, "^p", contrasts(df$p), mid_pi)
  p_max  <- total_effect(m1, "^p", contrasts(df$p), max_pi)
  pc_min <- total_effect(m1, "^p", contrasts(df$p), pi_min_gamma)
  pc_max <- total_effect(m1, "^p", contrasts(df$p), pi_max_gamma)
  
  pdf <- as.data.frame(cbind(rep("year", 4), seq(2003, 2018, 5), p_min, p_mid, 
                             p_max, pc_min, pc_max))
  colnames(pdf) <- c("apc", "index", "min", "mid", "max", "cmin", "cmax")
  
  c_min <- total_effect(m1, "^c", contrasts(df$c), min_gamma)
  c_mid <- total_effect(m1, "^c", contrasts(df$c), mid_gamma)
  c_max <- total_effect(m1, "^c", contrasts(df$c), max_gamma)
  
  cdf <- as.data.frame(cbind(rep("cohort", 16), seq(1924, 1999, 5), c_min, 
                             c_mid, c_max, c_min, c_max))
  colnames(cdf) <- c("apc", "index", "min", "mid", "max", "cmin", "cmax")
  
  df2 <- rbind(adf, pdf, cdf) |>
    mutate(index = as.numeric(index),
           min   = ifelse(as.numeric(min) > 0, as.numeric(min), 0),
           mid   = as.numeric(mid),
           max   = as.numeric(max),
           cmin  = ifelse(as.numeric(cmin) > 0, as.numeric(cmin), 0),
           cmax  = ifelse(as.numeric(cmax) > 0, as.numeric(cmax), 0))
  return(df2)
}