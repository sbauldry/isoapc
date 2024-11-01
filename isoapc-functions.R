### Purpose: Various helper functions for preliminary APC analyses
### Author:  S Bauldry 
### Date:    Oct 17, 2024

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
figUniAPC <- function(df, v, dv, lb, ub, tit) {
  mn <- df |>
    group_by(!!sym(v)) |>
    summarise(wm = sum(!!sym(dv)*rwt)/sum(rwt))
  f <- ggplot(mn, aes(x = !!sym(v), y = wm)) +
    geom_point() +
    geom_smooth(method = "gam", formula = y ~ s(x), se = F) +
    scale_y_continuous(name = "mean minutes", limits = c(lb, ub)) +
    ggtitle(tit) +
    theme_light() +
    theme(
      axis.text  = element_text(size = 10, family = "Times New Roman"),
      axis.title = element_text(size = 14, family = "Times New Roman"),
      plot.title = element_text(size = 16, family = "Times New Roman")
    )
  return(f)
}




### Function to fit APC model and return theta
apc_model_theta <- function(df, dv, pord) {
  contrasts(df$a)      <- contr.poly.weighted(df$a, width = 5)
  contrasts(df$p)      <- contr.poly.weighted(df$p, width = 5)
  contrasts(df$c)      <- contr.poly.weighted(df$c, width = 9)
  contrasts(df$p)[, 1] <- 0
  
  # smoothing above given order of polynomial
  for(i in seq(pord, 15, 1)) {
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


### Function to fit APC model and return non-linear deviations
apc_model_nld <- function(df, dv, pord) {
  contrasts(df$a)      <- contr.poly.weighted(df$a, width = 5)
  contrasts(df$p)      <- contr.poly.weighted(df$p, width = 5)
  contrasts(df$c)      <- contr.poly.weighted(df$c, width = 9)
  contrasts(df$p)[, 1] <- 0

  # smoothing above given order of polynomial  
  for(i in seq(pord, 15, 1)) {
    contrasts(df$c)[, i] <- 0
    if(i < 13) {
      contrasts(df$a)[, i] <- 0
    }
  }
  
  model  <- as.formula( paste0(dv, " ~ a + p + c + day + hol") )
  m1     <- lm(model, data = df, weights = rwt)
  
  a_nld <- deviations_w_intercept(m1, "^a", contrasts(df$a))
  p_nld <- deviations_w_intercept(m1, "^p", contrasts(df$p))
  c_nld <- deviations_w_intercept(m1, "^c", contrasts(df$c))
  
  a_nld_df <- data.frame( rep("age", 13), seq(15, 75, 5), a_nld )
  p_nld_df <- data.frame( rep("year", 4), seq(2003, 2018, 5), p_nld )
  c_nld_df <- data.frame( rep("cohort", 16), seq(1924, 1999, 5), c_nld )
  colnames(a_nld_df) <- c("apc", "index", "nld")
  colnames(p_nld_df) <- c("apc", "index", "nld")
  colnames(c_nld_df) <- c("apc", "index", "nld")
  nld <- rbind(a_nld_df, p_nld_df, c_nld_df)
  return(nld)
}


### Function for 2D APC plot
fig2D <- function(theta1, theta2, lim1, lim2, min_alpha, mid_alpha, max_alpha) {
  min_pi     <- theta1 - max_alpha
  mid_pi     <- theta1 - mid_alpha
  max_pi     <- theta1 - min_alpha

  f <- ggplot() +
    geom_abline(aes(intercept = theta1, slope = -1), linewidth = 1, color = "darkblue") +
    geom_point(aes(x = mid_pi, y = mid_alpha), size = 3) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_hline(yintercept = -(theta2 - theta1), linetype = "dotted") +
    geom_vline(xintercept = 0, linetype = "dotted") +
    geom_segment(aes(x = lim1, y = max_alpha, xend = min_pi, yend = max_alpha), linetype = "dashed", color = "darkblue") + 
    geom_segment(aes(x = lim1, y = min_alpha, xend = max_pi, yend = min_alpha), linetype = "dashed", color = "darkblue") + 
    geom_text(aes(x = lim1 + 0.1, y = max_alpha + 0.1, label = "maximum linear age effect"),  hjust = 0, vjust = 0, family = "Times New Roman") +
    geom_text(aes(x = lim1 + 0.1, y = min_alpha + 0.1, label = "minimum linear age effect"),  hjust = 0, vjust = 0, family = "Times New Roman") +
    scale_y_continuous("age linear effect", limits = c(lim1, lim2), sec.axis = sec_axis(~.+(theta2 - theta1), name = "cohort linear effect")) +
    scale_x_continuous("period linear effect", limits = c(lim1, lim2), expand = c(0, 0), sec.axis = dup_axis()) +
    theme_light() +
    theme(
      axis.text  = element_text(size = 10, family = "Times New Roman"),
      axis.title = element_text(size = 14, family = "Times New Roman"),
      plot.title = element_text(size = 16, family = "Times New Roman"))
  return(f)
}


### Function to calculate total effects
total_effects_range <- function(df, dv, pord, min_alpha, mid_alpha, max_alpha) {
  contrasts(df$a)      <- contr.poly.weighted(df$a, width = 5)
  contrasts(df$p)      <- contr.poly.weighted(df$p, width = 5)
  contrasts(df$c)      <- contr.poly.weighted(df$c, width = 9)
  contrasts(df$p)[, 1] <- 0
  
  for(i in seq(pord, 15, 1)) {
    contrasts(df$c)[, i] <- 0
    if(i < 13) {
      contrasts(df$a)[, i] <- 0
    }
  }
  
  model  <- as.formula( paste0(dv, " ~ a + p + c + day + hol") )
  m1     <- lm(model, data = df, weights = rwt)
  theta1 <- coef(m1)["a.L"]
  theta2 <- coef(m1)["c.L"]
  
  min_pi    <- theta1 - max_alpha
  mid_pi    <- theta1 - mid_alpha
  max_pi    <- theta1 - min_alpha
  min_gamma <- theta2 - max_pi
  mid_gamma <- theta2 - mid_pi
  max_gamma <- theta2 - min_pi
  
  a_min_te <- total_effect(m1, "^a", contrasts(df$a), min_alpha)
  a_mid_te <- total_effect(m1, "^a", contrasts(df$a), mid_alpha)
  a_max_te <- total_effect(m1, "^a", contrasts(df$a), max_alpha)
  adf      <- as.data.frame(cbind(rep("age", 13), seq(15, 75, 5), a_min_te, a_mid_te, a_max_te))
  colnames(adf) <- c("apc", "index", "min_te", "mid_te", "max_te")
  
  contrasts(df$p) <- contr.poly.weighted(df$p, width = 5)
  p_min_te <- total_effect(m1, "^p", contrasts(df$p), min_pi)
  p_mid_te <- total_effect(m1, "^p", contrasts(df$p), mid_pi)
  p_max_te <- total_effect(m1, "^p", contrasts(df$p), max_pi)
  pdf      <- as.data.frame(cbind(rep("year", 4), seq(2003, 2018, 5), p_min_te, p_mid_te, p_max_te))
  colnames(pdf) <- c("apc", "index", "min_te", "mid_te", "max_te")
  
  c_min_te <- total_effect(m1, "^c", contrasts(df$c), min_gamma)
  c_mid_te <- total_effect(m1, "^c", contrasts(df$c), mid_gamma)
  c_max_te <- total_effect(m1, "^c", contrasts(df$c), max_gamma)
  cdf <- as.data.frame(cbind(rep("cohort", 16), seq(1924, 1999, 5), c_min_te, c_mid_te, c_max_te))
  colnames(cdf) <- c("apc", "index", "min_te", "mid_te", "max_te")
  
  df2 <- rbind(adf, pdf, cdf) |>
    mutate(min_te = as.numeric(min_te),
           mid_te = as.numeric(mid_te),
           max_te = as.numeric(max_te))
  return(df2)
}


### Function for total effect plots
te_fig <- function(d, ll, ul, tit) {
  f <- ggplot(data = d, mapping = aes(x = index, group = 1)) +
    geom_line(mapping = aes(y = mid_te)) +
    geom_ribbon(mapping = aes(ymin = min_te, ymax = max_te), fill = "lightgreen", alpha = 0.3) +
    scale_y_continuous(name = "minutes alone", limits = c(ll, ul)) +
    ggtitle(tit) +
    theme_light() +
    theme(
      axis.text  = element_text(size = 10, family = "Times New Roman"),
      axis.title = element_text(size = 14, family = "Times New Roman"),
      plot.title = element_text(size = 16, family = "Times New Roman"))
  return(f)
}


