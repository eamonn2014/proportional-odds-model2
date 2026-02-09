# Panel A: Interim estimates. Red shows trials that stopped early. These are selected because they looked unusually good early on, so this subgroup tends to show more favourable estimates (winner’s curse).
# Panel B: Final estimates. Blue shows trials that reached the final analysis (they did not meet the early success rule). Because this is a selected subgroup, it can look less favourable on average — this is a selection/conditioning effect, not a flaw in the final statistical method.
# Panel C: Boxplots summarise the same selection effects. The percentages printed next to each box show the proportion of all trials in each category (e.g., stopped early vs reached final).
# Clinical takeaway: early-stopped results can look overly optimistic because of selection; trials that continue to final are a different (selected) subset. This does not mean the final OR is biased — it reflects the trial pathway and stopping rules.
rm(list=ls())

suppressPackageStartupMessages({
  library(MASS)   # polr
  library(dplyr)
})

# ============================================================
# 1) Core helpers
# ============================================================

ilogit <- function(z) 1 / (1 + exp(-z))

parse_probs <- function(txt) {
  as.numeric(trimws(unlist(strsplit(txt, ","))))
}

validate_probs <- function(p) {
  if (any(is.na(p)))          return("Probabilities must be numeric (comma-separated).")
  if (length(p) < 3)          return("Need at least 3 categories.")
  if (any(p <= 0))            return("All probabilities must be > 0.")
  if (abs(sum(p) - 1) > 1e-6) return(sprintf("Must sum to ~1 (got %.6f).", sum(p)))
  NULL
}

theta_from_control_pmf <- function(p_control) {
  cum_control <- cumsum(p_control)
  qlogis(cum_control[seq_len(length(p_control) - 1)])
}

pmf_from_beta <- function(theta, beta, x) {
  cdf <- ilogit(theta - beta * x)
  cdf <- c(cdf, 1)
  pmf <- diff(c(0, cdf))
  pmf[pmf < 0] <- 0
  pmf / sum(pmf)
}

# ============================================================
# 2) Fit log(COR) correctly (NO sign flip)
# ============================================================

fit_logCOR <- function(df) {
  fit <- try(MASS::polr(y ~ trt, data = df, Hess = TRUE), silent = TRUE)
  if (inherits(fit, "try-error") || is.null(fit$coefficients)) return(NULL)
  
  coef_names <- names(fit$coefficients)
  trt_name <- coef_names[grepl("^trt", coef_names)][1]
  if (is.na(trt_name)) return(NULL)
  
  logCOR_hat <- as.numeric(fit$coefficients[[trt_name]])
  
  v <- vcov(fit)
  if (!(trt_name %in% rownames(v))) return(NULL)
  
  se <- sqrt(v[trt_name, trt_name])
  if (!is.finite(se) || se < 1e-8) return(NULL)
  
  list(logCOR_hat = logCOR_hat, se = se)
}

safe_summ_cor <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n == 0) {
    return(c(N = 0L, Mean = NA_real_, Median = NA_real_, `2.5%` = NA_real_, `97.5%` = NA_real_))
  }
  q <- quantile(x, probs = c(0.025, 0.975), names = FALSE)
  c(N = n, Mean = mean(x), Median = median(x), `2.5%` = q[1], `97.5%` = q[2])
}

# ============================================================
# 3) Simulation (NI success LEFT tail, futility RIGHT side)
# ============================================================

simulate_obf_ordinal <- function(
    COR_true, COR_NI,
    n_total,
    futility_frac, info_frac,
    alpha1, alpha2, futility_p,
    p_control, seed = 1234, nSims = 500,
    show_progress = interactive()
) {
  msg <- validate_probs(p_control)
  if (!is.null(msg)) stop(msg)
  
  set.seed(seed)
  
  theta <- theta_from_control_pmf(p_control)
  K <- length(p_control) - 1
  
  beta_true <- log(COR_true)
  beta_NI   <- log(COR_NI)
  
  pi_control <- p_control
  pi_treat   <- pmf_from_beta(theta, beta_true, x = 1)
  
  split_n <- function(N) {
    nC <- floor(N / 2)
    list(nC = nC, nT = N - nC)
  }
  
  Z_fut_all  <- rep(NA_real_, nSims)
  COR_fut_all <- rep(NA_real_, nSims)
  Z1_all     <- rep(NA_real_, nSims)
  COR1_all   <- rep(NA_real_, nSims)
  Z2_all     <- rep(NA_real_, nSims)
  COR2_all   <- rep(NA_real_, nSims)
  
  stop_fut   <- logical(nSims)
  stop_ia    <- logical(nSims)
  stop_final <- logical(nSims)
  
  z_fut  <- qnorm(futility_p)
  zcrit1 <- qnorm(alpha1)
  zcrit2 <- qnorm(alpha2)
  
  n_fut <- round(futility_frac * n_total)
  n1    <- round(info_frac  * n_total)
  
  if (show_progress) pb <- txtProgressBar(min = 0, max = nSims, style = 3)
  
  for (i in seq_len(nSims)) {
    
    # Futility look
    s_f <- split_n(n_fut)
    yCf <- sample(0:K, s_f$nC, TRUE, pi_control)
    yTf <- sample(0:K, s_f$nT, TRUE, pi_treat)
    
    df_f <- data.frame(
      y   = factor(c(yCf, yTf), ordered = TRUE, levels = 0:K),
      trt = factor(rep(c("C", "T"), c(s_f$nC, s_f$nT)), levels = c("C", "T"))
    )
    
    res_f <- fit_logCOR(df_f)
    if (!is.null(res_f)) {
      Z_fut_all[i]   <- (res_f$logCOR_hat - beta_NI) / res_f$se
      COR_fut_all[i] <- exp(res_f$logCOR_hat)
      
      if (Z_fut_all[i] > z_fut) {
        stop_fut[i] <- TRUE
        if (show_progress) setTxtProgressBar(pb, i)
        next
      }
    }
    
    # Interim look
    n_add_ia <- n1 - n_fut
    s_add <- split_n(n_add_ia)
    yCa <- sample(0:K, s_add$nC, TRUE, pi_control)
    yTa <- sample(0:K, s_add$nT, TRUE, pi_treat)
    
    df_add <- data.frame(
      y   = factor(c(yCa, yTa), ordered = TRUE, levels = 0:K),
      trt = factor(rep(c("C", "T"), c(s_add$nC, s_add$nT)), levels = c("C", "T"))
    )
    
    df1 <- rbind(df_f, df_add)
    
    res1 <- fit_logCOR(df1)
    if (!is.null(res1)) {
      Z1_all[i]   <- (res1$logCOR_hat - beta_NI) / res1$se
      COR1_all[i] <- exp(res1$logCOR_hat)
      
      if (Z1_all[i] <= zcrit1) {
        stop_ia[i] <- TRUE
        if (show_progress) setTxtProgressBar(pb, i)
        next
      }
    }
    
    # Final look
    n_add_final <- n_total - n1
    s_final <- split_n(n_add_final)
    yCf2 <- sample(0:K, s_final$nC, TRUE, pi_control)
    yTf2 <- sample(0:K, s_final$nT, TRUE, pi_treat)
    
    df_final_add <- data.frame(
      y   = factor(c(yCf2, yTf2), ordered = TRUE, levels = 0:K),
      trt = factor(rep(c("C", "T"), c(s_final$nC, s_final$nT)), levels = c("C", "T"))
    )
    
    df2 <- rbind(df1, df_final_add)
    
    res2 <- fit_logCOR(df2)
    if (!is.null(res2)) {
      Z2_all[i]   <- (res2$logCOR_hat - beta_NI) / res2$se
      COR2_all[i] <- exp(res2$logCOR_hat)
      
      if (Z2_all[i] <= zcrit2) {
        stop_final[i] <- TRUE
      }
    }
    
    if (show_progress) setTxtProgressBar(pb, i)
  }
  
  if (show_progress) close(pb)
  
  list(
    Z_fut_all = Z_fut_all, COR_fut_all = COR_fut_all,
    Z1_all = Z1_all, COR1_all = COR1_all,
    Z2_all = Z2_all, COR2_all = COR2_all,
    stop_fut = stop_fut, stop_ia = stop_ia, stop_final = stop_final,
    nSims = nSims,
    z_fut = z_fut, zcrit1 = zcrit1, zcrit2 = zcrit2
  )
}

# ============================================================
# 4) Outputs
# ============================================================

 

sim_table <- function(sim) {
  
  # Helper to compute all the usual stats + min & max
  safe_summ_cor_ext <- function(x) {
    x <- x[is.finite(x)]
    n <- length(x)
    if (n == 0) {
      return(c(
        N = 0L, 
        Min = NA_real_, Max = NA_real_,
        Mean = NA_real_, Median = NA_real_, 
        `2.5%` = NA_real_, `97.5%` = NA_real_
      ))
    }
    
    q <- quantile(x, probs = c(0.025, 0.5, 0.975), names = FALSE)
    c(
      N     = n,
      Min   = min(x),
      Max   = max(x),
      Mean  = mean(x),
      Median = q[2],           # or median(x)
      `2.5%` = q[1],
      `97.5%`= q[3]
    )
  }
  
  fut <- safe_summ_cor_ext(sim$COR_fut_all[sim$stop_fut])
  ia  <- safe_summ_cor_ext(sim$COR1_all[sim$stop_ia])
  fin <- safe_summ_cor_ext(sim$COR2_all[sim$stop_final])
  
  out <- data.frame(
    Stage  = c("Futility stop", "IA success stop", "Final success stop"),
    N      = c(fut["N"], ia["N"], fin["N"]),
    Min    = c(fut["Min"], ia["Min"], fin["Min"]),
  #  Max    = c(fut["Max"], ia["Max"], fin["Max"]),
    Mean   = c(fut["Mean"], ia["Mean"], fin["Mean"]),
    Median = c(fut["Median"], ia["Median"], fin["Median"]),
    `2.5%` = c(fut["2.5%"], ia["2.5%"], fin["2.5%"]),
    `97.5%`= c(fut["97.5%"], ia["97.5%"], fin["97.5%"]),
  Max    = c(fut["Max"], ia["Max"], fin["Max"]),
    check.names = FALSE
  )
  
  # Round numeric columns nicely
  out |> mutate(across(c(Min, Max, Mean, Median, `2.5%`, `97.5%`), ~ round(.x, 3)))
}


 


# Your current selection_boxplot function (unchanged)

selection_boxplot <- function(sim, COR_true, COR_NI,
                              futility_frac = 0.50,
                              info_frac = 0.80,
                              main = "Selection effects (log scale)",
                              point_cex = 0.58,
                              point_alpha = 0.30,
                              jitter_height = 0.22,
                              box_width = 0.26,
                              ylim = c(-1.3, 1.3),
                              text_offset_factor = 0.13) {
  
  log_true <- log(COR_true)
  log_M    <- log(COR_NI)
  
  ok_f <- is.finite(sim$COR_fut_all)
  ok1  <- is.finite(sim$COR1_all)
  ok2  <- is.finite(sim$COR2_all)
  
  groups <- list(
    "All estimates"       = log(sim$COR_fut_all[ok_f]),
    "Futile"              = log(sim$COR_fut_all[ok_f & sim$stop_fut]),
    "IA all estimates"    = log(sim$COR1_all[ok1]),
    "IA efficacy success" = log(sim$COR1_all[ok1 & sim$stop_ia]),
    "Final all"           = log(sim$COR2_all[ok2]),
    "Final success"       = log(sim$COR2_all[ok2 & sim$stop_final])
  )
  
  keep <- sapply(groups, function(v) length(v) > 0)
  groups <- groups[keep]
  group_names <- names(groups)
  
  if (length(groups) == 0) {
    plot.new()
    text(0.5, 0.5, "No valid data for boxplot", cex = 1.3, col = "gray50")
    return(invisible(NULL))
  }
  
  op <- par(mar = c(9.5, 11, 5, 15))
  on.exit(par(op), add = TRUE)
  
  plot(0, type = "n", 
       xlim = ylim, ylim = c(0.4, length(groups) + 0.8),
       xlab = "", ylab = "", yaxt = "n",
       main = main, las = 1)
  
  axis(2, at = seq_along(groups), labels = group_names, las = 1, cex.axis = 0.95)
  
  # 1. Points first (Background)
  for (i in seq_along(groups)) {
    vals <- groups[[i]]
    n <- length(vals)
    if (n == 0) next
    
    set.seed(20250610 + i * 43)
    jitter_y <- runif(n, -jitter_height, jitter_height)
    
    grp_name <- group_names[i]
    col_pt <- if (grp_name == "Futile") {
      rgb(0.9, 0.1, 0.1, point_alpha)
    } else if (grp_name %in% c("IA efficacy success", "Final success")) {
      rgb(0.1, 0.7, 0.1, point_alpha)
    } else {
      rgb(0.2, 0.2, 0.2, point_alpha)
    }
    
    points(vals, i + jitter_y, pch = 19, cex = point_cex, col = col_pt)
  }
  
  # 2. Boxplots with fixed whiskers
  for (i in seq_along(groups)) {
    vals <- groups[[i]]
    if (length(vals) < 2) next # Need at least 2 points for a box
    
    q <- quantile(vals, probs = c(0.25, 0.5, 0.75))
    iqr <- q[3] - q[1]
    
    # Calculate whisker reach
    whisk_low  <- min(vals[vals >= q[1] - 1.5 * iqr])
    whisk_high <- max(vals[vals <= q[3] + 1.5 * iqr])
    
    # BOX BODY
    rect(q[1], i - box_width/2, q[3], i + box_width/2,
         col = rgb(0.88, 0.93, 1.0, 0.45), border = "steelblue", lwd = 1.4)
    
    # MEDIAN
    segments(q[2], i - box_width/2, q[2], i + box_width/2, lwd = 4.5, col = "royalblue3")
    
    # WHISKER LINES (The Fix: horizontal connection)
    segments(whisk_low,  i, q[1], i, lwd = 2.5, col = "midnightblue")
    segments(whisk_high, i, q[3], i, lwd = 2.5, col = "midnightblue")
    
    # WHISKER CAPS (Vertical bars at the ends)
    cap_size <- box_width / 3
    segments(whisk_low,  i - cap_size, whisk_low,  i + cap_size, lwd = 2.5, col = "midnightblue")
    segments(whisk_high, i - cap_size, whisk_high, i + cap_size, lwd = 2.5, col = "midnightblue")
  }
  
  # Reference lines and Labels (Logic remains the same)
  abline(v = c(log_true, log_M), lty = c(2, 3), col = c("darkgreen", "red"), lwd = 2.3)
  
  # ... [Rest of the labeling/footnote code remains unchanged] ...
  p_fut   <- mean(sim$stop_fut,    na.rm = TRUE)
  p_ia    <- mean(sim$stop_ia,     na.rm = TRUE)
  p_final <- mean(sim$stop_final,  na.rm = TRUE)
  
  props <- c("All estimates" = 1, "Futile" = p_fut, "IA all estimates" = 1 - p_fut,
             "IA efficacy success" = p_ia, "Final all" = 1 - p_fut - p_ia, "Final success" = p_final)[group_names]
  
  usr <- par("usr")
  x_text <- usr[2] + text_offset_factor * (usr[2] - usr[1])
  
  for (i in seq_along(groups)) {
    text(x_text, i, sprintf("%.1f%%", 100 * props[i]), adj = 0, cex = 1.05, font = 2, 
         col = if(group_names[i]=="Futile") "firebrick" else if(group_names[i] %in% c("IA efficacy success", "Final success")) "forestgreen" else "black", xpd = TRUE)
  }
  
  text(x_text, length(groups) + 0.75, "% of trials", adj = 0, cex = 1.0, 
       font = 2, col = "gray30", xpd = TRUE)
  
  power <- p_ia + p_final
  footnote_text <- sprintf("True COR = %.2f  •  Futility at ≈%.0f%%  •  IA at ≈%.0f%%  •  Power = %.1f%%",
                           COR_true, 100 * futility_frac, 100 * info_frac, 100 * power)
  text(usr[1] + (usr[2] - usr[1])/2, usr[3] - 0.32 * (usr[4] - usr[3]), footnote_text, adj = 0.5, cex = 0.90, col = "gray40", xpd = TRUE)
  
  invisible(NULL)
}
# ============================================================
# 5) Runner (simulation only)
# ============================================================

run_all <- function(
    N_total = 600,
    futility_frac = 0.50,
    info_frac = 0.80,
    futility_p = 0.70,
    alpha1 = 0.0122,
    alpha2 = 0.0214,
    M_margin = 1.6,
    COR_true = 1.0,
    p_control_txt = "0.04, 0.02, 0.45, 0.34, 0.15",
    nSims = 2000,
    seed = 1234,
    show_progress = interactive()
) {
  p_control <- parse_probs(p_control_txt)
  msg <- validate_probs(p_control)
  if (!is.null(msg)) stop(msg)
  
  sim <- simulate_obf_ordinal(
    COR_true = COR_true,
    COR_NI = M_margin,
    n_total = N_total,
    futility_frac = futility_frac,
    info_frac = info_frac,
    alpha1 = alpha1,
    alpha2 = alpha2,
    futility_p = futility_p,
    p_control = p_control,
    seed = seed,
    nSims = nSims,
    show_progress = show_progress
  )
  
  tbl_sim <- sim_table(sim)
  
  list(sim = sim, sim_summary = tbl_sim)
}

x <-1.3
res <- run_all( nSims = 1000, COR_true = x, M_margin = 1.6)

res$sim_summary

selection_boxplot(res$sim, COR_true = x, COR_NI = 1.6,
                  point_cex       = 0.58,
                  point_alpha     = 0.20,
                  jitter_height   = 0.3,     # increase to 0.28–0.35 if you want more vertical spread
                  box_width       = 0.3,
                  ylim            = c(-1.0, 1.4),
                  text_offset_factor = 0.13)  # push % further right if needed

