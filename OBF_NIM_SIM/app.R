############################################
# Ordinal NI Group Sequential Trial Simulator
# Using rms::rcs for conditional power curve + fit quality in title
############################################

library(shiny)
library(shinyWidgets)
library(MASS)
library(dplyr)
library(bslib)
library(rpact)
library(future)
library(future.apply)
library(progressr)
library(rms)

############################################################
# Helpers
############################################################

ilogit <- function(z) 1 / (1 + exp(-z))

parse_probs <- function(txt) as.numeric(trimws(unlist(strsplit(txt, ","))))

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

pmf_from_beta <- function(theta, beta, x = 1) {
  cdf <- ilogit(theta - beta * x)
  cdf <- c(cdf, 1)
  pmf <- diff(c(0, cdf))
  pmf[pmf < 0] <- 0
  pmf / sum(pmf)
}

conditional_power_ni_cp <- function(Z_i, I_i, I_T, theta_margin_adj, zcrit_T) {
  if (!is.finite(Z_i) || !is.finite(I_i) || !is.finite(I_T) ||
      !is.finite(theta_margin_adj) || !is.finite(zcrit_T)) {
    return(NA_real_)
  }
  if (I_T <= 0 || I_i < 0 || I_T < I_i) return(NA_real_)
  
  if (abs(I_T - I_i) < 1e-12) {
    return(as.numeric(Z_i <= zcrit_T))
  }
  
  I_rem <- I_T - I_i
  arg <- (zcrit_T * sqrt(I_T) - Z_i * sqrt(I_i) - theta_margin_adj * I_rem) / sqrt(I_rem)
  pnorm(arg)
}

expected_n_breakdown <- function(sim) {
  p_fut        <- mean(sim$stop_fut, na.rm = TRUE)
  p_ia_success <- mean(sim$stop_ia,  na.rm = TRUE)
  p_low_cp     <- mean(sim$stop_fut_low_cp, na.rm = TRUE)
  p_final      <- mean(!(sim$stop_fut | sim$stop_ia | sim$stop_fut_low_cp), na.rm = TRUE)
  
  df <- data.frame(
    Stage       = c("IA1 Futility stop", "IA2 success stop", "IA2 low-CP futility stop", "Final analysis"),
    Probability = c(p_fut, p_ia_success, p_low_cp, p_final),
    N_at_stage  = c(sim$n_at_fut, sim$n_at_ia, sim$n_at_ia, sim$n_total),
    check.names = FALSE
  )
  
  df$Contribution <- df$Probability * df$N_at_stage
  
  df |>
    mutate(
      N_at_stage   = as.integer(round(N_at_stage)),
      Probability  = round(Probability, 3),
      Contribution = round(Contribution, 1)
    )
}

############################################################
# Simulator
############################################################

simulate_obf_ordinal <- function(
    COR_true, COR_NI, n_total, futility_frac, info_frac,
    zcrit1, zcrit2, futility_p, p_control, cp_threshold = 0.2,
    seed = 1234, nSims = 1000, chunk_size = 100
) {
  msg <- validate_probs(p_control)
  if (!is.null(msg)) stop(msg)
  
  theta <- theta_from_control_pmf(p_control)
  K <- length(p_control) - 1L
  y_levels <- 0:K
  
  beta_true <- log(COR_true)
  beta_NI   <- log(COR_NI)
  
  pi_control <- p_control
  pi_treat   <- pmf_from_beta(theta, beta_true)
  
  split_n <- function(N) { nC <- floor(N / 2); c(nC = nC, nT = N - nC) }
  
  n_fut <- round(futility_frac * n_total)
  n1    <- round(info_frac     * n_total)
  
  s_fut <- split_n(n_fut)
  s_ia  <- split_n(n1)
  s_tot <- split_n(n_total)
  
  z_fut <- qnorm(futility_p)
  
  trt_fut <- factor(c(rep("C", s_fut["nC"]), rep("T", s_fut["nT"])), levels = c("C","T"))
  trt_ia  <- factor(c(rep("C", s_ia ["nC"]), rep("T", s_ia ["nT"])), levels = c("C","T"))
  trt_tot <- factor(c(rep("C", s_tot["nC"]), rep("T", s_tot["nT"])), levels = c("C","T"))
  
  fit_logCOR_vec <- function(y_int, trt_fac) {
    dat <- list(
      y   = factor(y_int, ordered = TRUE, levels = y_levels),
      trt = trt_fac
    )
    
    fit <- try(
      MASS::polr(y ~ trt, data = dat, Hess = TRUE, model = FALSE, method = "logistic"),
      silent = TRUE
    )
    if (inherits(fit, "try-error") || is.null(fit$coefficients)) return(NULL)
    if (!("trtT" %in% names(fit$coefficients))) return(NULL)
    
    logCOR_hat <- unname(fit$coefficients[["trtT"]])
    v <- vcov(fit)
    se <- sqrt(v["trtT","trtT"])
    if (!is.finite(se) || se < 1e-8) return(NULL)
    
    list(logCOR_hat = logCOR_hat, se = se)
  }
  
  one_sim <- function(i) {
    out <- rep(NA_real_, 15)
    out[7:10] <- 0
    
    yCf <- sample.int(K + 1, s_fut["nC"], TRUE, pi_control) - 1L
    yTf <- sample.int(K + 1, s_fut["nT"], TRUE, pi_treat)   - 1L
    y_fut <- c(yCf, yTf)
    
    fit_f <- fit_logCOR_vec(y_fut, trt_fut)
    
    I_fut <- theta_fut <- NA_real_
    if (!is.null(fit_f)) {
      out[11] <- fit_f$logCOR_hat
      out[1]  <- (fit_f$logCOR_hat - beta_NI) / fit_f$se
      out[2]  <- exp(fit_f$logCOR_hat)
      
      I_fut     <- 1 / (fit_f$se^2)
      theta_fut <- fit_f$logCOR_hat - beta_NI
      
      if (out[1] > z_fut) { out[7] <- 1; return(out) }
    }
    
    n_add_iaC <- s_ia["nC"] - s_fut["nC"]
    n_add_iaT <- s_ia["nT"] - s_fut["nT"]
    
    if (n_add_iaC > 0) yCf <- c(yCf, sample.int(K + 1, n_add_iaC, TRUE, pi_control) - 1L)
    if (n_add_iaT > 0) yTf <- c(yTf, sample.int(K + 1, n_add_iaT, TRUE, pi_treat)   - 1L)
    
    y_ia <- c(yCf, yTf)
    fit1 <- fit_logCOR_vec(y_ia, trt_ia)
    
    I_ia <- theta_ia <- NA_real_
    if (!is.null(fit1)) {
      out[12] <- fit1$logCOR_hat
      out[3]  <- (fit1$logCOR_hat - beta_NI) / fit1$se
      out[4]  <- exp(fit1$logCOR_hat)
      
      I_ia     <- 1 / (fit1$se^2)
      theta_ia <- fit1$logCOR_hat - beta_NI
      
      if (is.finite(out[1]) && is.finite(I_fut) && n_fut > 0) {
        I_ia_plan_from_fut <- I_fut * (n1 / n_fut)
        out[14] <- conditional_power_ni_cp(
          Z_i              = out[1],
          I_i              = I_fut,
          I_T              = I_ia_plan_from_fut,
          theta_margin_adj = theta_fut,
          zcrit_T          = zcrit1
        )
      }
      
      if (is.finite(I_ia) && n1 > 0) {
        I_final_plan <- I_ia * (n_total / n1)
        out[15] <- conditional_power_ni_cp(
          Z_i              = out[3],
          I_i              = I_ia,
          I_T              = I_final_plan,
          theta_margin_adj = theta_ia,
          zcrit_T          = zcrit2
        )
      }
      
      if (out[3] <= zcrit1) { out[8] <- 1; return(out) }
      
      if (is.finite(out[15]) && out[15] < cp_threshold) {
        out[10] <- 1; return(out)
      }
    }
    
    n_add_finC <- s_tot["nC"] - s_ia["nC"]
    n_add_finT <- s_tot["nT"] - s_ia["nT"]
    
    if (n_add_finC > 0) yCf <- c(yCf, sample.int(K + 1, n_add_finC, TRUE, pi_control) - 1L)
    if (n_add_finT > 0) yTf <- c(yTf, sample.int(K + 1, n_add_finT, TRUE, pi_treat)   - 1L)
    
    y_fin <- c(yCf, yTf)
    fit2 <- fit_logCOR_vec(y_fin, trt_tot)
    
    if (!is.null(fit2)) {
      out[13] <- fit2$logCOR_hat
      out[5]  <- (fit2$logCOR_hat - beta_NI) / fit2$se
      out[6]  <- exp(fit2$logCOR_hat)
      if (out[5] <= zcrit2) out[9] <- 1
    }
    
    out
  }
  
  idx <- split(seq_len(nSims), ceiling(seq_len(nSims) / chunk_size))
  
  set.seed(seed)
  
  res_mat <- progressr::withProgressShiny(
    message = "Running simulations...",
    #detail  = "Processing...",
    value   = 0,
    {
      p <- progressr::progressor(along = idx)
      
      chunks <- future.apply::future_lapply(
        idx,
        FUN = function(ii) {
          block <- vapply(ii, one_sim, FUN.VALUE = numeric(15))
          p()
          t(block)
        },
        future.seed = TRUE
      )
      
      do.call(rbind, chunks)
    }
  )
  
  stop_fut        <- as.logical(res_mat[, 7])
  stop_ia         <- as.logical(res_mat[, 8])
  stop_final      <- as.logical(res_mat[, 9])
  stop_fut_low_cp <- as.logical(res_mat[,10])
  
  logCOR_paths <- cbind(
    fut   = res_mat[,11],
    ia    = res_mat[,12],
    final = res_mat[,13]
  )
  
  list(
    Z_fut_all   = res_mat[,1],
    COR_fut_all = res_mat[,2],
    Z1_all      = res_mat[,3],
    COR1_all    = res_mat[,4],
    Z2_all      = res_mat[,5],
    COR2_all    = res_mat[,6],
    
    stop_fut        = stop_fut,
    stop_ia         = stop_ia,
    stop_final      = stop_final,
    stop_fut_low_cp = stop_fut_low_cp,
    
    logCOR_paths = logCOR_paths,
    
    CP_after_fut_to_ia_obs   = res_mat[,14],
    CP_after_ia_to_final_obs = res_mat[,15],
    
    n_at_fut = n_fut, n_at_ia = n1, n_total = n_total, nSims = nSims,
    z_fut = z_fut, zcrit1 = zcrit1, zcrit2 = zcrit2,
    
    rpact_design = NULL
  )
}

############################################################
# Summary table (UPDATED with MSE and RMSE)
############################################################



sim_table <- function(sim, COR_true) {
  
  safe_summ_ext <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(c(N=0L, Min=NA, Max=NA, Mean=NA, Median=NA, `2.5%`=NA, `97.5%`=NA))
    q <- quantile(x, c(0.025, 0.5, 0.975), names = FALSE)
    c(N=length(x), Min=min(x), Max=max(x), Mean=mean(x), Median=q[2], `2.5%`=q[1], `97.5%`=q[3])
  }
  
  # Helper: MSE & RMSE **on the log(COR) scale**
  safe_mse_rmse_log <- function(est_log, true_log) {
    valid <- is.finite(est_log)
    if (sum(valid) == 0) return(c(MSE_log = NA_real_, RMSE_log = NA_real_))
    err  <- est_log[valid] - true_log
    mse  <- mean(err^2)
    rmse <- sqrt(mse)
    c(MSE_log = round(mse, 4), RMSE_log = round(rmse, 4))
  }
  
  # Summaries of COR values (for display columns)
  fut      <- safe_summ_ext(sim$COR_fut_all[sim$stop_fut])
  ia_suc   <- safe_summ_ext(sim$COR1_all[sim$stop_ia])
  ia_lowcp <- safe_summ_ext(sim$COR1_all[sim$stop_fut_low_cp])
  fin      <- safe_summ_ext(sim$COR2_all[sim$stop_final])
  
  # MSE / RMSE on log scale using stored logCOR_paths
  log_true <- log(COR_true)
  
  mse_fut_log      <- safe_mse_rmse_log(sim$logCOR_paths[sim$stop_fut,     "fut"],   log_true)
  mse_ia_suc_log   <- safe_mse_rmse_log(sim$logCOR_paths[sim$stop_ia,      "ia"],   log_true)
  mse_ia_lowcp_log <- safe_mse_rmse_log(sim$logCOR_paths[sim$stop_fut_low_cp, "ia"], log_true)
  mse_fin_log      <- safe_mse_rmse_log(sim$logCOR_paths[sim$stop_final,   "final"], log_true)
  
  # Build the main display table
  df <- data.frame(
    Stage      = c("IA1 Futility stop", "IA2 success stop", "IA2 low-CP futility", "Final success stop"),
    N          = c(fut["N"], ia_suc["N"], ia_lowcp["N"], fin["N"]),
    Min        = c(fut["Min"], ia_suc["Min"], ia_lowcp["Min"], fin["Min"]),
    Mean       = c(fut["Mean"], ia_suc["Mean"], ia_lowcp["Mean"], fin["Mean"]),
    Median     = c(fut["Median"], ia_suc["Median"], ia_lowcp["Median"], fin["Median"]),
    `2.5%`     = c(fut["2.5%"], ia_suc["2.5%"], ia_lowcp["2.5%"], fin["2.5%"]),
    `97.5%`    = c(fut["97.5%"], ia_suc["97.5%"], ia_lowcp["97.5%"], fin["97.5%"]),
    Max        = c(fut["Max"], ia_suc["Max"], ia_lowcp["Max"], fin["Max"]),
    MSE_log    = c(mse_fut_log["MSE_log"],   mse_ia_suc_log["MSE_log"],   
                   mse_ia_lowcp_log["MSE_log"],   mse_fin_log["MSE_log"]),
    RMSE_log   = c(mse_fut_log["RMSE_log"],  mse_ia_suc_log["RMSE_log"],  
                   mse_ia_lowcp_log["RMSE_log"],  mse_fin_log["RMSE_log"]),
    check.names = FALSE
  ) |>
    mutate(
      N = as.integer(round(N)),
      across(c(Min, Mean, Median, `2.5%`, `97.5%`, Max, MSE_log, RMSE_log), ~ round(.x, 3))
    )
  
  # Return both the table and the success-stage RMSE values for footnote use
  list(
    table             = df,
    rmse_ia_success   = mse_ia_suc_log["RMSE_log"],
    rmse_final_success = mse_fin_log["RMSE_log"]
  )
}


#}

# sim_table <- function(sim, COR_true) {
#   
#   safe_summ_ext <- function(x) {
#     x <- x[is.finite(x)]
#     if (length(x) == 0) return(c(N=0L, Min=NA, Max=NA, Mean=NA, Median=NA, `2.5%`=NA, `97.5%`=NA))
#     q <- quantile(x, c(0.025, 0.5, 0.975), names = FALSE)
#     c(N=length(x), Min=min(x), Max=max(x), Mean=mean(x), Median=q[2], `2.5%`=q[1], `97.5%`=q[3])
#   }
#   
#   safe_mse_rmse <- function(est, true_val) {
#     valid <- is.finite(est)
#     if (sum(valid) == 0) return(c(MSE = NA_real_, RMSE = NA_real_))
#     err  <- est[valid] - true_val
#     mse  <- mean(err^2)
#     rmse <- sqrt(mse)
#     c(MSE = round(mse, 4), RMSE = round(rmse, 4))
#   }
#   
#   fut      <- safe_summ_ext(sim$COR_fut_all[sim$stop_fut])
#   ia_suc   <- safe_summ_ext(sim$COR1_all[sim$stop_ia])
#   ia_lowcp <- safe_summ_ext(sim$COR1_all[sim$stop_fut_low_cp])
#   fin      <- safe_summ_ext(sim$COR2_all[sim$stop_final])
#   
#   mse_fut      <- safe_mse_rmse(sim$COR_fut_all[sim$stop_fut],      COR_true)
#   mse_ia_suc   <- safe_mse_rmse(sim$COR1_all[sim$stop_ia],          COR_true)
#   mse_ia_lowcp <- safe_mse_rmse(sim$COR1_all[sim$stop_fut_low_cp],  COR_true)
#   mse_fin      <- safe_mse_rmse(sim$COR2_all[sim$stop_final],       COR_true)
#   
#   data.frame(
#     Stage      = c("IA1 Futility stop", "IA2 success stop", "IA2 low-CP futility", "Final success stop"),
#     N          = c(fut["N"], ia_suc["N"], ia_lowcp["N"], fin["N"]),
#     Min        = c(fut["Min"], ia_suc["Min"], ia_lowcp["Min"], fin["Min"]),
#     Mean       = c(fut["Mean"], ia_suc["Mean"], ia_lowcp["Mean"], fin["Mean"]),
#     Median     = c(fut["Median"], ia_suc["Median"], ia_lowcp["Median"], fin["Median"]),
#     `2.5%`     = c(fut["2.5%"], ia_suc["2.5%"], ia_lowcp["2.5%"], fin["2.5%"]),
#     `97.5%`    = c(fut["97.5%"], ia_suc["97.5%"], ia_lowcp["97.5%"], fin["97.5%"]),
#     Max        = c(fut["Max"], ia_suc["Max"], ia_lowcp["Max"], fin["Max"]),
#     MSE        = c(mse_fut["MSE"],   mse_ia_suc["MSE"],   mse_ia_lowcp["MSE"],   mse_fin["MSE"]),
#     RMSE       = c(mse_fut["RMSE"],  mse_ia_suc["RMSE"],  mse_ia_lowcp["RMSE"],  mse_fin["RMSE"]),
#     check.names = FALSE
#   ) |>
#     mutate(
#       N = as.integer(round(N)),
#       across(c(Min, Mean, Median, `2.5%`, `97.5%`, Max, MSE, RMSE), ~ round(.x, 3))
#     )
#}

############################################################
# Plots (selection_boxplot and cp_selection_boxplot unchanged)
############################################################

selection_boxplot <- function(sim, COR_true, COR_NI, futility_frac, info_frac,
                              show_traj_success = FALSE, show_traj_fail = FALSE,
                              use_cor_scale = FALSE,
                              xlim_log_low = -3, xlim_log_high = 4,
                              main = "Bias and Treatment Effect Estimates by Trial Stopping Stage") {
  
  log_min_nice <- xlim_log_low
  log_max_nice <- xlim_log_high
  
  if (use_cor_scale) {
    xlim_use     <- exp(c(log_min_nice, log_max_nice))
    x_transform  <- exp
    x_label_expr <- "Cumulative Odds Ratio (COR)"
  } else {
    xlim_use     <- c(log_min_nice, log_max_nice)
    x_transform  <- identity
    x_label_expr <- expression(log(Cumulative~Odds~Ratio))
  }
  
  p_fut        <- mean(sim$stop_fut, na.rm = TRUE)
  p_ia_success <- mean(sim$stop_ia,  na.rm = TRUE)
  p_low_cp     <- mean(sim$stop_fut_low_cp, na.rm = TRUE)
  p_final_suc  <- mean(sim$stop_final, na.rm = TRUE)
  
  empirical_power <- p_ia_success + p_final_suc
  p_reach_final <- mean(!(sim$stop_fut | sim$stop_ia | sim$stop_fut_low_cp), na.rm = TRUE)
  
  avg_n <- round(
    p_fut         * sim$n_at_fut +
      p_ia_success * sim$n_at_ia +
      p_low_cp     * sim$n_at_ia +
      p_reach_final* sim$n_total
  )
  
  n_sims_used <- sim$nSims
  k_success   <- round(empirical_power * n_sims_used)
  
  z <- qnorm(0.975)
  n <- n_sims_used
  center <- (k_success + z^2 / 2) / (n + z^2)
  margin <- z * sqrt((k_success/n * (1 - k_success/n) + z^2/(4*n)) / (n + z^2))
  ci_lower <- max(0, center - margin)
  ci_upper <- min(1, center + margin)
  
  power_text <- sprintf(
    "- Empirical Power: %.1f%% (%.1f%%–%.1f%%) (IA success + Final success) across %d simulations",
    empirical_power * 100, ci_lower * 100, ci_upper * 100, n_sims_used
  )
  
  groups <- list(
    "1 All @ futility"        = sim$logCOR_paths[, "fut"][is.finite(sim$logCOR_paths[, "fut"])],
    "1 Stopped futility"      = sim$logCOR_paths[, "fut"][sim$stop_fut & is.finite(sim$logCOR_paths[, "fut"])],
    "2 All @ interim"         = sim$logCOR_paths[, "ia"][is.finite(sim$logCOR_paths[, "ia"])],
    "2 Stopped IA success"    = sim$logCOR_paths[, "ia"][sim$stop_ia & is.finite(sim$logCOR_paths[, "ia"])],
    "2 Stopped low CP at IA2" = sim$logCOR_paths[, "ia"][sim$stop_fut_low_cp & is.finite(sim$logCOR_paths[, "ia"])],
    "3 All @ final"           = sim$logCOR_paths[, "final"][is.finite(sim$logCOR_paths[, "final"])],
    "3 final failed"          = sim$logCOR_paths[, "final"][!sim$stop_final & is.finite(sim$logCOR_paths[, "final"])],
    "3 Stopped final success" = sim$logCOR_paths[, "final"][sim$stop_final & is.finite(sim$logCOR_paths[, "final"])]
  )
  
  counts_actual <- sapply(groups, length)
  if (use_cor_scale) groups <- lapply(groups, exp)
  
  keep <- sapply(groups, function(g) length(g) > 0 && all(is.finite(g)))
  groups <- groups[keep]
  group_names <- names(groups)
  counts_actual <- counts_actual[keep]
  
  set.seed(202506)
  n_max <- nrow(sim$logCOR_paths)
  jitter_master <- runif(n_max, -0.30, 0.30)
  
  op <- par(no.readonly = TRUE)
  on.exit({ par(op); layout(1) }, add = TRUE)
  layout(matrix(c(1, 2), nrow = 1), widths = c(5.2, 2.2))
  
  par(mar = c(15, 12, 6, 2), xpd = FALSE)
  
  plot(
    0, type = "n",
    xlim = xlim_use,
    ylim = c(0.4, length(groups) + 0.9),
    xlab = "", ylab = "",
    yaxt = "n", las = 1,
    main = sprintf("%s\n(Expected Sample Size = %d)", main, avg_n)
  )
  
  footnote_cex <- 1.05
  mtext(x_label_expr, side = 1, line = 3.5, adj = 0.5, font = 2, cex = 1.2)
  
  mtext(
    sprintf(
      "- Futility Look: N = %d (%d%% info) | IA Success Look: N = %d (%d%% info) | Final: N = %d",
      sim$n_at_fut, round(futility_frac * 100),
      sim$n_at_ia,  round(info_frac * 100),
      sim$n_total
    ),
    side = 1, line = 7, adj = 0, cex = footnote_cex
  )
  
  mtext(
    sprintf(
      "- Vertical Lines: Green dashed = True Effect (%.2f) | Red dotted = NI Margin (%.2f)",
      if (use_cor_scale) COR_true else log(COR_true),
      if (use_cor_scale) COR_NI   else log(COR_NI)
    ),
    side = 1, line = 8.5, adj = 0, cex = footnote_cex
  )
  
  mtext(power_text, side = 1, line = 10, adj = 0, cex = footnote_cex)
  
  mtext(
    "- Expected Sample Size (ESS) is the average N across all sims, accounting for early stopping.",
    side = 1, line = 11.5, adj = 0, cex = footnote_cex
  )
  
  axis(2, at = seq_along(groups), labels = group_names, las = 1, cex.axis = 1.1)
  
  abline(v = x_transform(log(COR_true)), lty = 2, col = "darkgreen", lwd = 2.5)
  abline(v = x_transform(log(COR_NI)),   lty = 3, col = "red",       lwd = 2.5)
  
  abline(h = seq_along(groups), col = "gray92", lwd = 0.8)
  
  if (show_traj_success || show_traj_fail) {
    group_y_map <- setNames(seq_along(groups), names(groups))
    
    for (i in seq_len(nrow(sim$logCOR_paths))) {
      is_success <- isTRUE(sim$stop_ia[i]) || isTRUE(sim$stop_final[i])
      is_low_cp  <- isTRUE(sim$stop_fut_low_cp[i])
      
      if ((is_success && !show_traj_success) ||
          (!is_success && !show_traj_fail)) next
      
      path_x <- numeric(0)
      path_y <- numeric(0)
      j <- jitter_master[i]
      
      if (is.finite(sim$logCOR_paths[i, "fut"]) &&
          "1 All @ futility" %in% names(group_y_map)) {
        path_x <- c(path_x, sim$logCOR_paths[i, "fut"])
        path_y <- c(path_y, group_y_map["1 All @ futility"] + j)
      }
      
      if (is.finite(sim$logCOR_paths[i, "ia"])) {
        if (is_low_cp && "2 Stopped low CP at IA2" %in% names(group_y_map)) {
          g <- "2 Stopped low CP at IA2"
        } else if (sim$stop_ia[i] && "2 Stopped IA success" %in% names(group_y_map)) {
          g <- "2 Stopped IA success"
        } else if ("2 All @ interim" %in% names(group_y_map)) {
          g <- "2 All @ interim"
        } else {
          g <- NULL
        }
        
        if (!is.null(g)) {
          path_x <- c(path_x, sim$logCOR_paths[i, "ia"])
          path_y <- c(path_y, group_y_map[g] + j)
        }
      }
      
      if (is.finite(sim$logCOR_paths[i, "final"])) {
        if (sim$stop_final[i] && "3 Stopped final success" %in% names(group_y_map)) {
          g_final <- "3 Stopped final success"
        } else if (!sim$stop_final[i] && "3 final failed" %in% names(group_y_map)) {
          g_final <- "3 final failed"
        } else if ("3 All @ final" %in% names(group_y_map)) {
          g_final <- "3 All @ final"
        } else {
          g_final <- NULL
        }
        
        if (!is.null(g_final)) {
          path_x <- c(path_x, sim$logCOR_paths[i, "final"])
          path_y <- c(path_y, group_y_map[g_final] + j)
        }
      }
      
      if (length(path_x) >= 2) {
        col_line <- if (is_success) rgb(0, 0.7, 0, 0.18)
        else if (is_low_cp) rgb(0.6, 0, 0.8, 0.22)
        else rgb(0.9, 0.1, 0.1, 0.18)
        
        lines(x_transform(path_x), path_y, col = col_line, lwd = 1)
      }
    }
  }
  
  for (ii in seq_along(groups)) {
    nm <- group_names[ii]
    
    if      (grepl("1 All @ futility", nm))        idx <- which(is.finite(sim$logCOR_paths[, "fut"]))
    else if (grepl("1 Stopped futility", nm))      idx <- which(sim$stop_fut & is.finite(sim$logCOR_paths[, "fut"]))
    else if (grepl("2 All @ interim", nm))         idx <- which(is.finite(sim$logCOR_paths[, "ia"]))
    else if (grepl("2 Stopped IA success", nm))    idx <- which(sim$stop_ia & is.finite(sim$logCOR_paths[, "ia"]))
    else if (grepl("2 Stopped low CP at IA2", nm)) idx <- which(sim$stop_fut_low_cp & is.finite(sim$logCOR_paths[, "ia"]))
    else if (grepl("3 All @ final", nm))           idx <- which(is.finite(sim$logCOR_paths[, "final"]))
    else if (grepl("3 final failed", nm))          idx <- which(!sim$stop_final & is.finite(sim$logCOR_paths[, "final"]))
    else if (grepl("3 Stopped final success", nm)) idx <- which(sim$stop_final & is.finite(sim$logCOR_paths[, "final"]))
    else next
    
    vals <- groups[[ii]]
    jitter_y <- jitter_master[idx]
    
    col_p <- if (grepl("All @ futility", nm)) rgb(0.1, 0.4, 0.9, 0.25) else
      if (grepl("futility", nm))              rgb(0.9, 0.15,0.15,0.25) else
        if (grepl("failed", nm))                rgb(0.9, 0.15,0.15,0.25) else
          if (grepl("success", nm))               rgb(0.1, 0.65,0.1, 0.25) else
            if (grepl("low CP", nm))                rgb(0.6, 0.0, 0.8, 0.35) else
              rgb(0.3, 0.3, 0.3, 0.25)
    
    points(vals, ii + jitter_y, pch = 19, cex = 0.6, col = col_p)
    
    if (length(vals) >= 3) {
      q <- quantile(vals, c(0.25, 0.5, 0.75))
      iqr <- q[3] - q[1]
      
      w_upper <- max(vals[vals <= q[3] + 1.5 * iqr])
      w_lower <- min(vals[vals >= q[1] - 1.5 * iqr])
      
      segments(w_lower, ii, q[1], ii, col = "steelblue", lwd = 1.2)
      segments(q[3],   ii, w_upper, ii, col = "steelblue", lwd = 1.2)
      segments(w_lower, ii - 0.05, w_lower, ii + 0.05, col = "steelblue", lwd = 1.2)
      segments(w_upper, ii - 0.05, w_upper, ii + 0.05, col = "steelblue", lwd = 1.2)
      rect(q[1], ii - 0.14, q[3], ii + 0.14, col = rgb(0.88,0.93,1,0.5), border = "steelblue", lwd = 1.4)
      segments(q[2], ii - 0.14, q[2], ii + 0.14, lwd = 5, col = "royalblue3")
    }
  }
  
  par(mar = c(15, 1, 6, 1), xpd = FALSE)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0.4, length(groups) + 0.9))
  
  header_cex <- 1.15
  text(0.02, length(groups) + 0.75, "Trials (N)", adj = 0, font = 2, cex = header_cex)
  text(0.45, length(groups) + 0.75, "N/Trial",   adj = 0, font = 2, cex = header_cex)
  text(0.78, length(groups) + 0.75, "% Sims",    adj = 0, font = 2, cex = header_cex)
  
  data_cex <- 1.10
  
  n_per_group <- rep(NA, length(groups))
  for (k in seq_along(groups)) {
    nm <- names(groups)[k]
    if (grepl("futility", nm)) {
      n_per_group[k] <- sim$n_at_fut
    } else if (grepl("interim|IA", nm)) {
      n_per_group[k] <- sim$n_at_ia
    } else {
      n_per_group[k] <- sim$n_total
    }
  }
  
  props <- numeric(length(groups))
  for (k in seq_along(groups)) {
    nm <- names(groups)[k]
    if      (grepl("^1 All @ futility", nm))         props[k] <- 1
    else if (grepl("^1 Stopped futility", nm))       props[k] <- mean(sim$stop_fut, na.rm = TRUE)
    else if (grepl("^2 All @ interim", nm))          props[k] <- mean(!sim$stop_fut, na.rm = TRUE)
    else if (grepl("^2 Stopped IA success", nm))     props[k] <- mean(sim$stop_ia, na.rm = TRUE)
    else if (grepl("^2 Stopped low CP at IA2", nm))  props[k] <- mean(sim$stop_fut_low_cp, na.rm = TRUE)
    else if (grepl("^3 All @ final", nm))            props[k] <- mean(!(sim$stop_fut | sim$stop_ia | sim$stop_fut_low_cp), na.rm = TRUE)
    else if (grepl("^3 final failed", nm))           props[k] <- mean(!(sim$stop_fut | sim$stop_ia | sim$stop_fut_low_cp) & !sim$stop_final, na.rm = TRUE)
    else if (grepl("^3 Stopped final success", nm))  props[k] <- mean(sim$stop_final, na.rm = TRUE)
    else props[k] <- NA
  }
  
  for (ii in seq_along(groups)) {
    col_t <- if (grepl("All @ futility", names(groups)[ii])) "black" else
      if (grepl("futility", names(groups)[ii])) "firebrick" else
        if (grepl("failed", names(groups)[ii])) "firebrick" else
          if (grepl("success", names(groups)[ii])) "forestgreen" else
            if (grepl("low CP", names(groups)[ii])) "firebrick" else
              "gray30"
    
    text(0.02, ii, format(counts_actual[ii], big.mark=","), adj = 0, col = col_t, cex = data_cex)
    text(0.45, ii, sprintf("%d", n_per_group[ii]),           adj = 0, col = "gray30", cex = data_cex)
    text(0.78, ii, sprintf("%.1f%%", 100 * props[ii]),       adj = 0, font = 2, col = col_t, cex = data_cex)
  }
}

cp_selection_boxplot <- function(sim, COR_true, COR_NI,
                                 use_cor_scale = FALSE,
                                 xlim_log_low = -3, xlim_log_high = 4,
                                 main_prefix = "'2 All @ interim' stage– coloured by conditional probability") {
  
  log_min_nice <- xlim_log_low
  log_max_nice <- xlim_log_high
  
  if (use_cor_scale) {
    xlim_use     <- exp(c(log_min_nice, log_max_nice))
    x_transform  <- exp
    x_label_expr <- "Cumulative Odds Ratio (COR)"
  } else {
    xlim_use     <- c(log_min_nice, log_max_nice)
    x_transform  <- identity
    x_label_expr <- expression(log(Cumulative~Odds~Ratio))
  }
  
  logcor_ia <- sim$logCOR_paths[, "ia"][is.finite(sim$logCOR_paths[, "ia"])]
  vals <- if (use_cor_scale) exp(logcor_ia) else logcor_ia
  
  idx <- which(is.finite(sim$logCOR_paths[, "ia"]))
  cp_ia <- sim$CP_after_ia_to_final_obs[idx]
  
  n_points <- length(vals)
  n_at_ia  <- sim$n_at_ia
  
  main <- sprintf("%s\n(%d trials reaching IA2, N/trial = %d)", 
                  main_prefix, n_points, n_at_ia)
  
  if (n_points == 0) {
    plot(0, type = "n", main = "No finite estimates at interim", xaxt = "n", yaxt = "n")
    return(invisible(NULL))
  }
  
  colramp <- colorRampPalette(c("red", "yellow", "green"))(101)
  cp_norm <- pmax(0, pmin(1, cp_ia))
  col_idx <- round(cp_norm * 100) + 1
  point_colors <- colramp[col_idx]
  
  set.seed(202506)
  n_max <- nrow(sim$logCOR_paths)
  jitter_master <- runif(n_max, -0.30, 0.30)
  jitter_y <- jitter_master[idx]
  
  op <- par(no.readonly = TRUE)
  on.exit({ par(op); layout(1) }, add = TRUE)
  
  layout(matrix(c(1, 2), nrow = 1), widths = c(6, 1.4))
  
  par(mar = c(15, 12, 8, 2), xpd = FALSE)
  
  plot(
    0, type = "n",
    xlim = xlim_use,
    ylim = c(0.4, 1.9),
    xlab = "", ylab = "",
    yaxt = "n", las = 1,
    main = main
  )
  
  mtext(x_label_expr, side = 1, line = 3.5, adj = 0.5, font = 2, cex = 1.2)
  
  mtext(
    sprintf("- True effect: %.2f (green dashed) | NI margin: %.2f (red dotted)",
            if (use_cor_scale) COR_true else log(COR_true),
            if (use_cor_scale) COR_NI   else log(COR_NI)),
    side = 1, line = 8.5, adj = 0, cex = 1.0
  )
  
  axis(2, at = 1, labels = "2 All @ interim", las = 1, cex.axis = 1.1)
  
  abline(v = x_transform(log(COR_true)), lty = 2, col = "darkgreen", lwd = 2.5)
  abline(v = x_transform(log(COR_NI)),   lty = 3, col = "red",       lwd = 2.5)
  
  abline(h = 1, col = "gray92", lwd = 0.8)
  
  points(vals, 1 + jitter_y,
         pch = 19,
         cex = 1.0,
         col = adjustcolor(point_colors, alpha.f = 0.8))
  
  if (length(vals) >= 3) {
    q <- quantile(vals, c(0.25, 0.5, 0.75))
    iqr <- q[3] - q[1]
    w_upper <- max(vals[vals <= q[3] + 1.5 * iqr])
    w_lower <- min(vals[vals >= q[1] - 1.5 * iqr])
    
    segments(w_lower, 1, q[1], 1, col = "steelblue", lwd = 1.2)
    segments(q[3],   1, w_upper, 1, col = "steelblue", lwd = 1.2)
    segments(w_lower, 1 - 0.05, w_lower, 1 + 0.05, col = "steelblue", lwd = 1.2)
    segments(w_upper, 1 - 0.05, w_upper, 1 + 0.05, col = "steelblue", lwd = 1.2)
    rect(q[1], 1 - 0.14, q[3], 1 + 0.14, col = rgb(0.88,0.93,1,0.5), border = "steelblue", lwd = 1.4)
    segments(q[2], 1 - 0.14, q[2], 1 + 0.14, lwd = 5, col = "royalblue3")
  }
  
  par(mar = c(15, 1, 8, 4), xpd = TRUE)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  
  n_colors <- 100
  y_pos <- seq(0.15, 0.85, length.out = n_colors)
  for (i in seq_len(n_colors)) {
    rect(0.1, y_pos[i], 0.4, y_pos[i] + (0.7 / n_colors),
         col = colramp[i], border = NA)
  }
  
  text(0.5, 0.90, "", adj = 0.5, font = 2, cex = 1.1) #Conditional Power
  text(0.5, 0.10, "0", adj = 0.5, cex = 1.0)
  text(0.5, 0.90, "1", adj = 0.5, cex = 1.0)
  text(0.5, 0.50, "0.5", adj = 0.5, cex = 1.0)
  
  arrows(0.25, 0.12, 0.25, 0.05, length = 0.08, col = "gray40")
  arrows(0.25, 0.88, 0.25, 0.95, length = 0.08, col = "gray40")
  
  par(op)
  
}

############################################################
# UI
############################################################

ui <- page_sidebar(
  title = "Ordinal Endpoint, Non Inferiority, Group Sequential Trial Simulator v4.0",
  sidebar = sidebar(
    width = 350,
    
    actionButton(
      "run_btn", "Run Simulation",
      class = "btn-primary btn-lg",
      icon  = icon("play-circle"),
      width = "100%"
    ),
    
    hr(style = "margin: 1.2em 0;"),
    
    tags$div(
      style = "padding: 0 8px;",
      h5("Trial & Simulation Settings"),
      
      numericInput("n_sims",  "Number of simulations", value = 1000),
      numericInput("n_total", "Total sample size (1:1 rand.)", value = 600),
      textInput("p_control_txt", "Control probabilities",
                value = "0.04, 0.02, 0.45, 0.34, 0.15"),
      
      numericInput("COR_true", "True cumulative odds ratio (COR)", value = 1.0, step = 0.05),
      numericInput("COR_NI",   "Non inferiority margin COR",       value = 1.6, step = 0.1),
      
      numericInput("futility_p",    "Futility p-value threshold IA1", value = 0.70),
      numericInput("cp_threshold",  "CP futility threshold at IA2",   value = 0.1, min = 0, max = 1, step = 0.05),
      
      sliderInput("futility_frac", "Futility look fraction", min = 0.2, max = 0.7,  value = 0.5),
      sliderInput("info_frac",     "Interim look fraction",  min = 0.5, max = 0.95, value = 0.80),
      
      numericInput("chunk_size", "Chunk size (speed for large sims)", value = 100, min = 10, step = 10),
      
      numericInput("seed", "Random seed", value = 202506)
    ),
    
    hr(style = "margin: 1.5em 0;"),
    
    tags$div(
      style = "padding: 0 8px;",
      h5("Plot Options"),
      
      checkboxInput("use_cor_scale", "Display on COR scale (instead of log)", value = FALSE),
      sliderInput("xlim_log_low",  "X-axis lower limit (log scale)", min = -6, max = 0, value = -1, step = 0.1),
      sliderInput("xlim_log_high", "X-axis upper limit (log scale)", min = 0,  max = 7, value = 2,  step = 0.1),
      
      hr(style = "margin: 1.2em 0; border-top: 1px dashed #ccc;"),
      tags$strong("Show trajectories"),
      checkboxInput("show_traj_success", "Successful trajectories", value = FALSE),
      checkboxInput("show_traj_fail",    "Failed / futility trajectories", value = FALSE),
      
      hr(),
      
      downloadButton("download_plot", "Download Superb JPEG",
                     class = "btn-success", style = "width: 100%;")
    )
  ),
  
  card(
    tabsetPanel(
      tabPanel("Operating Characteristics",
               plotOutput("boxplot", height = "750px")),
      
      tabPanel("Cumulative Odds Ratio Distributions",
  verbatimTextOutput("status"),
  tableOutput("summary_table"),
  uiOutput("success_precision_note"),
  # ← Put the explanation here
  tags$div(
    style = "font-size: 0.9em; color: #555; margin-top: 10px; padding: 8px; border-left: 3px solid #ccc;",
    tags$p(
      tags$strong("Understanding RMSE_log (precision on the log(COR) scale):"),
      tags$br(),
      "RMSE_log tells you the typical size of error in log(estimated COR). ",
      "Because we exponentiate to get back to the COR scale, this error becomes a ",
      tags$strong("multiplicative factor"), " around the true value."
    ),
    tags$p(
      "Examples of what different RMSE_log values mean in practice:"
    ),
    tags$ul(
      style = "margin-left: 20px; margin-top: 6px; line-height: 1.5;",
      tags$li(
        tags$strong("RMSE_log = 0.10"), " → typical range: ×", tags$strong("0.90 – 1.11"),
        " (most estimates are within about ±10–11% of the true COR)"
      ),
      tags$li(
        tags$strong("RMSE_log = 0.25"), " → typical range: ×", tags$strong("0.78 – 1.28"),
        " (estimates usually fall within ±22–28% of the true COR — moderate/wide)"
      ),
      tags$li(
        tags$strong("RMSE_log = 0.50"), " → typical range: ×", tags$strong("0.61 – 1.65"),
        " (estimates can easily be 40% lower or 65% higher than true — very imprecise)"
      ),
      tags$li(
        tags$strong("RMSE_log = 0.75"), " → typical range: ×", tags$strong("0.47 – 2.12"),
        " (estimates may be half or more than double the true value — extremely wide uncertainty)"
      )
    ),
    tags$p(
      tags$small(
        "These intervals approximate a 68% range assuming roughly normal errors on the log scale. ",
        "Smaller RMSE_log = tighter, more reliable estimates."
      )
    )
  ),
  
  hr(),
  h5("rpact Design & Nominal P-values"),
  verbatimTextOutput("rpact_info")
),
      
      tabPanel("Expected Sample Size",
               tableOutput("ess_breakdown"),
               hr(),
               verbatimTextOutput("ess_total_note")),
      
      tabPanel("Conditional Probability",
               
               plotOutput("cp_boxplot", height = "600px"),
               
               hr(),
               h4("Conditional Probability vs log(COR) at IA2"),
               
               fluidRow(
                 column(4,
                        numericInput("cor_input",
                                     "Input COR for prediction (e.g. 0.6, 1.0, 1.6)",
                                     value = 1.0,
                                     min = 0.01, step = 0.05)
                 ),
                 column(4,
                        numericInput("spline_df",
                                     "Degrees of freedom for RCS (higher = more wiggly curve)",
                                     value = 4,
                                     min = 3, max = 8, step = 1)
                 ),
                 
                 textInput(
                   "custom_knots_logcor",
                   label = "Knot percentiles (of observed log(COR) at IA)",
                   value = "25,50,75,90,95,97,99",
                   placeholder = "e.g. 25,50,75,90,95,97,99 or -2,-1,0,1,2   (blank = auto quantiles)"
                 )
               ),
               
               plotOutput("cp_rcs_plot", height = "500px"),
               
               hr(),
               h5("Predicted CP for Input log(COR)"),
               verbatimTextOutput("cp_prediction")
      ),
      
      tabPanel("Wiki",
               div(style = "padding: 25px; max-width: 1000px;",
                   h2("Ordinal Non-Inferiority Group Sequential Trial Simulator"),
                   hr(),
                   h3("1) What This App Does"),
                   p("This application simulates a three-stage group sequential non-inferiority trial with an ordinal endpoint analyzed using a proportional odds model."),
                   tags$ul(
                     tags$li("Stage 1: Futility look (based on a user-defined p-value threshold)."),
                     tags$li("Stage 2: Interim efficacy look (O’Brien-Fleming alpha spending)."),
                     tags$li("Stage 3: Final analysis (if trial continues).")
                   ),
                   br(),
                   h4("Interpretation of the Cumulative Odds Ratio (COR)"),
                   tags$ul(
                     tags$li(strong("COR = 1"), " → No difference between treatment and control."),
                     tags$li(strong("COR > 1"), " → Treatment is worse than comparator."),
                     tags$li(strong("COR < 1"), " → Treatment is better than comparator.")
                   ),
                   p("Because this is a non-inferiority design, the red vertical line represents the NI margin. If the estimated log(COR) is sufficiently below that margin, non-inferiority is declared."),
                   br(),
                   h4("Futility at IA1 (P-value approach)"),
                   p("At the first look, the trial may stop for futility if the observed one-sided p-value is too large."),
                   br(),
                   h4("Conditional Power (CP) Futility at IA2"),
                   p("At the second look, the trial may stop for futility based on Conditional Power."),
                   hr(),
                   p(style = "font-style: italic; color: gray;",
                     "This simulator is for educational and design exploration purposes only.")
               ))
    )
  )
)

############################################################
# Server
############################################################

server <- function(input, output, session) {
  
  future::plan(future::multisession, workers = max(1, future::availableCores() - 1))
  onStop(function() future::plan(future::sequential))
  
  sim <- eventReactive(input$run_btn, {
    
    design <- getDesignGroupSequential(
      sided = 1,
      alpha = 0.025,
      informationRates = c(input$info_frac, 1),
      typeOfDesign = "asOF"
    )
    
    ans <- simulate_obf_ordinal(
      COR_true = input$COR_true,
      COR_NI   = input$COR_NI,
      n_total  = input$n_total,
      futility_frac = input$futility_frac,
      info_frac     = input$info_frac,
      zcrit1 = -design$criticalValues[1],
      zcrit2 = -design$criticalValues[2],
      futility_p = input$futility_p,
      p_control  = parse_probs(input$p_control_txt),
      cp_threshold = input$cp_threshold,
      seed  = input$seed,
      nSims = input$n_sims,
      chunk_size = input$chunk_size
    )
    
    ans$rpact_design <- design
    ans
    
  }, ignoreInit = TRUE)
  
  cp_rcs_fit <- reactive({
    s <- sim()
    req(s)
    
    df <- data.frame(
      logCOR = s$logCOR_paths[, "ia"],
      CP     = s$CP_after_ia_to_final_obs
    ) %>%
      filter(is.finite(logCOR) & is.finite(CP) & CP >= 0 & CP <= 1)
    
    shiny::validate(
      need(nrow(df) >= 30,
           paste0("Only ", nrow(df), " valid IA2 observations.\n",
                  "Need ≥30 for reliable spline fit.\n",
                  "Try: COR_true closer to NI margin, higher futility_p, more sims"))
    )
    
    df_val <- max(3L, min(8L, as.integer(input$spline_df %||% 4)))
    
    knots_text <- trimws(input$custom_knots_logcor %||% "")
    use_custom <- nzchar(knots_text)
    
    knots <- NULL
    knots_source <- "auto (quantiles)"
    actual_df <- df_val
    
    if (use_custom) {
      perc_vals <- suppressWarnings(as.numeric(unlist(strsplit(knots_text, "[, ]+"))))
      
      if (all(is.finite(perc_vals)) && all(perc_vals >= 0) && all(perc_vals <= 100) &&
          length(perc_vals) >= 3 && length(perc_vals) <= 7 &&
          length(unique(perc_vals)) == length(perc_vals)) {
        
        perc_sorted <- sort(perc_vals / 100)
        knots <- quantile(df$logCOR, probs = perc_sorted, names = FALSE)
        knots_source <- sprintf("custom percentiles: %s", paste(perc_vals, collapse = ", "))
        
      } else {
        knots <- as.numeric(unlist(strsplit(knots_text, "[, ]+")))
        
        if (any(!is.finite(knots)) || length(knots) < 2 ||
            any(diff(sort(knots)) <= 0)) {
          showNotification(
            "Invalid knots. Use comma-separated increasing numbers (direct log(COR)) or percentiles (0–100).",
            type = "error"
          )
          knots <- NULL
          use_custom <- FALSE
        } else {
          knots_source <- sprintf("custom log(COR): %s", paste(round(knots, 2), collapse = ", "))
        }
      }
      
      if (!is.null(knots)) {
        actual_df <- length(knots) - 1L
        if (actual_df < 2 || actual_df > 7) {
          showNotification(
            sprintf("With %d knots → df = %d (must be 2–7). Using automatic instead.", length(knots), actual_df),
            type = "warning"
          )
          knots <- NULL
          use_custom <- FALSE
          actual_df <- df_val
        }
      }
    }
    
    # ── Critical fix using <<- to place dd in global environment ─────────────
    dd <<- datadist(df)
    options(datadist = "dd")
    
    if (!is.null(knots) && use_custom) {
      fmla <- CP ~ rcs(logCOR, parms = knots)
    } else {
      nk <- df_val + 1L
      fmla <- as.formula(sprintf("CP ~ rcs(logCOR, %d)", nk))
    }
    
    fit <- tryCatch(
      ols(fmla, data = df, x = TRUE, y = TRUE),
      error = function(e) {
        showNotification(paste("Spline fit failed:", e$message), type = "error")
        NULL
      }
    )
    
    req(fit, "Restricted cubic spline did not converge")
    
    list(
      fit         = fit,
      data        = df,
      n_valid     = nrow(df),
      df          = actual_df,
      knots       = knots,
      knots_source = knots_source
    )
  })
  
  output$cp_rcs_plot <- renderPlot({
    input$spline_df          # trigger
    input$custom_knots_logcor # trigger
    
    obj <- cp_rcs_fit()
    req(obj, obj$fit, obj$data)
    
    df  <- obj$data
    fit <- obj$fit
    
    resids <- residuals(fit)
    mse    <- mean(resids^2, na.rm = TRUE)
    rse    <- sqrt(mse)
    
    rng <- range(df$logCOR, finite = TRUE)
    pad <- 0.08 * diff(rng)
    log_seq <- seq(rng[1] - pad, rng[2] + pad, length.out = 250)
    
    pred <- Predict(fit, logCOR = log_seq, conf.int = 0.95)
    pred <- pred[order(pred$logCOR), ]
    
    ok <- is.finite(pred$yhat) & is.finite(pred$lower) & is.finite(pred$upper)
    if (sum(ok) < 30) {
      plot(0, type = "n", main = "Prediction range invalid – check knots / data")
      return()
    }
    
    log_seq <- pred$logCOR[ok]
    fitv    <- pred$yhat[ok]
    lwr     <- pred$lower[ok]
    upr     <- pred$upper[ok]
    
    knot_txt <- if (!is.null(obj$knots)) {
      paste(round(obj$knots, 2), collapse = ", ")
    } else {
      "auto quantiles"
    }
    
    plot(
      range(log_seq), c(0, 1),
      type = "n",
      xlab = expression(log(Cumulative~Odds~Ratio)),
      ylab = "Conditional Power to final analysis",
      main = sprintf(
        "RCS fit to observed CP at IA2\n(n = %d reaching IA2, effective df = %d)\nKnots: %s\nResidual SE = %.4f   MSE = %.4f",
        obj$n_valid, obj$df, obj$knots_source, rse, mse
      ),
      las = 1, cex.main = 1.05, cex.lab = 1.1
    )
    
    polygon(c(log_seq, rev(log_seq)), c(lwr, rev(upr)),
            col = rgb(0.75, 0.88, 1, 0.4), border = NA)
    
    lines(log_seq, fitv, col = "blue3", lwd = 3)
    
    points(df$logCOR, df$CP, pch = 16, cex = 0.75,
           col = rgb(0,0,0, 0.3))
    
    abline(v = log(input$cor_input), col = "red2", lty = 2, lwd = 2.2)
    if (!is.null(obj$knots)) {
      abline(v = obj$knots, col = "darkorange", lty = 3, lwd = 1.4)
    }
    
    grid(col = "gray88", lty = "dotted")
    abline(h = seq(0,1,0.2), col = "gray90", lty = 3)
    
    legend("topright", inset = 0.02, bty = "n", cex = 0.95,
           legend = c("Fitted RCS", "95% pointwise CI", "Observed CP", "Input COR", "Knots (if custom)"),
           col = c("blue3", rgb(0.75,0.88,1,0.7), rgb(0,0,0,0.4), "red2", "darkorange"),
           lwd = c(3, NA, NA, 2, 1.4), lty = c(1, NA, NA, 2, 3),
           pch = c(NA, NA, 16, NA, NA), pt.bg = c(NA, rgb(0.75,0.88,1,0.7), NA, NA, NA))
  })
  
  output$cp_prediction <- renderPrint({
    obj <- cp_rcs_fit()
    req(obj, obj$fit)
    
    cor_input    <- input$cor_input
    logcor_input <- log(cor_input)
    
    if (!is.finite(logcor_input) || cor_input <= 0) {
      cat("Invalid COR value (must be positive)\n")
      return()
    }
    
    pred_df <- Predict(
      obj$fit,
      logCOR = logcor_input,
      conf.int = 0.95
    )
    
    if (nrow(pred_df) == 0 || !is.finite(pred_df$yhat)) {
      cat("Prediction failed (value out of range or model issue)\n")
      return()
    }
    
    est  <- pred_df$yhat
    lwr  <- pred_df$lower
    upr  <- pred_df$upper
    
    cat(sprintf("Input COR           = %.3f\n", cor_input))
    cat(sprintf("Corresponding log(COR) = %.3f\n", logcor_input))
    cat(sprintf("Predicted CP        : %.3f\n", est))
    cat(sprintf("95%% CI for CP       : %.3f – %.3f\n", max(0, lwr), min(1, upr)))
    cat(sprintf("(based on %d simulations reaching IA2)\n", obj$n_valid))
  })
  
  output$boxplot <- renderPlot({
    req(sim())
    selection_boxplot(
      sim(),
      COR_true = input$COR_true,
      COR_NI   = input$COR_NI,
      futility_frac = input$futility_frac,
      info_frac     = input$info_frac,
      show_traj_success = input$show_traj_success,
      show_traj_fail    = input$show_traj_fail,
      use_cor_scale     = input$use_cor_scale,
      xlim_log_low      = input$xlim_log_low,
      xlim_log_high     = input$xlim_log_high
    )
  })
  
  output$cp_boxplot <- renderPlot({
    req(sim())
    cp_selection_boxplot(
      sim(),
      COR_true = input$COR_true,
      COR_NI   = input$COR_NI,
      use_cor_scale     = input$use_cor_scale,
      xlim_log_low      = input$xlim_log_low,
      xlim_log_high     = input$xlim_log_high
    )
  })
  
  output$status <- renderText({
    if (is.null(sim())) "Click 'Run' to start." else "Complete."
  })
  
  output$summary_table <- renderTable({
    req(sim())
    sim_table(sim(), COR_true = input$COR_true)$table
  }, digits = 3)
  
  output$success_precision_note <- renderUI({
    req(sim())
    res <- sim_table(sim(), COR_true = input$COR_true)
    
    rmse_ia  <- res$rmse_ia_success
    rmse_fin <- res$rmse_final_success
    
    if (is.na(rmse_ia) && is.na(rmse_fin)) return(NULL)
    
    tags$div(
      style = "font-size: 0.92em; color: #444; margin-top: 12px; padding: 8px 12px; background: #f8f9fa; border-left: 4px solid #6c757d;",
      tags$p(
        tags$strong("Precision in successful arms (log-scale):")
      ),
      if (!is.na(rmse_ia)) {
        low  <- round(exp(-rmse_ia), 2)
        high <- round(exp(rmse_ia), 2)
        tags$p(
          sprintf("• IA2 success stop: RMSE_log = %.3f → typical range ≈ ×%.2f – ×%.2f", 
                  rmse_ia, low, high)
        )
      },
      if (!is.na(rmse_fin)) {
        low  <- round(exp(-rmse_fin), 2)
        high <- round(exp(rmse_fin), 2)
        tags$p(
          sprintf("• Final success stop: RMSE_log = %.3f → typical range ≈ ×%.2f – ×%.2f", 
                  rmse_fin, low, high)
        )
      },
      tags$small("(≈ 68% interval assuming roughly normal errors on the log scale)")
    )
  })
  
  output$ess_breakdown <- renderTable({
    req(sim())
    expected_n_breakdown(sim())
  }, digits = 4)
  
  output$ess_total_note <- renderPrint({
    req(sim())
    df <- expected_n_breakdown(sim())
    cat(sprintf("ESS (Average Sample Size) = %.1f\n", sum(df$Contribution)))
  })
  
  output$rpact_info <- renderPrint({
    req(sim())
    d <- sim()$rpact_design
    
    cat(sprintf(
      "Total One-sided Alpha: %.4f\nAlpha Spent Stage 1 (IA): %.4f\nAlpha Spent Stage 2 (Final): %.4f\n",
      d$alpha, d$alphaSpent[1], d$alphaSpent[2]
    ))
    
    cat("\n--- Nominal Stage-wise P-values (1-sided) ---\n")
    
    cat(sprintf(
      "Nominal p at Stage 1: %.6f\nNominal p at Stage 2: %.6f\n",
      pnorm(-d$criticalValues[1]),
      pnorm(-d$criticalValues[2])
    ))
  })
  
  output$download_plot <- downloadHandler(
    filename = function() {
      paste0("Simulation_Results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".jpg")
    },
    content = function(file) {
      req(sim())
      jpeg(file, width = 12, height = 10, units = "in", res = 300, quality = 100)
      selection_boxplot(
        sim(),
        COR_true = input$COR_true,
        COR_NI   = input$COR_NI,
        futility_frac = input$futility_frac,
        info_frac     = input$info_frac,
        show_traj_success = input$show_traj_success,
        show_traj_fail    = input$show_traj_fail,
        use_cor_scale     = input$use_cor_scale,
        xlim_log_low      = input$xlim_log_low,
        xlim_log_high     = input$xlim_log_high
      )
      dev.off()
    }
  )
}

shinyApp(ui, server)