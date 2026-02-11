# =============================================================================
#   Shiny App: Ordinal Non-Inferiority Trial Simulator with Winner's Curse
#   v1.14 — Integrated Empirical Power & Look Timings Footnotes
# =============================================================================

library(shiny)
library(shinyWidgets)
library(MASS)
library(dplyr)
library(bslib)
library(rpact)

# ─────────────────────────────────────────────────────────────────────────────
#   Helpers (Logic remains as provided)
# ─────────────────────────────────────────────────────────────────────────────

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
  qlogis(cum_control[seq_len(length(p_control)-1)])
}

pmf_from_beta <- function(theta, beta, x = 1) {
  cdf <- ilogit(theta - beta * x)
  cdf <- c(cdf, 1)
  pmf <- diff(c(0, cdf))
  pmf[pmf < 0] <- 0
  pmf / sum(pmf)
}

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

expected_n_breakdown <- function(sim) {
  p_fut   <- mean(sim$stop_fut, na.rm = TRUE)
  p_ia    <- mean(sim$stop_ia,  na.rm = TRUE)
  p_final <- mean(!(sim$stop_fut | sim$stop_ia), na.rm = TRUE)
  df <- data.frame(
    Stage       = c("Futility stop", "IA success stop", "Final analysis"),
    Probability = c(p_fut, p_ia, p_final),
    N_at_stage  = c(sim$n_at_fut, sim$n_at_ia, sim$n_total),
    check.names = FALSE
  )
  df$Contribution <- df$Probability * df$N_at_stage
  df |> mutate(Probability = round(Probability, 3), Contribution = round(Contribution, 1))
}

# ─────────────────────────────────────────────────────────────────────────────
#   Simulation Engine
# ─────────────────────────────────────────────────────────────────────────────

simulate_obf_ordinal <- function(
    COR_true, COR_NI, n_total, futility_frac, info_frac,
    zcrit1, zcrit2, futility_p, p_control,
    seed = 1234, nSims = 1000, show_progress = TRUE
) {
  msg <- validate_probs(p_control)
  if (!is.null(msg)) stop(msg)
  set.seed(seed)
  theta <- theta_from_control_pmf(p_control)
  K <- length(p_control) - 1
  beta_true <- log(COR_true); beta_NI <- log(COR_NI)
  pi_control <- p_control; pi_treat <- pmf_from_beta(theta, beta_true)
  split_n <- function(N) { nC <- floor(N/2); list(nC = nC, nT = N - nC) }
  logCOR_paths <- matrix(NA_real_, nrow = nSims, ncol = 3, dimnames = list(NULL, c("fut", "ia", "final")))
  res <- list(
    Z_fut_all = rep(NA_real_, nSims), COR_fut_all = rep(NA_real_, nSims),
    Z1_all    = rep(NA_real_, nSims), COR1_all    = rep(NA_real_, nSims),
    Z2_all    = rep(NA_real_, nSims), COR2_all    = rep(NA_real_, nSims),
    stop_fut  = logical(nSims), stop_ia = logical(nSims), stop_final = logical(nSims),
    logCOR_paths = logCOR_paths,
    n_at_fut  = round(futility_frac * n_total), n_at_ia = round(info_frac * n_total), n_total = n_total,
    nSims = nSims, z_fut = qnorm(futility_p), zcrit1 = zcrit1, zcrit2 = zcrit2
  )
  n_fut <- res$n_at_fut; n1 <- res$n_at_ia
  pb <- if(show_progress) shiny::Progress$new() else NULL
  if (!is.null(pb)) pb$set(message = "Running simulations...", value = 0)
  for (i in seq_len(nSims)) {
    s_f <- split_n(n_fut); yCf <- sample(0:K, s_f$nC, TRUE, pi_control); yTf <- sample(0:K, s_f$nT, TRUE, pi_treat)
    df_f <- data.frame(y = factor(c(yCf,yTf), ordered = TRUE, levels = 0:K), trt = factor(rep(c("C","T"), c(s_f$nC,s_f$nT)), levels = c("C","T")))
    fit_f <- fit_logCOR(df_f)
    if (!is.null(fit_f)) {
      res$Z_fut_all[i] <- (fit_f$logCOR_hat - beta_NI) / fit_f$se
      res$COR_fut_all[i] <- exp(fit_f$logCOR_hat); res$logCOR_paths[i, "fut"] <- fit_f$logCOR_hat
      if (res$Z_fut_all[i] > res$z_fut) { res$stop_fut[i] <- TRUE; if(!is.null(pb)) pb$inc(1/nSims); next }
    }
    n_add_ia <- n1 - n_fut; s_add <- split_n(n_add_ia)
    yCa <- sample(0:K, s_add$nC, TRUE, pi_control); yTa <- sample(0:K, s_add$nT, TRUE, pi_treat)
    df_add <- data.frame(y = factor(c(yCa,yTa), ordered = TRUE, levels = 0:K), trt = factor(rep(c("C","T"), c(s_add$nC,s_add$nT)), levels = c("C","T")))
    df1 <- rbind(df_f, df_add); fit1 <- fit_logCOR(df1)
    if (!is.null(fit1)) {
      res$Z1_all[i] <- (fit1$logCOR_hat - beta_NI) / fit1$se
      res$COR1_all[i] <- exp(fit1$logCOR_hat); res$logCOR_paths[i, "ia"] <- fit1$logCOR_hat
      if (res$Z1_all[i] <= res$zcrit1) { res$stop_ia[i] <- TRUE; if(!is.null(pb)) pb$inc(1/nSims); next }
    }
    n_add_final <- n_total - n1; s_final <- split_n(n_add_final)
    yCf2 <- sample(0:K, s_final$nC, TRUE, pi_control); yTf2 <- sample(0:K, s_final$nT, TRUE, pi_treat)
    df_final_add <- data.frame(y = factor(c(yCf2,yTf2), ordered = TRUE, levels = 0:K), trt = factor(rep(c("C","T"), c(s_final$nC,s_final$nT)), levels = c("C","T")))
    df2 <- rbind(df1, df_final_add); fit2 <- fit_logCOR(df2)
    if (!is.null(fit2)) {
      res$Z2_all[i] <- (fit2$logCOR_hat - beta_NI) / fit2$se
      res$COR2_all[i] <- exp(fit2$logCOR_hat); res$logCOR_paths[i, "final"] <- fit2$logCOR_hat
      if (res$Z2_all[i] <= res$zcrit2) res$stop_final[i] <- TRUE
    }
    if (!is.null(pb)) pb$inc(1/nSims)
  }
  if (!is.null(pb)) pb$close(); res
}

sim_table <- function(sim) {
  safe_summ_ext <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(c(N=0L, Min=NA, Max=NA, Mean=NA, Median=NA, `2.5%`=NA, `97.5%`=NA))
    q <- quantile(x, c(0.025,0.5,0.975), names=FALSE)
    c(N=length(x), Min=min(x), Max=max(x), Mean=mean(x), Median=q[2], `2.5%`=q[1], `97.5%`=q[3])
  }
  fut <- safe_summ_ext(sim$COR_fut_all[sim$stop_fut]); ia <- safe_summ_ext(sim$COR1_all[sim$stop_ia]); fin <- safe_summ_ext(sim$COR2_all[sim$stop_final])
  data.frame(Stage = c("Futility stop", "IA success stop", "Final success stop"), N = c(fut["N"], ia["N"], fin["N"]), Min = c(fut["Min"], ia["Min"], fin["Min"]), Mean = c(fut["Mean"], ia["Mean"], fin["Mean"]), Median = c(fut["Median"], ia["Median"], fin["Median"]), `2.5%` = c(fut["2.5%"], ia["2.5%"], fin["2.5%"]), `97.5%`= c(fut["97.5%"], ia["97.5%"], fin["97.5%"]), Max = c(fut["Max"], ia["Max"], fin["Max"]), check.names = FALSE) |> mutate(across(where(is.numeric), ~ round(.x, 3)))
}

# ─────────────────────────────────────────────────────────────────────────────
#   Plotting Function (Updated with Power Footnote)
# ─────────────────────────────────────────────────────────────────────────────

selection_boxplot <- function(sim, COR_true, COR_NI, futility_frac, info_frac, 
                              show_traj_success = FALSE, show_traj_fail = FALSE, 
                              use_cor_scale = FALSE, 
                              xlim_log_low = -3, xlim_log_high = 4, 
                              main = "Winner's curse & selection bias") {
  
  point_cex <- 0.6
  point_alpha <- 0.25
  jitter_height <- 0.30
  box_width <- 0.28
  
  if (use_cor_scale) {
    xlim_use <- c(exp(xlim_log_low), exp(xlim_log_high))
    x_transform <- exp
    x_label_expr <- "Cumulative Odds Ratio (COR)"
  } else {
    xlim_use <- c(xlim_log_low, xlim_log_high)
    x_transform <- identity
    x_label_expr <- expression(log(Cumulative~Odds~Ratio))
  }
  
  # ────────────────────────────────────────────────────────────────
  # Power and average sample size
  # ────────────────────────────────────────────────────────────────
  p_fut <- mean(sim$stop_fut, na.rm = TRUE)
  p_ia  <- mean(sim$stop_ia, na.rm = TRUE)
  p_final_suc <- mean(sim$stop_final, na.rm = TRUE)
  empirical_power <- p_ia + p_final_suc
  
  p_reach_final <- mean(!(sim$stop_fut | sim$stop_ia), na.rm = TRUE)
  avg_n <- round(p_fut * sim$n_at_fut + p_ia * sim$n_at_ia + p_reach_final * sim$n_total)
  
  # ────────────────────────────────────────────────────────────────
  # Monte Carlo 95% CI for empirical power (Wilson score interval)
  # ────────────────────────────────────────────────────────────────
  n_sims_used <- sim$nSims
  k_success <- round(empirical_power * n_sims_used)  # integer successes
  z <- qnorm(0.975)                                  # ≈1.96
  n <- n_sims_used
  center <- (k_success + z^2 / 2) / (n + z^2)
  margin <- z * sqrt( (k_success/n * (1 - k_success/n) + z^2/(4*n)) / (n + z^2) )
  ci_lower <- max(0, center - margin)
  ci_upper <- min(1, center + margin)
  
  power_text <- sprintf("Empirical Power: %.1f%% (%.1f%%–%.1f%%) (IA success + Final success) across %d simulations",
                        empirical_power * 100,
                        ci_lower * 100,
                        ci_upper * 100,
                        n_sims_used)
  
  # ────────────────────────────────────────────────────────────────
  # Define groups
  # ────────────────────────────────────────────────────────────────
  groups <- list(
    "All @ futility"       = sim$logCOR_paths[,"fut"][is.finite(sim$logCOR_paths[,"fut"])],
    "Stopped futility"     = sim$logCOR_paths[,"fut"][sim$stop_fut & is.finite(sim$logCOR_paths[,"fut"])],
    "All @ interim"        = sim$logCOR_paths[,"ia"][is.finite(sim$logCOR_paths[,"ia"])],
    "Stopped IA success"   = sim$logCOR_paths[,"ia"][sim$stop_ia & is.finite(sim$logCOR_paths[,"ia"])],
    "All @ final"          = sim$logCOR_paths[,"final"][is.finite(sim$logCOR_paths[,"final"])],
    "Stopped final success"= sim$logCOR_paths[,"final"][sim$stop_final & is.finite(sim$logCOR_paths[,"final"])]
  )
  
  counts_actual <- sapply(groups, length)
  
  if (use_cor_scale) groups <- lapply(groups, exp)
  
  keep <- sapply(groups, function(g) length(g) > 0 && all(is.finite(g)))
  groups <- groups[keep]
  group_names <- names(groups)
  counts_actual <- counts_actual[keep]
  
  if (length(groups) == 0) {
    plot.new()
    text(0.5, 0.5, "No valid data", cex = 1.4, col = "gray")
    return(invisible(NULL))
  }
  
  # ────────────────────────────────────────────────────────────────
  # Consistent jitter vector (per original simulation row)
  # ────────────────────────────────────────────────────────────────
  set.seed(202506)
  n_max <- nrow(sim$logCOR_paths)
  jitter_master <- runif(n_max, -jitter_height, jitter_height)
  
  # ────────────────────────────────────────────────────────────────
  # Plot setup
  # ────────────────────────────────────────────────────────────────
  op <- par(mar = c(14, 12, 6, 28))
  on.exit(par(op))
  
  plot(0, type = "n", 
       xlim = xlim_use, 
       ylim = c(0.4, length(groups) + 0.9), 
       xlab = "", ylab = "", yaxt = "n", las = 1,
       main = sprintf("%s\n(avg sample size = %d)", main, avg_n))
  
  # Footnotes
  mtext(x_label_expr, side = 1, line = 3, adj = 0.5, font = 2, cex = 1.1)
  
  mtext(sprintf("Futility Look: N = %d (%d%% info) | IA Success Look: N = %d (%d%% info) | Final Analysis: N = %d", 
                sim$n_at_fut, round(futility_frac*100), 
                sim$n_at_ia, round(info_frac*100), sim$n_total), 
        side = 1, line = 6, adj = 0, cex = 0.9, col = "black")
  
  mtext(sprintf("Vertical Lines: Green dashed = True Effect (%.2f) | Red dotted = NI Margin (%.2f)", 
                if(use_cor_scale) COR_true else log(COR_true), 
                if(use_cor_scale) COR_NI else log(COR_NI)), 
        side = 1, line = 7.5, adj = 0, cex = 0.9, col = "darkblue")
  
  mtext(power_text, 
        side = 1, line = 9, adj = 0, cex = 1.0, font = 2, col = "forestgreen")
  
  mtext("Note: ESS (Expected Sample Size) is the average N across all sims, accounting for early stopping. Power CI is Monte Carlo 95% Wilson interval.",
        side = 1, line = 10.5, adj = 0, cex = 0.85, col = "gray30")
  
  # ────────────────────────────────────────────────────────────────
  # Trajectories: connect through jittered positions (as in your last request)
  # ────────────────────────────────────────────────────────────────
  if (show_traj_success || show_traj_fail) {
    
    group_y_map <- setNames(seq_along(groups), group_names)
    
    for (i in seq_len(nrow(sim$logCOR_paths))) {
      is_success <- sim$stop_ia[i] || sim$stop_final[i]
      
      if (is_success && !show_traj_success) next
      if (!is_success && !show_traj_fail)    next
      
      has_fut   <- is.finite(sim$logCOR_paths[i, "fut"])
      has_ia    <- is.finite(sim$logCOR_paths[i, "ia"])
      has_final <- is.finite(sim$logCOR_paths[i, "final"])
      
      if (sum(c(has_fut, has_ia, has_final)) < 2) next
      
      path_x <- numeric(0)
      path_y <- numeric(0)
      
      if (has_fut) {
        group_name <- ifelse(sim$stop_fut[i], "Stopped futility", "All @ futility")
        if (group_name %in% names(groups)) {
          idx_match <- which(sim$logCOR_paths[,"fut"] == sim$logCOR_paths[i, "fut"] & is.finite(sim$logCOR_paths[,"fut"]))
          if (length(idx_match) > 0) {
            jitter_this <- jitter_master[idx_match[1]]
            path_x <- c(path_x, sim$logCOR_paths[i, "fut"])
            path_y <- c(path_y, group_y_map[group_name] + jitter_this)
          }
        }
      }
      
      if (has_ia) {
        group_name <- ifelse(sim$stop_ia[i], "Stopped IA success", "All @ interim")
        if (group_name %in% names(groups)) {
          idx_match <- which(sim$logCOR_paths[,"ia"] == sim$logCOR_paths[i, "ia"] & is.finite(sim$logCOR_paths[,"ia"]))
          if (length(idx_match) > 0) {
            jitter_this <- jitter_master[idx_match[1]]
            path_x <- c(path_x, sim$logCOR_paths[i, "ia"])
            path_y <- c(path_y, group_y_map[group_name] + jitter_this)
          }
        }
      }
      
      if (has_final) {
        group_name <- ifelse(sim$stop_final[i], "Stopped final success", "All @ final")
        if (group_name %in% names(groups)) {
          idx_match <- which(sim$logCOR_paths[,"final"] == sim$logCOR_paths[i, "final"] & is.finite(sim$logCOR_paths[,"final"]))
          if (length(idx_match) > 0) {
            jitter_this <- jitter_master[idx_match[1]]
            path_x <- c(path_x, sim$logCOR_paths[i, "final"])
            path_y <- c(path_y, group_y_map[group_name] + jitter_this)
          }
        }
      }
      
      if (length(path_x) < 2) next
      
      if (use_cor_scale) path_x <- exp(path_x)
      
      lines(path_x, path_y,
            col = if(is_success) rgb(0.1, 0.65, 0.1, 0.18) else rgb(0.9, 0.2, 0.2, 0.18),
            lwd = 0.8, lty = 2)
    }
  }
  
  # ────────────────────────────────────────────────────────────────
  # Axes & reference lines
  # ────────────────────────────────────────────────────────────────
  axis(2, at = seq_along(groups), labels = group_names, las = 1, cex.axis = 0.95)
  abline(v = x_transform(log(COR_true)),  lty = 2, col = "darkgreen", lwd = 2.5)
  abline(v = x_transform(log(COR_NI)),    lty = 3, col = "red",       lwd = 2.5)
  
  abline(h = seq_along(groups), col = "gray92", lwd = 0.8)
  
  # ────────────────────────────────────────────────────────────────
  # Points & boxes
  # ────────────────────────────────────────────────────────────────
  for (i in seq_along(groups)) {
    name <- group_names[i]
    
    if (grepl("All @ futility", name)) {
      idx <- which(is.finite(sim$logCOR_paths[,"fut"]))
    } else if (grepl("Stopped futility", name)) {
      idx <- which(sim$stop_fut & is.finite(sim$logCOR_paths[,"fut"]))
    } else if (grepl("All @ interim", name)) {
      idx <- which(is.finite(sim$logCOR_paths[,"ia"]))
    } else if (grepl("Stopped IA success", name)) {
      idx <- which(sim$stop_ia & is.finite(sim$logCOR_paths[,"ia"]))
    } else if (grepl("All @ final", name)) {
      idx <- which(is.finite(sim$logCOR_paths[,"final"]))
    } else if (grepl("Stopped final success", name)) {
      idx <- which(sim$stop_final & is.finite(sim$logCOR_paths[,"final"]))
    } else next
    
    vals <- groups[[i]]
    n <- length(vals)
    if (n == 0) next
    
    jitter_y <- jitter_master[idx]
    
    if (length(jitter_y) != n) {
      warning(sprintf("Jitter mismatch in %s", name))
      next
    }
    
    col_point <- if (grepl("All @ futility", name))       rgb(0.1,0.4,0.9,point_alpha) else
      if (grepl("futility", name))             rgb(0.9,0.15,0.15,point_alpha) else
        if (grepl("success", name))              rgb(0.1,0.65,0.1,point_alpha) else
          rgb(0.3,0.3,0.3,point_alpha)
    
    points(vals, i + jitter_y, pch = 19, cex = point_cex, col = col_point)
    
    if (n >= 3) {
      q <- quantile(vals, c(0.25, 0.5, 0.75))
      iqr <- q[3] - q[1]
      lower_whisk <- min(vals[vals >= q[1] - 1.5*iqr])
      upper_whisk <- max(vals[vals <= q[3] + 1.5*iqr])
      
      rect(q[1], i - box_width/2, q[3], i + box_width/2, 
           col = rgb(0.88,0.93,1,0.5), border = "steelblue", lwd = 1.4)
      segments(q[2], i - box_width/2, q[2], i + box_width/2, lwd = 5, col = "royalblue3")
      segments(lower_whisk, i, q[1], i, lwd = 2.4, col = "midnightblue")
      segments(q[3], i, upper_whisk, i, lwd = 2.4, col = "midnightblue")
    }
  }
  
  # ────────────────────────────────────────────────────────────────
  # Side table
  # ────────────────────────────────────────────────────────────────
  usr <- par("usr")
  x_range <- usr[2] - usr[1]
  x_text_actual <- usr[2] + 0.05 * x_range
  x_text_n      <- usr[2] + 0.16 * x_range
  x_text_pct    <- usr[2] + 0.27 * x_range
  
  props_all <- c(1, p_fut, 1 - p_fut, p_ia, 1 - p_fut - p_ia, p_final_suc)
  n_all     <- c(sim$n_at_fut, sim$n_at_fut, sim$n_at_ia, sim$n_at_ia, sim$n_total, sim$n_total)
  
  props <- props_all[keep]
  n_col <- n_all[keep]
  
  for (i in seq_along(groups)) {
    col_text <- if (grepl("All @ futility", group_names[i])) "black" else
      if (grepl("futility", group_names[i]))   "firebrick" else
        if (grepl("success", group_names[i]))    "forestgreen" else
          "gray30"
    
    text(x_text_actual, i, format(counts_actual[i], big.mark=","), adj = 0, cex = 1.05, col = col_text, xpd = TRUE)
    text(x_text_n,      i, sprintf("%d", n_col[i]), adj = 0, cex = 1.05, col = "gray30", xpd = TRUE)
    text(x_text_pct,    i, sprintf("%.1f%%", 100 * props[i]), adj = 0, cex = 1.05, font = 2, col = col_text, xpd = TRUE)
  }
  
  text(x_text_actual, length(groups) + 0.75, "No. of trials", adj = 0, cex = 1.0, font = 2, xpd = TRUE)
  text(x_text_n,      length(groups) + 0.75, "N per trial",  adj = 0, cex = 1.0, font = 2, xpd = TRUE)
  text(x_text_pct,    length(groups) + 0.75, "% of sims",    adj = 0, cex = 1.0, font = 2, xpd = TRUE)
  
  invisible(NULL)
}

# ─────────────────────────────────────────────────────────────────────────────
#   UI & Server (Logic remains as provided)
# ─────────────────────────────────────────────────────────────────────────────

# ui <- page_sidebar(
#   title = "Ordinal NI Trial Simulator + Winner's Curse v1.14",
#   sidebar = sidebar(
#     h4("Simulation Settings"),
#     numericInput("n_total", "Total sample size", value = 600),
#     numericInput("n_sims", "Number of simulations", value = 1000),
#     numericInput("seed", "Random seed", value = 202506),
#     numericInput("COR_true", "True COR", value = 1.3, step = 0.05),
#     numericInput("COR_NI", "NI margin (M)", value = 1.6, step = 0.1),
#     sliderInput("futility_frac", "Futility look fraction", min = 0.2, max = 0.7, value = 0.5),
#     sliderInput("info_frac", "Interim look fraction", min = 0.5, max = 0.95, value = 0.80),
#     numericInput("futility_p", "Futility p-threshold", value = 0.70),
#     textInput("p_control_txt", "Control probabilities", value = "0.04, 0.02, 0.45, 0.34, 0.15"),
#     checkboxInput("show_traj_success", "Show successful trajectories", value = FALSE),
#     checkboxInput("show_traj_fail", "Show failed trajectories", value = FALSE),
#     hr(), checkboxInput("use_cor_scale", "Display on COR scale", value = FALSE),
#     sliderInput("xlim_log_low", "X lower (log)", min = -6, max = 0, value = -1, step = 0.5),
#     sliderInput("xlim_log_high", "X upper (log)", min = 0, max = 7, value = 2, step = 0.5),
#     actionButton("run_btn", "Run Simulation", class = "btn-primary", icon = icon("play"))
#   ),
#   card(tabsetPanel(
#     tabPanel("Summary Table", verbatimTextOutput("status"), tableOutput("summary_table"), hr(), h5("rpact Design & Nominal P-values"), verbatimTextOutput("rpact_info")),
#     tabPanel("Expected Sample Size", tableOutput("ess_breakdown"), hr(), verbatimTextOutput("ess_total_note")),
#     tabPanel("Winner's Curse Plot", plotOutput("boxplot", height = "750px"))
#   ))
# )


# ui <- page_sidebar(
#   title = "Ordinal NI Trial Simulator + Winner's Curse v1.14",
#   
#   sidebar = sidebar(
#     # ────────────────────────────────────────────────────────────────
#     # Run button at the very top
#     # ────────────────────────────────────────────────────────────────
#     actionButton("run_btn", "Run Simulation", 
#                  class = "btn-primary btn-lg", 
#                  icon = icon("play"), 
#                  width = "100%"),
#     
#     hr(),
#     
#     h4("Simulation Settings"),
#     numericInput("n_total", "Total sample size", value = 600),
#     numericInput("n_sims", "Number of simulations", value = 1000),
#     numericInput("seed", "Random seed", value = 202506),
#     numericInput("COR_true", "True COR", value = 1.3, step = 0.05),
#     numericInput("COR_NI", "NI margin (M)", value = 1.6, step = 0.1),
#     sliderInput("futility_frac", "Futility look fraction", min = 0.2, max = 0.7, value = 0.5),
#     sliderInput("info_frac", "Interim look fraction", min = 0.5, max = 0.95, value = 0.80),
#     numericInput("futility_p", "Futility p-threshold", value = 0.70),
#     textInput("p_control_txt", "Control probabilities", value = "0.04, 0.02, 0.45, 0.34, 0.15"),
#     
#     hr(),
#     
#     # Other display options in the middle
#     checkboxInput("use_cor_scale", "Display on COR scale", value = FALSE),
#     sliderInput("xlim_log_low", "X lower (log)", min = -6, max = 0, value = -1, step = 0.5),
#     sliderInput("xlim_log_high", "X upper (log)", min = 0, max = 7, value = 2, step = 0.5),
#     
#     # ────────────────────────────────────────────────────────────────
#     # Trajectory checkboxes moved to the bottom
#     # ────────────────────────────────────────────────────────────────
#     hr(),
#     h5("Show trajectories"),
#     checkboxInput("show_traj_success", "Show successful trajectories", value = FALSE),
#     checkboxInput("show_traj_fail", "Show failed trajectories", value = FALSE)
#   ),
#   
#   card(tabsetPanel(
#     tabPanel("Summary Table", 
#              verbatimTextOutput("status"), 
#              tableOutput("summary_table"), 
#              hr(), 
#              h5("rpact Design & Nominal P-values"), 
#              verbatimTextOutput("rpact_info")),
#     
#     tabPanel("Expected Sample Size", 
#              tableOutput("ess_breakdown"), 
#              hr(), 
#              verbatimTextOutput("ess_total_note")),
#     
#     tabPanel("Winner's Curse Plot", 
#              plotOutput("boxplot", height = "750px"))
#   ))
# )
# ui <- page_sidebar(
#   
#   title = "Ordinal NI Trial Simulator + Winner's Curse v1.14",
#   
#   sidebar = sidebar(
#     width = 350,   # optional: make sidebar a bit wider for comfort
#     
#     # ────────────────────────────────────────────────────────────────
#     # Primary action — prominent at the top
#     # ────────────────────────────────────────────────────────────────
#     actionButton(
#       "run_btn", 
#       "Run Simulation", 
#       class = "btn-primary btn-lg", 
#       icon = icon("play-circle"), 
#       width = "100%"
#     ),
#     
#     hr(style = "margin: 1.2em 0;"),
#     
#     # ────────────────────────────────────────────────────────────────
#     # Core simulation parameters
#     # ────────────────────────────────────────────────────────────────
#     tags$div(
#       style = "padding: 0 8px;",
#       h5("Trial & Simulation Settings", style = "margin-bottom: 1em; color: #444;"),
#       
#       numericInput("n_total",     "Total sample size",          value = 600,      width = "100%"),
#       numericInput("n_sims",      "Number of simulations",      value = 1000,     width = "100%"),
#       numericInput("seed",        "Random seed",                value = 202506,   width = "100%"),
#       
#       hr(style = "margin: 1em 0; border-top: 1px dashed #ccc;"),
#       
#       numericInput("COR_true",    "True COR",                   value = 1.3,      step = 0.05, width = "100%"),
#       numericInput("COR_NI",      "NI margin (M)",              value = 1.6,      step = 0.1,  width = "100%"),
#       
#       sliderInput("futility_frac", "Futility look fraction",   min = 0.2, max = 0.7, value = 0.5, width = "100%"),
#       sliderInput("info_frac",     "Interim look fraction",     min = 0.5, max = 0.95, value = 0.80, width = "100%"),
#       numericInput("futility_p",   "Futility p-threshold",      value = 0.70,     width = "100%"),
#       
#       textInput("p_control_txt",   "Control probabilities",     value = "0.04, 0.02, 0.45, 0.34, 0.15", width = "100%")
#     ),
#     
#     hr(style = "margin: 1.5em 0;"),
#     
#     # ────────────────────────────────────────────────────────────────
#     # Display & visualization options — lower priority
#     # ────────────────────────────────────────────────────────────────
#     tags$div(
#       style = "padding: 0 8px;",
#       h5("Plot Options", style = "margin-bottom: 1em; color: #444;"),
#       
#       checkboxInput("use_cor_scale", "Display on COR scale (instead of log)", value = FALSE),
#       
#       sliderInput("xlim_log_low",  "X-axis lower limit (log scale)",  min = -6, max = 0, value = -1, step = 0.5),
#       sliderInput("xlim_log_high", "X-axis upper limit (log scale)",  min = 0,  max = 7, value = 2,  step = 0.5),
#       
#       hr(style = "margin: 1.2em 0; border-top: 1px dashed #ccc;"),
#       
#       tags$strong("Show trajectories"),
#       checkboxInput("show_traj_success", "Successful trajectories", value = FALSE),
#       checkboxInput("show_traj_fail",    "Failed / futility trajectories", value = FALSE)
#     )
#   ),
#   
#   # ────────────────────────────────────────────────────────────────
#   # Main content remains unchanged
#   # ────────────────────────────────────────────────────────────────
#   card(
#     tabsetPanel(
#       tabPanel("Summary Table", 
#                verbatimTextOutput("status"), 
#                tableOutput("summary_table"), 
#                hr(), 
#                h5("rpact Design & Nominal P-values"), 
#                verbatimTextOutput("rpact_info")),
#       
#       tabPanel("Expected Sample Size", 
#                tableOutput("ess_breakdown"), 
#                hr(), 
#                verbatimTextOutput("ess_total_note")),
#       
#       tabPanel("Winner's Curse Plot", 
#                plotOutput("boxplot", height = "750px"))
#     )
#   )
# )

ui <- page_sidebar(
  
  title = "Ordinal NI Trial Simulator + Winner's Curse v1.14",
  
  sidebar = sidebar(
    width = 350,
    
    # Primary action — prominent at the top
    actionButton(
      "run_btn", 
      "Run Simulation", 
      class = "btn-primary btn-lg", 
      icon = icon("play-circle"), 
      width = "100%"
    ),
    
    hr(style = "margin: 1.2em 0;"),
    
    # Core simulation parameters
    tags$div(
      style = "padding: 0 8px;",
      h5("Trial & Simulation Settings", style = "margin-bottom: 1em; color: #444;"),
      
      numericInput("n_total",     "Total sample size",          value = 600,      width = "100%"),
      numericInput("n_sims",      "Number of simulations",      value = 1000,     width = "100%"),
      numericInput("seed",        "Random seed",                value = 202506,   width = "100%"),
      
      hr(style = "margin: 1em 0; border-top: 1px dashed #ccc;"),
      
      numericInput("COR_true",    "True COR",                   value = 1.3,      step = 0.05, width = "100%"),
      numericInput("COR_NI",      "NI margin (M)",              value = 1.6,      step = 0.1,  width = "100%"),
      
      sliderInput("futility_frac", "Futility look fraction",   min = 0.2, max = 0.7, value = 0.5, width = "100%"),
      sliderInput("info_frac",     "Interim look fraction",     min = 0.5, max = 0.95, value = 0.80, width = "100%"),
      numericInput("futility_p",   "Futility p-threshold",      value = 0.70,     width = "100%"),
      
      textInput("p_control_txt",   "Control probabilities",     value = "0.04, 0.02, 0.45, 0.34, 0.15", width = "100%")
    ),
    
    hr(style = "margin: 1.5em 0;"),
    
    # Display & visualization options
    tags$div(
      style = "padding: 0 8px;",
      h5("Plot Options", style = "margin-bottom: 1em; color: #444;"),
      
      checkboxInput("use_cor_scale", "Display on COR scale (instead of log)", value = FALSE),
      
      sliderInput("xlim_log_low",  "X-axis lower limit (log scale)",  min = -6, max = 0, value = -1, step = 0.5),
      sliderInput("xlim_log_high", "X-axis upper limit (log scale)",  min = 0,  max = 7, value = 2,  step = 0.5),
      
      hr(style = "margin: 1.2em 0; border-top: 1px dashed #ccc;"),
      
      tags$strong("Show trajectories"),
      checkboxInput("show_traj_success", "Successful trajectories", value = FALSE),
      checkboxInput("show_traj_fail",    "Failed / futility trajectories", value = FALSE)
    )
  ),
  
  # ────────────────────────────────────────────────────────────────
  # Main content — added new tab at the end
  # ────────────────────────────────────────────────────────────────
  card(
    tabsetPanel(
      tabPanel("Summary Table", 
               verbatimTextOutput("status"), 
               tableOutput("summary_table"), 
               hr(), 
               h5("rpact Design & Nominal P-values"), 
               verbatimTextOutput("rpact_info")),
      
      tabPanel("Expected Sample Size", 
               tableOutput("ess_breakdown"), 
               hr(), 
               verbatimTextOutput("ess_total_note")),
      
      tabPanel("Winner's Curse Plot", 
               plotOutput("boxplot", height = "750px")),
      
      # New tab — last position
      tabPanel("About & References",
               h4("What this app does", style = "margin-top: 1.2em;"),
               p("This Shiny application simulates a **group-sequential non-inferiority trial** with an ordinal primary endpoint using the cumulative odds ratio (COR) as the effect measure."),
               p("It implements a simple **two-look design** with one interim analysis (for early success) and a futility look (for early stopping for futility)."),
               p("The simulation illustrates **selection bias / winner's curse** by comparing the distribution of estimated CORs at different stages across many simulated trials:"),
               tags$ul(
                 tags$li("All trials at futility look vs. those stopped for futility"),
                 tags$li("All trials at interim vs. those stopped for early success"),
                 tags$li("All trials at final analysis vs. those that succeeded at the final look")
               ),
               p("The plot helps visualize how early stopping rules inflate effect estimates in the subsets that cross the boundary."),
               
               h4("What is the winner's curse?", style = "margin-top: 2em;"),
               p("The **winner's curse** (also called selection bias or regression to the mean after selection) occurs when you only look at trials / estimates that passed a statistical threshold (e.g. crossed efficacy boundary, avoided futility)."),
               p("Because of random variation, trials that just barely meet the stopping criterion tend to overestimate the true treatment effect — especially when the true effect is close to the boundary."),
               p("In adaptive / group-sequential designs this phenomenon is particularly relevant:"),
               tags$ul(
                 tags$li("Early efficacy stopping → inflated effect estimates in the 'success' subgroup"),
                 tags$li("Futility stopping → the remaining trials at later looks tend to look better than they really are"),
                 tags$li("The plot shows this visually by comparing 'all trials' vs 'selected trials' at each look.")
               ),
               
               h4("References & further reading", style = "margin-top: 2em;"),
               tags$ul(
                 tags$li("Proschan MA, Lan KKG, Wittes JT (2006). ", em("Statistical Monitoring of Clinical Trials: A Unified Approach.")),
                 tags$li("Bauer P, Bretz F, Dragalin V, König F, Wassmer G (2016). Twenty-five years of confirmatory adaptive designs: opportunities and pitfalls. ", em("Statistics in Medicine"), " 35(3):325–347."),
                 tags$li("Whitehead J (1986). The evaluation of a sequential trial when a positive treatment effect is anticipated. ", em("Applied Statistics"), " 35:155–163."),
                 tags$li("rpact package documentation: ", a("https://www.rpact.org", href = "https://www.rpact.org", target = "_blank")),
                 tags$li("The statistical properties of RCTs and a proposal for shrinkage: ", a("https://doi.org/10.1002/sim.9173", href = "https://doi.org/10.1002/sim.9173", target = "_blank")),
                 tags$li("The performance of that shrinkage estimator to address the winner's curse: ", a("https://doi.org/10.1002/sim.9992", href = "https://doi.org/10.1002/sim.9992", target = "_blank"))
                 
               ),
               
               p(style = "margin-top: 2em; font-style: italic; color: #555;",
                 "App version 1.14 – February 2026 • Built with Shiny, rpact & polr")
      )
    )
  )
)

server <- function(input, output, session) {
  sim_result <- eventReactive(input$run_btn, {
    req(input$run_btn)
    p_control <- parse_probs(input$p_control_txt)
    design <- getDesignGroupSequential(sided = 1, alpha = 0.025, informationRates = c(input$info_frac, 1), typeOfDesign = "asOF")
    zcrit1 <- -design$criticalValues[1]; zcrit2 <- -design$criticalValues[2]
    sim <- simulate_obf_ordinal(input$COR_true, input$COR_NI, input$n_total, input$futility_frac, input$info_frac, zcrit1, zcrit2, input$futility_p, p_control, input$seed, input$n_sims)
    sim$rpact_design <- design; sim
  })
  output$status <- renderText({ if (is.null(sim_result())) "Click 'Run' to start." else "Complete." })
  output$summary_table <- renderTable({ req(sim_result()); sim_table(sim_result()) }, digits = 3)
  output$ess_breakdown <- renderTable({ req(sim_result()); expected_n_breakdown(sim_result()) })
  output$ess_total_note <- renderPrint({ req(sim_result()); df <- expected_n_breakdown(sim_result()); cat(sprintf("ESS (Average Sample Size) = %.1f\n", sum(df$Contribution))) })
  output$rpact_info <- renderPrint({ 
    req(sim_result()); d <- sim_result()$rpact_design
    cat(sprintf("Total One-sided Alpha: %.4f\n", d$alpha))
    cat(sprintf("Alpha Spent Stage 1 (IA): %.4f\n", d$alphaSpent[1]))
    cat(sprintf("Alpha Spent Stage 2 (Final): %.4f\n", d$alphaSpent[2]))
    cat("\n--- Nominal Stage-wise P-values (1-sided) ---\n")
    cat(sprintf("Nominal p at Stage 1: %.6f\n", pnorm(-d$criticalValues[1])))
    cat(sprintf("Nominal p at Stage 2: %.6f\n", pnorm(-d$criticalValues[2])))
    cat("\n--- Critical Boundaries (Lower-tail Z) ---\n")
    cat(sprintf("Z-Boundary IA:    %.4f\n", -d$criticalValues[1]))
    cat(sprintf("Z-Boundary Final: %.4f\n", -d$criticalValues[2]))
  })
  output$boxplot <- renderPlot({ req(sim_result()); selection_boxplot(sim_result(), input$COR_true, input$COR_NI, input$futility_frac, input$info_frac, input$show_traj_success, input$show_traj_fail, input$use_cor_scale, input$xlim_log_low, input$xlim_log_high) })
}

shinyApp(ui, server)