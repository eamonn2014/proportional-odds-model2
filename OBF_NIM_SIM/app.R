# =============================================================================
#   Shiny App: Ordinal Non-Inferiority Trial Simulator with Winner's Curse
#   v1.8 — efficacy boundaries from rpact (one-sided alpha = 0.025) using IA slider
#          removed manual alpha inputs; prints rpact info under summary table
#          + TAB: Expected Sample Size Breakdown (stage probs × N contributions)
#          + Plot: add "N per sim" column (planned N at that stage) left of "% of sims"
# =============================================================================

library(shiny)
library(shinyWidgets)
library(MASS)
library(dplyr)
library(bslib)
library(rpact)

# ─────────────────────────────────────────────────────────────────────────────
#   Helpers
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

# Expected sample size breakdown table
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
  
  df |>
    mutate(
      Probability  = round(Probability, 3),
      Contribution = round(Contribution, 1)
    )
}

# ─────────────────────────────────────────────────────────────────────────────
#   Simulation (efficacy boundaries provided as Z, derived from rpact in server)
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
  beta_true <- log(COR_true)
  beta_NI   <- log(COR_NI)
  
  pi_control <- p_control
  pi_treat   <- pmf_from_beta(theta, beta_true)
  
  split_n <- function(N) {
    nC <- floor(N/2); list(nC = nC, nT = N - nC)
  }
  
  logCOR_paths <- matrix(NA_real_, nrow = nSims, ncol = 3,
                         dimnames = list(NULL, c("fut", "ia", "final")))
  
  res <- list(
    Z_fut_all = rep(NA_real_, nSims), COR_fut_all = rep(NA_real_, nSims),
    Z1_all    = rep(NA_real_, nSims), COR1_all    = rep(NA_real_, nSims),
    Z2_all    = rep(NA_real_, nSims), COR2_all    = rep(NA_real_, nSims),
    stop_fut  = logical(nSims), stop_ia = logical(nSims), stop_final = logical(nSims),
    logCOR_paths = logCOR_paths,
    n_at_fut  = round(futility_frac * n_total),
    n_at_ia   = round(info_frac  * n_total),
    n_total   = n_total,
    nSims = nSims,
    z_fut  = qnorm(futility_p),
    zcrit1 = zcrit1,
    zcrit2 = zcrit2
  )
  
  n_fut <- res$n_at_fut
  n1    <- res$n_at_ia
  
  pb <- if(show_progress) shiny::Progress$new() else NULL
  if (!is.null(pb)) pb$set(message = "Running simulations...", value = 0)
  
  for (i in seq_len(nSims)) {
    
    # ---- Futility look ----
    s_f  <- split_n(n_fut)
    yCf  <- sample(0:K, s_f$nC, TRUE, pi_control)
    yTf  <- sample(0:K, s_f$nT, TRUE, pi_treat)
    df_f <- data.frame(
      y   = factor(c(yCf,yTf), ordered = TRUE, levels = 0:K),
      trt = factor(rep(c("C","T"), c(s_f$nC,s_f$nT)), levels = c("C","T"))
    )
    fit_f <- fit_logCOR(df_f)
    if (!is.null(fit_f)) {
      res$Z_fut_all[i]   <- (fit_f$logCOR_hat - beta_NI) / fit_f$se
      res$COR_fut_all[i] <- exp(fit_f$logCOR_hat)
      res$logCOR_paths[i, "fut"] <- fit_f$logCOR_hat
      
      # futility stop rule: stop if Z > z_fut
      if (res$Z_fut_all[i] > res$z_fut) {
        res$stop_fut[i] <- TRUE
        if(!is.null(pb)) pb$inc(1/nSims)
        next
      }
    }
    
    # ---- Interim efficacy look ----
    n_add_ia <- n1 - n_fut
    s_add    <- split_n(n_add_ia)
    yCa      <- sample(0:K, s_add$nC, TRUE, pi_control)
    yTa      <- sample(0:K, s_add$nT, TRUE, pi_treat)
    df_add   <- data.frame(
      y   = factor(c(yCa,yTa), ordered = TRUE, levels = 0:K),
      trt = factor(rep(c("C","T"), c(s_add$nC,s_add$nT)), levels = c("C","T"))
    )
    df1  <- rbind(df_f, df_add)
    fit1 <- fit_logCOR(df1)
    if (!is.null(fit1)) {
      res$Z1_all[i]   <- (fit1$logCOR_hat - beta_NI) / fit1$se
      res$COR1_all[i] <- exp(fit1$logCOR_hat)
      res$logCOR_paths[i, "ia"] <- fit1$logCOR_hat
      
      # success at IA if Z <= zcrit1 (lower-tail bound)
      if (res$Z1_all[i] <= res$zcrit1) {
        res$stop_ia[i] <- TRUE
        if(!is.null(pb)) pb$inc(1/nSims)
        next
      }
    }
    
    # ---- Final look ----
    n_add_final <- n_total - n1
    s_final     <- split_n(n_add_final)
    yCf2        <- sample(0:K, s_final$nC, TRUE, pi_control)
    yTf2        <- sample(0:K, s_final$nT, TRUE, pi_treat)
    df_final_add <- data.frame(
      y   = factor(c(yCf2,yTf2), ordered = TRUE, levels = 0:K),
      trt = factor(rep(c("C","T"), c(s_final$nC,s_final$nT)), levels = c("C","T"))
    )
    df2  <- rbind(df1, df_final_add)
    fit2 <- fit_logCOR(df2)
    if (!is.null(fit2)) {
      res$Z2_all[i]   <- (fit2$logCOR_hat - beta_NI) / fit2$se
      res$COR2_all[i] <- exp(fit2$logCOR_hat)
      res$logCOR_paths[i, "final"] <- fit2$logCOR_hat
      
      # success at final if Z <= zcrit2
      if (res$Z2_all[i] <= res$zcrit2) res$stop_final[i] <- TRUE
    }
    
    if (!is.null(pb)) pb$inc(1/nSims)
  }
  
  if (!is.null(pb)) pb$close()
  res
}

sim_table <- function(sim) {
  safe_summ_ext <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(c(N=0L, Min=NA, Max=NA, Mean=NA, Median=NA, `2.5%`=NA, `97.5%`=NA))
    q <- quantile(x, c(0.025,0.5,0.975), names=FALSE)
    c(N=length(x), Min=min(x), Max=max(x), Mean=mean(x), Median=q[2], `2.5%`=q[1], `97.5%`=q[3])
  }
  
  fut <- safe_summ_ext(sim$COR_fut_all[sim$stop_fut])
  ia  <- safe_summ_ext(sim$COR1_all[sim$stop_ia])
  fin <- safe_summ_ext(sim$COR2_all[sim$stop_final])
  
  data.frame(
    Stage  = c("Futility stop", "IA success stop", "Final success stop"),
    N      = c(fut["N"], ia["N"], fin["N"]),
    Min    = c(fut["Min"], ia["Min"], fin["Min"]),
    Mean   = c(fut["Mean"], ia["Mean"], fin["Mean"]),
    Median = c(fut["Median"], ia["Median"], fin["Median"]),
    `2.5%` = c(fut["2.5%"], ia["2.5%"], fin["2.5%"]),
    `97.5%`= c(fut["97.5%"], ia["97.5%"], fin["97.5%"]),
    Max    = c(fut["Max"], ia["Max"], fin["Max"]),
    check.names = FALSE
  ) |>
    mutate(across(where(is.numeric), ~ round(.x, 3)))
}

selection_boxplot <- function(
    sim,
    COR_true,
    COR_NI,
    futility_frac,
    info_frac,
    show_traj_success = FALSE,
    show_traj_fail    = FALSE,
    use_cor_scale     = FALSE,
    xlim_log_low      = -3,
    xlim_log_high     = 4,
    main = "Winner's curse & selection bias"
) {
  point_cex     <- 0.6
  point_alpha   <- 0.25
  jitter_height <- 0.30
  box_width     <- 0.28
  
  # ── Scale setup ───────────────────────────────────────────────────────────
  if (use_cor_scale) {
    xlim_use     <- c(exp(xlim_log_low), exp(xlim_log_high))
    x_transform  <- exp
    x_label_expr <- "Cumulative Odds Ratio (COR)"
    ref_text     <- c(sprintf("True COR = %.2f", COR_true),
                      sprintf("NI margin = %.2f", COR_NI))
    scale_note   <- "Scale: Odds Ratio (COR)"
  } else {
    xlim_use     <- c(xlim_log_low, xlim_log_high)
    x_transform  <- identity
    x_label_expr <- expression(log(Cumulative~Odds~Ratio))
    ref_text     <- c(sprintf("True log(COR) = %.2f", log(COR_true)),
                      sprintf("log(NI margin) = %.2f", log(COR_NI)))
    scale_note   <- "Scale: log(Odds Ratio)"
  }
  
  # ── Avg sample size (expected N) ──────────────────────────────────────────
  p_fut   <- mean(sim$stop_fut, na.rm = TRUE)
  p_ia    <- mean(sim$stop_ia,  na.rm = TRUE)
  p_reach_final <- mean(!(sim$stop_fut | sim$stop_ia), na.rm = TRUE)
  
  avg_n <- round(
    p_fut * sim$n_at_fut +
      p_ia * sim$n_at_ia +
      p_reach_final * sim$n_total
  )
  
  # ── Groups ────────────────────────────────────────────────────────────────
  groups <- list(
    "All @ futility"            = sim$logCOR_paths[,"fut"][is.finite(sim$logCOR_paths[,"fut"])],
    "Stopped futility"          = sim$logCOR_paths[,"fut"][sim$stop_fut & is.finite(sim$logCOR_paths[,"fut"])],
    "All @ interim"             = sim$logCOR_paths[,"ia"][is.finite(sim$logCOR_paths[,"ia"])],
    "Stopped IA success"        = sim$logCOR_paths[,"ia"][sim$stop_ia & is.finite(sim$logCOR_paths[,"ia"])],
    "All @ final"               = sim$logCOR_paths[,"final"][is.finite(sim$logCOR_paths[,"final"])],
    "Stopped final success"     = sim$logCOR_paths[,"final"][sim$stop_final & is.finite(sim$logCOR_paths[,"final"])]
  )
  
  if (use_cor_scale) groups <- lapply(groups, exp)
  
  keep <- sapply(groups, function(g) length(g) > 0 && all(is.finite(g)))
  groups <- groups[keep]
  group_names <- names(groups)
  
  if (length(groups) == 0) {
    plot.new()
    text(0.5, 0.5, "No valid data", cex = 1.4, col = "gray")
    return(invisible(NULL))
  }
  
  # Increased right margin for two text columns (N + %)
  op <- par(mar = c(13, 12, 6, 22))
  on.exit(par(op))
  
  plot(0, type = "n",
       xlim = xlim_use,
       ylim = c(0.4, length(groups) + 0.9),
       xlab = "", ylab = "",
       yaxt = "n", las = 1,
       main = sprintf("%s\n(avg sample size = %d)", main, avg_n))
  
  # X-label moved left
  mtext(x_label_expr, side = 1, line = 5.5, adj = 0, cex = 1.1, xpd = TRUE)
  
  axis(2, at = seq_along(groups), labels = group_names, las = 1, cex.axis = 0.95)
  
  # Reference lines
  abline(v = x_transform(log(COR_true)), lty = 2, col = "darkgreen", lwd = 2.5)
  abline(v = x_transform(log(COR_NI)),   lty = 3, col = "red",       lwd = 2.5)
  
  text(x_transform(log(COR_true)), length(groups) + 0.45, ref_text[1],
       col = "darkgreen", cex = 0.9, pos = 4, xpd = TRUE)
  text(x_transform(log(COR_NI)),   length(groups) + 0.25, ref_text[2],
       col = "red", cex = 0.9, pos = 4, xpd = TRUE)
  
  # ── Trajectories — optional via checkboxes ────────────────────────────────
  if (show_traj_success || show_traj_fail) {
    
    idx_have_fut <- which(is.finite(sim$logCOR_paths[,"fut"]))
    
    if (length(idx_have_fut) > 0) {
      for (i in idx_have_fut) {
        
        is_success     <- sim$stop_final[i]
        is_interim_suc <- sim$stop_ia[i]
        
        # Treat interim success as "success" category (green only)
        if (is_success || is_interim_suc) {
          if (!show_traj_success) next
          col_rgb <- rgb(0.1, 0.65, 0.1, 0.18)  # green
        } else {
          if (!show_traj_fail) next
          col_rgb <- rgb(0.9, 0.2, 0.2, 0.18)   # red
        }
        
        has_fut   <- is.finite(sim$logCOR_paths[i,"fut"])
        has_ia    <- is.finite(sim$logCOR_paths[i,"ia"])
        has_final <- is.finite(sim$logCOR_paths[i,"final"])
        
        if (!has_fut) next
        
        y_pos <- numeric(0)
        vals  <- numeric(0)
        
        if (has_fut)   { y_pos <- c(y_pos, 1); vals <- c(vals, sim$logCOR_paths[i,"fut"]) }
        if (has_ia)    { y_pos <- c(y_pos, 3); vals <- c(vals, sim$logCOR_paths[i,"ia"])  }
        if (has_final) { y_pos <- c(y_pos, 5); vals <- c(vals, sim$logCOR_paths[i,"final"]) }
        
        if (use_cor_scale) vals <- exp(vals)
        
        if (length(vals) >= 1) {
          lines(vals, y_pos, col = col_rgb, lwd = 0.9, lty = 2)
        }
      }
    }
  }
  
  # ── Points & boxes ────────────────────────────────────────────────────────
  for (i in seq_along(groups)) {
    vals <- groups[[i]]
    n <- length(vals)
    if (n == 0) next
    
    set.seed(2025 + i * 100)
    jitter_y <- runif(n, -jitter_height, jitter_height)
    
    col_pt <- if (group_names[i] == "All @ futility") rgb(0.1,0.4,0.9,point_alpha) else
      if (grepl("futility", group_names[i])) rgb(0.9,0.15,0.15,point_alpha) else
        if (grepl("success", group_names[i])) rgb(0.1,0.65,0.1,point_alpha) else
          rgb(0.3,0.3,0.3,point_alpha)
    
    points(vals, i + jitter_y, pch = 19, cex = point_cex, col = col_pt)
    
    if (n >= 3) {
      q <- quantile(vals, c(0.25,0.5,0.75))
      iqr <- q[3] - q[1]
      lower_whisk <- min(vals[vals >= q[1] - 1.5*iqr])
      upper_whisk <- max(vals[vals <= q[3] + 1.5*iqr])
      
      rect(q[1], i - box_width/2, q[3], i + box_width/2,
           col = rgb(0.88,0.93,1,0.5), border = "steelblue", lwd = 1.4)
      segments(q[2], i - box_width/2, q[2], i + box_width/2, lwd = 5, col = "royalblue3")
      segments(lower_whisk, i, q[1], i, lwd = 2.4, col = "midnightblue")
      segments(q[3], i, upper_whisk, i, lwd = 2.4, col = "midnightblue")
      cap_len <- box_width * 0.4
      segments(lower_whisk, i - cap_len, lower_whisk, i + cap_len, lwd = 2.4, col = "midnightblue")
      segments(upper_whisk, i - cap_len, upper_whisk, i + cap_len, lwd = 2.4, col = "midnightblue")
    }
  }
  
  # ── Two right-side columns: N per sim (planned N at that stage) + % of sims ─
  p_fut   <- mean(sim$stop_fut,   na.rm = TRUE)
  p_ia    <- mean(sim$stop_ia,    na.rm = TRUE)
  p_final <- mean(sim$stop_final, na.rm = TRUE)
  
  # Probabilities aligned with the original 6 groups (then filtered by keep):
  # 1 All@futility         -> 1
  # 2 Stopped futility     -> p_fut
  # 3 All@interim          -> 1 - p_fut
  # 4 Stopped IA success   -> p_ia
  # 5 All@final            -> 1 - p_fut - p_ia
  # 6 Stopped final success-> p_final
  props_all <- c(
    1,
    p_fut,
    1 - p_fut,
    p_ia,
    1 - p_fut - p_ia,
    p_final
  )
  props <- props_all[keep]
  
  # Planned N used in each simulation for each of the 6 groups (then filtered by keep)
  n_all <- c(
    sim$n_at_fut,  # All @ futility
    sim$n_at_fut,  # Stopped futility
    sim$n_at_ia,   # All @ interim
    sim$n_at_ia,   # Stopped IA success
    sim$n_total,   # All @ final
    sim$n_total    # Stopped final success
  )
  n_col <- n_all[keep]
  
  usr <- par("usr")
  x_range <- usr[2] - usr[1]
  
  # right column (%)
  x_text_pct <- usr[2] + 0.14 * x_range
  # left column (N)
  x_text_n   <- usr[2] + 0.05 * x_range
  
  for (i in seq_along(groups)) {
    
    col_text <- if (group_names[i] == "All @ futility") "black" else
      if (grepl("futility", group_names[i])) "firebrick" else
        if (grepl("success", group_names[i])) "forestgreen" else "gray30"
    
    # N column (planned N at this stage)
    text(x_text_n, i, sprintf("%d", as.integer(round(n_col[i]))),
         adj = 0, cex = 1.05, font = 2, col = "gray30", xpd = TRUE)
    
    # % column
    text(x_text_pct, i, sprintf("%.1f%%", 100 * props[i]),
         adj = 0, cex = 1.05, font = 2, col = col_text, xpd = TRUE)
  }
  
  # Column headers
  text(x_text_n,   length(groups) + 0.75, "N per sim", adj = 0, cex = 1.05, font = 2,
       col = "gray40", xpd = TRUE)
  text(x_text_pct, length(groups) + 0.75, "% of sims", adj = 0, cex = 1.05, font = 2,
       col = "gray40", xpd = TRUE)
  
  # ── Footnotes ─────────────────────────────────────────────────────────────
  power <- 100 * (mean(sim$stop_ia, na.rm = TRUE) + mean(sim$stop_final, na.rm = TRUE))
  mtext(sprintf("True COR = %.2f • Futility @ %.0f%% • IA @ %.0f%% • Power ≈ %.1f%% • %s",
                COR_true, 100*futility_frac, 100*info_frac, power, scale_note),
        side = 1, line = 2.5, cex = 0.85, col = "gray50")
  
  mtext("Winner's curse: success groups are systematically too optimistic (shifted left)",
        side = 1, line = 4.5, cex = 0.85, col = "gray50")
  
  invisible(NULL)
}

# ─────────────────────────────────────────────────────────────────────────────
#   UI
# ─────────────────────────────────────────────────────────────────────────────

ui <- page_sidebar(
  title = "Ordinal Non-Inferiority Trial Simulator + Winner's Curse v1.8 (rpact α-spending)",
  
  sidebar = sidebar(
    h4("Simulation Settings"),
    
    numericInput("n_total", "Total sample size", value = 600, min = 200, step = 50),
    numericInput("n_sims",  "Number of simulations", value = 1000, min = 100, max = 10000, step = 100),
    numericInput("seed",    "Random seed", value = 202506, min = 1),
    
    numericInput("COR_true", "True COR (treatment effect)", value = 1.3, min = 0.5, step = 0.05),
    numericInput("COR_NI",   "Non-inferiority margin (M)", value = 1.6, min = 1.0, step = 0.1),
    
    sliderInput("futility_frac", "Futility look fraction", min = 0.2, max = 0.7, value = 0.5, step = 0.05),
    sliderInput("info_frac",     "Interim look fraction",  min = 0.5, max = 0.95, value = 0.80, step = 0.05),
    
    numericInput("futility_p", "Futility p-value threshold", value = 0.70, min = 0.5, max = 0.95, step = 0.05),
    
    helpText("Efficacy boundaries are derived from rpact using one-sided α = 0.025 and the IA fraction above."),
    
    textInput("p_control_txt", "Control probabilities (comma sep)",
              value = "0.04, 0.02, 0.45, 0.34, 0.15"),
    
    checkboxInput("show_traj_success",
                  "Show trajectories – successful trials (green dashed)",
                  value = FALSE),
    
    checkboxInput("show_traj_fail",
                  "Show trajectories – failed/continued trials (red dashed)",
                  value = FALSE),
    
    hr(),
    h5("Plot options"),
    checkboxInput("use_cor_scale",
                  "Display on odds ratio (COR) scale instead of log",
                  value = FALSE),
    sliderInput("xlim_log_low",  "X lower limit (log scale)",  min = -6, max = 0, value = -3, step = 0.25),
    sliderInput("xlim_log_high", "X upper limit (log scale)",  min = 0,  max = 7, value = 4,  step = 0.25),
    
    actionButton("run_btn", "Run Simulation", class = "btn-primary btn-lg", icon = icon("play")),
    
    hr(),
    helpText("Simulation may take 10–90 s depending on nSims.")
  ),
  
  card(
    card_header("Results"),
    tabsetPanel(
      tabPanel(
        "Summary Table",
        verbatimTextOutput("status"),
        tableOutput("summary_table"),
        tags$hr(),
        verbatimTextOutput("rpact_info")
      ),
      
      tabPanel(
        "Expected Sample Size Breakdown",
        p("Decomposition of the expected (average) sample size by stopping stage."),
        tableOutput("ess_breakdown"),
        tags$hr(),
        verbatimTextOutput("ess_total_note")
      ),
      
      tabPanel("Winner's Curse Plot", plotOutput("boxplot", height = "750px")),
      
      tabPanel("About",
               h4("About"),
               p("Ordinal NI trial simulator with group-sequential stopping."),
               p("Efficacy boundaries are computed via rpact O'Brien–Fleming type alpha-spending at the selected interim fraction."),
               p("Winner's curse / selection bias may occur when stopping early for efficacy."))
    )
  )
)

# ─────────────────────────────────────────────────────────────────────────────
#   Server
# ─────────────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {
  
  sim_result <- eventReactive(input$run_btn, {
    req(input$run_btn)
    output$status <- renderText("Running simulations...")
    
    p_control <- parse_probs(input$p_control_txt)
    msg <- validate_probs(p_control)
    if (!is.null(msg)) {
      showNotification(msg, type = "error", duration = 10)
      return(NULL)
    }
    
    # --- rpact efficacy boundaries (one-sided alpha=0.025) ---
    # Analyses at informationRates = c(info_frac, 1)
    # typeOfDesign = "asOF" gives an O'Brien–Fleming type alpha-spending approximation
    design <- rpact::getDesignGroupSequential(
      sided = 1,
      alpha = 0.025,
      informationRates = c(input$info_frac, 1),
      typeOfDesign = "asOF"
    )
    
    # rpact criticalValues are upper-tail; your NI success uses lower-tail (Z <= bound)
    zcrit1 <- -design$criticalValues[1]
    zcrit2 <- -design$criticalValues[2]
    
    sim <- simulate_obf_ordinal(
      COR_true      = input$COR_true,
      COR_NI        = input$COR_NI,
      n_total       = input$n_total,
      futility_frac = input$futility_frac,
      info_frac     = input$info_frac,
      zcrit1        = zcrit1,
      zcrit2        = zcrit2,
      futility_p    = input$futility_p,
      p_control     = p_control,
      seed          = input$seed,
      nSims         = input$n_sims,
      show_progress = TRUE
    )
    
    sim$rpact_design <- design
    sim
  }, ignoreNULL = TRUE)
  
  output$status <- renderText({
    if (is.null(sim_result())) "Click 'Run Simulation' to start." else "Simulation complete."
  })
  
  output$summary_table <- renderTable({
    req(sim_result())
    sim_table(sim_result())
  }, digits = 3, align = "l")
  
  output$rpact_info <- renderPrint({
    req(sim_result())
    d <- sim_result()$rpact_design
    
    cat("rpact group-sequential design (efficacy)\n")
    cat("--------------------------------------\n")
    cat(sprintf("Type: %s | sided: %d | overall one-sided alpha: %.3f\n",
                d$typeOfDesign, d$sided, d$alpha))
    cat(sprintf("Information rates: %s\n",
                paste0(round(d$informationRates, 3), collapse = ", ")))
    
    cat(sprintf("rpact critical values (upper-tail Z): %s\n",
                paste0(round(d$criticalValues, 4), collapse = ", ")))
    
    lower_bounds <- -d$criticalValues
    cat(sprintf("=> NI success boundaries used here (lower-tail Z): %s\n",
                paste0(round(lower_bounds, 4), collapse = ", ")))
    
    stage_alpha_lower <- pnorm(lower_bounds)
    cat(sprintf("Stage-wise alpha (lower-tail): %s\n",
                paste0(signif(stage_alpha_lower, 4), collapse = ", ")))
    
    if (!is.null(d$alphaSpent)) {
      cat(sprintf("Cumulative alpha spent (rpact alphaSpent): %s\n",
                  paste0(signif(d$alphaSpent, 4), collapse = ", ")))
    }
  })
  
  output$ess_breakdown <- renderTable({
    req(sim_result())
    expected_n_breakdown(sim_result())
  }, align = "l")
  
  output$ess_total_note <- renderPrint({
    req(sim_result())
    df <- expected_n_breakdown(sim_result())
    
    cat(sprintf("Expected (average) sample size = %.1f\n\n", sum(df$Contribution)))
    cat("Note:\n")
    cat("- Probabilities are unconditional.\n")
    cat("- Contribution = Probability × Sample size at that stage.\n")
    cat("- The sum of contributions equals the expected sample size.\n")
  })
  
  output$boxplot <- renderPlot({
    req(sim_result())
    selection_boxplot(
      sim               = sim_result(),
      COR_true          = input$COR_true,
      COR_NI            = input$COR_NI,
      futility_frac     = input$futility_frac,
      info_frac         = input$info_frac,
      show_traj_success = input$show_traj_success,
      show_traj_fail    = input$show_traj_fail,
      use_cor_scale     = input$use_cor_scale,
      xlim_log_low      = input$xlim_log_low,
      xlim_log_high     = input$xlim_log_high
    )
  })
}

# ─────────────────────────────────────────────────────────────────────────────
#   Run the app
# ─────────────────────────────────────────────────────────────────────────────

shinyApp(ui, server)