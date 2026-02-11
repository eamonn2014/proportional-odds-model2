# =============================================================================
#   Shiny App: Ordinal Non-Inferiority Trial Simulator with Winner's Curse
#   v3.3 — Added user-adjustable CP futility threshold at IA2
# =============================================================================

# [R] Core Shiny framework for reactive web apps in R.
library(shiny)

# [R] Adds nicer UI widgets (e.g., styled buttons, sliders).
library(shinyWidgets)

# [R] MASS provides polr() for proportional-odds ordinal regression.
library(MASS)

# [R] dplyr for data manipulation verbs (mutate, across, etc.).
library(dplyr)

# [R] bslib gives Bootstrap theming support for Shiny UIs (page_sidebar, cards).
library(bslib)

# [STAT] rpact provides group sequential design utilities, including O'Brien-Fleming
#        alpha-spending critical values under various design families.
library(rpact)

# [R] future / future.apply allow parallelized simulation across CPU cores.
library(future)
library(future.apply)

# [R] progressr supports progress reporting that works with Shiny + parallel code.
library(progressr)

# ─────────────────────────────────────────────────────────────────────────────
#   Helpers
# ─────────────────────────────────────────────────────────────────────────────

# [R] Inverse logit; used to map linear predictors to probabilities.
# [STAT] In a cumulative logit/proportional odds model, P(Y <= k | x) = logit^{-1}(theta_k - beta x).
ilogit <- function(z) 1 / (1 + exp(-z))

# [R] Parse a comma-separated probability vector from a textbox string.
#     trimws() removes whitespace; strsplit() tokenizes; as.numeric coerces.
parse_probs <- function(txt) as.numeric(trimws(unlist(strsplit(txt, ","))))

# [R] Defensive input validation (fail fast with a clear message).
# [STAT] Ensures a proper probability mass function (pmf): positive probs, sum to 1,
#        and at least 3 categories (ordinal needs >= 3 for meaningful structure).
validate_probs <- function(p) {
  if (any(is.na(p)))           return("Probabilities must be numeric (comma-separated).")
  if (length(p) < 3)           return("Need at least 3 categories.")
  if (any(p <= 0))             return("All probabilities must be > 0.")
  if (abs(sum(p) - 1) > 1e-6)  return(sprintf("Must sum to ~1 (got %.6f).", sum(p)))
  NULL
}

# [STAT] Convert a *control* pmf into proportional-odds "cutpoints" (thresholds) theta.
#        If x=0 denotes control, then:
#           P(Y <= k | x=0) = ilogit(theta_k)
#        Hence theta_k = qlogis( P(Y <= k | control) ) for k = 1,...,K (excluding last category).
# [R] cumsum builds the cumulative distribution; qlogis is the logit quantile function.
theta_from_control_pmf <- function(p_control) {
  cum_control <- cumsum(p_control)
  qlogis(cum_control[seq_len(length(p_control) - 1)])
}

# [STAT] Given thresholds theta and a treatment effect beta in a proportional-odds model,
#        derive the treated pmf by:
#           CDF_k(x) = P(Y <= k | x) = ilogit(theta_k - beta*x)
#        then pmf = differences of the CDF (with 0 and 1 boundaries).
#        Here, x=1 represents treatment arm, x=0 control.
# [R] diff(c(0,cdf)) converts cumulative probabilities to category probabilities.
pmf_from_beta <- function(theta, beta, x = 1) {
  cdf <- ilogit(theta - beta * x)
  cdf <- c(cdf, 1)                # [STAT] Force CDF(K) = 1 for the top category boundary.
  pmf <- diff(c(0, cdf))          # [STAT] pmf_j = P(Y=j) = CDF_j - CDF_{j-1}
  pmf[pmf < 0] <- 0               # [R] Numerical guard (can happen due to floating error).
  pmf / sum(pmf)                  # [STAT] Renormalize so pmf sums to 1.
}

# [STAT] Fit proportional-odds model via polr() and extract log cumulative odds ratio (logCOR)
#        for trt (T vs C).
# [R] try(..., silent=TRUE) prevents a single failure from stopping the simulation batch.
# [R] Hess=TRUE requests Hessian for vcov() standard errors.
fit_logCOR <- function(df) {
  fit <- try(MASS::polr(y ~ trt, data = df, Hess = TRUE), silent = TRUE)
  if (inherits(fit, "try-error") || is.null(fit$coefficients)) return(NULL)
  
  # [R] polr() stores regression coefficients in fit$coefficients; thresholds elsewhere.
  coef_names <- names(fit$coefficients)
  
  # [R] Identify the treatment coefficient by name (e.g., "trtT").
  trt_name <- coef_names[grepl("^trt", coef_names)][1]
  if (is.na(trt_name)) return(NULL)
  
  # [STAT] Under proportional odds, exp(beta_hat) is the cumulative odds ratio (COR).
  logCOR_hat <- as.numeric(fit$coefficients[[trt_name]])
  
  # [R] vcov(fit) uses the Hessian to approximate var-cov matrix of coefficients.
  v <- vcov(fit)
  se <- sqrt(v[trt_name, trt_name])
  
  # [R] Protect against degenerate SEs (near 0) or numerical issues.
  if (!is.finite(se) || se < 1e-8) return(NULL)
  list(logCOR_hat = logCOR_hat, se = se)
}

# [STAT] Conditional power for one-sided non-inferiority (NI) with information fraction scaling.
#        Convention used here:
#          - Lower Z is better (more evidence treatment is NOT worse than NI margin).
#          - Test statistic at look i: Z_i = (beta_hat_i - beta_NI) / SE_i
#        Success criterion at next look uses zcrit_next (a negative boundary) and we compute:
#          CP = P(Z_next <= zcrit_next | current Z_obs, assumed theta_assumed)
#
# [STAT] The information scaling is approximated by:
#          Z_next ≈ delta_next + sqrt(info_current/info_next) * Z_obs + Normal(0,1)
#        where delta_next = (theta_assumed - beta_NI) * sqrt(info_next)
#
# [R] pnorm(zcrit_next, mean=..., sd=1) yields P(Z <= boundary) since we want lower-tail.
conditional_power_ni <- function(z_obs, info_current, info_next, theta_assumed, zcrit_next, beta_NI) {
  delta <- (theta_assumed - beta_NI) * sqrt(info_next)
  mean_final_z <- delta + sqrt(info_current / info_next) * z_obs
  pnorm(zcrit_next, mean = mean_final_z, sd = 1)
}

# [STAT] Expected sample size (ESS) decomposition:
#        E[N] = sum_s P(stop at stage s) * N_s.
# [R] Produces a tidy data frame with stage probabilities and contributions.
expected_n_breakdown <- function(sim) {
  p_fut        <- mean(sim$stop_fut, na.rm = TRUE)
  p_ia_success <- mean(sim$stop_ia,  na.rm = TRUE)
  p_low_cp     <- mean(sim$stop_fut_low_cp, na.rm = TRUE)
  
  # [STAT] Reaching final = not stopped by any earlier stopping rule.
  p_final      <- mean(!(sim$stop_fut | sim$stop_ia | sim$stop_fut_low_cp), na.rm = TRUE)
  
  df <- data.frame(
    Stage       = c("Futility stop", "IA success stop", "IA low-CP futility stop", "Final analysis"),
    Probability = c(p_fut, p_ia_success, p_low_cp, p_final),
    N_at_stage  = c(sim$n_at_fut, sim$n_at_ia, sim$n_at_ia, sim$n_total),
    check.names = FALSE
  )
  
  # [STAT] Contribution is the term P(stage)*N(stage) that sums to ESS.
  df$Contribution <- df$Probability * df$N_at_stage
  
  # [R] Base pipe `|>` (R >= 4.1) used here; mutate rounds for presentation only.
  df |> mutate(Probability = round(Probability, 3), Contribution = round(Contribution, 1))
}

# ─────────────────────────────────────────────────────────────────────────────
#   Simulation Engine
# ─────────────────────────────────────────────────────────────────────────────

# [STAT] simulate_obf_ordinal:
#        Simulates a 3-look design:
#          IA1: Futility (p-value threshold rule)
#          IA2: Efficacy using O'Brien-Fleming boundary (via rpact critical values)
#          IA2: Additional CP-based futility override (user-defined CP threshold)
#          Final: Efficacy boundary
#
# [STAT] Endpoint is ordinal; analysis uses proportional odds model -> COR (cumulative odds ratio).
#        COR > 1 implies treatment worse (per your Wiki). So "success" means COR sufficiently < NI margin.
#
# [R] Parallelization: future_lapply runs `one_sim` over nSims in parallel.
simulate_obf_ordinal <- function(COR_true, COR_NI, n_total, futility_frac, info_frac,
                                 zcrit1, zcrit2, futility_p, p_control, cp_threshold = 0.2,
                                 seed = 1234, nSims = 1000, workers = NULL) {
  msg <- validate_probs(p_control)
  if (!is.null(msg)) stop(msg)
  
  # [STAT] Compute ordinal cutpoints consistent with the specified control pmf.
  theta <- theta_from_control_pmf(p_control)
  
  # [STAT] K = number of cumulative logits = (#categories - 1).
  K <- length(p_control) - 1
  
  # [STAT] True effect and NI margin on log scale (since model coefficient is log(COR)).
  beta_true <- log(COR_true)
  beta_NI   <- log(COR_NI)
  
  # [STAT] Control pmf is as specified; treatment pmf is derived under proportional odds with beta_true.
  pi_control <- p_control
  pi_treat   <- pmf_from_beta(theta, beta_true)
  
  # [R] Utility for 1:1 randomization that handles odd N deterministically.
  split_n <- function(N) { nC <- floor(N / 2); list(nC = nC, nT = N - nC) }
  
  # [STAT] Planned sample sizes at each look using information-rate approximation:
  #        - IA1 futility look at futility_frac * n_total
  #        - IA2 efficacy look at info_frac * n_total
  n_fut <- round(futility_frac * n_total)
  n1    <- round(info_frac     * n_total)
  
  # [STAT] Futility threshold expressed on Z scale.
  #        If futility_p = 0.70, z_fut ~ 0.524. Rule below stops when Z_fut > z_fut,
  #        i.e., not sufficiently negative (not enough evidence for NI).
  z_fut <- qnorm(futility_p)
  
  # [R] Single simulation replication; returns a list so we can vapply() later.
  one_sim <- function(i) {
    
    # [R] Pre-allocate outputs for speed and to keep shape consistent even on failures.
    out <- list(
      Z_fut = NA_real_, COR_fut = NA_real_,
      Z1    = NA_real_, COR1    = NA_real_,
      Z2    = NA_real_, COR2    = NA_real_,
      
      stop_fut = FALSE, stop_ia = FALSE, stop_final = FALSE, stop_fut_low_cp = FALSE,
      
      logCOR_fut   = NA_real_,
      logCOR_ia    = NA_real_,
      logCOR_final = NA_real_,
      
      CP_after_fut_to_ia_obs   = NA_real_,
      CP_after_ia_to_final_obs = NA_real_
    )
    
    # ---------------------------
    # Stage 1: Futility look (IA1)
    # ---------------------------
    
    s_f <- split_n(n_fut)
    
    # [R] Sample ordinal outcomes as integers 0..K by drawing category indices 1..(K+1) then subtract 1.
    # [STAT] This implements multinomial sampling from the specified pmf (control/treat).
    yCf <- sample.int(K + 1, s_f$nC, TRUE, pi_control) - 1
    yTf <- sample.int(K + 1, s_f$nT, TRUE, pi_treat)   - 1
    
    # [R] Build an ordered factor outcome (polr expects ordered factor for ordinal regression).
    df_f <- data.frame(
      y   = factor(c(yCf, yTf), ordered = TRUE, levels = 0:K),
      trt = factor(rep(c("C", "T"), c(s_f$nC, s_f$nT)), levels = c("C", "T"))
    )
    
    fit_f <- fit_logCOR(df_f)
    
    if (!is.null(fit_f)) {
      # [STAT] NI Z-statistic: (beta_hat - beta_NI)/SE
      #        - If beta_hat < beta_NI, Z is negative => evidence that COR < COR_NI.
      out$Z_fut <- (fit_f$logCOR_hat - beta_NI) / fit_f$se
      
      # [STAT] Store effect on COR and log(COR) scales.
      out$COR_fut    <- exp(fit_f$logCOR_hat)
      out$logCOR_fut <- fit_f$logCOR_hat
      
      # [STAT] Futility stopping:
      #        stop when Z_fut is *too large* (i.e., not sufficiently negative).
      #        This matches a "lack of promise" rule.
      if (out$Z_fut > z_fut) {
        out$stop_fut <- TRUE
        return(out)
      }
    }
    
    # -----------------------------------------
    # Stage 2: Interim efficacy look (IA2 / IA)
    # -----------------------------------------
    
    # [R] Incremental sampling: add only the extra participants to reach n1,
    #     then rbind() to combine with prior stage.
    n_add_ia <- n1 - n_fut
    s_add    <- split_n(n_add_ia)
    
    yCa <- sample.int(K + 1, s_add$nC, TRUE, pi_control) - 1
    yTa <- sample.int(K + 1, s_add$nT, TRUE, pi_treat)   - 1
    
    df_add <- data.frame(
      y   = factor(c(yCa, yTa), ordered = TRUE, levels = 0:K),
      trt = factor(rep(c("C", "T"), c(s_add$nC, s_add$nT)), levels = c("C", "T"))
    )
    
    # [R] Combine stage-1 and added interim data.
    df1  <- rbind(df_f, df_add)
    fit1 <- fit_logCOR(df1)
    
    if (!is.null(fit1)) {
      out$Z1 <- (fit1$logCOR_hat - beta_NI) / fit1$se
      out$COR1 <- exp(fit1$logCOR_hat)
      out$logCOR_ia <- fit1$logCOR_hat
      
      # [STAT] Approximate information fractions with sample-size fractions.
      #        This is common in simple simulations when variance scales ~ 1/N.
      info_fut <- n_fut / n_total
      info_ia  <- n1    / n_total
      
      # [STAT] "Observed CP" from futility look to interim look:
      #        Uses theta_assumed = observed logCOR at futility (winner's-curse prone).
      out$CP_after_fut_to_ia_obs <- conditional_power_ni(
        z_obs         = out$Z_fut,
        info_current  = info_fut,
        info_next     = info_ia,
        theta_assumed = out$logCOR_fut,
        zcrit_next    = zcrit1,
        beta_NI       = beta_NI
      )
      
      # [STAT] Interim efficacy stopping rule:
      #        Stop for success if Z1 <= zcrit1 (boundary is negative because "lower is better").
      if (out$Z1 <= zcrit1) {
        out$stop_ia <- TRUE
        return(out)
      }
      
      # --------------------------------------------
      # IA2 futility override using Conditional Power
      # --------------------------------------------
      
      info_final <- 1.0
      
      # [STAT] Conditional power from interim to final:
      #        Uses theta_assumed = observed interim logCOR. This is a *plug-in* approach.
      #        It tends to be optimistic when early estimates are inflated (winner's curse),
      #        but here it is intentionally part of the educational design.
      out$CP_after_ia_to_final_obs <- conditional_power_ni(
        z_obs         = out$Z1,
        info_current  = info_ia,
        info_next     = info_final,
        theta_assumed = out$logCOR_ia,
        zcrit_next    = zcrit2,
        beta_NI       = beta_NI
      )
      
      # [STAT] CP-based futility: stop if the probability of ultimate success is below threshold.
      #        This is a common operational futility rule in adaptive / group sequential designs.
      if (out$CP_after_ia_to_final_obs < cp_threshold) {
        out$stop_fut_low_cp <- TRUE
        return(out)
      }
    }
    
    # -------------------------
    # Stage 3: Final analysis
    # -------------------------
    
    n_add_final <- n_total - n1
    s_final     <- split_n(n_add_final)
    
    yCf2 <- sample.int(K + 1, s_final$nC, TRUE, pi_control) - 1
    yTf2 <- sample.int(K + 1, s_final$nT, TRUE, pi_treat)   - 1
    
    df_final_add <- data.frame(
      y   = factor(c(yCf2, yTf2), ordered = TRUE, levels = 0:K),
      trt = factor(rep(c("C", "T"), c(s_final$nC, s_final$nT)), levels = c("C", "T"))
    )
    
    df2  <- rbind(df1, df_final_add)
    fit2 <- fit_logCOR(df2)
    
    if (!is.null(fit2)) {
      out$Z2 <- (fit2$logCOR_hat - beta_NI) / fit2$se
      out$COR2 <- exp(fit2$logCOR_hat)
      out$logCOR_final <- fit2$logCOR_hat
      
      # [STAT] Final success criterion:
      #        stop_final=TRUE means "declared NI success at final" (your naming is consistent with plot labels).
      if (out$Z2 <= zcrit2) out$stop_final <- TRUE
    }
    
    out
  }
  
  # [R] Choose parallel workers: default = all cores minus one (keep UI responsive).
  if (is.null(workers)) workers <- max(1, future::availableCores() - 1)
  
  # [R] Save & restore prior plan (good citizen: don’t permanently change user’s future plan).
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  
  # [R] multisession uses separate R processes (portable across platforms).
  future::plan(multisession, workers = workers)
  
  # [R] Set seed for reproducibility of the *main* process; future.seed handles worker RNG.
  set.seed(seed)
  
  # [R] progressr::withProgressShiny wires progress updates into Shiny UI.
  # [R] future_lapply parallelizes over sims; p() increments progress each finished sim.
  sims <- progressr::withProgressShiny(
    message = "Running simulations...",
    detail  = "Processing...",
    value   = 0,
    {
      p <- progressr::progressor(along = seq_len(nSims))
      future.apply::future_lapply(seq_len(nSims), FUN = function(i) {
        res <- one_sim(i)
        p()
        res
      }, future.seed = TRUE)
    }
  )
  
  # [R] Extract per-sim paths for plotting trajectories:
  #     fut, ia, final are the *logCOR* values at each look (may be NA if model failed / stage not reached).
  # [R] vapply is type-safe and faster than sapply.
  logCOR_paths <- cbind(
    fut   = vapply(sims, `[[`, numeric(1), "logCOR_fut"),
    ia    = vapply(sims, `[[`, numeric(1), "logCOR_ia"),
    final = vapply(sims, `[[`, numeric(1), "logCOR_final")
  )
  
  # [R] Return a list-of-vectors structure that is easy for downstream plotting and summarizing.
  list(
    Z_fut_all   = vapply(sims, `[[`, numeric(1),  "Z_fut"),
    COR_fut_all = vapply(sims, `[[`, numeric(1),  "COR_fut"),
    Z1_all      = vapply(sims, `[[`, numeric(1),  "Z1"),
    COR1_all    = vapply(sims, `[[`, numeric(1),  "COR1"),
    Z2_all      = vapply(sims, `[[`, numeric(1),  "Z2"),
    COR2_all    = vapply(sims, `[[`, numeric(1),  "COR2"),
    
    stop_fut        = vapply(sims, `[[`, logical(1), "stop_fut"),
    stop_ia         = vapply(sims, `[[`, logical(1), "stop_ia"),
    stop_final      = vapply(sims, `[[`, logical(1), "stop_final"),
    stop_fut_low_cp = vapply(sims, `[[`, logical(1), "stop_fut_low_cp"),
    
    logCOR_paths = logCOR_paths,
    
    CP_after_fut_to_ia_obs   = vapply(sims, `[[`, numeric(1), "CP_after_fut_to_ia_obs"),
    CP_after_ia_to_final_obs = vapply(sims, `[[`, numeric(1), "CP_after_ia_to_final_obs"),
    
    # [R] Store design constants for downstream display.
    n_at_fut = n_fut, n_at_ia = n1, n_total = n_total, nSims = nSims,
    z_fut = z_fut, zcrit1 = zcrit1, zcrit2 = zcrit2,
    
    # [R] Placeholder; server fills this with rpact design object for printing.
    rpact_design = NULL
  )
}

# [R] sim_table creates stage-wise summaries of COR estimates under each stopping condition.
# [STAT] This table illustrates selection effects: early stopping can distort effect estimate distributions (winner's curse).
sim_table <- function(sim) {
  
  # [R] Helper that safely summarizes numeric vectors, guarding against all-NA.
  # [STAT] Reports central tendency + 95% empirical interval (2.5% and 97.5% quantiles).
  safe_summ_ext <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(c(N=0L, Min=NA, Max=NA, Mean=NA, Median=NA, `2.5%`=NA, `97.5%`=NA))
    q <- quantile(x, c(0.025, 0.5, 0.975), names = FALSE)
    c(N=length(x), Min=min(x), Max=max(x), Mean=mean(x), Median=q[2], `2.5%`=q[1], `97.5%`=q[3])
  }
  
  # [R] Subset estimates by which stopping rule fired.
  fut      <- safe_summ_ext(sim$COR_fut_all[sim$stop_fut])
  ia_suc   <- safe_summ_ext(sim$COR1_all[sim$stop_ia])
  ia_lowcp <- safe_summ_ext(sim$COR1_all[sim$stop_fut_low_cp])
  fin      <- safe_summ_ext(sim$COR2_all[sim$stop_final])
  
  # [R] Build a data.frame and round numeric columns for display.
  data.frame(
    Stage = c("Futility stop", "IA success stop", "IA low-CP futility", "Final success stop"),
    N     = c(fut["N"], ia_suc["N"], ia_lowcp["N"], fin["N"]),
    Min   = c(fut["Min"], ia_suc["Min"], ia_lowcp["Min"], fin["Min"]),
    Mean  = c(fut["Mean"], ia_suc["Mean"], ia_lowcp["Mean"], fin["Mean"]),
    Median= c(fut["Median"], ia_suc["Median"], ia_lowcp["Median"], fin["Median"]),
    `2.5%` = c(fut["2.5%"], ia_suc["2.5%"], ia_lowcp["2.5%"], fin["2.5%"]),
    `97.5%`= c(fut["97.5%"], ia_suc["97.5%"], ia_lowcp["97.5%"], fin["97.5%"]),
    Max   = c(fut["Max"], ia_suc["Max"], ia_lowcp["Max"], fin["Max"]),
    check.names = FALSE
  ) |>
    mutate(
      N = as.integer(round(N)),
      across(-c(Stage, N), ~ round(as.numeric(.x), 3))
    )
}

# ─────────────────────────────────────────────────────────────────────────────
#   Plotting Function (unchanged from previous version)
# ─────────────────────────────────────────────────────────────────────────────

# [R] selection_boxplot does “two-panel base R plotting”:
#     - Left panel: strip/box style view of estimates by stage group + optional trajectories
#     - Right panel: aligned "table" showing counts, N/Trial, % sims
# [STAT] Intended to visualize winner’s curse / selection bias by comparing distributions
#        of estimates conditional on stopping stage.
selection_boxplot <- function(sim, COR_true, COR_NI, futility_frac, info_frac,
                              show_traj_success = FALSE, show_traj_fail = FALSE,
                              use_cor_scale = FALSE,
                              xlim_log_low = -3, xlim_log_high = 4,
                              main = "Bias and Treatment Effect Estimates by Trial Stopping Stage") {
  
  # [R] User controls for axis limits on log scale; used either directly (log scale)
  #     or exponentiated (COR scale).
  log_min_nice <- xlim_log_low
  log_max_nice <- xlim_log_high
  
  if (use_cor_scale) {
    # [R] If plotting on COR scale, transform log values with exp() and axis limits accordingly.
    xlim_use     <- exp(c(log_min_nice, log_max_nice))
    x_transform  <- exp
    x_label_expr <- "Cumulative Odds Ratio (COR)"
  } else {
    xlim_use     <- c(log_min_nice, log_max_nice)
    x_transform  <- identity
    x_label_expr <- expression(log(Cumulative~Odds~Ratio))
  }
  
  # [STAT] Empirical operating characteristics from simulation:
  p_fut        <- mean(sim$stop_fut, na.rm = TRUE)
  p_ia_success <- mean(sim$stop_ia,  na.rm = TRUE)
  p_low_cp     <- mean(sim$stop_fut_low_cp, na.rm = TRUE)
  p_final_suc  <- mean(sim$stop_final, na.rm = TRUE)
  
  # [STAT] Power (for NI success) = success at interim + success at final.
  empirical_power <- p_ia_success + p_final_suc
  
  # [STAT] Probability of reaching final sample size.
  p_reach_final <- mean(!(sim$stop_fut | sim$stop_ia | sim$stop_fut_low_cp), na.rm = TRUE)
  
  # [STAT] Expected sample size (ESS) computed via stage probabilities.
  avg_n <- round(
    p_fut        * sim$n_at_fut +
      p_ia_success * sim$n_at_ia +
      p_low_cp     * sim$n_at_ia +
      p_reach_final* sim$n_total
  )
  
  # [STAT] Binomial CI for power via Wilson interval (approx) using z=1.96 for 95%.
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
  
  # [R] Group definitions are vectors of *logCOR* values by stage and stopping condition.
  # [STAT] Conditioning on stopping stage produces different sampling distributions (selection bias).
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
  
  # [R] Track counts before filtering empty groups.
  counts_actual <- sapply(groups, length)
  
  # [R] If plotting on COR scale, exponentiate all logCOR values.
  if (use_cor_scale) groups <- lapply(groups, exp)
  
  # [R] Remove empty or non-finite groups to prevent plotting errors.
  keep <- sapply(groups, function(g) length(g) > 0 && all(is.finite(g)))
  groups <- groups[keep]
  group_names <- names(groups)
  counts_actual <- counts_actual[keep]
  
  # [R] Precompute jitter for every simulated trial so:
  #     - points align with trajectories
  #     - jitter stays consistent across groups
  set.seed(202506)
  n_max <- nrow(sim$logCOR_paths)
  jitter_master <- runif(n_max, -0.30, 0.30)
  
  # [R] Save/restore graphics params; use layout for two panels.
  op <- par(no.readonly = TRUE)
  on.exit({ par(op); layout(1) }, add = TRUE)
  layout(matrix(c(1, 2), nrow = 1), widths = c(5.2, 2.2))
  
  # ---------------- Left panel: points + custom box overlays ----------------
  
  par(mar = c(15, 12, 6, 2), xpd = FALSE)
  
  # [R] Empty plot as canvas; y-axis is categorical group index.
  plot(
    0, type = "n",
    xlim = xlim_use,
    ylim = c(0.4, length(groups) + 0.9),
    xlab = "", ylab = "",
    yaxt = "n", las = 1,
    main = sprintf("%s\n(Expected Sample Size = %d)", main, avg_n)
  )
  
  footnote_cex <- 1.05
  
  # [R] Axis/annotation text; mtext places labels in margins.
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
  
  # [STAT] Vertical reference lines:
  #        - True effect (green dashed)
  #        - NI margin (red dotted)
  #        On log scale: log(COR). On COR scale: COR.
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
  
  # [R] Group labels on y axis.
  axis(2, at = seq_along(groups), labels = group_names, las = 1, cex.axis = 1.1)
  
  # [STAT] Use x_transform so lines are correct on both scales.
  abline(v = x_transform(log(COR_true)), lty = 2, col = "darkgreen", lwd = 2.5)
  abline(v = x_transform(log(COR_NI)),   lty = 3, col = "red",       lwd = 2.5)
  
  # [R] Light horizontal guides.
  abline(h = seq_along(groups), col = "gray92", lwd = 0.8)
  
  # --------------------
  # Trajectories – centered + fainter
  # --------------------
  # [R] Draw per-trial line segments across looks (fut -> ia -> final),
  #     using the SAME jitter as the points.
  # [STAT] Highlights regression-to-the-mean and selection effects (winner’s curse).
  if (show_traj_success || show_traj_fail) {
    
    group_y_map <- setNames(seq_along(groups), names(groups))
    
    for (i in seq_len(nrow(sim$logCOR_paths))) {
      
      # [STAT] Define "success" as NI declared at interim or final.
      is_success <- isTRUE(sim$stop_ia[i]) || isTRUE(sim$stop_final[i])
      
      # [STAT] Identify the CP-futility subgroup (distinct from stage-1 futility).
      is_low_cp  <- isTRUE(sim$stop_fut_low_cp[i])
      
      # [R] Respect UI toggles.
      if ((is_success && !show_traj_success) ||
          (!is_success && !show_traj_fail)) next
      
      path_x <- numeric(0)
      path_y <- numeric(0)
      
      # [R] Reuse deterministic jitter for this trial index.
      j <- jitter_master[i]
      
      # --- Futility ---
      if (is.finite(sim$logCOR_paths[i, "fut"]) &&
          "1 All @ futility" %in% names(group_y_map)) {
        
        path_x <- c(path_x, sim$logCOR_paths[i, "fut"])
        path_y <- c(path_y, group_y_map["1 All @ futility"] + j)
      }
      
      # --- Interim ---
      if (is.finite(sim$logCOR_paths[i, "ia"])) {
        
        # [R] Decide which interim group row the point belongs to for this trial.
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
      
      # --- Final ---
      if (is.finite(sim$logCOR_paths[i, "final"])) {
        
        # [R] Decide which final group row the point belongs to for this trial.
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
      
      # --- Draw line ---
      if (length(path_x) >= 2) {
        
        # [R] Color codes: green for success, purple for low-CP stop, red for non-success.
        col_line <- if (is_success) rgb(0,   0.7, 0,   0.18)  # faint green
        else if (is_low_cp) rgb(0.6, 0,   0.8, 0.22)         # faint purple
        else rgb(0.9, 0.1, 0.1, 0.18)                        # faint red
        
        lines(x_transform(path_x), path_y, col = col_line, lwd = 1)
      }
    }
  }
  
  # [R] Overlay points and a custom boxplot-like summary per group.
  # [STAT] The box/whisker shows IQR and typical range; median as thick line.
  for (ii in seq_along(groups)) {
    nm <- group_names[ii]
    
    # [R] Identify which trials correspond to this group for jitter alignment.
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
    
    # [R] Point colors by group type.
    col_p <- if (grepl("All @ futility", nm)) rgb(0.1, 0.4, 0.9, 0.25) else
      if (grepl("futility", nm))              rgb(0.9, 0.15,0.15,0.25) else
        if (grepl("failed", nm))              rgb(0.9, 0.15,0.15,0.25) else
          if (grepl("success", nm))           rgb(0.1, 0.65,0.1, 0.25) else
            if (grepl("low CP", nm))          rgb(0.6, 0.0, 0.8, 0.35) else
              rgb(0.3, 0.3, 0.3, 0.25)
    
    points(vals, ii + jitter_y, pch = 19, cex = 0.6, col = col_p)
    
    # [STAT] Draw IQR box and whiskers manually (to keep layout control).
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
  
  # ---------------- Right panel: aligned summary "table" ----------------
  
  par(mar = c(15, 1, 6, 1), xpd = FALSE)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0.4, length(groups) + 0.9))
  
  header_cex <- 1.15
  text(0.02, length(groups) + 0.75, "Trials (N)", adj = 0, font = 2, cex = header_cex)
  text(0.45, length(groups) + 0.75, "N/Trial",   adj = 0, font = 2, cex = header_cex)
  text(0.78, length(groups) + 0.75, "% Sims",    adj = 0, font = 2, cex = header_cex)
  
  data_cex <- 1.10
  
  # [R] Determine N per trial at each group row based on stage.
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
  
  # [STAT] Robust % of simulations per row:
  #        Each "All @" row is a *conditioning event* (e.g., all who reached interim),
  #        whereas "Stopped ..." rows represent stopping probabilities.
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
    
    # [R] Color of text mirrors point colors for interpretability.
    col_t <- if (grepl("All @ futility", names(groups)[ii])) "black" else
      if (grepl("futility", names(groups)[ii])) "firebrick" else
        if (grepl("failed", names(groups)[ii])) "firebrick" else
          if (grepl("success", names(groups)[ii])) "forestgreen" else
            if (grepl("low CP", names(groups)[ii])) "firebrick" else
              "gray30"
    
    text(0.02, ii, format(counts_actual[ii], big.mark=","), adj = 0, col = col_t,  cex = data_cex)
    text(0.45, ii, sprintf("%d", n_per_group[ii]),           adj = 0, col = "gray30", cex = data_cex)
    text(0.78, ii, sprintf("%.1f%%", 100 * props[ii]),       adj = 0, font = 2, col = col_t, cex = data_cex)
  }
}

# ─────────────────────────────────────────────────────────────────────────────
#   UI
# ─────────────────────────────────────────────────────────────────────────────

# [R] page_sidebar is a bslib layout: sidebar + main content card.
# [R] Inputs live in sidebar; outputs live in tabs in main panel.
ui <- page_sidebar(
  title = "Ordinal Endpoint, Non-Inferiority, Group Sequential Design, Trial Simulator v3.3",
  sidebar = sidebar(
    width = 350,
    
    # [R] actionButton triggers computation explicitly (better than auto-reacting on every input change).
    actionButton("run_btn", "Run Simulation",
                 class = "btn-primary btn-lg", icon = icon("play-circle"), width = "100%"),
    
    hr(style = "margin: 1.2em 0;"),
    
    tags$div(
      style = "padding: 0 8px;",
      h5("Trial & Simulation Settings"),
      
      # [R] numericInput/textInput/sliderInput define reactive inputs accessible as input$...
      numericInput("n_sims", "Number of simulations", value = 1000),
      numericInput("n_total", "Total sample size (1:1 rand.)", value = 600),
      textInput("p_control_txt", "Control probabilities", value = "0.04, 0.02, 0.45, 0.34, 0.15"),
      
      # [STAT] COR_true defines data-generating treatment effect; COR_NI is NI margin (worse-than limit).
      numericInput("COR_true", "True cumulative odds ratio (COR)", value = 1.0, step = 0.05),
      numericInput("COR_NI", "Non inferiority margin COR", value = 1.6, step = 0.1),
      
      # [STAT] Futility uses a p-value threshold at IA1; CP futility threshold at IA2 is additional rule.
      numericInput("futility_p", "Futility p-value threshold IA1", value = 0.70),
      numericInput("cp_threshold", "CP futility threshold at IA2", value = 0.2, min = 0, max = 1, step = 0.05),
      
      # [STAT] These determine look timings (information/sample-size fractions).
      sliderInput("futility_frac", "Futility look fraction", min = 0.2, max = 0.7, value = 0.5),
      sliderInput("info_frac", "Interim look fraction", min = 0.5, max = 0.95, value = 0.80),
      
      # [R] Seed makes simulations reproducible for a given parameter set.
      numericInput("seed", "Random seed", value = 202506)
    ),
    
    hr(style = "margin: 1.5em 0;"),
    
    tags$div(
      style = "padding: 0 8px;",
      h5("Plot Options"),
      
      checkboxInput("use_cor_scale", "Display on COR scale (instead of log)", value = FALSE),
      sliderInput("xlim_log_low",  "X-axis lower limit (log scale)", min = -6, max = 0, value = -1, step = 0.5),
      sliderInput("xlim_log_high", "X-axis upper limit (log scale)", min = 0,  max = 7, value = 2,  step = 0.5),
      
      hr(style = "margin: 1.2em 0; border-top: 1px dashed #ccc;"),
      tags$strong("Show trajectories"),
      checkboxInput("show_traj_success", "Successful trajectories", value = FALSE),
      checkboxInput("show_traj_fail",    "Failed / futility trajectories", value = FALSE),
      
      hr(),
      
      # [R] downloadButton triggers server-side downloadHandler to write a JPEG.
      downloadButton("download_plot", "Download Superb JPEG",
                     class = "btn-success", style = "width: 100%;")
    )
  ),
  
  # [R] card() holds the tabsetPanel of outputs.
  card(
    tabsetPanel(
      tabPanel("Operating Characteristics Plot",
               plotOutput("boxplot", height = "750px")),
      
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
      
      # [R] A static Wiki tab explaining the app; no server logic.
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
                   p("At the first look, the trial may stop for futility if the observed one-sided p-value is too large (i.e., insufficient evidence that treatment is non-inferior)."),
                   p("This approach evaluates the observed test statistic directly against a user-defined futility threshold."),
                   br(),
                   h4("Conditional Power (CP) Futility at IA2"),
                   p("At the second look, the trial may stop for futility based on Conditional Power."),
                   p("Conditional Power estimates the probability of ultimately achieving success at the final analysis, given:"),
                   tags$ul(
                     tags$li("The data observed so far, and"),
                     tags$li("An assumed true treatment effect.")
                   ),
                   p("In this simulator, the assumed future treatment effect is the observed interim COR."),
                   p(strong("Important:"), " This means the observed interim estimate is used to model future behavior of the trial. If the resulting conditional power falls below the user-defined threshold (e.g., 20%), the trial stops early for futility."),
                   hr(),
                   h3("2) Winner's Curse in Group Sequential Trials"),
                   p("Group sequential designs are vulnerable to selection bias, often referred to as the "),
                   strong("Winner's Curse."),
                   br(), br(),
                   p("When a trial stops early for success, it typically does so because the observed treatment effect is more extreme than average."),
                   p("As a result, the observed COR at the time of stopping tends to overestimate the true underlying treatment effect."),
                   br(),
                   h4("Why This Happens"),
                   tags$ul(
                     tags$li("Stopping rules select for extreme results."),
                     tags$li("Random variability inflates early large effects."),
                     tags$li("Trials that continue tend to regress toward the true effect.")
                   ),
                   p("Therefore, early-stopped trials frequently report treatment effects that are larger than the truth."),
                   br(),
                   h4("Implications"),
                   tags$ul(
                     tags$li("Overestimation of treatment benefit."),
                     tags$li("Replication studies often show smaller effects."),
                     tags$li("Regulatory and HTA interpretation should consider potential inflation.")
                   ),
                   br(),
                   h4("Selected References"),
                   tags$ul(
                     tags$li("Proschan MA, Lan KKG, Wittes JT. Statistical Monitoring of Clinical Trials."),
                     tags$li("Pocock SJ. Group Sequential Methods in the Design and Analysis of Clinical Trials."),
                     tags$li("Jennison C, Turnbull BW. Group Sequential Methods with Applications to Clinical Trials."),
                     tags$li("Schönbrodt FD, Wagenmakers EJ (2018). Sequential hypothesis testing can inflate effect sizes.")
                   ),
                   br(),
                   hr(),
                   p(style = "font-style: italic; color: gray;",
                     "This simulator is for educational and design exploration purposes only. It demonstrates how stopping rules influence power, bias, and expected sample size.")
               ))
    )
  )
)

# ─────────────────────────────────────────────────────────────────────────────
#   Server
# ─────────────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {
  
  # [R] reactiveVal stores mutable reactive state. Think: a single reactive "slot".
  #     When sim_store(sim) is called, dependent outputs re-render.
  sim_store <- reactiveVal(NULL)
  
  # [R] observeEvent runs code when input$run_btn changes (clicked).
  #     ignoreInit=TRUE prevents running on app start.
  observeEvent(input$run_btn, {
    
    # [STAT] Group sequential design via rpact:
    #        sided=1 => one-sided alpha
    #        alpha=0.025 => typical one-sided 2.5%
    #        informationRates=c(info_frac,1) => two looks: interim at info_frac and final at 1
    #        typeOfDesign="asOF" => O'Brien-Fleming-like alpha spending
    design <- getDesignGroupSequential(
      sided = 1,
      alpha = 0.025,
      informationRates = c(input$info_frac, 1),
      typeOfDesign = "asOF"
    )
    
    # [STAT] rpact criticalValues are for an upper-tail Z by default; your NI test is "lower is better".
    #        Hence you negate them so zcrit1 and zcrit2 are negative thresholds.
    sim <- simulate_obf_ordinal(
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
      nSims = input$n_sims
    )
    
    # [R] Store design object for later printing in UI.
    sim$rpact_design <- design
    
    # [R] Update reactive state.
    sim_store(sim)
  }, ignoreInit = TRUE)
  
  # [R] renderPlot re-runs whenever sim_store() or relevant inputs change.
  #     req() blocks until sim_store is non-NULL.
  output$boxplot <- renderPlot({
    req(sim_store())
    selection_boxplot(
      sim_store(),
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
  
  # [R] downloadHandler creates a file on demand; `content` is called with a temp filepath.
  output$download_plot <- downloadHandler(
    filename = function() {
      paste0("Simulation_Results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".jpg")
    },
    content = function(file) {
      req(sim_store())
      
      # [R] Open a JPEG device; any plot commands write into it until dev.off().
      jpeg(file, width = 12, height = 10, units = "in", res = 300, quality = 100)
      
      selection_boxplot(
        sim_store(),
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
  
  # [R] Simple status indicator for the UI.
  output$status <- renderText({
    if (is.null(sim_store())) "Click 'Run' to start." else "Complete."
  })
  
  # [R] renderTable shows a data.frame in UI. digits=3 controls printing formatting.
  output$summary_table <- renderTable({
    req(sim_store())
    sim_table(sim_store())
  }, digits = 3)
  
  output$ess_breakdown <- renderTable({
    req(sim_store())
    expected_n_breakdown(sim_store())
  }, digits = 4)
  
  output$ess_total_note <- renderPrint({
    req(sim_store())
    df <- expected_n_breakdown(sim_store())
    # [STAT] ESS is sum of probability-weighted N contributions across stages.
    cat(sprintf("ESS (Average Sample Size) = %.1f\n", sum(df$Contribution)))
  })
  
  output$rpact_info <- renderPrint({
    req(sim_store())
    d <- sim_store()$rpact_design
    
    # [STAT] alphaSpent indicates cumulative alpha used at each analysis in the spending function.
    cat(sprintf(
      "Total One-sided Alpha: %.4f\nAlpha Spent Stage 1 (IA): %.4f\nAlpha Spent Stage 2 (Final): %.4f\n",
      d$alpha, d$alphaSpent[1], d$alphaSpent[2]
    ))
    
    cat("\n--- Nominal Stage-wise P-values (1-sided) ---\n")
    
    # [STAT] Convert Z boundary to one-sided nominal p:
    #        For lower-tail test, p = Phi(z_boundary). You print pnorm(-criticalValue)
    #        because rpact returns upper-tail criticalValue by default.
    cat(sprintf(
      "Nominal p at Stage 1: %.6f\nNominal p at Stage 2: %.6f\n",
      pnorm(-d$criticalValues[1]),
      pnorm(-d$criticalValues[2])
    ))
  })
}

# [R] Launch the application (UI + server).
shinyApp(ui, server)