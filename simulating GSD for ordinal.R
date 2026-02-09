# rm(list=ls())
# 
# suppressPackageStartupMessages({
#   library(MASS)
#   library(dplyr)
#   library(ggplot2) # Added for the new plot
# })
# 
# 
# # ============================================================
# # 1) Core helpers
# # ============================================================
# ilogit <- function(z) 1 / (1 + exp(-z))
# parse_probs <- function(txt) as.numeric(trimws(unlist(strsplit(txt, ","))))
# 
# validate_probs <- function(p) {
#   if (any(is.na(p))) return("Probabilities must be numeric.")
#   if (length(p) < 3) return("Need at least 3 categories.")
#   if (abs(sum(p) - 1) > 1e-6) return(sprintf("Must sum to ~1 (got %.6f).", sum(p)))
#   NULL
# }
# 
# theta_from_control_pmf <- function(p_control) {
#   cum_control <- cumsum(p_control)
#   qlogis(cum_control[seq_len(length(p_control) - 1)])
# }
# 
# pmf_from_beta <- function(theta, beta, x) {
#   cdf <- ilogit(theta - beta * x)
#   cdf <- c(cdf, 1)
#   pmf <- diff(c(0, cdf))
#   pmf[pmf < 0] <- 0
#   pmf / sum(pmf)
# }
# 
# # ============================================================
# # 2) Updated Fit Function (Collects P-values)
# # ============================================================
# fit_logCOR <- function(df, beta_NI = 0) {
#   fit <- try(MASS::polr(y ~ trt, data = df, Hess = TRUE), silent = TRUE)
#   if (inherits(fit, "try-error") || is.null(fit$coefficients)) return(NULL)
#   
#   trt_name <- names(fit$coefficients)[grepl("^trt", names(fit$coefficients))][1]
#   logCOR_hat <- as.numeric(fit$coefficients[[trt_name]])
#   se <- sqrt(vcov(fit)[trt_name, trt_name])
#   
#   # Calculate p-value relative to the NI Margin (beta_NI)
#   z_stat <- (logCOR_hat - beta_NI) / se
#   p_val  <- pnorm(z_stat) # One-sided for NI
#   
#   list(logCOR_hat = logCOR_hat, se = se, p_val = p_val)
# }
# 
# # ============================================================
# # 3) Updated Simulation (Stores P-values)
# # ============================================================
# simulate_obf_ordinal <- function(COR_true, COR_NI, n_total, futility_frac, info_frac,
#                                  alpha1, alpha2, futility_p, p_control, 
#                                  seed = 1234, nSims = 500) {
#   set.seed(seed)
#   theta <- theta_from_control_pmf(p_control)
#   beta_true <- log(COR_true); beta_NI <- log(COR_NI)
#   pi_control <- p_control
#   pi_treat <- pmf_from_beta(theta, beta_true, x = 1)
#   
#   # Initialize storage
#   res_storage <- data.frame(
#     sim = rep(1:nSims, each=3),
#     stage = rep(c("Futility", "IA", "Final"), nSims),
#     logCOR = NA_real_,
#     pval = NA_real_,
#     stopped = FALSE
#   )
#   
#   n_fut <- round(futility_frac * n_total)
#   n1    <- round(info_frac * n_total)
#   z_fut <- qnorm(futility_p); zcrit1 <- qnorm(alpha1); zcrit2 <- qnorm(alpha2)
#   
#   for (i in 1:nSims) {
#     # Data Generation helpers
#     gen_data <- function(n) {
#       nC <- floor(n/2); nT <- n - nC
#       data.frame(
#         y = factor(c(sample(0:(length(p_control)-1), nC, TRUE, pi_control),
#                      sample(0:(length(p_control)-1), nT, TRUE, pi_treat)), ordered=T),
#         trt = factor(rep(c("C", "T"), c(nC, nT)), levels=c("C", "T"))
#       )
#     }
#     
#     # --- Futility ---
#     df_f <- gen_data(n_fut)
#     fit_f <- fit_logCOR(df_f, beta_NI)
#     idx_f <- (i-1)*3 + 1
#     if(!is.null(fit_f)) {
#       res_storage$logCOR[idx_f] <- fit_f$logCOR_hat
#       res_storage$pval[idx_f]   <- fit_f$p_val
#       if ((fit_f$logCOR_hat - beta_NI)/fit_f$se > z_fut) {
#         res_storage$stopped[idx_f] <- TRUE; next
#       }
#     }
#     
#     # --- IA ---
#     df_ia <- rbind(df_f, gen_data(n1 - n_fut))
#     fit_ia <- fit_logCOR(df_ia, beta_NI)
#     idx_ia <- (i-1)*3 + 2
#     if(!is.null(fit_ia)) {
#       res_storage$logCOR[idx_ia] <- fit_ia$logCOR_hat
#       res_storage$pval[idx_ia]   <- fit_ia$p_val
#       if ((fit_ia$logCOR_hat - beta_NI)/fit_ia$se <= zcrit1) {
#         res_storage$stopped[idx_ia] <- TRUE; next
#       }
#     }
#     
#     # --- Final ---
#     df_fin <- rbind(df_ia, gen_data(n_total - n1))
#     fit_fin <- fit_logCOR(df_fin, beta_NI)
#     idx_fin <- (i-1)*3 + 3
#     if(!is.null(fit_fin)) {
#       res_storage$logCOR[idx_fin] <- fit_fin$logCOR_hat
#       res_storage$pval[idx_fin]   <- fit_fin$p_val
#       if ((fit_fin$logCOR_hat - beta_NI)/fit_fin$se <= zcrit2) res_storage$stopped[idx_fin] <- TRUE
#     }
#   }
#   return(res_storage)
# }
#  
# 
# plot_pvals_vs_or <- function(sim_data, alpha1 = 0.0122, alpha2 = 0.0214, M_margin = 1.6) {
#   
#   # 1. Filter and prepare data (Remove Futility stage for clarity of Success thresholds)
#   plot_df <- sim_data %>%
#     filter(!is.na(logCOR), !is.na(pval)) %>%
#     mutate(OR = exp(logCOR))
#   
#   # 2. Calculate thresholds based on median SEs found in the data
#   # IA and Final have different SEs because they have different sample sizes
#   ia_se  <- median(abs(log(plot_df$OR/M_margin) / qnorm(plot_df$pval))[sim_data$stage == "IA"], na.rm=TRUE)
#   fin_se <- median(abs(log(plot_df$OR/M_margin) / qnorm(plot_df$pval))[sim_data$stage == "Final"], na.rm=TRUE)
#   
#   crit_OR_IA    <- exp(qnorm(alpha1) * ia_se + log(M_margin))
#   crit_OR_Final <- exp(qnorm(alpha2) * fin_se + log(M_margin))
#   
#   # Print to console
#   cat(sprintf("\n--- Decision Thresholds (One Panel) ---\n"))
#   cat(sprintf("At IA (N_ia):    Must see OR <= %.3f\n", crit_OR_IA))
#   cat(sprintf("At Final (N_tot): Must see OR <= %.3f\n", crit_OR_Final))
#   
#   # 3. Create the Plot
#   ggplot(plot_df, aes(x = OR, y = pval)) +
#     # Background points
#     geom_point(aes(color = stage, alpha = stopped), size = 1.5) +
#     
#     # Threshold Lines
#     geom_hline(yintercept = alpha1, linetype = "dashed", color = "blue", linewidth = 0.7) +
#     geom_hline(yintercept = alpha2, linetype = "dotdash", color = "darkgreen", linewidth = 0.7) +
#     geom_vline(xintercept = M_margin, linetype = "dotted", color = "red3", linewidth = 1) +
#     
#     # Left Label (IA)
#     annotate("label", x = -Inf, y = alpha1, 
#              label = paste("IA: OR <=", round(crit_OR_IA, 2)), 
#              fill = "blue", color = "white", fontface = "bold", size = 3.5, hjust = -0.1) +
#     
#     # Right Label (Final)
#     annotate("label", x = Inf, y = alpha2, 
#              label = paste("Final: OR <=", round(crit_OR_Final, 2)), 
#              fill = "darkgreen", color = "white", fontface = "bold", size = 3.5, hjust = 1.1) +
#     
#     # Formatting
#     scale_y_log10(breaks = c(1, 0.1, alpha2, alpha1, 0.001),
#                   labels = function(x) round(x, 4)) +
#     scale_alpha_manual(values = c(0.2, 0.8), guide = "none") + # Dim non-stopping points
#     scale_color_brewer(palette = "Set1", name = "Trial Stage") +
#     theme_minimal() +
#     theme(panel.grid.minor = element_blank(),
#           legend.position = "top") +
#     labs(title = "P-value vs. OR: Combined Success Thresholds",
#          subtitle = paste("Red Line = NI Margin (", M_margin, ")"),
#          x = "Estimated Odds Ratio (OR)",
#          y = "P-value (Log Scale)")
# }
# 
# # ============================================================
# # 5) Run
# # ============================================================
# results <- simulate_obf_ordinal(
#   COR_true = 1.1, COR_NI = 1.6, n_total = 600,
#   futility_frac = 0.5, info_frac = 0.8,
#   alpha1 = 0.0122, alpha2 = 0.0214, futility_p = 0.70,
#   p_control = parse_probs("0.04, 0.02, 0.45, 0.34, 0.15"),
#   nSims = 1000 # Reduced for speed
# )
# 
# plot_pvals_vs_or(results)


# rm(list=ls())
# 
# suppressPackageStartupMessages({
#   library(MASS)
#   library(dplyr)
#   library(ggplot2)
# })
# 
# # ============================================================
# # 1) Core helpers (Restored and Complete)
# # ============================================================
# ilogit <- function(z) 1 / (1 + exp(-z))
# 
# parse_probs <- function(txt) {
#   as.numeric(trimws(unlist(strsplit(txt, ","))))
# }
# 
# validate_probs <- function(p) {
#   if (any(is.na(p)))          return("Probabilities must be numeric.")
#   if (length(p) < 3)          return("Need at least 3 categories.")
#   if (abs(sum(p) - 1) > 1e-6) return(sprintf("Must sum to ~1 (got %.6f).", sum(p)))
#   NULL
# }
# 
# # The missing function:
# theta_from_control_pmf <- function(p_control) {
#   cum_control <- cumsum(p_control)
#   # Take the logit of cumulative probabilities except the last one (which is 1)
#   qlogis(cum_control[seq_len(length(p_control) - 1)])
# }
# 
# pmf_from_beta <- function(theta, beta, x) {
#   cdf <- ilogit(theta - beta * x)
#   cdf <- c(cdf, 1)
#   pmf <- diff(c(0, cdf))
#   pmf[pmf < 0] <- 0
#   pmf / sum(pmf)
# }
# 
# # ============================================================
# # 2) Fit Function (Collects P-values)
# # ============================================================
# fit_logCOR <- function(df, beta_NI = 0) {
#   fit <- try(MASS::polr(y ~ trt, data = df, Hess = TRUE), silent = TRUE)
#   if (inherits(fit, "try-error") || is.null(fit$coefficients)) return(NULL)
#   
#   trt_name <- names(fit$coefficients)[grepl("^trt", names(fit$coefficients))][1]
#   logCOR_hat <- as.numeric(fit$coefficients[[trt_name]])
#   se <- sqrt(vcov(fit)[trt_name, trt_name])
#   
#   # Calculate p-value relative to the NI Margin (beta_NI)
#   z_stat <- (logCOR_hat - beta_NI) / se
#   p_val  <- pnorm(z_stat) 
#   
#   list(logCOR_hat = logCOR_hat, se = se, p_val = p_val)
# }
# 
# # ============================================================
# # 3) Simulation Logic
# # ============================================================
# simulate_obf_ordinal <- function(COR_true, COR_NI, n_total, futility_frac, info_frac,
#                                  alpha1, alpha2, futility_p, p_control, 
#                                  seed = 1234, nSims = 500) {
#   set.seed(seed)
#   theta <- theta_from_control_pmf(p_control)
#   beta_true <- log(COR_true)
#   beta_NI   <- log(COR_NI)
#   pi_control <- p_control
#   pi_treat   <- pmf_from_beta(theta, beta_true, x = 1)
#   
#   res_storage <- data.frame(
#     sim = rep(1:nSims, each=3),
#     stage = rep(c("Futility", "IA", "Final"), nSims),
#     logCOR = NA_real_,
#     pval = NA_real_,
#     stopped = FALSE
#   )
#   
#   n_fut <- round(futility_frac * n_total)
#   n1    <- round(info_frac * n_total)
#   z_fut <- qnorm(futility_p)
#   zcrit1 <- qnorm(alpha1)
#   zcrit2 <- qnorm(alpha2)
#   
#   for (i in 1:nSims) {
#     gen_data <- function(n) {
#       nC <- floor(n/2); nT <- n - nC
#       data.frame(
#         y = factor(c(sample(0:(length(p_control)-1), nC, TRUE, pi_control),
#                      sample(0:(length(p_control)-1), nT, TRUE, pi_treat)), ordered=T),
#         trt = factor(rep(c("C", "T"), c(nC, nT)), levels=c("C", "T"))
#       )
#     }
#     
#     # --- Futility Look ---
#     df_f <- gen_data(n_fut); fit_f <- fit_logCOR(df_f, beta_NI)
#     idx_f <- (i-1)*3 + 1
#     if(!is.null(fit_f)) {
#       res_storage$logCOR[idx_f] <- fit_f$logCOR_hat
#       res_storage$pval[idx_f]   <- fit_f$p_val
#       if ((fit_f$logCOR_hat - beta_NI)/fit_f$se > z_fut) {
#         res_storage$stopped[idx_f] <- TRUE; next
#       }
#     }
#     
#     # --- IA Look ---
#     df_ia <- rbind(df_f, gen_data(n1 - n_fut)); fit_ia <- fit_logCOR(df_ia, beta_NI)
#     idx_ia <- (i-1)*3 + 2
#     if(!is.null(fit_ia)) {
#       res_storage$logCOR[idx_ia] <- fit_ia$logCOR_hat
#       res_storage$pval[idx_ia]   <- fit_ia$p_val
#       if ((fit_ia$logCOR_hat - beta_NI)/fit_ia$se <= zcrit1) {
#         res_storage$stopped[idx_ia] <- TRUE; next
#       }
#     }
#     
#     # --- Final Look ---
#     df_fin <- rbind(df_ia, gen_data(n_total - n1)); fit_fin <- fit_logCOR(df_fin, beta_NI)
#     idx_fin <- (i-1)*3 + 3
#     if(!is.null(fit_fin)) {
#       res_storage$logCOR[idx_fin] <- fit_fin$logCOR_hat
#       res_storage$pval[idx_fin]   <- fit_fin$p_val
#       if ((fit_fin$logCOR_hat - beta_NI)/fit_fin$se <= zcrit2) res_storage$stopped[idx_fin] <- TRUE
#     }
#   }
#   return(res_storage)
# }
# 
# # ============================================================
# # 4) Success Plot
# # ============================================================
# plot_pvals_vs_or <- function(sim_data, alpha1 = 0.0122, alpha2 = 0.0214, M_margin = 1.6) {
#   plot_df <- sim_data %>% filter(!is.na(logCOR), !is.na(pval)) %>% mutate(OR = exp(logCOR))
#   
#   # Calculate critical ORs based on the median SE found in simulation
#   ia_se  <- median(abs(log(plot_df$OR/M_margin) / qnorm(plot_df$pval))[plot_df$stage == "IA"], na.rm=TRUE)
#   fin_se <- median(abs(log(plot_df$OR/M_margin) / qnorm(plot_df$pval))[plot_df$stage == "Final"], na.rm=TRUE)
#   
#   crit_OR_IA    <- exp(qnorm(alpha1) * ia_se + log(M_margin))
#   crit_OR_Final <- exp(qnorm(alpha2) * fin_se + log(M_margin))
#   
#   ggplot(plot_df, aes(x = OR, y = pval)) +
#     geom_point(aes(color = stage, alpha = stopped), size = 1.5) +
#     geom_hline(yintercept = c(alpha1, alpha2), color = c("blue", "darkgreen"), linetype = "dashed") +
#     geom_vline(xintercept = M_margin, linetype = "dotted", color = "red3") +
#     annotate("label", x = -Inf, y = alpha1, label = paste("IA: OR <=", round(crit_OR_IA, 2)), fill = "blue", color = "white", hjust = -0.1) +
#     annotate("label", x = Inf, y = alpha2, label = paste("Final: OR <=", round(crit_OR_Final, 2)), fill = "darkgreen", color = "white", hjust = 1.1) +
#     scale_y_log10(breaks = c(1, 0.1, alpha2, alpha1, 0.001), labels = function(x) round(x, 4)) +
#     scale_alpha_manual(values = c(0.2, 0.8), guide = "none") +
#     theme_minimal() + 
#     labs(title = "Success Thresholds", x = "Estimated OR", y = "P-value (Log Scale)")
# }
# 
# # ============================================================
# # 5) Futility Plot
# # ============================================================
# plot_futility_vs_or <- function(sim_data, futility_p = 0.70, M_margin = 1.6) {
#   fut_df <- sim_data %>% filter(stage == "Futility", !is.na(logCOR)) %>% mutate(OR = exp(logCOR))
#   fut_se <- median(abs(log(fut_df$OR/M_margin) / qnorm(fut_df$pval)), na.rm=TRUE)
#   crit_OR_fut <- exp(qnorm(futility_p) * fut_se + log(M_margin))
#   
#   ggplot(fut_df, aes(x = OR, y = pval)) +
#     geom_point(aes(color = stopped), alpha = 0.6) +
#     geom_hline(yintercept = futility_p, linetype = "dashed", color = "firebrick") +
#     geom_vline(xintercept = M_margin, linetype = "dotted") +
#     annotate("label", x = Inf, y = futility_p, label = paste("Stop Futility: OR >=", round(crit_OR_fut, 2)), fill = "firebrick", color = "white", hjust = 1.1) +
#     scale_color_manual(values = c("gray60", "firebrick")) +
#     theme_minimal() + labs(title = "Futility Threshold", x = "Estimated OR", y = "P-value")
# }
# 
# # ============================================================
# # 6) Main Execution
# # ============================================================
# set.seed(1234) 
# results <- simulate_obf_ordinal(
#   COR_true = 1.1, COR_NI = 1.6, n_total = 600,
#   futility_frac = 0.5, info_frac = 0.8,
#   alpha1 = 0.0122, alpha2 = 0.0214, futility_p = 0.70,
#   p_control = parse_probs("0.04, 0.02, 0.45, 0.34, 0.15"),
#   nSims = 1000
# )
# 
# # Render both plots
# plot_pvals_vs_or(results)
# plot_futility_vs_or(results)

rm(list=ls())

suppressPackageStartupMessages({
  library(MASS)
  library(dplyr)
  library(ggplot2)
})

# ============================================================
# 1) Core helpers
# ============================================================
ilogit <- function(z) 1 / (1 + exp(-z))
parse_probs <- function(txt) as.numeric(trimws(unlist(strsplit(txt, ","))))

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
# 2) Fit Function
# ============================================================
fit_logCOR <- function(df, beta_NI = 0) {
  fit <- try(MASS::polr(y ~ trt, data = df, Hess = TRUE), silent = TRUE)
  if (inherits(fit, "try-error") || is.null(fit$coefficients)) return(NULL)
  trt_name <- names(fit$coefficients)[grepl("^trt", names(fit$coefficients))][1]
  logCOR_hat <- as.numeric(fit$coefficients[[trt_name]])
  se <- sqrt(vcov(fit)[trt_name, trt_name])
  z_stat <- (logCOR_hat - beta_NI) / se
  p_val  <- pnorm(z_stat) 
  list(logCOR_hat = logCOR_hat, se = se, p_val = p_val)
}

# ============================================================
# 3) Simulation (Fixed Data Generation)
# ============================================================
simulate_obf_ordinal <- function(COR_true, COR_NI, n_total, futility_frac, info_frac,
                                 alpha1, alpha2, futility_p, p_control, 
                                 seed = 1234, nSims = 500) {
  set.seed(seed)
  theta <- theta_from_control_pmf(p_control)
  K <- length(p_control) - 1 # Categories 0 to K
  beta_true <- log(COR_true); beta_NI <- log(COR_NI)
  
  pi_control <- p_control
  pi_treat   <- pmf_from_beta(theta, beta_true, x = 1)
  
  res_storage <- data.frame(
    sim = rep(1:nSims, each=3),
    stage = rep(c("Futility", "IA", "Final"), nSims),
    OR = NA_real_, pval = NA_real_, stopped = FALSE
  )
  
  n_fut <- round(futility_frac * n_total)
  n1    <- round(info_frac * n_total)
  z_fut <- qnorm(futility_p); zcrit1 <- qnorm(alpha1); zcrit2 <- qnorm(alpha2)
  
  for (i in 1:nSims) {
    # Helper to generate data correctly for both arms
    gen_data <- function(N) {
      nC <- floor(N / 2); nT <- N - nC
      yC <- sample(0:K, nC, TRUE, pi_control)
      yT <- sample(0:K, nT, TRUE, pi_treat)
      data.frame(
        y = factor(c(yC, yT), ordered = TRUE, levels = 0:K),
        trt = factor(rep(c("C", "T"), c(nC, nT)), levels = c("C", "T"))
      )
    }
    
    # --- Futility Look ---
    df_f <- gen_data(n_fut); fit_f <- fit_logCOR(df_f, beta_NI)
    idx_f <- (i-1)*3 + 1
    if(!is.null(fit_f)) {
      res_storage$OR[idx_f] <- exp(fit_f$logCOR_hat)
      res_storage$pval[idx_f] <- fit_f$p_val
      if ((fit_f$logCOR_hat - beta_NI)/fit_f$se > z_fut) {
        res_storage$stopped[idx_f] <- TRUE; next
      }
    }
    # --- IA Look ---
    df_ia <- rbind(df_f, gen_data(n1 - n_fut)); fit_ia <- fit_logCOR(df_ia, beta_NI)
    idx_ia <- (i-1)*3 + 2
    if(!is.null(fit_ia)) {
      res_storage$OR[idx_ia] <- exp(fit_ia$logCOR_hat)
      res_storage$pval[idx_ia] <- fit_ia$p_val
      if ((fit_ia$logCOR_hat - beta_NI)/fit_ia$se <= zcrit1) {
        res_storage$stopped[idx_ia] <- TRUE; next
      }
    }
    # --- Final Look ---
    df_fin <- rbind(df_ia, gen_data(n_total - n1)); fit_fin <- fit_logCOR(df_fin, beta_NI)
    idx_fin <- (i-1)*3 + 3
    if(!is.null(fit_fin)) {
      res_storage$OR[idx_fin] <- exp(fit_fin$logCOR_hat)
      res_storage$pval[idx_fin] <- fit_fin$p_val
      if ((fit_fin$logCOR_hat - beta_NI)/fit_fin$se <= zcrit2) res_storage$stopped[idx_fin] <- TRUE
    }
  }
  return(res_storage)
}

# ============================================================
# 4) Plotting
# ============================================================
 

plot_all_stages <- function(sim_data, alpha1, alpha2, futility_p, M_margin) {
  plot_df <- sim_data %>% 
    filter(!is.na(OR), !is.na(pval)) %>%
    mutate(
      # Detailed status for specific coloring
      status = case_when(
        stage == "IA"    & pval <= alpha1 ~ "IA Success",
        stage == "Final" & pval <= alpha2 ~ "Final Success",
        pval >= futility_p                ~ "Futility",
        TRUE                              ~ "Inconclusive"
      )
    )
  
  # Calculate OR thresholds
  calc_crit <- function(stg, p_val) {
    sub <- plot_df %>% filter(stage == stg)
    if(nrow(sub) == 0) return(NA)
    se <- median(abs(log(sub$OR/M_margin) / qnorm(sub$pval)), na.rm=TRUE)
    exp(qnorm(p_val) * se + log(M_margin))
  }
  
  c_fut <- calc_crit("Futility", futility_p)
  c_ia  <- calc_crit("IA", alpha1)
  c_fin <- calc_crit("Final", alpha2)
  
  # Label y-position
  y_lab <- 0.00025
  
  ggplot(plot_df, aes(x = OR, y = pval)) +
    # 1. Decision Points with specific color mapping
    geom_point(aes(color = status, alpha = stopped), size = 1.3) +
    
    # 2. Vertical Threshold Lines
    geom_vline(xintercept = c_ia, color = "dodgerblue3", linetype = "solid", linewidth = 0.8) +
    geom_vline(xintercept = c_fin, color = "forestgreen", linetype = "solid", linewidth = 0.8) +
    geom_vline(xintercept = c_fut, color = "firebrick", linetype = "solid", linewidth = 0.8) +
    
    # 3. Horizontal P-value Lines
    geom_hline(yintercept = c(alpha1, alpha2, futility_p), 
               color = c("dodgerblue3", "forestgreen", "firebrick"), 
               linetype = "dashed", alpha = 0.3) +
    
    # 4. NI Margin Reference
    geom_vline(xintercept = M_margin, linetype = "dotted", color = "black", linewidth = 1) +
    
    # 5. Adjusted Labels
    # IA Success: Moved to the LEFT (hjust = 1.1)
    annotate("label", x = c_ia, y = y_lab, 
             label = paste0("IA SUCCESS\nOR <= ", round(c_ia, 2), "\np=", alpha1), 
             fill = "dodgerblue3", color = "white", size = 2.8, fontface = "bold", hjust = 1.1) +
    
    # Final Success: Moved to the RIGHT (hjust = -0.1)
    annotate("label", x = c_fin, y = y_lab, 
             label = paste0("FINAL SUCCESS\nOR <= ", round(c_fin, 2), "\np=", alpha2), 
             fill = "forestgreen", color = "white", size = 2.8, fontface = "bold", hjust = -0.1) +
    
    # Futility: Centered or slightly offset
    annotate("label", x = c_fut, y = y_lab, 
             label = paste0("FUTILITY\nOR >= ", round(c_fut, 2), "\np=", futility_p), 
             fill = "firebrick", color = "white", size = 2.8, fontface = "bold", hjust = 0.5) +
    
    # Formatting
    scale_y_log10(limits = c(0.0001, 1.1), breaks = c(1, futility_p, 0.1, alpha2, alpha1, 0.001)) +
    scale_color_manual(values = c(
      "IA Success"    = "dodgerblue3", 
      "Final Success" = "forestgreen", 
      "Futility"      = "firebrick", 
      "Inconclusive"  = "black"
    )) +
    scale_alpha_manual(values = c(0.2, 0.9), guide = "none") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(), legend.position = "top") +
    labs(title = "Trial Decision Architecture: Stage-Specific Success",
         subtitle = "Blue = IA Success | Green = Final Success | Red = Futility | Black = Continued",
         x = "Estimated Odds Ratio (OR)", y = "P-value (Log Scale)")
}
# ============================================================
# 5) RUN
# ============================================================
res <- simulate_obf_ordinal(COR_true = 1., COR_NI = 1.6, n_total = 600, 
                            futility_frac = 0.5, info_frac = 0.8, 
                            alpha1 = 0.0122, alpha2 = 0.0214, futility_p = 0.7,
                            p_control = parse_probs("0.04, 0.02, 0.45, 0.34, 0.15"), 
                            nSims = 500, seed=123)

plot_all_stages(res, 0.0122, 0.0214, 0.7, 1.6)