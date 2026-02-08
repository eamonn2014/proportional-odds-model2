# Shiny app: Classic OBF NI (Ordinal PO) — Theory vs Monte-Carlo + Winner’s Curse
# Endpoint: ordinal 0..K (BEST → WORST), COR < 1 favors treatment
#
# NOTES:
# - Theoretical "borderline_COR" is a DESIGN THRESHOLD (depends on N, alpha, margin, control PMF),
#   so it will NOT change when COR_true changes. That is expected.
# - FIX: theoretical borderline calculation is consistent with simulation definition:
#       logCOR_hat = - polr(trt coefficient).
# - Boxplot proportions are DRAWN on the plot (Panel C), not in axis labels.
# - Interpretation text is added under each plot, with clarified language:
#   the final estimator is not biased; any shifts are due to selection/conditioning.

suppressPackageStartupMessages({
  library(shiny)
  library(MASS)   # polr
})

# ────────────────────────────────────────
# Core helpers
# ────────────────────────────────────────

ilogit <- function(z) 1 / (1 + exp(-z))

parse_probs <- function(txt) {
  as.numeric(trimws(unlist(strsplit(txt, ","))))
}

validate_probs <- function(p) {
  if (any(is.na(p)))               return("Probabilities must be numeric (comma-separated).")
  if (length(p) < 3)               return("Need at least 3 categories.")
  if (any(p <= 0))                 return("All probabilities must be > 0.")
  if (abs(sum(p) - 1) > 1e-6)      return(sprintf("Must sum to ~1 (got %.6f).", sum(p)))
  NULL
}

theta_from_control_pmf <- function(p_control) {
  cum_control <- cumsum(p_control)
  qlogis(cum_control[seq_len(length(p_control) - 1)])
}

pmf_from_beta <- function(theta, beta, x) {
  # PO model with cumulative logits: logit(P(Y<=k | x)) = theta_k - beta * x
  cdf <- ilogit(theta - beta * x)
  cdf <- c(cdf, 1)
  pmf <- diff(c(0, cdf))
  pmf[pmf < 0] <- 0
  pmf / sum(pmf)
}

build_weighted_data <- function(N_total, theta, beta) {
  n_c <- floor(N_total / 2)
  n_t <- N_total - n_c
  pmf_c <- pmf_from_beta(theta, beta, 0)
  pmf_t <- pmf_from_beta(theta, beta, 1)
  y_levels <- seq_along(pmf_c) - 1
  data.frame(
    y   = factor(rep(y_levels, 2), ordered = TRUE, levels = y_levels),
    trt = rep(c(0, 1), each = length(y_levels)),
    w   = c(n_c * pmf_c, n_t * pmf_t)
  ) |> transform(w = pmax(w, 1e-6))
}

se_beta_polr <- function(N_total, theta, beta) {
  df <- build_weighted_data(N_total, theta, beta)
  fit <- suppressWarnings(polr(y ~ trt, data = df, weights = w, Hess = TRUE))
  sqrt(vcov(fit)["trt", "trt"])
}

# ────────────────────────────────────────
# Theoretical borderline (one look)
# ────────────────────────────────────────
#
# GOAL (unchanged): Find the effect (COR) such that the asymptotic Z-statistic hits the boundary:
#   Z = (logCOR - log(M_margin))/SE_logCOR  = z_alpha
#
# IMPORTANT: In this app, simulation uses:
#   logCOR_hat = - (polr coefficient for treatment)
# so logCOR = -beta_polr
#
# FIX: Solve on logCOR scale, but compute SE by passing beta_polr = -logCOR into se_beta_polr().

borderline_beta_one_look <- function(
    N_total, alpha, M_margin, theta,
    beta_lower = -6,
    beta_upper = log(M_margin) - 1e-8
) {
  logM <- log(M_margin)
  z_a  <- qnorm(alpha)
  
  logCOR_lower <- beta_lower
  logCOR_upper <- beta_upper
  
  f <- function(logCOR) {
    beta_polr <- -logCOR
    se_logCOR <- se_beta_polr(N_total, theta, beta_polr)
    (logCOR - logM) / se_logCOR - z_a
  }
  
  fL <- f(logCOR_lower)
  fU <- f(logCOR_upper)
  
  if (anyNA(c(fL, fU))) return(NA_real_)
  if (!is.finite(fL) || !is.finite(fU)) return(NA_real_)
  
  if (fL > 0) return(NA_real_)
  if (fU < 0) return(logCOR_upper)
  
  uniroot(f, lower = logCOR_lower, upper = logCOR_upper)$root
}

# ────────────────────────────────────────
# Monte-Carlo simulation (Classic OBF)
# ────────────────────────────────────────

fit_logCOR <- function(df) {
  fit <- try(MASS::polr(y ~ trt, data = df, Hess = TRUE), silent = TRUE)
  if (inherits(fit, "try-error") || is.null(fit$coefficients)) return(NULL)
  
  coef_names <- names(fit$coefficients)
  trt_name <- coef_names[grepl("^trt", coef_names)][1]
  if (is.na(trt_name)) return(NULL)
  
  # Simulation convention: logCOR_hat = - polr(trt coefficient)
  logCOR_hat <- - as.numeric(fit$coefficients[[trt_name]])
  
  v <- vcov(fit)
  if (!(trt_name %in% rownames(v))) return(NULL)
  
  se <- sqrt(v[trt_name, trt_name])
  if (!is.finite(se) || se < 1e-8) return(NULL)
  
  list(logCOR_hat = logCOR_hat, se = se)
}

simulate_obf_ordinal <- function(
    COR_true, COR_NI,
    n1_total, n_total,
    alpha1, alpha2,
    p_control, seed = 1, nSims = 2000
) {
  set.seed(seed)
  eps <- 1e-12
  
  theta <- theta_from_control_pmf(p_control)
  K <- length(p_control) - 1
  
  beta_true <- log(COR_true)
  beta_NI   <- log(COR_NI)
  
  pi_control <- p_control
  pi_treat   <- diff(c(0, ilogit(theta - beta_true), 1))
  
  split_n <- function(N) {
    nC <- floor(N / 2)
    list(nC = nC, nT = N - nC)
  }
  
  Z1_all <- OR1_all <- Z2_all <- OR2_all <- rep(NA_real_, nSims)
  stop_stage <- rep(0L, nSims)
  OR_at_stop <- numeric()
  stop_stage_success <- integer()
  
  stop1 <- stop2 <- 0L
  zcrit1 <- qnorm(alpha1)
  zcrit2 <- qnorm(alpha2)
  
  for (i in seq_len(nSims)) {
    # ---- Stage 1 ----
    s1 <- split_n(n1_total)
    yC1 <- sample(0:K, s1$nC, TRUE, pi_control)
    yT1 <- sample(0:K, s1$nT, TRUE, pi_treat)
    
    df1 <- data.frame(
      y   = factor(c(yC1, yT1), ordered = TRUE, levels = 0:K),
      trt = factor(rep(c("C", "T"), c(s1$nC, s1$nT)), levels = c("C", "T"))
    )
    
    res1 <- fit_logCOR(df1)
    p1 <- 1
    
    if (!is.null(res1)) {
      Z1_all[i]  <- (res1$logCOR_hat - beta_NI) / res1$se
      OR1_all[i] <- exp(res1$logCOR_hat)
      p1 <- pnorm(Z1_all[i])
      p1 <- min(max(p1, eps), 1 - eps)
      
      if (p1 <= alpha1) {
        stop1 <- stop1 + 1L
        stop_stage[i] <- 1L
        if (is.finite(OR1_all[i])) {
          OR_at_stop <- c(OR_at_stop, OR1_all[i])
          stop_stage_success <- c(stop_stage_success, 1L)
        }
        next
      }
    }
    
    # ---- Stage 2 ----
    n2_total <- n_total - n1_total
    s2 <- split_n(n2_total)
    yC2 <- sample(0:K, s2$nC, TRUE, pi_control)
    yT2 <- sample(0:K, s2$nT, TRUE, pi_treat)
    
    df2 <- data.frame(
      y   = factor(c(yC2, yT2), ordered = TRUE, levels = 0:K),
      trt = factor(rep(c("C", "T"), c(s2$nC, s2$nT)), levels = c("C", "T"))
    )
    
    df12 <- rbind(df1, df2)
    res12 <- fit_logCOR(df12)
    
    if (!is.null(res12)) {
      z2 <- (res12$logCOR_hat - beta_NI) / res12$se
      p2 <- pnorm(z2)
      p2 <- min(max(p2, eps), 1 - eps)
      OR2 <- exp(res12$logCOR_hat)
      
      Z2_all[i]  <- z2
      OR2_all[i] <- OR2
      
      if (p2 <= alpha2) {
        stop2 <- stop2 + 1L
        stop_stage[i] <- 2L
        if (is.finite(OR2)) {
          OR_at_stop <- c(OR_at_stop, OR2)
          stop_stage_success <- c(stop_stage_success, 2L)
        }
      }
    }
  }
  
  list(
    stop1 = stop1, stop2 = stop2, nSims = nSims,
    OR_at_stop = OR_at_stop, stop_stage_success = stop_stage_success,
    Z1_all = Z1_all, OR1_all = OR1_all,
    Z2_all = Z2_all, OR2_all = OR2_all,
    zcrit1 = zcrit1, zcrit2 = zcrit2
  )
}

summ_or <- function(x) {
  if (length(x) == 0) return(c(N=0L, mean=NA_real_, median=NA_real_, q025=NA_real_, q975=NA_real_))
  c(
    N      = length(x),
    mean   = mean(x),
    median = median(x),
    q025   = unname(quantile(x, 0.025)),
    q975   = unname(quantile(x, 0.975))
  )
}

# ────────────────────────────────────────
# UI
# ────────────────────────────────────────

ui <- fluidPage(
  titlePanel("OBF Non-Inferiority – Ordinal PO: Theory vs Monte-Carlo + Winner’s Curse"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("N_total", "Total N (final)", 600, min = 50, step = 10),
      numericInput("info_frac", "Info fraction at IA", 0.80, 0.1, 0.95, 0.05),
      numericInput("alpha1", "α₁ (one-sided IA)", 0.0122, 1e-6, 0.2, 0.001),
      numericInput("alpha2", "α₂ (one-sided final)", 0.0214, 1e-6, 0.2, 0.001),
      
      numericInput("M_margin", "NI margin (COR)", 1.6, 1.01, step = 0.05),
      numericInput("COR_true", "True COR (simulation)", 1.0, 0.2, step = 0.05),
      
      textInput("p_control", "Control probs (best→worst, comma sep)",
                "0.04, 0.02, 0.45, 0.34, 0.15"),
      
      tags$hr(),
      checkboxInput("run_grid", "Run grid over proportions", FALSE),
      numericInput("prop_start", "Start prop", 0.30, 0.05, 1, 0.05),
      numericInput("prop_end",   "End prop",   1.00, 0.05, 1, 0.05),
      numericInput("prop_step",  "Step",       0.05, 0.01, 0.5, 0.01),
      
      tags$hr(),
      numericInput("nSims", "# Sims per N", 500, 100, step = 500),
      numericInput("seed", "Seed", 1234, 1),
      
      tags$hr(),
      actionButton("run", "Run Simulation", class = "btn-primary", width = "100%"),
      tags$br(), tags$br(),
      helpText("Early stopping selects overly optimistic estimates (winner’s curse). Reaching final selects the opposite selection effect.")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Interim (Z1 & OR1)",
                 plotOutput("plot_compare", height = "520px"),
                 tags$hr(),
                 tags$div(style="font-size: 15px; line-height: 1.4;",
                          tags$b("Interpretation:"),
                          tags$br(), tags$br(),
                          tags$ul(
                            tags$li("The left panel shows interim (early) results across many simulated trials. Further left means the treatment looks more favourable relative to the non-inferiority margin."),
                            tags$li("The red dashed line is the pre-set early stopping rule for success. Trials left of this line would stop early and declare non-inferiority at interim."),
                            tags$li("The right panel shows the estimated treatment effect (OR) at interim. Blue is all trials; red is only those that stopped early."),
                            tags$li("Because early stopping requires an unusually strong early signal, the trials that stop early tend to look more favourable than average. This is the expected ‘winner’s curse’: the subgroup that stops early often has an over-optimistic estimate compared with the true underlying effect.")
                          )
                 )
        ),
        
        tabPanel("Final (Z2 & OR2)",
                 plotOutput("plot_final", height = "520px"),
                 tags$hr(),
                 tags$div(style="font-size: 15px; line-height: 1.4;",
                          tags$b("Interpretation:"),
                          tags$br(), tags$br(),
                          tags$ul(
                            tags$li("These results are for trials that did not stop early at interim and therefore continued to the planned final analysis."),
                            tags$li("The red dashed line is the final decision boundary for non-inferiority at the end of the study."),
                            tags$li("On the right, blue shows the estimated treatment effect among trials that reached final; red shows the subset that succeeded at final."),
                            tags$li(tags$b("Important clarification: "),
                                    "the final statistical method is not ‘biased’. Any shift you see is because we are looking at a selected subgroup (only trials that did not stop early). Trials that continue tend to be those without a strong early signal, so their distribution can look less favourable on average — but the final estimator remains valid.")
                          )
                 )
        ),
        
        tabPanel("Selection Bias Overview",
                 plotOutput("plot_overview", height = "560px"),
                 tags$hr(),
                 tags$div(style="font-size: 15px; line-height: 1.4;",
                          tags$b("Interpretation (how stopping rules affect what you ‘see’):"),
                          tags$br(), tags$br(),
                          tags$ul(
                            tags$li("Panel A: Interim estimates. Red shows trials that stopped early. These are selected because they looked unusually good early on, so this subgroup tends to show more favourable estimates (winner’s curse)."),
                            tags$li("Panel B: Final estimates. Blue shows trials that reached the final analysis (they did not meet the early success rule). Because this is a selected subgroup, it can look less favourable on average — this is a selection/conditioning effect, not a flaw in the final statistical method."),
                            tags$li("Panel C: Boxplots summarise the same selection effects. The percentages printed next to each box show the proportion of all trials in each category (e.g., stopped early vs reached final)."),
                            tags$li(tags$b("Clinical takeaway: "),
                                    "early-stopped results can look overly optimistic because of selection; trials that continue to final are a different (selected) subset. This does not mean the final OR is biased — it reflects the trial pathway and stopping rules.")
                          )
                 )
        ),
        
        tabPanel("Tables",
                 h4("Theoretical borderline COR"),
                 tableOutput("tbl_theory"),
                 tags$hr(),
                 h4("Observed OR at stopping"),
                 tableOutput("tbl_sim")
        )
      )
    )
  )
)

# ────────────────────────────────────────
# Server
# ────────────────────────────────────────

server <- function(input, output, session) {
  
  invalid_probs_msg <- reactive({
    p <- parse_probs(input$p_control)
    validate_probs(p)
  })
  
  results <- eventReactive(input$run, {
    p <- parse_probs(input$p_control)
    if (!is.null(validate_probs(p))) return(NULL)
    
    props <- if (input$run_grid) seq(input$prop_start, input$prop_end, by = input$prop_step) else 1
    N_eff <- round(input$N_total * props)
    N1    <- round(input$info_frac * N_eff)
    
    theta <- theta_from_control_pmf(p)
    
    theory <- data.frame(
      prop = props, N_eff = N_eff, N_IA = N1,
      borderline_COR_IA = NA_real_, borderline_COR_final = NA_real_
    )
    
    sim_out <- data.frame(
      prop = rep(props, each = 2),
      N_eff = rep(N_eff, each = 2),
      stage = rep(c("Stop @ IA", "Stop @ Final"), length(props)),
      N = NA_real_, mean = NA_real_, median = NA_real_, q025 = NA_real_, q975 = NA_real_,
      Pr_stop = NA_real_
    )
    
    sim_for_plot <- NULL
    N_eff_plot <- max(N_eff)
    N1_plot    <- N1[N_eff == N_eff_plot][1]
    
    withProgress(message = "Running simulations...", value = 0, {
      for (i in seq_along(N_eff)) {
        incProgress(1 / length(N_eff), detail = paste("N =", N_eff[i]))
        
        # Theory: borderline effect hitting boundary (DESIGN threshold)
        b1 <- borderline_beta_one_look(N1[i],    input$alpha1, input$M_margin, theta)
        b2 <- borderline_beta_one_look(N_eff[i], input$alpha2, input$M_margin, theta)
        
        theory$borderline_COR_IA[i]    <- if (is.na(b1)) NA_real_ else exp(b1)
        theory$borderline_COR_final[i] <- if (is.na(b2)) NA_real_ else exp(b2)
        
        # Simulation
        sim <- simulate_obf_ordinal(
          COR_true  = input$COR_true,
          COR_NI    = input$M_margin,
          n1_total  = N1[i],
          n_total   = N_eff[i],
          alpha1    = input$alpha1,
          alpha2    = input$alpha2,
          p_control = p,
          seed      = input$seed + i,
          nSims     = input$nSims
        )
        
        if (N_eff[i] == N_eff_plot) {
          sim_for_plot <- sim
          sim_for_plot$N_eff_plot <- N_eff_plot
          sim_for_plot$N1_plot    <- N1_plot
        }
        
        OR <- sim$OR_at_stop
        st <- sim$stop_stage_success
        
        s1 <- summ_or(OR[st == 1])
        s2 <- summ_or(OR[st == 2])
        
        idx1 <- (i-1)*2 + 1
        idx2 <- (i-1)*2 + 2
        
        sim_out[idx1, c("N","mean","median","q025","q975")] <- as.numeric(s1)
        sim_out[idx2, c("N","mean","median","q025","q975")] <- as.numeric(s2)
        
        sim_out$Pr_stop[idx1] <- sim$stop1 / sim$nSims
        sim_out$Pr_stop[idx2] <- sim$stop2 / sim$nSims
      }
    })
    
    list(
      theory = theory,
      sim = sim_out,
      sim_for_plot = sim_for_plot,
      M_margin = input$M_margin,
      COR_true = input$COR_true
    )
  })
  
  # ── Plot: Interim ────────────────────────────────────────────────────────
  
  output$plot_compare <- renderPlot({
    msg <- invalid_probs_msg()
    if (!is.null(msg)) {
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
      text(0.5, 0.5, msg, cex = 1.4, col = "darkred")
      return()
    }
    
    sim <- req(results())$sim_for_plot
    
    ok <- is.finite(sim$Z1_all) & is.finite(sim$OR1_all)
    Z1  <- sim$Z1_all[ok]
    OR1 <- sim$OR1_all[ok]
    
    cond <- Z1 <= sim$zcrit1
    OR1_cond <- OR1[cond]
    has_cond <- sum(cond) > 10
    
    op <- par(mfrow = c(1,2), mar = c(4,4,3,1))
    
    hist(Z1, breaks = 50, freq = FALSE, col = "gray90", border = "white",
         main = sprintf("Interim Z1 (N1 = %d of %d)", sim$N1_plot, sim$N_eff_plot),
         xlab = "Z1 = (loĝCOR – log M) / SE")
    abline(v = sim$zcrit1, col = "red3", lwd = 2.5, lty = 2)
    if (length(Z1) > 5) lines(density(Z1), lwd = 2)
    
    if (length(OR1) < 5) {
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
      text(0.5, 0.5, "Too few valid OR1 values", cex = 1.3, col = "gray50")
    } else {
      d_all <- density(OR1, na.rm = TRUE)
      
      if (has_cond) {
        d_cond <- density(OR1_cond, na.rm = TRUE)
        ymax <- max(d_all$y, d_cond$y, na.rm = TRUE) * 1.1
      } else {
        d_cond <- NULL
        ymax <- max(d_all$y, na.rm = TRUE) * 1.1
      }
      
      xlim <- unname(quantile(OR1, c(0.005, 0.995), na.rm = TRUE))
      xlim <- c(xlim[1] * 0.95, xlim[2] * 1.05)
      
      plot(d_all, lwd = 2, col = "#1f77b4",
           xlim = xlim, ylim = c(0, ymax),
           main = "Interim OR1: all vs stop@IA",
           xlab = "OR1", ylab = "Density")
      
      if (has_cond && !is.null(d_cond)) lines(d_cond, lwd = 2, col = "#d62728")
      
      legend("topright", inset = 0.02,
             legend = c(sprintf("All (n=%d)", length(OR1)),
                        if (has_cond) sprintf("Stopped @ IA (n=%d)", sum(cond))),
             col = c("#1f77b4", if (has_cond) "#d62728"),
             lwd = 2, bty = "n")
    }
    
    par(op)
  })
  
  # ── Plot: Final ──────────────────────────────────────────────────────────
  
  output$plot_final <- renderPlot({
    msg <- invalid_probs_msg()
    if (!is.null(msg)) {
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
      text(0.5, 0.5, msg, cex = 1.4, col = "darkred")
      return()
    }
    
    sim <- req(results())$sim_for_plot
    
    ok <- is.finite(sim$Z2_all) & is.finite(sim$OR2_all)
    Z2  <- sim$Z2_all[ok]
    OR2 <- sim$OR2_all[ok]
    
    cond <- Z2 <= sim$zcrit2
    OR2_cond <- OR2[cond]
    has_cond <- sum(cond) > 10
    
    n_final <- length(Z2)
    n_stop  <- sum(cond)
    
    op <- par(mfrow = c(1,2), mar = c(4,4,3,1))
    
    hist(Z2, breaks = 50, freq = FALSE, col = "gray90", border = "white",
         main = sprintf("Final Z2 (reached: %d/%d; stop: %d)", n_final, sim$nSims, n_stop),
         xlab = "Z2 (pooled)")
    abline(v = sim$zcrit2, col = "red3", lwd = 2.5, lty = 2)
    if (length(Z2) > 5) lines(density(Z2), lwd = 2)
    
    if (length(OR2) < 5) {
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
      text(0.5, 0.5, "Too few valid OR2 values", cex = 1.3, col = "gray50")
    } else {
      d_all <- density(OR2, na.rm = TRUE)
      
      if (has_cond) {
        d_cond <- density(OR2_cond, na.rm = TRUE)
        ymax <- max(d_all$y, d_cond$y, na.rm = TRUE) * 1.1
      } else {
        d_cond <- NULL
        ymax <- max(d_all$y, na.rm = TRUE) * 1.1
      }
      
      xlim <- unname(quantile(OR2, c(0.005, 0.995), na.rm = TRUE))
      xlim <- c(xlim[1] * 0.95, xlim[2] * 1.05)
      
      plot(d_all, lwd = 2, col = "#1f77b4",
           xlim = xlim, ylim = c(0, ymax),
           main = sprintf("Final OR2 (reached: %d/%d)", n_final, sim$nSims),
           xlab = "OR2", ylab = "Density")
      
      if (has_cond && !is.null(d_cond)) lines(d_cond, lwd = 2, col = "#d62728")
      
      legend("topright", inset = 0.02,
             legend = c(sprintf("Reached final (n=%d)", n_final),
                        if (has_cond) sprintf("Stopped @ final (n=%d)", n_stop)),
             col = c("#1f77b4", if (has_cond) "#d62728"),
             lwd = 2, bty = "n")
    }
    
    par(op)
  })
  
  # ── Plot: Overview ───────────────────────────────────────────────────────
  
  output$plot_overview <- renderPlot({
    msg <- invalid_probs_msg()
    if (!is.null(msg)) {
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
      text(0.5, 0.5, msg, cex = 1.4, col = "darkred")
      return()
    }
    
    res <- req(results())
    sim <- req(res$sim_for_plot)
    
    log_true <- log(res$COR_true)
    log_M    <- log(res$M_margin)
    
    ok1 <- is.finite(sim$OR1_all) & is.finite(sim$Z1_all)
    L1_all  <- log(sim$OR1_all[ok1])
    L1_stop <- log(sim$OR1_all[ok1 & (sim$Z1_all <= sim$zcrit1)])
    has_stop1 <- length(L1_stop) > 5
    
    ok2 <- is.finite(sim$OR2_all) & is.finite(sim$Z2_all)
    L2_reach <- log(sim$OR2_all[ok2])
    L2_stop  <- log(sim$OR2_all[ok2 & (sim$Z2_all <= sim$zcrit2)])
    has_stop2 <- length(L2_stop) > 5
    
    op <- par(mfrow = c(1,3), mar = c(4,4,3,1))
    
    # A) Interim log(COR)
    if (length(L1_all) < 5) {
      plot(1, type = "n", axes = FALSE, main = "A) Interim log(COR)")
      text(0.5, 0.5, "Too few points", cex = 1.2, col = "gray50")
    } else {
      d_all <- density(L1_all, na.rm = TRUE)
      if (has_stop1) {
        d_stop <- density(L1_stop, na.rm = TRUE)
        ymax <- max(d_all$y, d_stop$y) * 1.1
      } else {
        d_stop <- NULL
        ymax <- max(d_all$y) * 1.1
      }
      
      plot(d_all, lwd = 2, col = "#1f77b4", ylim = c(0, ymax),
           main = "A) Interim log(COR)", xlab = "log(OR1)", ylab = "Density")
      if (has_stop1) lines(d_stop, lwd = 2, col = "#d62728")
      abline(v = c(log_true, log_M), lty = c(2,3), col = c("black","gray50"), lwd = 2)
    }
    
    # B) Final log(COR)
    if (length(L2_reach) < 5) {
      plot(1, type = "n", axes = FALSE, main = "B) Final log(COR)")
      text(0.5, 0.5, "Too few points", cex = 1.2, col = "gray50")
    } else {
      d_all <- density(L2_reach, na.rm = TRUE)
      if (has_stop2) {
        d_stop <- density(L2_stop, na.rm = TRUE)
        ymax <- max(d_all$y, d_stop$y) * 1.1
      } else {
        d_stop <- NULL
        ymax <- max(d_all$y) * 1.1
      }
      
      plot(d_all, lwd = 2, col = "#1f77b4", ylim = c(0, ymax),
           main = "B) Final log(COR)", xlab = "log(OR2)", ylab = "Density")
      if (has_stop2) lines(d_stop, lwd = 2, col = "#d62728")
      abline(v = c(log_true, log_M), lty = c(2,3), col = c("black","gray50"), lwd = 2)
    }
    
    # C) Boxplot summary + proportions DRAWN ON THE PLOT (no margin changes)
    groups <- list(
      "IA all"        = L1_all,
      "IA stopped"    = L1_stop,
      "Final reached" = L2_reach,
      "Final stopped" = L2_stop
    )
    groups <- groups[sapply(groups, function(v) length(v) > 1)]
    
    if (length(groups) > 0) {
      boxplot(groups, horizontal = TRUE, las = 1, col = "gray92",
              main = "C) Selection effects (log scale)",
              xlab = "log(COR)")
      abline(v = c(log_true, log_M), lty = c(2,3), col = c("black","gray50"), lwd = 2)
      
      # proportions computed from counters (proportion of ALL trials)
      p_stop1 <- sim$stop1 / sim$nSims
      p_reach <- (sim$nSims - sim$stop1) / sim$nSims
      p_stop2 <- sim$stop2 / sim$nSims
      
      # Map to displayed groups (in case any were dropped due to too few points)
      props <- c(
        "IA all"        = 1,
        "IA stopped"    = p_stop1,
        "Final reached" = p_reach,
        "Final stopped" = p_stop2
      )[names(groups)]
      
      # Place percents near the right edge of panel C, aligned to each box
      usr <- par("usr")
      x_text <- usr[2] - 0.03 * (usr[2] - usr[1])
      y_pos  <- seq_along(groups)
      
      text(x = x_text, y = max(y_pos) + 0.6, labels = "% trials", cex = 0.85, adj = 1, col = "gray30")
      text(x = x_text, y = y_pos, labels = sprintf("%.1f%%", 100 * props), cex = 0.9, adj = 1, col = "black")
    } else {
      plot(1, type = "n", axes = FALSE, main = "C) Selection effects")
      text(0.5, 0.5, "No valid data", cex = 1.2, col = "gray50")
    }
    
    par(op)
  })
  
  # ── Tables ───────────────────────────────────────────────────────────────
  
  output$tbl_theory <- renderTable({
    msg <- invalid_probs_msg()
    if (!is.null(msg)) {
      data.frame(Note = msg)
    } else {
      req(results())$theory |>
        transform(
          borderline_COR_IA    = round(borderline_COR_IA, 3),
          borderline_COR_final = round(borderline_COR_final, 3)
        )
    }
  }, digits = 3)
  
  output$tbl_sim <- renderTable({
    msg <- invalid_probs_msg()
    if (!is.null(msg)) {
      data.frame(Note = msg)
    } else {
      req(results())$sim |>
        transform(
          mean    = round(mean,   3),
          median  = round(median, 3),
          q025    = round(q025,   3),
          q975    = round(q975,   3),
          Pr_stop = round(Pr_stop, 3)
        )
    }
  }, digits = 3)
}

shinyApp(ui, server)
