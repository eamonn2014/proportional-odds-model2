suppressPackageStartupMessages({
  library(shiny)
  library(MASS)     # polr
  library(dplyr)    # table helpers
  library(tidyr)    # unnest_wider
})

# ────────────────────────────────────────
# Core helpers (unchanged)
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

fit_logCOR <- function(df) {
  fit <- try(MASS::polr(y ~ trt, data = df, Hess = TRUE), silent = TRUE)
  if (inherits(fit, "try-error") || is.null(fit$coefficients)) return(NULL)
  
  coef_names <- names(fit$coefficients)
  trt_name <- coef_names[grepl("^trt", coef_names)][1]
  if (is.na(trt_name)) return(NULL)
  
  logCOR_hat <- - as.numeric(fit$coefficients[[trt_name]])
  
  v <- vcov(fit)
  if (!(trt_name %in% rownames(v))) return(NULL)
  
  se <- sqrt(v[trt_name, trt_name])
  if (!is.finite(se) || se < 1e-8) return(NULL)
  
  list(logCOR_hat = logCOR_hat, se = se)
}

# ────────────────────────────────────────
# Simulation
# ────────────────────────────────────────

simulate_obf_ordinal <- function(
    COR_true, COR_NI,
    n_total,
    futility_frac, info_frac,
    alpha1, alpha2, futility_p,
    p_control, seed = 1234, nSims = 500
) {
  set.seed(seed)
  
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
  
  Z_fut_all  <- rep(NA_real_, nSims)
  OR_fut_all <- rep(NA_real_, nSims)
  Z1_all     <- rep(NA_real_, nSims)
  OR1_all    <- rep(NA_real_, nSims)
  Z2_all     <- rep(NA_real_, nSims)
  OR2_all    <- rep(NA_real_, nSims)
  
  stop_fut   <- logical(nSims)
  stop_ia    <- logical(nSims)
  stop_final <- logical(nSims)
  
  z_fut  <- qnorm(futility_p)
  zcrit1 <- qnorm(alpha1)
  zcrit2 <- qnorm(alpha2)
  
  n_fut <- round(futility_frac * n_total)
  n1    <- round(info_frac  * n_total)
  
  for (i in seq_len(nSims)) {
    incProgress(1/nSims, detail = sprintf("Sim %d/%d", i, nSims))
    
    s_f <- split_n(n_fut)
    yCf <- sample(0:K, s_f$nC, TRUE, pi_control)
    yTf <- sample(0:K, s_f$nT, TRUE, pi_treat)
    df_f <- data.frame(
      y   = factor(c(yCf, yTf), ordered = TRUE, levels = 0:K),
      trt = factor(rep(c("C", "T"), c(s_f$nC, s_f$nT)), levels = c("C", "T"))
    )
    
    res_f <- fit_logCOR(df_f)
    if (!is.null(res_f)) {
      Z_fut_all[i]  <- (res_f$logCOR_hat - beta_NI) / res_f$se
      OR_fut_all[i] <- exp(res_f$logCOR_hat)
      if (Z_fut_all[i] > z_fut) {
        stop_fut[i] <- TRUE
        next
      }
    }
    
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
      Z1_all[i]  <- (res1$logCOR_hat - beta_NI) / res1$se
      OR1_all[i] <- exp(res1$logCOR_hat)
      if (Z1_all[i] <= zcrit1) {
        stop_ia[i] <- TRUE
        next
      }
    }
    
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
      Z2_all[i]  <- (res2$logCOR_hat - beta_NI) / res2$se
      OR2_all[i] <- exp(res2$logCOR_hat)
      if (Z2_all[i] <= zcrit2) {
        stop_final[i] <- TRUE
      }
    }
  }
  
  list(
    Z_fut_all = Z_fut_all, OR_fut_all = OR_fut_all,
    Z1_all = Z1_all, OR1_all = OR1_all,
    Z2_all = Z2_all, OR2_all = OR2_all,
    stop_fut = stop_fut, stop_ia = stop_ia, stop_final = stop_final,
    nSims = nSims,
    z_fut = z_fut, zcrit1 = zcrit1, zcrit2 = zcrit2
  )
}

safe_summ_or <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n == 0) {
    return(c(N = 0L, Mean = NA_real_, Median = NA_real_, `2.5%` = NA_real_, `97.5%` = NA_real_))
  }
  q <- quantile(x, probs = c(0.025, 0.975), names = FALSE)
  c(
    N      = n,
    Mean   = mean(x),
    Median = median(x),
    `2.5%` = q[1],
    `97.5%`= q[2]
  )
}

# ────────────────────────────────────────
# UI
# ────────────────────────────────────────

ui <- fluidPage(
  titlePanel("OBF Non-Inferiority – Ordinal PO: Theory vs Monte-Carlo + Winner’s Curse"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("N_total", "Total N (final)", 600, min = 200, step = 50),
      numericInput("futility_frac", "Futility IA info fraction", 0.50, 0.1, 0.95, 0.05),
      numericInput("info_frac", "Efficacy IA info fraction", 0.80, 0.1, 0.95, 0.05),
      numericInput("futility_p", "Futility p-value threshold", 0.70, 0.1, 0.95, 0.05),
      numericInput("alpha1", "α₁ (one-sided IA success)", 0.0122, 1e-6, 0.05, 0.001),
      numericInput("alpha2", "α₂ (one-sided final)", 0.0214, 1e-6, 0.05, 0.001),
      
      numericInput("M_margin", "NI margin (COR)", 1.6, 1.01, step = 0.05),
      numericInput("COR_true", "True COR (simulation)", 1.0, 0.2, step = 0.05),
      
      textInput("p_control", "Control probs (best→worst, comma sep)",
                "0.04, 0.02, 0.45, 0.34, 0.15"),
      
      tags$hr(),
      numericInput("nSims", "# Simulations", 500, 100, step = 100),
      numericInput("seed", "Seed", 1234, 1),
      
      tags$hr(),
      actionButton("run", "Run Simulation", class = "btn-primary", width = "100%")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Futility (Z_f & OR_f)",
                 plotOutput("plot_futility", height = "520px"),
                 tags$hr(),
                 uiOutput("interp_futility")
        ),
        
        tabPanel("Interim (Z1 & OR1)",
                 plotOutput("plot_interim", height = "520px"),
                 tags$hr(),
                 uiOutput("interp_interim")
        ),
        
        tabPanel("Final (Z2 & OR2)",
                 plotOutput("plot_final", height = "520px"),
                 tags$hr(),
                 uiOutput("interp_final")
        ),
        
        tabPanel("Selection Bias Overview",
                 plotOutput("plot_overview", height = "560px"),
                 tags$hr(),
                 uiOutput("interp_overview")
        ),
        
        tabPanel("Selection Effects Boxplot",
                 plotOutput("plot_boxplot", height = "560px"),
                 tags$hr(),
                 uiOutput("interp_boxplot")
        ),
        
        tabPanel("Tables",
                 h4("Theoretical Borderline COR (design threshold for success)"),
                 tableOutput("tbl_theory"),
                 tags$hr(),
                 h4("Simulated OR at Stopping (95% percentile interval)"),
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
  
  sim_res <- eventReactive(input$run, {
    msg <- invalid_probs_msg()
    if (!is.null(msg)) {
      showNotification(msg, type = "error")
      return(NULL)
    }
    
    withProgress(message = "Running Monte-Carlo simulations...", value = 0, {
      simulate_obf_ordinal(
        COR_true   = input$COR_true,
        COR_NI     = input$M_margin,
        n_total    = input$N_total,
        futility_frac = input$futility_frac,
        info_frac  = input$info_frac,
        alpha1     = input$alpha1,
        alpha2     = input$alpha2,
        futility_p = input$futility_p,
        p_control  = parse_probs(input$p_control),
        seed       = input$seed,
        nSims      = input$nSims
      )
    })
  })
  
  # Theory table
  output$tbl_theory <- renderTable({
    req(input$N_total, input$info_frac, input$M_margin, input$p_control)
    p <- parse_probs(input$p_control)
    if (!is.null(validate_probs(p))) return(data.frame(Note = validate_probs(p)))
    
    theta <- theta_from_control_pmf(p)
    N_total <- input$N_total
    N_ia    <- round(input$info_frac * N_total)
    
    b_ia    <- borderline_beta_one_look(N_ia,    input$alpha1, input$M_margin, theta)
    b_final <- borderline_beta_one_look(N_total, input$alpha2, input$M_margin, theta)
    
    data.frame(
      Look           = c("Efficacy IA", "Final"),
      N              = c(N_ia, N_total),
      Alpha          = c(input$alpha1, input$alpha2),
      Borderline_COR = round(exp(c(b_ia, b_final)), 3)
    )
  }, digits = 3)
  
  # Simulation table – row-by-row construction to avoid any mismatch
  output$tbl_sim <- renderTable({
    req(sim_res())
    sim <- sim_res()
    
    fut <- safe_summ_or(sim$OR_fut_all[sim$stop_fut])
    ia  <- safe_summ_or(sim$OR1_all[sim$stop_ia])
    fin <- safe_summ_or(sim$OR2_all[sim$stop_final])
    
    data.frame(
      Stage  = c("Futility stop", "IA success stop", "Final success stop"),
      N      = c(fut["N"], ia["N"], fin["N"]),
      Mean   = c(fut["Mean"], ia["Mean"], fin["Mean"]),
      Median = c(fut["Median"], ia["Median"], fin["Median"]),
      `2.5%` = c(fut["2.5%"], ia["2.5%"], fin["2.5%"]),
      `97.5%`= c(fut["97.5%"], ia["97.5%"], fin["97.5%"]),
      check.names = FALSE
    ) |>
      dplyr::mutate(
        across(c(Mean, Median, `2.5%`, `97.5%`), ~round(.x, 3))
      )
  }, digits = 3, na = "—")
  
  # Clinician-friendly interpretation text
  output$interp_futility <- renderUI({
    tags$div(style="font-size: 15px; line-height: 1.4;",
             tags$b("What this means for clinicians – Futility look:"),
             tags$br(), tags$br(),
             "This is an early safety check. If the treatment looks much worse than expected at this point, the trial stops for 'futility' (no point continuing).",
             tags$br(), tags$br(),
             "Red curve = trials that stopped here. These tend to show higher (worse) odds ratios because we only stop when things look bad early.",
             tags$br(), tags$br(),
             "Blue curve = all trials at this stage. If red is shifted right of blue, it means the futility rule is working as intended.")
  })
  
  output$interp_interim <- renderUI({
    tags$div(style="font-size: 15px; line-height: 1.4;",
             tags$b("What this means for clinicians – Interim efficacy check:"),
             tags$br(), tags$br(),
             "This is an early chance to declare the treatment non-inferior (good enough) if the data look very promising.",
             tags$br(), tags$br(),
             "Red curve = trials that stopped early and declared non-inferiority. These often look **too good** (lower odds ratios) because we only stop when the result is unusually strong at this point.",
             tags$br(), tags$br(),
             "This is called 'winner’s curse' — early positive results tend to overestimate the true benefit. Blue curve = all trials that reached this stage.",
             tags$br(), tags$br(),
             "Clinically: be cautious interpreting very impressive early results — they are selected for optimism.")
  })
  
  output$interp_final <- renderUI({
    tags$div(style="font-size: 15px; line-height: 1.4;",
             tags$b("What this means for clinicians – Final analysis:"),
             tags$br(), tags$br(),
             "These are the trials that did not stop early (neither for futility nor for strong early success).",
             tags$br(), tags$br(),
             "Red curve = those that met the non-inferiority criterion at the end.",
             tags$br(), tags$br(),
             "Because trials with very strong early signals were already stopped, the remaining ones can look less impressive on average — this is normal conditioning, not bias in the final estimate.")
  })
  
  output$interp_overview <- renderUI({
    tags$div(style="font-size: 15px; line-height: 1.4;",
             tags$b("What clinicians should understand from these graphs:"),
             tags$br(), tags$br(),
             "A) Futility: stopped trials look worse (higher OR) — this is expected and protects patients.",
             tags$br(), tags$br(),
             "B) Interim success: stopped trials look better (lower OR) — this is the 'winner’s curse'. Early positive results often look overly impressive.",
             tags$br(), tags$br(),
             "C) Final: results among trials that continued to the end — may appear more modest because strong early signals were already counted separately.",
             tags$br(), tags$br(),
             "Bottom line: early positive stopping can make results look too good; final results are valid but reflect a different (selected) set of trials.")
  })
  
  output$interp_boxplot <- renderUI({
    tags$div(style="font-size: 15px; line-height: 1.4;",
             tags$b("Boxplot summary for clinicians:"),
             tags$br(), tags$br(),
             "Boxes show the spread of estimated odds ratios (log scale) at each stage.",
             tags$br(), tags$br(),
             "'All' = all trials that reached that point.",
             tags$br(), tags$br(),
             "'Stopped' = trials that stopped at that point (selected for strong/weak signal).",
             tags$br(), tags$br(),
             "Percentages = how many trials fell into each group out of all simulations.",
             tags$br(), tags$br(),
             "Green dashed line = true underlying effect. Red dashed line = non-inferiority margin.")
  })
  
  # ── Plots with fallbacks ────────────────────────────────────────────────
  
  output$plot_futility <- renderPlot({
    req(sim_res())
    sim <- sim_res()
    
    ok <- is.finite(sim$Z_fut_all) & is.finite(sim$OR_fut_all)
    Z   <- sim$Z_fut_all[ok]
    OR  <- sim$OR_fut_all[ok]
    
    cond <- Z > sim$z_fut
    OR_cond <- OR[cond]
    
    op <- par(mfrow = c(1,2), mar = c(4,4,3,1))
    
    if (length(Z) >= 10) {
      hist(Z, breaks = 50, freq = FALSE, col = "gray90", main = "Futility Z", xlab = "Z_fut")
      abline(v = sim$z_fut, col = "red3", lwd = 2, lty = 2)
    } else {
      plot.new()
      text(0.5, 0.5, "Too few trials reached futility look", cex = 1.3, col = "gray50")
    }
    
    if (length(OR) >= 10) {
      d_all <- density(OR, na.rm = TRUE)
      d_cond <- if (sum(cond) >= 5) density(OR_cond, na.rm = TRUE) else NULL
      ymax <- max(c(d_all$y, if (!is.null(d_cond)) d_cond$y else 0)) * 1.15
      xquant <- quantile(OR, c(0.005, 0.995), na.rm = TRUE)
      xlim <- c(xquant[1] * 0.9, xquant[2] * 1.1)
      plot(d_all, lwd = 2, col = "#1f77b4", ylim = c(0, ymax), xlim = xlim,
           main = "OR at Futility", xlab = "OR", ylab = "Density")
      if (!is.null(d_cond)) lines(d_cond, lwd = 2, col = "#d62728")
    } else {
      plot.new()
      text(0.5, 0.5, "Too few valid OR values at futility", cex = 1.3, col = "gray50")
    }
    
    par(op)
  })
  
  output$plot_interim <- renderPlot({
    req(sim_res())
    sim <- sim_res()
    
    ok <- is.finite(sim$Z1_all) & is.finite(sim$OR1_all)
    Z   <- sim$Z1_all[ok]
    OR  <- sim$OR1_all[ok]
    
    cond <- Z <= sim$zcrit1
    OR_cond <- OR[cond]
    
    op <- par(mfrow = c(1,2), mar = c(4,4,3,1))
    
    if (length(Z) >= 10) {
      hist(Z, breaks = 50, freq = FALSE, col = "gray90", main = "Interim Z", xlab = "Z1")
      abline(v = sim$zcrit1, col = "red3", lwd = 2, lty = 2)
    } else {
      plot.new()
      text(0.5, 0.5, "Too few trials reached interim", cex = 1.3, col = "gray50")
    }
    
    if (length(OR) >= 10) {
      d_all <- density(OR, na.rm = TRUE)
      d_cond <- if (sum(cond) >= 5) density(OR_cond, na.rm = TRUE) else NULL
      ymax <- max(c(d_all$y, if (!is.null(d_cond)) d_cond$y else 0)) * 1.15
      xquant <- quantile(OR, c(0.005, 0.995), na.rm = TRUE)
      xlim <- c(xquant[1] * 0.9, xquant[2] * 1.1)
      plot(d_all, lwd = 2, col = "#1f77b4", ylim = c(0, ymax), xlim = xlim,
           main = "OR at Interim", xlab = "OR", ylab = "Density")
      if (!is.null(d_cond)) lines(d_cond, lwd = 2, col = "#d62728")
    } else {
      plot.new()
      text(0.5, 0.5, "Too few valid OR values at interim", cex = 1.3, col = "gray50")
    }
    
    par(op)
  })
  
  output$plot_final <- renderPlot({
    req(sim_res())
    sim <- sim_res()
    
    ok <- is.finite(sim$Z2_all) & is.finite(sim$OR2_all)
    Z   <- sim$Z2_all[ok]
    OR  <- sim$OR2_all[ok]
    
    cond <- Z <= sim$zcrit2
    OR_cond <- OR[cond]
    
    op <- par(mfrow = c(1,2), mar = c(4,4,3,1))
    
    if (length(Z) >= 10) {
      hist(Z, breaks = 50, freq = FALSE, col = "gray90", main = "Final Z", xlab = "Z2")
      abline(v = sim$zcrit2, col = "red3", lwd = 2, lty = 2)
    } else {
      plot.new()
      text(0.5, 0.5, "Too few trials reached final", cex = 1.3, col = "gray50")
    }
    
    if (length(OR) >= 10) {
      d_all <- density(OR, na.rm = TRUE)
      d_cond <- if (sum(cond) >= 5) density(OR_cond, na.rm = TRUE) else NULL
      ymax <- max(c(d_all$y, if (!is.null(d_cond)) d_cond$y else 0)) * 1.15
      xquant <- quantile(OR, c(0.005, 0.995), na.rm = TRUE)
      xlim <- c(xquant[1] * 0.9, xquant[2] * 1.1)
      plot(d_all, lwd = 2, col = "#1f77b4", ylim = c(0, ymax), xlim = xlim,
           main = "OR at Final", xlab = "OR", ylab = "Density")
      if (!is.null(d_cond)) lines(d_cond, lwd = 2, col = "#d62728")
    } else {
      plot.new()
      text(0.5, 0.5, "Too few valid OR values at final", cex = 1.3, col = "gray50")
    }
    
    par(op)
  })
  
  output$plot_overview <- renderPlot({
    req(sim_res())
    sim <- sim_res()
    
    log_true <- log(input$COR_true)
    log_M    <- log(input$M_margin)
    
    ok_f <- is.finite(sim$Z_fut_all) & is.finite(sim$OR_fut_all)
    Lf_all  <- log(sim$OR_fut_all[ok_f])
    cond_f  <- sim$Z_fut_all[ok_f] > sim$z_fut
    Lf_stop <- Lf_all[cond_f]
    
    ok1 <- is.finite(sim$Z1_all) & is.finite(sim$OR1_all)
    L1_all  <- log(sim$OR1_all[ok1])
    cond1   <- sim$Z1_all[ok1] <= sim$zcrit1
    L1_stop <- L1_all[cond1]
    
    ok2 <- is.finite(sim$Z2_all) & is.finite(sim$OR2_all)
    L2_all  <- log(sim$OR2_all[ok2])
    cond2   <- sim$Z2_all[ok2] <= sim$zcrit2
    L2_stop <- L2_all[cond2]
    
    op <- par(mfrow = c(1,3), mar = c(4,4,3,1))
    
    if (length(Lf_all) >= 10) {
      d_all <- density(Lf_all, na.rm = TRUE)
      d_stop <- if (sum(cond_f) >= 5) density(Lf_stop, na.rm = TRUE) else NULL
      ymax <- max(c(d_all$y, if (!is.null(d_stop)) d_stop$y else 0)) * 1.15
      plot(d_all, lwd = 2, col = "#1f77b4", ylim = c(0, ymax),
           main = "A) Futility log(OR)", xlab = "log(OR)", ylab = "Density")
      if (!is.null(d_stop)) lines(d_stop, lwd = 2, col = "#d62728")
      abline(v = c(log_true, log_M), lty = c(2,3), col = c("black","gray50"), lwd = 2)
    } else {
      plot.new()
      text(0.5, 0.5, "Too few trials reached futility", cex = 1.3, col = "gray50")
    }
    
    if (length(L1_all) >= 10) {
      d_all <- density(L1_all, na.rm = TRUE)
      d_stop <- if (sum(cond1) >= 5) density(L1_stop, na.rm = TRUE) else NULL
      ymax <- max(c(d_all$y, if (!is.null(d_stop)) d_stop$y else 0)) * 1.15
      plot(d_all, lwd = 2, col = "#1f77b4", ylim = c(0, ymax),
           main = "B) Interim log(OR)", xlab = "log(OR)", ylab = "Density")
      if (!is.null(d_stop)) lines(d_stop, lwd = 2, col = "#d62728")
      abline(v = c(log_true, log_M), lty = c(2,3), col = c("black","gray50"), lwd = 2)
    } else {
      plot.new()
      text(0.5, 0.5, "Too few trials reached interim", cex = 1.3, col = "gray50")
    }
    
    if (length(L2_all) >= 10) {
      d_all <- density(L2_all, na.rm = TRUE)
      d_stop <- if (sum(cond2) >= 5) density(L2_stop, na.rm = TRUE) else NULL
      ymax <- max(c(d_all$y, if (!is.null(d_stop)) d_stop$y else 0)) * 1.15
      plot(d_all, lwd = 2, col = "#1f77b4", ylim = c(0, ymax),
           main = "C) Final log(OR)", xlab = "log(OR)", ylab = "Density")
      if (!is.null(d_stop)) lines(d_stop, lwd = 2, col = "#d62728")
      abline(v = c(log_true, log_M), lty = c(2,3), col = c("black","gray50"), lwd = 2)
    } else {
      plot.new()
      text(0.5, 0.5, "Too few trials reached final", cex = 1.3, col = "gray50")
    }
    
    par(op)
  })
  
  # Boxplot – labels moved right with extra margin and xpd
  output$plot_boxplot <- renderPlot({
    req(sim_res())
    sim <- sim_res()
    
    log_true <- log(input$COR_true)
    log_M    <- log(input$M_margin)
    
    ok_f <- is.finite(sim$OR_fut_all)
    ok1  <- is.finite(sim$OR1_all)
    ok2  <- is.finite(sim$OR2_all)
    
    groups <- list(
      "Fut all"      = log(sim$OR_fut_all[ok_f]),
      "Fut stopped"  = log(sim$OR_fut_all[ok_f & sim$stop_fut]),
      "IA all"       = log(sim$OR1_all[ok1]),
      "IA stopped"   = log(sim$OR1_all[ok1 & sim$stop_ia]),
      "Final all"    = log(sim$OR2_all[ok2]),
      "Final stopped"= log(sim$OR2_all[ok2 & sim$stop_final])
    )
    groups <- groups[sapply(groups, function(v) length(v) > 1)]
    
    if (length(groups) == 0) {
      plot.new()
      text(0.5, 0.5, "No valid data for boxplot", cex = 1.4, col = "gray50")
      return()
    }
    
    # Extra right margin for labels
    op <- par(mar = c(5, 4, 4, 10))  # right margin = 10 lines
    
    boxplot(groups, horizontal = TRUE, las = 1, col = "gray92",
            main = "Selection effects (log scale)", xlab = "log(COR)")
    abline(v = c(log_true, log_M), lty = c(2,3), col = c("darkgreen","red"), lwd = 2)
    
    p_fut   <- mean(sim$stop_fut)
    p_ia    <- mean(sim$stop_ia)
    p_final <- mean(sim$stop_final)
    
    props <- c(
      "Fut all"      = 1,
      "Fut stopped"  = p_fut,
      "IA all"       = 1 - p_fut,
      "IA stopped"   = p_ia,
      "Final all"    = 1 - p_fut - p_ia,
      "Final stopped"= p_final
    )[names(groups)]
    
    usr <- par("usr")
    x_text <- usr[2] + 0.08 * (usr[2] - usr[1])   # slightly closer than before
    
    text(x_text, seq_along(groups), sprintf("%.1f%%", 100 * props), 
         adj = 0, cex = 0.85, font = 2, col = "black", xpd = TRUE)
    
    text(x_text, max(seq_along(groups)) + 0.8, "% trials", 
         adj = 0, cex = 0.8, col = "gray30", font = 2, xpd = TRUE)
    
    par(op)
  })
}

shinyApp(ui, server)