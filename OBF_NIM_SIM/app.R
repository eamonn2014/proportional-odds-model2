# =============================================================================
#   Shiny App: Ordinal Non-Inferiority Trial Simulator with Winner's Curse
# =============================================================================

library(shiny)
library(shinyWidgets)
library(MASS)     # polr
library(dplyr)
library(bslib)    # nicer theme

# ─────────────────────────────────────────────────────────────────────────────
#   Helpers (same as before — copied here for standalone app)
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

simulate_obf_ordinal <- function(
    COR_true, COR_NI, n_total, futility_frac, info_frac,
    alpha1, alpha2, futility_p, p_control,
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
  
  res <- list(
    Z_fut_all = rep(NA_real_, nSims), COR_fut_all = rep(NA_real_, nSims),
    Z1_all    = rep(NA_real_, nSims), COR1_all    = rep(NA_real_, nSims),
    Z2_all    = rep(NA_real_, nSims), COR2_all    = rep(NA_real_, nSims),
    stop_fut  = logical(nSims), stop_ia = logical(nSims), stop_final = logical(nSims),
    nSims = nSims, z_fut = qnorm(futility_p), zcrit1 = qnorm(alpha1), zcrit2 = qnorm(alpha2)
  )
  
  n_fut <- round(futility_frac * n_total)
  n1    <- round(info_frac  * n_total)
  
  pb <- if(show_progress) shiny::Progress$new() else NULL
  if (!is.null(pb)) pb$set(message = "Running simulations...", value = 0)
  
  for (i in seq_len(nSims)) {
    # Futility
    s_f <- split_n(n_fut)
    yCf <- sample(0:K, s_f$nC, TRUE, pi_control)
    yTf <- sample(0:K, s_f$nT, TRUE, pi_treat)
    df_f <- data.frame(
      y = factor(c(yCf,yTf), ordered=TRUE, levels=0:K),
      trt = factor(rep(c("C","T"), c(s_f$nC,s_f$nT)), levels=c("C","T"))
    )
    fit_f <- fit_logCOR(df_f)
    if (!is.null(fit_f)) {
      res$Z_fut_all[i]   <- (fit_f$logCOR_hat - beta_NI) / fit_f$se
      res$COR_fut_all[i] <- exp(fit_f$logCOR_hat)
      if (res$Z_fut_all[i] > res$z_fut) { res$stop_fut[i] <- TRUE; if(!is.null(pb)) pb$inc(1/nSims); next }
    }
    
    # Interim
    n_add_ia <- n1 - n_fut
    s_add <- split_n(n_add_ia)
    yCa <- sample(0:K, s_add$nC, TRUE, pi_control)
    yTa <- sample(0:K, s_add$nT, TRUE, pi_treat)
    df_add <- data.frame(
      y = factor(c(yCa,yTa), ordered=TRUE, levels=0:K),
      trt = factor(rep(c("C","T"), c(s_add$nC,s_add$nT)), levels=c("C","T"))
    )
    df1 <- rbind(df_f, df_add)
    fit1 <- fit_logCOR(df1)
    if (!is.null(fit1)) {
      res$Z1_all[i]   <- (fit1$logCOR_hat - beta_NI) / fit1$se
      res$COR1_all[i] <- exp(fit1$logCOR_hat)
      if (res$Z1_all[i] <= res$zcrit1) { res$stop_ia[i] <- TRUE; if(!is.null(pb)) pb$inc(1/nSims); next }
    }
    
    # Final
    n_add_final <- n_total - n1
    s_final <- split_n(n_add_final)
    yCf2 <- sample(0:K, s_final$nC, TRUE, pi_control)
    yTf2 <- sample(0:K, s_final$nT, TRUE, pi_treat)
    df_final_add <- data.frame(
      y = factor(c(yCf2,yTf2), ordered=TRUE, levels=0:K),
      trt = factor(rep(c("C","T"), c(s_final$nC,s_final$nT)), levels=c("C","T"))
    )
    df2 <- rbind(df1, df_final_add)
    fit2 <- fit_logCOR(df2)
    if (!is.null(fit2)) {
      res$Z2_all[i]   <- (fit2$logCOR_hat - beta_NI) / fit2$se
      res$COR2_all[i] <- exp(fit2$logCOR_hat)
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
    Stage = c("Futility stop", "IA success stop", "Final success stop"),
    N     = c(fut["N"], ia["N"], fin["N"]),
    Min   = c(fut["Min"], ia["Min"], fin["Min"]),
    Max   = c(fut["Max"], ia["Max"], fin["Max"]),
    Mean  = c(fut["Mean"], ia["Mean"], fin["Mean"]),
    Median= c(fut["Median"], ia["Median"], fin["Median"]),
    `2.5%` = c(fut["2.5%"], ia["2.5%"], fin["2.5%"]),
    `97.5%`= c(fut["97.5%"], ia["97.5%"], fin["97.5%"]),
    check.names = FALSE
  ) |>
    mutate(across(where(is.numeric), ~ round(.x, 3)))
}

# Boxplot function (slightly simplified for Shiny)
selection_boxplot <- function(
    sim,
    COR_true,
    COR_NI,
    main = "Winner's curse & selection bias (log COR scale)",
    point_cex     = 0.6,
    point_alpha   = 0.25,
    jitter_height = 0.30,
    box_width     = 0.28,
    ylim          = c(-1.5, 2.0)   # adjust as needed for your data range
) {
  log_true <- log(COR_true)
  log_M    <- log(COR_NI)
  
  # Prepare groups (same as before)
  groups <- list(
    "All @ futility"            = log(sim$COR_fut_all[is.finite(sim$COR_fut_all)]),
    "Stopped futility"          = log(sim$COR_fut_all[sim$stop_fut & is.finite(sim$COR_fut_all)]),
    "All @ interim"             = log(sim$COR1_all[is.finite(sim$COR1_all)]),
    "Stopped IA success"        = log(sim$COR1_all[sim$stop_ia & is.finite(sim$COR1_all)]),
    "All @ final"               = log(sim$COR2_all[is.finite(sim$COR2_all)]),
    "Stopped final success"     = log(sim$COR2_all[sim$stop_final & is.finite(sim$COR2_all)])
  )
  
  keep <- sapply(groups, function(g) length(g) > 0 && all(is.finite(g)))
  groups <- groups[keep]
  group_names <- names(groups)
  if (length(groups) == 0) {
    plot.new()
    text(0.5, 0.5, "No valid data", cex = 1.4, col = "gray")
    return(invisible(NULL))
  }
  
  # Plot setup – note x and y swapped visually for horizontal layout
  op <- par(mar = c(10, 12, 5, 16))
  on.exit(par(op))
  
  plot(0, type = "n",
       xlim = ylim,   # log(COR) horizontal
       ylim = c(0.4, length(groups) + 0.9),
       xlab = "", ylab = "",
       yaxt = "n", las = 1, main = main)
  
  axis(2, at = seq_along(groups), labels = group_names, las = 1, cex.axis = 0.95)
  
  # ── Draw points + custom horizontal boxes ────────────────────────────────
  for (i in seq_along(groups)) {
    vals <- groups[[i]]
    n <- length(vals)
    if (n == 0) next
    
    # Jittered points
    set.seed(2025 + i * 100)
    jitter_y <- runif(n, -jitter_height, jitter_height)
    col_pt <- if (grepl("futility", group_names[i])) {
      rgb(0.9, 0.15, 0.15, point_alpha)
    } else if (grepl("success", group_names[i])) {
      rgb(0.1, 0.65, 0.1, point_alpha)
    } else {
      rgb(0.3, 0.3, 0.3, point_alpha)
    }
    points(vals, i + jitter_y, pch = 19, cex = point_cex, col = col_pt)
    
    # ── Custom box (horizontal) ─────────────────────────────────────────────
    if (n >= 3) {
      q <- quantile(vals, probs = c(0.25, 0.5, 0.75))
      iqr <- q[3] - q[1]
      
      lower_whisk <- min(vals[vals >= q[1] - 1.5 * iqr])
      upper_whisk <- max(vals[vals <= q[3] + 1.5 * iqr])
      
      # Box body
      rect(q[1], i - box_width/2, q[3], i + box_width/2,
           col = rgb(0.88, 0.93, 1.0, 0.5), border = "steelblue", lwd = 1.4)
      
      # Median line (thick)
      segments(q[2], i - box_width/2, q[2], i + box_width/2,
               lwd = 5, col = "royalblue3")
      
      # Whiskers (horizontal lines)
      segments(lower_whisk, i, q[1], i, lwd = 2.4, col = "midnightblue")
      segments(q[3], i, upper_whisk, i, lwd = 2.4, col = "midnightblue")
      
      # Whisker caps (short vertical ticks)
      cap_len <- box_width * 0.4
      segments(lower_whisk, i - cap_len, lower_whisk, i + cap_len, lwd = 2.4, col = "midnightblue")
      segments(upper_whisk, i - cap_len, upper_whisk, i + cap_len, lwd = 2.4, col = "midnightblue")
    }
  }
  
  # ── Reference lines ───────────────────────────────────────────────────────
  abline(v = c(log_true, log_M), lty = c(2, 3), col = c("darkgreen", "red"), lwd = 2.5)
  
  # ── Percentage labels on right ───────────────────────────────────────────
  p_fut   <- mean(sim$stop_fut,   na.rm = TRUE)
  p_ia    <- mean(sim$stop_ia,    na.rm = TRUE)
  p_final <- mean(sim$stop_final, na.rm = TRUE)
  
  props <- c(1, p_fut, 1 - p_fut, p_ia, 1 - p_fut - p_ia, p_final)[keep]
  
  usr <- par("usr")
  x_text <- usr[2] + 0.14 * (usr[2] - usr[1])
  
  for (i in seq_along(groups)) {
    text(x_text, i, sprintf("%.1f%%", 100 * props[i]), adj = 0, cex = 1.05, font = 2,
         col = if (grepl("futility", group_names[i])) "firebrick" else
           if (grepl("success",   group_names[i])) "forestgreen" else "gray30",
         xpd = TRUE)
  }
  
  text(x_text, length(groups) + 0.75, "% of sims", adj = 0, cex = 1.1, font = 2,
       col = "gray40", xpd = TRUE)
  
  # Footer note
  power <- p_ia + p_final
  mtext(sprintf("True COR = %.2f  •  Futility @ %.0f%%  •  IA @ %.0f%%  •  Power ≈ %.1f%%",
                COR_true, 100 * 0.5, 100 * 0.8, 100 * power),
        side = 1, line = 3.2, cex = 0.9, col = "gray50")
  
  mtext("Winner's curse: success groups (green) are systematically too optimistic → shifted left",
        side = 1, line = 5, cex = 0.9, col = "gray50")
  
  invisible(NULL)
}
# ─────────────────────────────────────────────────────────────────────────────
#   Shiny UI
# ─────────────────────────────────────────────────────────────────────────────

ui <- page_sidebar(
  title = "Ordinal Non-Inferiority Trial Simulator + Winner's Curse",
  
  sidebar = sidebar(
    h4("Simulation Settings"),
    
    numericInput("n_total",       "Total sample size",           value = 600,   min = 200, step = 50),
    numericInput("n_sims",        "Number of simulations",       value = 1000,  min = 100, max = 10000, step = 100),
    numericInput("seed",          "Random seed",                 value = 202506, min = 1),
    
    numericInput("COR_true",      "True COR (treatment effect)", value = 1.3,   min = 0.5, step = 0.05),
    numericInput("COR_NI",        "Non-inferiority margin (M)",  value = 1.6,   min = 1.0, step = 0.1),
    
    sliderInput("futility_frac",  "Futility look fraction",      min = 0.2, max = 0.7, value = 0.5, step = 0.05),
    sliderInput("info_frac",      "Interim look fraction",        min = 0.5, max = 0.95, value = 0.80, step = 0.05),
    
    numericInput("futility_p",    "Futility p-value threshold",  value = 0.70, min = 0.5, max = 0.95, step = 0.05),
    numericInput("alpha1",        "Interim alpha (one-sided)",   value = 0.0122, min = 0.001, step = 0.001),
    numericInput("alpha2",        "Final alpha (one-sided)",     value = 0.0214, min = 0.001, step = 0.001),
    
    textInput("p_control_txt",    "Control probabilities (comma sep)",
              value = "0.04, 0.02, 0.45, 0.34, 0.15"),
    
    actionButton("run_btn", "Run Simulation", class = "btn-primary btn-lg", icon = icon("play")),
    
    hr(),
    helpText("Simulation may take 10–90 seconds depending on nSims.")
  ),
  
  card(
    card_header("Results"),
    tabsetPanel(
      tabPanel("Summary Table",
               verbatimTextOutput("status"),
               tableOutput("summary_table")
      ),
      tabPanel("Winner's Curse Plot",
               plotOutput("boxplot", height = "700px")
      ),
      tabPanel("About / Winner's Curse",
               h4("What does this app do?"),
               p("This simulator lets you explore operating characteristics of an ordinal outcome non-inferiority trial with group-sequential interim analyses (futility at ~50% and efficacy at ~80% information by default)."),
               p("It shows how often the trial stops early for futility or overwhelming efficacy, what the observed cumulative odds ratios look like at each stage, and — most importantly — how selection at interim looks creates the classical ", strong("winner's curse / conditional selection bias"), "."),
               
               h4("Winner's Curse & Selection Bias in Group-Sequential Trials"),
               p("When a trial can stop early for overwhelming efficacy, the observed treatment effect among trials/arms that ", em("stop early for success"), " is systematically ", strong("overestimated"), " (looks too optimistic / too good to be true)."),
               p("This happens because only paths where — by chance — the interim estimate is unusually large cross the early efficacy boundary. This is a form of ", em("conditional selection bias"), " or ", strong("winner's curse"), ". The green groups in the boxplot are typically shifted toward more positive values (in log-COR scale — i.e. leftward in your current horizontal layout)."),
               
               p("Conversely, trials that continue all the way to the final analysis are a ", em("selected"), " subset — they avoided both futility stopping ", em("and"), " early efficacy stopping. Their observed effects tend to be closer to (or even slightly underestimate) the truth on average."),
               
               p(strong("Important clarification:"), " The final estimate from a trial that reaches the last analysis ", em("is not biased"), " in the classical sense — it correctly reflects the treatment effect ", em("conditional on having followed the pre-specified trial pathway and stopping rules"), ". The apparent optimism is concentrated in the early-stopped successes."),
               
               h5("Note on O'Brien-Fleming boundaries (commonly used here)"),
               p("The efficacy stopping rules in designs like this often follow (or approximate) ", strong("O'Brien-Fleming"), " boundaries. These do ", em("not"), " apply explicit weighting to test statistics or p-values at different stages. Instead, they simply use very conservative (very low) nominal p-value thresholds early on (e.g., p << 0.001 at the first look), becoming progressively less stringent until the final look is close to the usual p ≈ 0.025 (one-sided)."),
               p("This approach spends almost no Type I error early (preserving most for the final analysis) and requires an extremely strong signal to stop early — which helps limit over-optimism in early-stopped results while still allowing stopping when evidence is overwhelming. No per-stage re-weighting of p-values is involved; the adjustment is purely through changing the critical value (or nominal significance level) at each pre-specified look."),
               
               p("In practice this phenomenon can potentially lead to:"),
               tags$ul(
                 tags$li("Overly optimistic power calculations for future trials"),
                 tags$li("Inflated effect size estimates in meta-analyses if early-stopped trials dominate"),
                 tags$li("Overconfidence when selecting promising doses / arms / populations")
               ),
               p(style = "font-style: italic; color: #666;",
                 "Note: While estimates are generally trustworthy unconditionally (i.e., averaging across all trials, with typically small overall bias under conservative designs like O'Brien-Fleming), the conditional overestimation in early successes can still impact these downstream applications if unadjusted estimates are used naively."),
               
               p("Mitigation options include:"),
               tags$ul(
                 tags$li("Very conservative (spend little alpha early) efficacy boundaries — like O'Brien-Fleming"),
                 tags$li("Bias-adjusted point estimates / confidence intervals"),
                 tags$li("Bayesian shrinkage methods"),
                 tags$li("Simply reporting the design and cautioning about interpretation of early successes")
               ),
               
               h5("Selected references"),
               tags$ul(
                 tags$li("O'Brien & Fleming (1979). A multiple testing procedure for clinical trials. Biometrics."),
                 tags$li("Bassler et al. (2010). Systematic reviewers neglect bias that results from trials stopped early for benefit. J Clin Epidemiol."),
                 tags$li("Pocock & Hughes (1990). Estimation issues in clinical trials with sequential monitoring. Stat Med."),
                 tags$li("Whitehead (1986). On the bias of maximum likelihood estimation following group sequential tests. Biometrika."),
                 tags$li("Good reviews: Proschan, Lan & Wittes (2006). Group sequential methods with applications to clinical trials.")
               ),
               
               p(style = "color: #555; font-size: 0.95em; margin-top: 2em;",
                 "The boxplot deliberately shows both unconditional (all simulations at a look) and conditional (stopped) distributions so you can see the selection effect visually.")
      )
    )
  )
)

# ─────────────────────────────────────────────────────────────────────────────
#   Server
# ─────────────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {
  
  sim_result <- eventReactive(input$run_btn, {
    req(input$run_btn)
    
    output$status <- renderText("Running simulations... please wait.")
    
    p_control <- parse_probs(input$p_control_txt)
    msg <- validate_probs(p_control)
    if (!is.null(msg)) {
      showNotification(msg, type = "error", duration = 10)
      return(NULL)
    }
    
    simulate_obf_ordinal(
      COR_true      = input$COR_true,
      COR_NI        = input$COR_NI,
      n_total       = input$n_total,
      futility_frac = input$futility_frac,
      info_frac     = input$info_frac,
      alpha1        = input$alpha1,
      alpha2        = input$alpha2,
      futility_p    = input$futility_p,
      p_control     = p_control,
      seed          = input$seed,
      nSims         = input$n_sims,
      show_progress = TRUE
    )
  }, ignoreNULL = TRUE)
  
  output$status <- renderText({
    if (is.null(sim_result())) {
      "Click 'Run Simulation' to start."
    } else {
      "Simulation complete."
    }
  })
  
  output$summary_table <- renderTable({
    req(sim_result())
    sim_table(sim_result())
  }, digits = 3, align = "l")
  
  output$boxplot <- renderPlot({
    req(sim_result())
    selection_boxplot(
      sim_result(),
      COR_true = input$COR_true,
      COR_NI   = input$COR_NI
    )
  })
}

# Run the app
shinyApp(ui, server)