#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Rshiny ideas from on https://gallery.shinyapps.io/multi_regression/
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list=ls()) 
# library(ggplot2)
library(shiny) 
library(shinyWidgets)
library(shinythemes)  # more funky looking apps
library(DT)
library(shinyalert)
# library(gtools)
library(Hmisc)
# library(scales)
# library(LaplacesDemon)
# library(bayesboot)
# library(boot)
library(reshape)
library(rms)
require(ordinal)
require(ggplot2)
require(tidyverse)
#options(mc.cores = parallel::detectCores())
#rstan_options(auto_write = TRUE)
options(max.print=1000000)    
fig.width <- 400
fig.height <- 300
fig.width1 <- 1380
fig.height1 <- 700
fig.width2 <- 1400
fig.height2 <- 300
fig.width3 <- 1400  
fig.height3 <- 600
fig.width4 <- 1380
fig.height4 <- 450
fig.width5 <- 1380
fig.height5 <- 225
fig.width6 <- 400
fig.height6 <- 550
fig.width7 <- 600
fig.height7 <- 600
fig.width9 <- 1380
fig.height9 <- 500
p0 <- function(x) {formatC(x, format="f", digits=1)}
p1 <- function(x) {formatC(x, format="f", digits=1)}
p2 <- function(x) {formatC(x, format="f", digits=2)}
p3 <- function(x) {formatC(x, format="f", digits=3)}
p5 <- function(x) {formatC(x, format="f", digits=5)}

## convenience functions
logit <- function(p) log(1/(1/p-1))
expit <- function(x) 1/(1/exp(x) + 1)
inv_logit <- function(logit) exp(logit) / (1 + exp(logit))

options(width=200)
set.seed(12345) # reproducible

is.even <- function(x){ x %% 2 == 0 } # function to id. odd maybe useful
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ui <- fluidPage(theme = shinytheme("journal"), #https://www.rdocumentation.org/packages/shinythemes/versions/1.1.2
                # paper
                useShinyalert(),  # Set up shinyalert
                setBackgroundColor(
                  color = c( "#2171B5", "#F7FBFF"), 
                  gradient = "linear",
                  direction = "bottom"
                ),
                
                h2("The Proportional Odds Model"), 
                
                h4("The proportional odds model is a recommended approach for modelliing an ordinal response. 
                Patient reported outcomes are often reported using an ordinal reponse. Often they are analysed using a linear model treating the outcome as continuous. 
                This is not a correct assumption. For example the outcomes are strict integers there are no in between levels.
                The scales too have a distinct bottom level and top level, a likert scale may have 5 levels, so there is no level 6. 
                A linear model may predict responses above or below the only levels possible!
         "), 
                
                h3("  "), 
                
                
                sidebarLayout(
                  
                  sidebarPanel( width=3 ,
                                
                                tags$style(type="text/css", ".span8 .well { background-color: #00FFFF; }"),
                                
                                
                                actionButton(inputId='ab1', label="R code",   icon = icon("th"),   
                                             onclick ="window.open('https://raw.githubusercontent.com/eamonn2014/Bayesian_bootstrap/master/Bayesian_bootstrap/app.R', '_blank')"),    
                                actionButton("resample", "simulate a new sample"),
                                br(),  
                                tags$style(".well {background-color:#b6aebd ;}"), ##ABB0B4AF
                                
                                h4("Instructions: The first input below is the number of total patients randomised 1:1 to treatment vrs placebo. 
                                     The next input is the number of ordinal levels in the response. This is followed by the 
                                     treatment proportional odds ratio. The last input is the proportional odds ratio for the baseline version of the outcome.
                                     The distribution of the baseline version of the outcome can be specified by selecting a Beta distribution that approximates that 
                                   which is expected."),
                                div(
                                  
                                  tags$head(
                                    tags$style(HTML('#ab1{background-color:orange}'))
                                  ),
                                  
                                  tags$head(
                                    tags$style(HTML('#resample{background-color:orange}'))
                                  ),
                                  
                                  textInput('n', 
                                            div(h5(tags$span(style="color:blue", "Total sample size"))), "100"),
                                  
                                  tags$hr(),
                                  textInput('dist', 
                                            div(h5(tags$span(style="color:blue", "Approximate the distribution of the baseline version of the response by specifying
                                                             Beta shape parameters"))), "2,1"),
                                  
                                  textInput('levels', 
                                            div(h5(tags$span(style="color:blue", "Number of ordinal categories in response"))), "15"),
                                  tags$hr(), 
                                  textInput('or1', 
                                            div(h5(tags$span(style="color:blue", "Treatment odds ratio"))), "1.2"),
                                 
                                  
                                  textInput('or2', 
                                            #   div(h5("Number of samples for Correlation (tab 2)")), "10"),
                                            div(h5(tags$span(style="color:blue", "Odds ratio of effect of baseline version of outcome"))), "1"),
                                  
                                  #  textInput('n2y2', 
                                  # #      div(h5("Enter the true correlation (tab 2)")), ".8"),
                                  # div(h5(tags$span(style="color:blue", "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"))), "0.8"),
                                  # tags$hr(),
                                  
                                  div(h5("References:")),  
                                  tags$a(href = "https://en.wikipedia.org/wiki/Bootstrapping_%28statistics%29", tags$span(style="color:blue", "[1] PRO"),),   
                                  div(p(" ")),
                                  tags$a(href = "https://projecteuclid.org/download/pdf_1/euclid.aos/1176345338",  tags$span(style="color:blue", "[2] PO"),),   
                                  div(p(" ")),
                                  tags$a(href = "https://projecteuclid.org/download/pdf_1/euclid.aos/1176344552", tags$span(style="color:blue", "[3] Krushke"),),
                                  div(p(" ")),
                                  tags$a(href = "https://blogs.sas.com/content/iml/2017/09/20/fishers-transformation-correlation.html", tags$span(style="color:blue", "[4] xxxxxx"),),  
                                  div(p(" ")),
                                  tags$a(href = "https://en.wikisource.org/wiki/Keats;_poems_published_in_1820/Bards_of_Passion_and_of_Mirth", tags$span(style="color:blue", "xxxxxx"),),  
                                  div(p(" ")),
                                  tags$hr()
                                )
                                
                                
                  ),
                  
                  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~tab panels
                  mainPanel(width=9,
                            
                            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                            #    tabsetPanel(type = "tabs", 
                            navbarPage(       
                              #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
                              tags$style(HTML("
                            .navbar-default .navbar-brand {color: orange;}
                            .navbar-default .navbar-brand:hover {color: blue;}
                            .navbar { background-color: #b6aebd;}
                            .navbar-default .navbar-nav > li > a {color:black;}
                            .navbar-default .navbar-nav > .active > a,
                            .navbar-default .navbar-nav > .active > a:focus,
                            .navbar-default .navbar-nav > .active > a:hover {color: pink;background-color: purple;}
                            .navbar-default .navbar-nav > li > a:hover {color: black;background-color:yellow;text-decoration:underline;}
                            .navbar-default .navbar-nav > li > a[data-value='t1'] {color: red;background-color: pink;}
                            .navbar-default .navbar-nav > li > a[data-value='t2'] {color: blue;background-color: lightblue;}
                            .navbar-default .navbar-nav > li > a[data-value='t3'] {color: green;background-color: lightgreen;}
                   ")),
                              
                              #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                              tabPanel("1 Baseline version of response", value=7, 
                                       h4("The distribution of the baseline version of the response variable is specified here.
                                          By selecting a beta distribution using the shape parameters the
                                          expected baseline counts in categories can be approximated. The default is Beta(2,1)."),
                                       
                                       #    h4(paste("Figure 1. Bayesian and frequentist bootstrap distributions, estimating one sample mean")), 
                                       #   div(plotOutput("diff", width=fig.width4, height=fig.height4)),       
                                       
                                       # fluidRow(
                                       #     column(width = 6, offset = 0, style='padding:1px;',
                                       #            h4("Proprtional odds model"), 
                                       #            div( verbatimTextOutput("reg.summary2") )
                                       #     )
                                       #     ),
                                       
                                       
                                       ###############
                                       
                                       
                                       fluidRow(
                                         column(width = 6, offset = 0, style='padding:1px;',
                                                #h4("Proportional odds model"), 
                                                # div( verbatimTextOutput("reg.summary2") )
                                                div(plotOutput("beta",  width=fig.width7, height=fig.height7)),
                                         ) ,
                                         
                                         fluidRow(
                                           column(width = 5, offset = 0, style='padding:1px;',
                                                  # h4("Proportional odds ratio summaries. Do we recover the input odds ratios..."),
                                                  # div( verbatimTextOutput("reg.summary3")),
                                                  div(plotOutput("reg.plotx",  width=fig.width7, height=fig.height7)) 
                                                  #  h4(htmlOutput("textWithNumber",) ),
                                           ))),
                                       
                                       
                                       #  h4(htmlOutput("textWithNumber",) ),
                                       
                                       # h4(htmlOutput("textWithNumber1",) ),
                              ) ,
                              #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                              
                              
                              tabPanel("2 Proportional odds model", value=7, 
                                       h4("  (when all other variables are set to zero)"),
                                       
                                       #    h4(paste("Figure 1. Bayesian and frequentist bootstrap distributions, estimating one sample mean")), 
                                       #   div(plotOutput("diff", width=fig.width4, height=fig.height4)),       
                                       
                                       # fluidRow(
                                       #     column(width = 6, offset = 0, style='padding:1px;',
                                       #            h4("Proprtional odds model"), 
                                       #            div( verbatimTextOutput("reg.summary2") )
                                       #     )
                                       #     ),
                                       
                                       
                                       ###############
                                       
                                       
                                       fluidRow(
                                         column(width = 6, offset = 0, style='padding:1px;',
                                                h4("Proportional odds model lrm function"), 
                                                div( verbatimTextOutput("reg.summary2") )
                                         ) ,
                                         
                                         fluidRow(
                                           column(width = 6, offset = 0, style='padding:1px;',
                                                  h4("Proportional odds model orm function"), 
                                                  div( verbatimTextOutput("reg.summary1")),
                                                  h4("Proportional odds ratio summaries. Do we recover the input odds ratios..."),
                                                  div( verbatimTextOutput("reg.summary3")),
                                                  
                                                  h4(htmlOutput("textWithNumber",) ),
                                           ))),
                                       
                                       
                                       #  h4(htmlOutput("textWithNumber",) ),
                                       
                                       # h4(htmlOutput("textWithNumber1",) ),
                              ) ,
                              
                              
                              ################
                              
                              
                              
                              
                              
                              #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                              tabPanel("3 Barplots of outcome", value=3, 
                                       h4("xxxxxxxxxxxxxxxxxxxxxx."),
                                       h4(paste("Figure 3. xxxxxxxxxxxxxxxxxx")),  
                                       div(plotOutput("reg.plot99", width=fig.width1, height=fig.height1)),
                                       
                                       fluidRow(
                                         column(width = 7, offset = 0, style='padding:1px;',
                                                h4("xxxxxxxxxxxxxxxxxxxxxn"), 
                                                # div( verbatimTextOutput("reg.summary4"))
                                         )),
                                       
                                       
                              ),
                              
                              #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                              tabPanel("2 Barplot", value=3, 
                                       h4("xxxxxxxxxxxxxxxxxx."),
                                       #h4(paste("Figure 2. xxxxxxxxxxxxxxxxxx")), 
                                       
                                       div(plotOutput("diff",  width=fig.width5, height=fig.height5)),
                                       div(plotOutput("diff2", width=fig.width5, height=fig.height5)),     
                                       div(plotOutput("diff3", width=fig.width5, height=fig.height5)),  
                                       h4(paste("Figure 1. Barplots of counts in ordinal variable, baseline, outcome separately for treated and placebo patients (not very informative unless 
                         you select a large number of patients and large treatment odds ratio) ")), 
                                       fluidRow(
                                         column(width = 7, offset = 0, style='padding:1px;',
                                                h4("xxxxxxxxxxxxxxx"), 
                                                # div( verbatimTextOutput("reg.summary3"))
                                         )),
                                       
                                       
                                       
                              ),
                              #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                              tabPanel("3 Barplot", value=3, 
                                       h4("xxxxxxxxxxxxxxxxxxxxxx."),
                                       h4(paste("Figure 3. xxxxxxxxxxxxxxxxxx")),  
                                       div(plotOutput("reg.plot", width=fig.width1, height=fig.height1)),
                                       
                                       fluidRow(
                                         column(width = 7, offset = 0, style='padding:1px;',
                                                h4("xxxxxxxxxxxxxxxxxxxxxn"), 
                                                # div( verbatimTextOutput("reg.summary4"))
                                         )),
                                       
                                       
                              ),
                              #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                              tabPanel("4 Linear model", value=3, 
                                       h4("ANCOVA model"),
                                       
                                       fluidRow(
                                         column(
                                           DT::dataTableOutput("tablex"),width = 6)
                                       ),
                                       
                                       fluidRow(
                                         column(width = 6, offset = 0, style='padding:1px;',
                                                #h4("ANCOVA model"), 
                                                div( verbatimTextOutput("reg.summary4") )
                                         ) ,
                                         
                                         fluidRow(
                                           column(width = 5, offset = 0, style='padding:1px;',
                                                  #   h4("Proportional odds ratio summaries. Do we recover the input odds ratios..."),
                                                  div( verbatimTextOutput("reg.summary5")),
                                                  
                                                  # h4(htmlOutput("textWithNumber",) ),
                                           ))),
                                       
                                       
                                       
                              ) ,
                              
                              
                              
                              
                              
                              #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~NEW
                              tabPanel("5 Assumptions", fluid = TRUE, width = 4,
                                       
                                       h4(("Upload your own data for correlation analysis. Requires 2 columns of numeric data. Select 'Header' 
                         if your data columns have names. 
                         The top two radio button options are to help load,
                                 the bottom option is to either print the top six rows of the data or show all the data.
                               ")) ,
                                       
                                       h4(("Here are two example data sets (download either file and click 'Browse...' to locate and upload for the analysis):")) ,
                                       
                                       tags$a(href = "https://raw.githubusercontent.com/eamonn2014/Bayesian_bootstrap/master/icreamsales", tags$span(style="color:blue", "Example 1 Ice cream sales and temperature N=12, has a header"),), 
                                       div(p(" ")),
                                       
                                       tags$a(href = "https://raw.githubusercontent.com/eamonn2014/Bayesian_bootstrap/master/height.selfesteem", tags$span(style="color:blue", "Example 2 height and self esteem N=20, no header"),), 
                                       div(p(" ")),
                                       
                                       sidebarLayout(
                                         
                                         # Sidebar panel for inputs ----
                                         sidebarPanel(
                                           
                                           # Input: Select a file ----
                                           fileInput("file1", "Choose CSV File",
                                                     multiple = TRUE,
                                                     accept = c("text/csv",
                                                                "text/comma-separated-values,text/plain",
                                                                ".csv")),
                                           
                                           # Horizontal line ----
                                           tags$hr(),
                                           
                                           # Input: Checkbox if file has header ----
                                           checkboxInput("header", "Header", TRUE),
                                           
                                           # Input: Select separator ----
                                           radioButtons("sep", "Separator",
                                                        choices = c(Comma = ",",
                                                                    Semicolon = ";",
                                                                    Tab = "\t",
                                                                    Whitespace = ""),
                                                        selected = ""),
                                           
                                           # Input: Select quotes ----
                                           radioButtons("quote", "Quote",
                                                        choices = c(None = "",
                                                                    "Double Quote" = '"',
                                                                    "Single Quote" = "'"),
                                                        selected = ''),
                                           
                                           # Horizontal line ----
                                           tags$hr(),
                                           
                                           # Input: Select number of rows to display ----
                                           radioButtons("disp", "List all the data or first 6 rows only",
                                                        choices = c(Head = "head",
                                                                    All = "all"),
                                                        selected = "head"),
                                           
                                           # Horizontal line ----
                                           # tags$hr(),
                                           
                                           # Input: Select number of rows to display ----
                                           # radioButtons("what", "Output",
                                           #              choices = c(Analysis = "Analysis",
                                           #                          Plot = "plot"),
                                           #              selected = "Analysis")
                                           
                                         ),
                                         
                                         # Main panel for displaying outputs ----
                                         mainPanel(
                                           
                                           # Output: Data file ----
                                           h4(paste("Figure 4. Bayesian and frequentist bootstrap distributions of user's data, estimating correlation")),  
                                           div(plotOutput("contents2", width=fig.width6, height=fig.height6)),
                                           #div(verbatimTextOutput("contents2")),
                                           #plotOutput("plotx"),
                                           tags$hr(),
                                           h4("Correlation and 95% confidence interval from R cor.test function"), 
                                           #div( verbatimTextOutput("reg.summary5")),
                                           tags$hr(),
                                           h4("Print the data listing"),
                                           tableOutput("contents") 
                                           
                                           
                                         ),
                                         
                                       )
                              ) ,
                              
                              
                              #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                              # tabPanel("6 inference", value=3, 
                              #          h4("xxxxxxxxxxxxxxxxxxxxxx."),
                              #          h4(paste("Figure 3. xxxxxxxxxxxxxxxxxx")),  
                              #          #div(plotOutput("preds", width=fig.width1, height=fig.height1)),
                              #          
                              #          fluidRow(
                              #            column(width = 7, offset = 0, style='padding:1px;',
                              #                   h4("xxxxxxxxxxxxxxxxxxxxxn"), 
                              #                   # div( verbatimTextOutput("reg.summary4"))
                              #            )),
                              #          
                              #          
                              # ),
                              
                              tabPanel("4 Take home messages", 
                                       
                                       div(plotOutput("preds", width=fig.width9, height=fig.height9)),
                                       
                                       fluidRow(
                                    
                                         textInput('base', 
                                        div(h5(tags$span(style="color:blue", 
                                                         "Enter a patient's baseline category and see their predicted probabilities for response outcomes"))), "1")
                                         
                                         
                                       ),
                                       
                                       #h4(htmlOutput("textWithNumber4",) ) ,
                                       width = 30 ),
                              #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                              
                              tabPanel("5 pred probs", 
                                      
                                      div(plotOutput("predicts", width=fig.width9, height=fig.height9)),
                                      
                                      fluidRow(
                                        
                                        textInput('group', 
                                                  div(h5(tags$span(style="color:blue", 
                                                                   "select treatment group: 0 for placebo, 1 for treatment, 2 for both"))), "1"),
                                        
                                        textInput('rcat', 
                                                  div(h5(tags$span(style="color:blue", 
                                                                   "Response category"))), "1,2,3,4,5")
                                        
                                        
                                      ),
                                      
                                      
                                      #h4(htmlOutput("textWithNumber4",) ) ,
                                      width = 30 )
                              
                              
                              #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                              
                              # textInput('levels', 
                              #           div(h5(tags$span(style="color:blue", "Number of ordinal categories in response"))), "15"),
                              # tags$hr(), 
                              
                              
                              
                              #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   END NEW   
                            )
                            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  )
                ) 
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~end tab panels 
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

server <- shinyServer(function(input, output   ) {
  
  shinyalert("Welcome! \nExplore the Proportional odds model!",
             "It's in progress!",
             type = "info")
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This is where a new sample is instigated 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  random.sample <- reactive({
    
    foo <- input$resample
    
    #sample sizes
    dis <- as.numeric(unlist(strsplit(input$dist,",")))
    
    #sample sizes
    trt <- as.numeric(unlist(strsplit(input$n,",")))
    # mean and sD
    ctr <- as.numeric(unlist(strsplit(input$levels,",")))
    
    #sample size for correlation
    n1y1 <- log(as.numeric(unlist(strsplit(input$or1,","))))   # user enter odds , need log for the maths
    # R
    n2y2 <- log(as.numeric(unlist(strsplit(input$or2,","))))    # user enter odds , need log for the maths
    
    
    base<- as.numeric(unlist(strsplit(input$base,",")))
    
   group <- (as.numeric(unlist(strsplit(input$group,","))))    
    # R
    rcat <- (as.numeric(unlist(strsplit(input$rcat,","))))     
    
    
    return(list(  
      n=trt[1],  
      lev=ctr[1],
      or1=n1y1[1], 
      or2=n2y2[1],
      shape1=dis[1], 
      shape2=dis[2],
      base=base[1],
      group=group[1],
      rcat=rcat[1]
      
      
      
    ))
    
  })
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # tab 1 simulate po model data and analyse
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mcmc <- reactive({
    
    sample <- random.sample()
    
    n    <- sample$n
    levz <- sample$lev
    b1  <- sample$or1
    b2  <- sample$or2
    shape1  <- sample$shape1
    shape2  <- sample$shape2
    #base  <- sample$base
    group  <- sample$group
    rcat  <- sample$rcat
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Parameters 
    
    # treatment assignment 1:1 
    treatment <- 1*(runif(n)<0.5)  
    
    # set up baseline distribution of the version of outcome, even distribution
    # prbs <- sample(1:1,levz,replace=T)
    #  prbs <- prbs/sum(prbs)                  # all equal 
    # baseline, select a baseline for everyone
    # baseline <- sample(1:(length(prbs)), n, prob=prbs, replace=TRUE)
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    prbs <- rbeta(10000,shape1,shape2)
    x <- hist(prbs, breaks= seq(min(prbs), max(prbs), length.out=levz+1), plot=FALSE)  ##
    ##str(x)
    prbs <- x$counts  # get 30 counts from beta dist.
    prbs <- prbs/sum(prbs)
    
    
    # baseline
    baseline <- sample(1:(length(prbs)), n, prob=prbs, replace=TRUE)
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Steps to generate a dataset
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ turn probs into logits 
    #  1. using the true probability of each category, calculate the levels-1 cummulative logits
    glevz <- levz-1
    
    b0 <- rep(NA, glevz)
    
    for (i in 1:glevz) {
      
      b0[i] <-  logit(1-sum(prbs[1:i])) 
      
    }
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 2. create an empty matrix columns are the categories and rows samples
    l0 <- matrix(1, nrow = n, ncol = levz)  # make a space for all patients and levels
    
    # 3. for each patient calculate linear predictor, not including step 1 yet
    lin <- b1 * treatment +  b2 * baseline   # make a linear combination
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
    # 4. combine 1 and 3 into the matrix 2.
    for (i in 2:levz) {
      
      l0[,i] <- inv_logit(b0[i-1] + lin)
      
    }   
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 5. make another matrix preparing for calculating the probabilites of membership
    fi<- matrix(1, nrow = n, ncol = levz)  
    
    # 6. Do the calculations, on the columns  
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for (i in 2:levz) {
      
      fi[,i-1] <-    l0[,i-1] - l0[,i]   # 1-2; 2-3; 3-4; 4-5; 
      
    } 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # add in the last column which need no manipulation...prob in highest level.
    fi[,levz] <- l0[, levz]
    
    # check .. should all sum to prob =1
    apply(fi, 1, sum)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # for each patient sample 1 level based on the probabilities associated with each level
    y <- c()
    for (i in 1:n) {
      y[i] <- sample(
        x = c(1:levz), 
        size = 1, 
        prob = c(fi[i,])
      )
    }
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # put the data together
    dat <- data.frame(treatment, baseline, y = factor(y))
  # dat$baseline <- factor(dat$baseline)
    # use harrell's po function analyse the data
    d <<- datadist(dat)
    options(datadist="d") 
    f1 <- lrm(y ~treatment + baseline, data=dat )  # would prefer orm here but coeffs dont show up !
    f2 <- orm(y ~treatment + baseline, data=dat ) 
    sf1 <- summary(f1, antilog=TRUE, verbose=FALSE)
    
    #f1 <- f1$coefficients 
    
    return(list(res= f1  , sf1=sf1 , dat=dat, f2=f2)) 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  })
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # non cummulative predicted probabilities plot run the analysis again
  # not efficient I know
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  output$preds <- renderPlot({
    
    sample <- random.sample()
    n    <- sample$n
    levz <- sample$lev
    base <- sample$base
    
    datx <- mcmc()$dat
    
    
    
    
  # I can get non cummulative probabilites using clm
  # Don't know how to do it with rms? yet
 
  Res.clm <- clm(y ~treatment + baseline, data=datx)
  
  # levz <-   length(unique(Res.clm$model$baseline))
  
  # summary(Res.clm)
  
  newdat <- data.frame(
    baseline =   (rep(1:levz)),
    treatment = rep(0:1, each = levz)
    )
  # newdat <- data.frame(
  #   baseline =  factor(sort(unique(Res.clm$model$baseline))),
  #   treatment = rep(0:1, each = levz)
  # )
 
  
  newdat <- cbind(newdat, predict(Res.clm, newdata=newdat, se.fit=TRUE,
                                  interval=TRUE, type="prob"))
  
  
  
  A<- cbind( newdat[,1:2], newdat[,grep("^fit", colnames(newdat)) ])
  B<- cbind( newdat[,1:2], newdat[,grep("^lwr", colnames(newdat)) ])
  C<- cbind( newdat[,1:2], newdat[,grep("^upr", colnames(newdat)) ])

  lA <- melt(data = A, id.vars = c("baseline","treatment") )
  lB <- melt(data = B, id.vars = c("baseline","treatment") )
  lC <- melt(data = C, id.vars = c("baseline","treatment") )

  lA$variable <-  gsub(".*\\.","", lA$variable)
  lB$variable <-   gsub(".*\\.","", lB$variable)
  lC$variable <-   gsub(".*\\.","", lC$variable)

  l <- cbind(lA,lB,lC)

  l <- l[,c(1:4,8,12)]
  
  names(l) <- c("baseline","treatment","response","estimate","lower", "upper")
  
  pd <- position_dodge(0.5) # move them .05 to the left and right
  
  # l$treatment <- factor(l$treatment)
  # br1 <- length(unique(l$baseline))
  # 
  # ggplot(l, aes(baseline,estimate, color=treatment)) +
  #   geom_point(aes(shape=treatment),size=4, position=pd) + 
  #   scale_color_manual(name="treatment",values=c("coral","steelblue")) + 
  #   theme_bw() + 
  #   scale_x_continuous("baseline", breaks=1:br1, labels=1:br1) + 
  #   scale_y_continuous("Probability")   + 
  #   geom_errorbar(aes(ymin=lower,ymax=upper),width=0.1,position=pd)
  # 
  
  l$response <- as.numeric(l$response)
  l$treatment <- factor( l$treatment)
  l$treatment <- ifelse( l$treatment %in% 0,"Placebo","Treatment")
  br1 <- length(unique(l$baseline))
  
  lx <- l[l$baseline %in% base,]     # user input selects this
  
 gp <- ggplot(lx, aes(response,estimate, color=treatment)) +
    geom_point(aes(shape=treatment),size=4, position=pd) + 
    scale_color_manual(name="treatment",values=c("coral","steelblue")) + 
    theme_bw() + 
    scale_x_continuous( breaks=1:br1, labels=1:br1) +   
  #  scale_y_continuous("Probability")   + 
   geom_line(position=pd, linetype = "dashed")+
    geom_errorbar(aes(ymin=lower,ymax=upper),width=0.2,position=pd) +
 
  theme(panel.background=element_blank(),
                  plot.title=element_text(), plot.margin = unit(c(5.5,12,5.5,5.5), "pt"), 
                  legend.text=element_text(size=12),
                  #legend.title=element_text(size=14),
                  legend.title=element_blank(),
                  axis.text.x = element_text(size=10),
                  axis.text.y = element_text(size=10),
                  axis.line.x = element_line(color="black"),
                  axis.line.y = element_line(color="black"),
                  axis.title.y=element_text(size=16),  
                  axis.title.x=element_text(size=16),  
                  axis.title = element_text(size = 20) , 
           plot.caption=element_text(hjust = 0, size = 7) ) +
   
    labs(title=paste0(c("xxxxxxxxxxxxxxxx"), collapse=" "), 
        x = "Response category",
        y = "Predicted probability",
          subtitle =c("xxxxxxxxxxxxxx"),
          caption = "")  
    # guides(fill=guide_legend(title="Treatment"))
    # 
  
  print(gp)
  
  })
  
   
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~not used
  
  output$predicts <- renderPlot({   
    
      sample <- random.sample()
    # 
      f    <- mcmc()$res
    # 
     levz <- sample$lev
    rcat <- sample$rcat
     group <- sample$group
     
     require(reshape)
     
     newdat <- data.frame(
       baseline = rep(1:levz),
       treatment = rep(0:1, each = levz))
     
     
     xx <- predict(f, newdat, type="fitted.ind")    #
     
     mm <- melt(data.frame(xx))
     
     mm <- cbind(newdat,mm )
     
     mm$variable <-  gsub(".*\\.","", mm$variable)
     
     mm <- plyr::arrange(mm,treatment,variable,baseline )
      
     mm$flag <- rep(seq_along( rle(mm$variable)$values ), times = rle(mm$variable)$lengths )
     
      
    mm <-  mm[mm$treatment %in% group,]
     
    mm <-  mm[mm$variable %in% rcat,]
     
     
     
     
    gpp <- ggplot(mm, aes(baseline, value, group=factor(flag))) +
       geom_line(aes(color=factor(flag))) +
     
     scale_x_continuous( breaks=1:levz, labels=1:levz) +  
 
       theme(panel.background=element_blank(),
             plot.title=element_text(), plot.margin = unit(c(5.5,12,5.5,5.5), "pt"), 
             legend.text=element_text(size=12),
             #legend.title=element_text(size=14),
             legend.title=element_blank(),
             axis.text.x = element_text(size=10),
             axis.text.y = element_text(size=10),
             axis.line.x = element_line(color="black"),
             axis.line.y = element_line(color="black"),
             axis.title.y=element_text(size=16),  
             axis.title.x=element_text(size=16),  
             axis.title = element_text(size = 20) , 
             plot.caption=element_text(hjust = 0, size = 7) ) +
       
       labs(title=paste0(c("xxxxxxxxxxxxxxxx"), collapse=" "), 
            x = "Baseline category",
            y = "Predicted probability",
            subtitle =c("xxxxxxxxxxxxxx"),
            caption = "")  
     # guides(fill=guide_legend(title="Treatment"))
     # 
     
     print(gpp)
    
    
    
    
    
    
    # 
    # # d <- datadist(newdat)
    # # options(datadist="d")
    # 
    # # L <- predict(f, newdata=newdat, se.fit=TRUE)          #omitted kint= so use 1st intercept
    # plogis(with(L, linear.predictors + 1.96*cbind(-se.fit,se.fit)))
    # predict(f, type="fitted.ind")#[1:5,]   #gets Prob(better) and all others
    # d1 <- newdat
    # predict(f, d1, type="fitted")        # Prob(Y>=j) for new observation
    # predict(f, d1, type="fitted.ind")    # Prob(Y=j)
    # predict(f, d1, type='mean', codes=TRUE) # predicts mean(y) using codes 1,2,3
    # m <- Mean(f, codes=TRUE)
    # lp <- predict(f, d1)
    # m(lp)
    # # Can use function m as an argument to Predict or nomogram to
    # # get predicted means instead of log odds or probabilities
    # dd <- datadist(baseline,treatment); options(datadist='dd')
    # m
    # plot(Predict(f, x1, fun=m), ylab='Predicted Mean')
    # # Note: Run f through bootcov with coef.reps=TRUE to get proper confidence
    # # limits for predicted means from the prop. odds model
    # options(datadist=NULL)
    # 
    # 
    # return(list(   sf1=sf1 , dat=dat)) 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  })
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # tab 2 barplot
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  output$diff <- renderPlot({
    
    sample <- random.sample()
    levz <- sample$lev
    
    dat <- mcmc()$dat
    d <- dat
    
    ym <- max(c(as.vector(table(dat$y)), as.vector(table(dat$baseline))))
    
    b1 <- barplot(table(d$baseline), las=1, main = ("Baseline distribution of outcome response all patients" ), 
                  xlab="Ordinal category", ylab = "Count", col=rainbow(50), ylim =c(0,ym))
    
    
  })
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # tab 2 barplot
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  output$diff2 <- renderPlot({
    
    sample <- random.sample()
    levz <- sample$lev
    
    dat <- mcmc()$dat
    d <- dat
    
    ym <- max(c(as.vector(table(dat$y)), as.vector(table(dat$baseline))))
    
    d <- d[d$treatment %in% 1,]
    
    treat <- barplot(table(d$y), las=1, main = ("Distribution of outcome responses treated patients only" ), 
                     xlab="Ordinal category", ylab = "Count", col=rainbow(50), ylim =c(0,ym))
    
  })
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # tab 2 barplot
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  output$diff3 <- renderPlot({
    
    sample <- random.sample()
    levz <- sample$lev
    
    dat <- mcmc()$dat
    d <- dat
    
    ym <- max(c(as.vector(table(dat$y)), as.vector(table(dat$baseline))))
    
    d <- d[d$treatment %in% 0,]
    
    treat <- barplot(table(d$y), las=1, main = ("Distribution of outcome responses placebo patients only" ), 
                     xlab="Ordinal category", ylab = "Count", col=rainbow(50), ylim =c(0,ym))
    
  })
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   ggplot barplot on side
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
  
  output$reg.plot <- renderPlot({         
    
    # Get the current regression data
    #  sample <- random.sample()
    # levz <- sample$lev
    #n   <- sample$n
    
    dat <- mcmc()$dat
    f <-   dat
    f <-   as.data.frame(table(f$y))
    
    f$Percentage <- round(f$Freq / sum(f$Freq)*100,1)
    
    z <- f              # data set for plot
    variable <- "Freq"  # variable of interest
    pN <- sum(f$Freq)   
    pN <- format(pN, big.mark=","
                 ,scientific=FALSE)
    roundUp <- function(x) 10^ceiling(log10(x))/2
    gupper <- roundUp((max(f$Freq)))  # plot upper limit
    gupper <- ceiling((max(f$Freq)))  # plot upper limit
    glower <- 0                       # plot lower limit
    gstep <- 5                        # grid steps
    
    # text for plot
    ylabel <- "Counts" 
    
    z$N <- z$Freq
    
    Gplotx <- function(data,  l1,l2,l3 ) {
      
      mlimit=l1
      
      p1 <- ggplot(data = data, aes(x =  Var1, y = N, fill = Var1)) + 
        
        geom_bar(stat = "identity", width =0.7) 
      p1 <- p1 + ggtitle( paste("Horizontal bar plot with counts and percentages, N =",pN), ) +
        
        
        theme(plot.title = element_text(size = 20, face = "bold")) +
        
        coord_flip()
      
      p1 <- p1 + ylab(ylabel ) + 
        
        xlab("Ordinal categories") +
        
        guides(fill=guide_legend(title=paste0("(",2,"-digit - ICD9 code)")), size = 14) 
      
      p1 <- p1 + geom_text(aes(label=paste0(format(N, big.mark=","
                                                   ,scientific=FALSE)," (",Percentage,"%)")),position = "stack", 
                           hjust=-0.2, size = 4.2, check_overlap = F)
      
      p1 <- p1 + scale_y_continuous(limits = c(0, mlimit)) 
      
      p1 <- p1 + theme(panel.background=element_blank(),
                       plot.title=element_text(), plot.margin = unit(c(5.5,12,5.5,5.5), "pt"), 
                       legend.text=element_text(size=12),
                       legend.title=element_text(size=14),
                       axis.text.x = element_text(size=13),
                       axis.text.y = element_text(size=15),
                       axis.line.x = element_line(color="black"),
                       axis.line.y = element_line(color="black"),
                       axis.title = element_text(size = 20) , 
                       plot.caption=element_text(hjust = 0, size = 7))
      
      g <- p1 + theme(legend.position="none") 
      
    }
    
    gx <- Gplotx(data = z,   l1=gupper,l2=glower,l3=gstep ) 
    
    print(gx)
    
  })
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # beta dist plot 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
  
  output$beta <- renderPlot({        
    
    sample <- random.sample()
    
    shape1. <- sample$shape1
    shape2. <- sample$shape2
    
    x_values <- seq(0,1, length.out = 1000)
    
    # require(ggplot2)
    # require(tidyverse)
    data.frame(x_values) %>%
      ggplot(aes(x_values))+
      stat_function(fun=dbeta, args=list(shape1=shape1.,shape2=shape2.)) +
      
      labs(title=paste0(c("Beta distribution, shape 1 =", shape1.,", shape 2 =", shape2.,""), collapse=" "), 
           x = "Latent underlying distribution of baseline version of response ",
           y = "Degree of belief",
           #subtitle =paste0(c("Note probabilites", prob," are equivalent to log odds: -4,-2, 0 ,2, 4 "), collapse=", "),
           caption = "") +
      guides(fill=FALSE) +
      theme_bw() +
      #  theme(legend.justification=c(1,0), legend.position=c(.96,.6)) +
      # scale_x_continuous("log odds", breaks=xs, labels=xs, limits=c(x1,x2)) +
      theme(legend.position="none") +
      theme(#panel.background=element_blank(),
        # axis.text.y=element_blank(),
        # axis.ticks.y=element_blank(),
        # https://stackoverflow.com/questions/46482846/ggplot2-x-axis-extreme-right-tick-label-clipped-after-insetting-legend
        # stop axis being clipped
        
        
        
        
        
        plot.title=element_text(size = 20, face = "bold"), plot.margin = unit(c(5.5,12,5.5,5.5), "pt"), 
        legend.text=element_text(size=12),
        legend.title=element_text(size=14),
        axis.text.x = element_text(size=13),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        axis.title = element_text(size = 20) , 
        plot.caption=element_text(hjust = 0, size = 7)
        
        
      )
    
    
    
  })
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   ggplot barplot of beta distribution
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
  
  output$reg.plotx <- renderPlot({         
    
    # Get the current regression data
    sample <- random.sample()
    levz <- sample$lev
    n   <- sample$n
    
    dat <- mcmc()$dat
    
    f <-   dat
    f <-   as.data.frame(table(f$baseline))
    
    f$Percentage <- round(f$Freq / sum(f$Freq)*100,1)
    
    
    z <- f              # data set for plot
    variable <- "Freq"  # variable of interest
    pN <- sum(f$Freq)   
    pN <- format(pN, big.mark=","
                 ,scientific=FALSE)
    roundUp <- function(x) 10^ceiling(log10(x))/2
    gupper <- roundUp((max(f$Freq)))  # plot upper limit
    gupper <- ceiling((max(f$Freq)))  # plot upper limit
    glower <- 0                       # plot lower limit
    gstep <- 5                        # grid steps
    
    # text for plot
    ylabel <- "Counts" 
    
    z$N <- z$Freq
    
    Gplotx <- function(data,  l1,l2,l3 ) {
      
      mlimit=l1
      
      p1 <- ggplot(data = data, aes(x =  Var1, y = N, fill = Var1)) + 
        
        geom_bar(stat = "identity", width =0.7) 
      
      p1 <- p1 + ggtitle( paste("Theorized dist. of baseline version of response, N =",pN), ) +
        theme(plot.title = element_text(size = 20, face = "bold")) #+
      
      #  coord_flip()
      
      p1 <- p1 + ylab(ylabel ) + 
        
        xlab("Ordinal categories") +
        
        guides(fill=guide_legend(title=paste0("(",2,"-digit - ICD9 code)")), size = 14) 
      
      p1 <- p1 + geom_text(aes(label=paste0(format(N, big.mark=","
                                                   ,scientific=FALSE)," (",Percentage,"%)")),position = "stack", 
                           vjust=-1.0,  hjust=.5, size = 3.1, check_overlap = F)
      
      p1 <- p1 + scale_y_continuous(limits = c(0, mlimit)) 
      
      p1 <- p1 + theme(panel.background=element_blank(),
                       plot.title=element_text(), plot.margin = unit(c(5.5,12,5.5,5.5), "pt"), 
                       legend.text=element_text(size=12),
                       legend.title=element_text(size=14),
                       axis.text.x = element_text(size=13),
                       axis.text.y = element_text(size=15),
                       axis.line.x = element_line(color="black"),
                       axis.line.y = element_line(color="black"),
                       axis.title = element_text(size = 20) , 
                       plot.caption=element_text(hjust = 0, size = 7))
      
      g <- p1 + theme(legend.position="none") 
      
    }
    
    gx <- Gplotx(data = z,   l1=gupper,l2=glower,l3=gstep ) 
    
    print(gx)
    
  })
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # text 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # side by side ggplot
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
  
  output$reg.plot99 <- renderPlot({         
    
    # Get the current regression data
    sample <- random.sample()
    levz <- sample$lev
    n   <- sample$n
    
    dat <- mcmc()$dat
    
    f <-   dat
    f <-   as.data.frame(table(f$y, f$treatment))
    
  #  f$Percentage <- round(f$Freq / sum(f$Freq)*100,1)
    
    
    library(dplyr)
    res <- group_by(f, Var2) %>% mutate(percent = 100*Freq/sum(Freq))
    
    res$Percentage <- round(res$percent,1)
    
    f <- res             # data set for plot
    variable <- "Freq"  # variable of interest
    pN <- sum(f$Freq)   
    pN <- format(pN, big.mark=","
                 ,scientific=FALSE)
    
    roundUp <- function(x) 10^ceiling(log10(x))/2
    gupper <- roundUp((max(f$Freq)))  # plot upper limit
    gupper <- ceiling((max(f$Freq)))*1.15  # plot upper limit
    glower <- 0                       # plot lower limit
    gstep <- 5                        # grid steps
    
    # text for plot
    ylabel <- "Counts" 
    
    f$N <- f$Freq
    
    z <- f

    NN <- tapply(z$Freq, z$Var2, sum)
    
    
    lab1 <- paste0("Placebo N = ",as.vector(NN[1]),"")
    lab2 <- paste0("Treatment N = ",as.vector(NN[2]),"")
    z$Var2 <- factor(z$Var2 , levels = c("0", "1"),
                      labels = c(lab1, lab2)
    )
    
    
    
    
    
    Gplotx <- function(data,  l1,l2,l3 ) {
      
      mlimit=l1
      
      p1 <- ggplot(data = data, aes(x =  Var1, y = N, fill = Var1)) + 
        
        geom_bar(stat = "identity", width =0.7) 
      
      p1 <- p1 + ggtitle( paste("Observed responses at follow up in trial arms, N =",pN), ) +
        theme(plot.title = element_text(size = 20, face = "bold")) #+
      
    
      p1 <- p1 + ylab(ylabel ) + 
        
        coord_flip() +
        
        xlab("Ordinal categories") +
        
        guides(fill=guide_legend(title=paste0("(",2,"-digit - ICD9 code)")), size = 14) 
    
      
      p1 <- p1 + geom_text(aes(label=paste0(format(N, big.mark=","
                                                   ,scientific=FALSE)," (",Percentage,"%)")),position = "stack", 
                           hjust=-0.2, size = 4.2, check_overlap = F)
      
      
      p1 <- p1 + labs(
              caption = "- Percentages calculated with respect to randomized group" 
            ) 
            
      
      p1 <- p1 + scale_y_continuous(limits = c(0, mlimit)) 
      
      p1 <- p1 + theme(panel.background=element_blank(),
                       plot.title=element_text(), plot.margin = unit(c(5.5,12,5.5,5.5), "pt"), 
                       legend.text=element_text(size=12),
                       legend.title=element_text(size=14),
                       axis.text.x = element_text(size=13),
                       axis.text.y = element_text(size=15),
                       axis.line.x = element_line(color="black"),
                       axis.line.y = element_line(color="black"),
                       axis.title = element_text(size = 20) , 
                       plot.caption=element_text(hjust = 0, size = 13),  #left align
                       strip.text = element_text(size=20)
                       )
      
      g <- p1 + theme(legend.position="none") +
        
        facet_wrap(Var2~.)
      
    }
    
    gx <- Gplotx(data = z,   l1=gupper,l2=glower,l3=gstep ) 
    
    print(gx)
    
  })
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # text 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
  
  output$textWithNumber <- renderText({ 
    
    A <- mcmc()$res     
    
    
    sample <- random.sample()
    levz <- sample$lev
    
    f <- A$coefficients
    x <-length(f) -2
    
    HTML(paste0( "Let's interpret the output on the left. The coefficient alongside y>=2 is "
                 , tags$span(style="color:red", p2( f     [1][[1]]) ) ,
                 " this is the log odds of having a response in categories 2 and above, so convert this to a probability "
                 , tags$span(style="color:red", p3(expit(A$coefficients[1][[1]]) )) , 
                 " and subtract from one to give the probability of being in the lowest category "
                 , tags$span(style="color:red", p3(1-  expit(f[1][[1]]) )) ,".",
                 br(), br(),  
                 
                 " The coefficient alongside y>=",levz," is "
                 , tags$span(style="color:red", p2( f     [x][[1]]) ) ,
                 " this is the log odds of having a response in the top category only, converting this to a probability gives "
                 , tags$span(style="color:red", p3(expit(f[x][[1]]) )) , 
                 "   "
                 , tags$span(style="color:red",  ) ,
                 
                 ""))    
    
  })
  
  
  output$textWithNumber1 <- renderText({ 
    
    A <- mcmc()$res     
    
    
  })
  
  
  
  
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # not used replaced by diff, ggplot 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # output$diffx <- renderPlot({         
  
  # z <- mcmc()$A
  # mu1 <- mcmc()$mu1
  # 
  # q <- quantile(z,c(.025, 0.25, 0.5, 0.75, 0.975))
  # par(bg = 'lightgoldenrodyellow') 
  # par(mfrow=c(1,3))
  # plot(density(z),
  #      xlab="Frequentist Bootstrap, Mean estimate",
  #      ylab="Density",
  #      main="",
  #      ylim=c(0,max(density(z)$y)),
  #      frame.plot=FALSE,cex.lab=1.5,lwd=3,yaxt="no")
  # abline(v=q[1], col="blue") #95% credible interval
  # abline(v=q[5], col="blue")
  # abline(v=q[3], col="red", lty='dashed')
  # abline(v=mu1, col="black", lty='dashed')
  # 
  # z <- (mcmc()$B)
  # 
  # q <- quantile(z,c(.025, 0.25, 0.5, 0.75, 0.975))
  # 
  # plot(density(z), #log="x",
  #      xlab="Bayesian Bootstrap, Mean estimate",
  #      ylab="Density",
  #      main="",
  #      ylim=c(0, max(density(z)$y)),##
  #      frame.plot=FALSE,cex.lab=1.5,lwd=3,yaxt="no")
  # abline(v=q[1], col="blue") #95% credible interval
  # abline(v=q[5], col="blue")
  # abline(v=q[3], col="red", lty='dashed')
  # abline(v=mu1, col="black", lty='dashed')
  # 
  # 
  # z <- mcmc()$C
  # 
  # q <- quantile(z,c(.025, 0.25, 0.5, 0.75, 0.975))
  # 
  # plot(density(z), #  log="x",
  #      xlab="Bayesian Bootstrap 2, Mean estimate",
  #      ylab="Density",
  #      main="",
  #      ylim=c(0, max(density(z)$y)),
  #      frame.plot=FALSE,cex.lab=1.5,lwd=3,yaxt="no")
  # abline(v=q[1], col="blue") #95% credible interval
  # abline(v=q[5], col="blue")
  # abline(v=q[3], col="red", lty='dashed')
  # abline(v=mu1, col="black", lty='dashed')
  # 
  # #caption=("xxxx")
  # 
  # par(mfrow=c(1,1))
  
  #})
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$reg.summary1 <- renderPrint({
    
    return( (mcmc()$f2 ))
    
  })
  output$reg.summary2 <- renderPrint({
    
    return( (mcmc()$res ))
    
  })
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # correlation from simulation
  output$reg.summary3 <- renderPrint({
    
    return(print(mcmc()$sf1, digits=4))
    
  })
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # correlation Efron data
  output$reg.summary4 <- renderPrint({
    
    return(print(lmx()$linear, digits=4))
    
  })
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # correlation user data
  output$reg.summary5 <- renderPrint({
    
    return(print(lmx()$an, digits=4))
    
  })
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # run the correlation analysis for tab 2
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  lmx <- reactive({
    
    #sample <- random.sample()
    
    dat <- mcmc()$dat
    
    dat$y <- as.numeric(as.character(dat$y))
    
    f <- lm(y ~treatment + baseline, data=dat)
    
    linear <- summary(f)
    
    
    d <<- datadist(dat)
    options(datadist="d")
    linear <- ols(y ~treatment + (baseline), data=dat)
    an <- anova(linear)
    
    return(list(linear=linear , an=an)) 
    
    #return(list(linear=linear ))
    
  })
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # NOT USED REPLACED BY FACET PLOT!!!!!!
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  #output$diff2 <- renderPlot({         
  # 
  #         f <- cor1()$f
  #         b  <- cor1()$b
  #         BB <-cor1()$BB
  #         xx1 <- cor1()$xx1
  #         
  #         ff <- qlogis(f)
  #         bb <- qlogis(b)
  #         BBB <- qlogis(BB)
  #         xx <- qlogis(xx1)
  #         
  #         all <- c(f,b,BB,xx1)
  #         minx <- min(all)
  #         maxx <- max(all)
  #         bz <- c(minx, maxx)
  #          
  # 
  #         a <- x <- c(0.001, 0.003,0.01, 0.05,seq(.1,0.9,0.1),.95,.98,.99,0.995,0.999, 0.9995,0.9999) 
  #         
  #         q <- qlogis(x)   
  #         
  #         lims <- unique(a[sapply(bz,function(x) which.min(abs(x-a)))])
  #         
  #         indx <- which(x %in% lims)
  #          i = indx[1]-1
  #          if (i %in% 0) {i=1} 
  #          j=  indx[2]+1
  #          limitz = c(q[i], q[j])
  #          breakz = q[i:j]
  #          labelz = x[i:j]
  #         
  #         require(ggplot2)
  #         x1 <- xlab("")
  #         est <- quantile(f, c(.025,.5,.975)) 
  #         ff <- as.data.frame(ff)
  #         
  #         pL1<- ggplot(data = ff, aes(x = ff)) + x1+   
  #             geom_histogram(bins = 100, fill = rainbow(100))+
  #             scale_x_continuous(limits =limitz,
  #                                 breaks= breakz,  # this is where the values go
  #                                 labels= labelz)   + 
  #           
  #           
  #             labs(title = paste("Frequentist bootstrap: Median",p3(est[2][[1]]),", 95%CI ("
  #                                , p3(est[1][[1]]) ,", ",  p3(est[3][[1]]) ,")") ) +
  #             theme_bw()  
  #         pL1 <- pL1 + theme(axis.line=element_blank(),
  #                            #axis.text.x=element_blank(),
  #                            #axis.text.y=element_blank(),
  #                            #axis.ticks=element_blank(),
  #                            #axis.title.x=element_blank(),
  #                            axis.text=element_text(size=14),
  #                            axis.title=element_text(size=12,face="bold"),
  #                            #axis.title.y=element_blank(),
  #                            # legend.position="none",
  #                            panel.background=element_blank(),
  #                            panel.border=element_blank(),
  #                            #panel.grid.major=element_blank(),
  #                            #panel.grid.minor=element_blank(),
  #                            # plot.background=element_blank())
  #                            #plot.margin = unit(c(1,1,1,1), "cm")
  #                            plot.title = element_text(size = 14)
  #                            
  #         )
  #         
  #         est <- quantile(b, c(.025,.5,.975))  
  #         bb <- as.data.frame(bb)
  #         pL2<- ggplot(data = bb, aes(x = bb)) +x1 +
  #             geom_histogram(bins = 100, fill = rainbow(100))+
  #           scale_x_continuous(limits =limitz,
  #                              breaks= breakz,  # this is where the values go
  #                              labels= labelz)   + 
  #             labs(title = paste("Bayesian bootstrap: Median",p3(est[2][[1]]),", 95%CI (", p3(est[1][[1]]) ,", ",  p3(est[3][[1]]) ,")") ) +
  #             theme_bw()  
  #         pL2 <- pL2 + theme(axis.line=element_blank(),
  #                            #axis.text.x=element_blank(),
  #                            #axis.text.y=element_blank(),
  #                            #axis.ticks=element_blank(),
  #                            #axis.title.x=element_blank(),
  #                            axis.text=element_text(size=14),
  #                            axis.title=element_text(size=12,face="bold"),
  #                            #axis.title.y=element_blank(),
  #                            # legend.position="none",
  #                            panel.background=element_blank(),
  #                            panel.border=element_blank(),
  #                            #panel.grid.major=element_blank(),
  #                            #panel.grid.minor=element_blank(),
  #                            # plot.background=element_blank())
  #                            #plot.margin = unit(c(1,1,1,1), "cm")
  #                            plot.title = element_text(size = 14)
  #         )
  #         
  #         est <- quantile(BB, c(.025,.5,.975))  
  #         BBB<- as.data.frame(BBB)
  #         pL3<- ggplot(data = BBB, aes(x = BBB)) +x1+ 
  #             geom_histogram(bins = 100, fill = rainbow(100))+
  #           scale_x_continuous(limits =limitz,
  #                              breaks= breakz,  # this is where the values go
  #                              labels= labelz)   + 
  #             labs(title = paste("LaplaceDemon bootstrap: Median",p3(est[2][[1]]),", 95%CI (", p3(est[1][[1]]) ,", ",  p3(est[3][[1]]) ,")") ) +
  #             theme_bw()  
  #         pL3 <- pL3 + theme(axis.line=element_blank(),
  #                            #axis.text.x=element_blank(),
  #                            #axis.text.y=element_blank(),
  #                            #axis.ticks=element_blank(),
  #                            #axis.title.x=element_blank(),
  #                            axis.text=element_text(size=14),
  #                            axis.title=element_text(size=12,face="bold"),
  #                            #axis.title.y=element_blank(),
  #                            # legend.position="none",
  #                            panel.background=element_blank(),
  #                            panel.border=element_blank(),
  #                            #panel.grid.major=element_blank(),
  #                            #panel.grid.minor=element_blank(),
  #                            # plot.background=element_blank())
  #                            #plot.margin = unit(c(1,1,1,1), "cm")
  #                            plot.title = element_text(size = 14)
  #         )
  #         
  #         
  #         x1 <- xlab(" ")
  #         est <- quantile(xx1, c(.025,.5,.975))  
  #         xx<- as.data.frame(xx)
  #         pL4<- ggplot(data = xx, aes(x = xx)) +x1+ 
  #             geom_histogram(bins = 100, fill = rainbow(100))+
  #           scale_x_continuous(limits =limitz,
  #                              breaks= breakz,  # this is where the values go
  #                              labels= labelz)   + 
  #             labs(title = paste("bayesboot bootstrap: Median",p3(est[2][[1]]),", 95%CI (", p3(est[1][[1]]) ,", ",  p3(est[3][[1]]) ,")") ) +
  #             theme_bw()  
  #         pL4 <- pL4 + theme(axis.line=element_blank(),
  #                            #axis.text.x=element_blank(),
  #                            #axis.text.y=element_blank(),
  #                            #axis.ticks=element_blank(),
  #                            #axis.title.x=element_blank(),
  #                            axis.text=element_text(size=14),
  #                            axis.title=element_text(size=12,face="bold"),
  #                            #axis.title.y=element_blank(),
  #                            # legend.position="none",
  #                            panel.background=element_blank(),
  #                            panel.border=element_blank(),
  #                            #panel.grid.major=element_blank(),
  #                            #panel.grid.minor=element_blank(),
  #                            # plot.background=element_blank())
  #                            #plot.margin = unit(c(1,1,1,1), "cm")
  #                            plot.title = element_text(size = 14)
  #         )
  # 
  #         gridExtra::grid.arrange(pL1,  pL2, pL3,  pL4, nrow=2) 
  
  #   })
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # correlation plot tab 2
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #output$diff3 <- renderPlot({         
  
  #       f <- cor1()$f
  #       b  <- cor1()$b
  #       BB <-cor1()$BB
  #       xx1 <- cor1()$xx1
  #       
  #       sample <- random.sample()
  #    
  #       r    <- sample$n2
  # 
  #     ff <-  (f)
  #     bb <- (b)
  #     BBB <- (BB)
  #     xx <- (xx1)
  #     foo <- cbind(ff,bb,BBB,xx)
  #     foo1 <- reshape::melt(foo)
  # 
  #     
  #     levels(foo1$X2)
  #     
  #   Ce <-  est <- quantile(f, c(.025,.5,.975)) 
  #   C <-  paste("Frequentist : Median",p3(est[2][[1]]),", 95%CI (", p3(est[1][[1]]) ,", ",  p3(est[3][[1]]) ,")")  
  # 
  #  Ae <-   est <- quantile(b, c(.025,.5,.975))  
  #   A <-  paste("Bayesian : Median",p3(est[2][[1]]),", 95%CI (", p3(est[1][[1]]) ,", ",  p3(est[3][[1]]) ,")") 
  # 
  #   Be<-  est <- quantile(BB, c(.025,.5,.975))  
  #   B <-  paste("LaplaceDemon : Median",p3(est[2][[1]]),", 95%CI (", p3(est[1][[1]]) ,", ",  p3(est[3][[1]]) ,")") 
  # 
  #  De<-    est <- quantile(xx1, c(.025,.5,.975))  
  #  D <-   paste("bayesboot : Median",p3(est[2][[1]]),", 95%CI (", p3(est[1][[1]]) ,", ",  p3(est[3][[1]]) ,")")  
  # 
  #  # make a dataset to add lines to ggplot facets
  #  
  #  dummy2 <- data.frame(X2=c(paste("",A),
  #                         paste("",B),
  #                         paste("",C),
  #                         paste("",D)
  #  ),  
  #                       q1 =  c(Ae[1], Be[1], Ce[1], De[1]),
  #                       q50 = c(Ae[2], Be[2], Ce[2], De[2]),
  #                       q3 =  c(Ae[3], Be[3], Ce[3], De[3])
  #                       )
  # 
  #  
  #      levels(foo1$X2) <- c(paste("",A),
  #                           paste("",B),
  #                           paste("",C),
  #                           paste("",D)
  #                           )
  #      
  #     p <- sort(c(-.9, .9,-.99,.99 ,.999,.9999,1,-.95, .95,.8,-.8,seq(-.6,.6,0.3))  )
  #     
  #    g0 <- ggplot(data=foo1, aes(x = value)) +#
  #       geom_vline(data = dummy2, aes(xintercept = q1,  colour="red", linetype = "dotdash")) +
  #       geom_vline(data = dummy2, aes(xintercept = q50, colour="red", linetype = "dotdash")) +
  #       geom_vline(data = dummy2, aes(xintercept = q3,  colour="red", linetype = "dotdash")) +
  #      geom_vline( aes(xintercept = r,  colour="black", linetype = "dash")) +
  #      
  #      geom_histogram( aes(y = ..density..), bins=100, colour="black" , fill=rainbow(400))+     ylab("")+
  #       geom_density(alpha = 0.1, fill = "red") +
  #       facet_wrap(X2~ .) 
  #     
  # 
  #       g0 <- g0  + scale_x_continuous(trans = atanh_trans()  ,
  #                                     breaks= p, 
  #                                     xlab("Correlation"),
  #                                     oob=discard) +
  #         scale_y_continuous(breaks = NULL) +
  # 
  #   theme_bw()  
  # 
  #   g0 <- g0 + theme(#axis.line=element_blank(),
  #                      #axis.text.x=element_blank(),
  #                      #axis.text.y=element_blank(),
  #                      #axis.ticks=element_blank(),
  #                      #axis.title.x=element_blank(),
  #                      axis.text=element_text(size=12),
  #                      axis.title=element_text(size=12,face="bold"),
  #                      #axis.title.y=element_blank(),
  #                      legend.position="none",
  #                      #anel.background=element_blank(),
  #                      #panel.grid.major=element_blank(),
  #                      #panel.grid.minor=element_blank(),
  #                      # plot.background=element_blank())
  #                      #plot.margin = unit(c(1,1,1,1), "cm")
  #                      plot.title = element_text(size = 16),
  #                    strip.text.x = element_text(size = 16, colour = "black", angle = 0),
  #                    strip.background = element_rect(fill="ivory"),
  #                    panel.border = element_blank(),
  #                    panel.grid.major = element_blank(), 
  #                    panel.grid.minor = element_blank(),
  #                    panel.background = element_blank(), 
  #                    axis.line = element_line(colour = "black")
  #                    )
  #   
  #   g0 <- g0 + labs(
  #     caption = "The dotted lines indicate the median, 2.5 and 97.5 percentiles, the red line is the true population value" 
  #   ) 
  #   
  #     print(g0)
  #     
  # })
  #     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #     # ruben's data
  #     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  #   ruben <- reactive({
  #       
  #       dye <- c(1.15, 1.7, 1.42, 1.38, 2.8, 4.7, 4.8, 1.41, 3.9)
  #       efp <- c(1.38, 1.72, 1.59, 1.47, 1.66, 3.45, 3.87, 1.31, 3.75)
  #       data.set <- data.frame(dye,efp)
  #       
  #       z1 <- cor.test(dye,efp)
  #       xx <- as.integer(as.vector(z1$parameter)+2) 
  #       z1 <- c( (xx), unlist(z1$estimate), unlist(z1$conf.int)[1:2])
  #       names(z1) <- c("N","Estimate","Lower","Upper")
  # 
  #       return(list(data.set=data.set, z2=z1)) 
  #       
  #     })
  # 
  #   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   # ruben's data analysis correlation plot tab 3
  #   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #     output$diff4 <- renderPlot({      
  #       
  #       sample <- random.sample()
  #       
  #       data.set <- ruben()$data.set
  #       
  #       reps <- sims <- sample$sims
  #       
  #       len <- length(data.set$dye)
  #       
  #         sboot <- function() {
  #           cor(data.set[sample(1:len, replace=T),])[1,2]
  #         }
  #         
  #         bboot <- function() {
  #           cov.wt(data.set, diff(c(0,sort(runif(len-1)),1)), cor=T)$cor[1,2]
  #         }
  #         
  #        A <- data.set$dye
  #        B <- data.set$efp
  #         
  #       X <- matrix(c(A, B), len, 2)
  #       colnames(X) <- c("dye","efp")
  #       BB <- BayesianBootstrap(X=X, n=sims,
  #                               Method=function(x,w) cov.wt(x, w, cor=TRUE)$cor[1,2]) 
  #        
  #       # Using the weighted correlation (corr) from the boot package.
  #       b4 <- bayesboot(data.set, corr, R = sims, use.weights = TRUE)
  #        
  #       
  #       f<-replicate(sims, sboot())
  #       b<-replicate(sims, bboot())
  #       BB <- unlist(BB)
  #       xx1 <- b4$V1
  #         
  #       ff <-  (f)
  #       bb <- (b)
  #       BBB <- (BB)
  #       xx <- (xx1)
  #       foo <- cbind(ff,bb,BBB,xx)
  #       foo1 <- reshape::melt(foo)
  #       
  #       levels(foo1$X2)
  #       
  #       Ce <-  est <- quantile(f, c(.025,.5,.975)) 
  #       C <-  paste("Frequentist : Median",p3(est[2][[1]]),", 95%CI (", p3(est[1][[1]]) ,", ",  p3(est[3][[1]]) ,")")  
  #       
  #       Ae <-   est <- quantile(b, c(.025,.5,.975))  
  #       A <-  paste("Bayesian : Median",p3(est[2][[1]]),", 95%CI (", p3(est[1][[1]]) ,", ",  p3(est[3][[1]]) ,")") 
  #       
  #       Be<-  est <- quantile(BB, c(.025,.5,.975))  
  #       B <-  paste("LaplaceDemon : Median",p3(est[2][[1]]),", 95%CI (", p3(est[1][[1]]) ,", ",  p3(est[3][[1]]) ,")") 
  #       
  #       De<-    est <- quantile(xx1, c(.025,.5,.975))  
  #       D <-   paste("bayesboot : Median",p3(est[2][[1]]),", 95%CI (", p3(est[1][[1]]) ,", ",  p3(est[3][[1]]) ,")")  
  #       
  #       # make a dataset to add lines to ggplot facets
  #       dummy2 <- data.frame(X2=c(paste("",A),
  #                                 paste("",B),
  #                                 paste("",C),
  #                                 paste("",D)
  #       ),  
  #       q1 =  c(Ae[1], Be[1], Ce[1], De[1]),
  #       q50 = c(Ae[2], Be[2], Ce[2], De[2]),
  #       q3 =  c(Ae[3], Be[3], Ce[3], De[3])
  #       )
  #       
  #       levels(foo1$X2) <- c(paste("",A),
  #                            paste("",B),
  #                            paste("",C),
  #                            paste("",D)
  #       )
  #       
  #       p <- c(-.9, -.99,.99 , -.95, .95, .999, .9999,c(-.8,-.4,0,.5,.8))  
  #       g0 <- ggplot(data=foo1, aes(x = value)) +#
  #         geom_vline(data = dummy2, aes(xintercept = q1,  colour="red", linetype = "dotdash")) +
  #         geom_vline(data = dummy2, aes(xintercept = q50, colour="red", linetype = "dotdash")) +
  #         geom_vline(data = dummy2, aes(xintercept = q3,  colour="red", linetype = "dotdash")) +
  #         geom_histogram(aes(y = ..density..), bins=100, colour="black" , fill=rainbow(400))+  ylab("")+
  #         geom_density(alpha = 0.1, fill = "red") +
  #         facet_wrap(X2~ .) 
  #       
  #        
  #       
  #       g0 <- g0  + scale_x_continuous(trans = atanh_trans()  ,
  #                                      breaks= p, xlab("Correlation"),
  #                                      oob=discard) +
  #         scale_y_continuous(breaks = NULL) +
  #         
  #         theme_bw()  
  #       
  #       g0 <- g0 + theme(#axis.line=element_blank(),
  #         #axis.text.x=element_blank(),
  #         #axis.text.y=element_blank(),
  #         #axis.ticks=element_blank(),
  #         #axis.title.x=element_blank(),
  #         axis.text=element_text(size=12),
  #         axis.title=element_text(size=12,face="bold"),
  #         #axis.title.y=element_blank(),
  #         legend.position="none",
  #         #anel.background=element_blank(),
  #         #panel.grid.major=element_blank(),
  #         #panel.grid.minor=element_blank(),
  #         # plot.background=element_blank())
  #         #plot.margin = unit(c(1,1,1,1), "cm")
  #         plot.title = element_text(size = 16),
  #         strip.text.x = element_text(size = 16, colour = "black", angle = 0),
  #         strip.background = element_rect(fill="ivory"),
  #         panel.border = element_blank(),
  #         panel.grid.major = element_blank(), 
  #         panel.grid.minor = element_blank(),
  #         panel.background = element_blank(), 
  #         axis.line = element_line(colour = "black")
  #       )
  #       
  #       g0 <- g0 + labs(
  #         caption = "The lines indicate the median, 2.5 and 97.5 percentiles" 
  #       ) 
  #       print(g0)
  #       
  #       
  #     })
  #     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #     # loading in user data, see my PBE bio equiv app
  #     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  #     output$contents <- renderTable({
  #       
  #       # input$file1 will be NULL initially. After the user selects
  #       # and uploads a file, head of that data file by default,
  #       # or all rows if selected, will be shown.
  #       
  #       df<- NULL
  #       req(input$file1)
  #       df <- read.csv(input$file1$datapath,
  #                      header = input$header,
  #                      sep = input$sep,
  #                      quote = input$quote)
  #       
  #       df <- as.data.frame(df)
  #       
  #       if(input$disp == "head") {
  #         
  #         return(head(df))
  #       }
  #       else {
  #         
  #         return(df)
  #       } 
  #       
  #     })
  #     
  #     
  #     
  #    usercor <- reactive({ 
  #     
  #      sample <- random.sample()
  #      reps <- sims <- sample$sims
  #      
  #      df<-NULL
  #      req(input$file1)
  #      df <- read.csv(input$file1$datapath,
  #                     header = input$header,
  #                     sep = input$sep,
  #                     quote = input$quote)
  #      
  #      df<- as.data.frame(df)
  #      
  #      
  #      names(df) <-c("A","B")
  #      
  #       z1 <- cor.test(df$A,df$B)
  #       xx <- as.integer(as.vector(z1$parameter)+2) 
  #       z1 <- c( (xx), unlist(z1$estimate), unlist(z1$conf.int)[1:2])
  #    
  #      names(z1) <- c("N","Estimate","Lower","Upper")
  #     
  #      return(list( z1=z1)) 
  #     
  #     })
  #     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #     # analyse user data
  #     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #     output$contents2 <- renderPlot({   
  #       
  #       sample <- random.sample()
  #       reps <- sims <- sample$sims
  #       
  #         df<-NULL
  #         req(input$file1)
  #         df <- read.csv(input$file1$datapath,
  #                        header = input$header,
  #                        sep = input$sep,
  #                        quote = input$quote)
  #         
  #         df<- as.data.frame(df)
  #         
  #         
  #         names(df) <-c("A","B")
  #         
  #          data.set <- df
  #         
  #         len <- length(data.set$A)
  #         
  #         sboot <- function() {
  #           cor(data.set[sample(1:len, replace=T),])[1,2]
  #         }
  #         
  #         bboot <- function() {
  #           cov.wt(data.set, diff(c(0,sort(runif(len-1)),1)), cor=T)$cor[1,2]
  #         }
  #         
  #         A <- data.set$A
  #         B <- data.set$B
  #         
  #         X <- matrix(c(A, B), len, 2)
  #         colnames(X) <- c("A","B")
  #         BB <- BayesianBootstrap(X=X, n=sims,
  #                                 Method=function(x,w) cov.wt(x, w, cor=TRUE)$cor[1,2]) 
  #         
  #         # Using the weighted correlation (corr) from the boot package.
  #         b4 <- bayesboot(data.set, corr, R = sims, use.weights = TRUE)
  #         
  #         
  #         f<-replicate(sims, sboot())
  #         b<-replicate(sims, bboot())
  #         BB <- unlist(BB)
  #         xx1 <- b4$V1
  #         
  #         ff <-  (f)
  #         bb <- (b)
  #         BBB <- (BB)
  #         xx <- (xx1)
  #         foo <- cbind(ff,bb,BBB,xx)
  #         foo1 <- reshape::melt(foo)
  #         
  #         levels(foo1$X2)
  #         
  #         Ce <-  est <- quantile(f, c(.025,.5,.975)) 
  #         C <-  paste("Frequentist : Median",p3(est[2][[1]]),", 95%CI (", p3(est[1][[1]]) ,", ",  p3(est[3][[1]]) ,")")  
  #         
  #         Ae <-   est <- quantile(b, c(.025,.5,.975))  
  #         A <-  paste("Bayesian : Median",p3(est[2][[1]]),", 95%CI (", p3(est[1][[1]]) ,", ",  p3(est[3][[1]]) ,")") 
  #         
  #         Be<-  est <- quantile(BB, c(.025,.5,.975))  
  #         B <-  paste("LaplaceDemon : Median",p3(est[2][[1]]),", 95%CI (", p3(est[1][[1]]) ,", ",  p3(est[3][[1]]) ,")") 
  #         
  #         De<-    est <- quantile(xx1, c(.025,.5,.975))  
  #         D <-   paste("bayesboot : Median",p3(est[2][[1]]),", 95%CI (", p3(est[1][[1]]) ,", ",  p3(est[3][[1]]) ,")")  
  #         
  #         # make a dataset to add lines to ggplot facets
  #         
  #         dummy2 <- data.frame(X2=c(paste("",A),
  #                                   paste("",B),
  #                                   paste("",C),
  #                                   paste("",D)
  #         ),  
  #         q1 =  c(Ae[1], Be[1], Ce[1], De[1]),
  #         q50 = c(Ae[2], Be[2], Ce[2], De[2]),
  #         q3 =  c(Ae[3], Be[3], Ce[3], De[3])
  #         )
  #         
  #         levels(foo1$X2) <- c(paste("",A),
  #                              paste("",B),
  #                              paste("",C),
  #                              paste("",D)
  #         )
  #         
  #         p <- c(-.9, -.99,.99 , -.95, .95, .999, .9999,c(-.8,-.4,0,.5,.8))  
  #         g0 <- ggplot(data=foo1, aes(x = value)) +#
  #           geom_vline(data = dummy2, aes(xintercept = q1,  colour="red", linetype = "dotdash")) +
  #           geom_vline(data = dummy2, aes(xintercept = q50, colour="red", linetype = "dotdash")) +
  #           geom_vline(data = dummy2, aes(xintercept = q3,  colour="red", linetype = "dotdash")) +
  #           geom_histogram(aes(y = ..density..), bins=100, colour="black" , fill=rainbow(400))+  ylab("")+
  #           geom_density(alpha = 0.1, fill = "red") +
  #           facet_wrap(X2~ .) 
  #         
  #         
  #         
  #         g0 <- g0  + scale_x_continuous(trans = atanh_trans()  ,
  #                                        breaks= p, xlab("Correlation"),
  #                                        oob=discard) +
  #           scale_y_continuous(breaks = NULL) +
  #           
  #           theme_bw()  
  #         
  #         g0 <- g0 + theme(#axis.line=element_blank(),
  #           #axis.text.x=element_blank(),
  #           #axis.text.y=element_blank(),
  #           #axis.ticks=element_blank(),
  #           #axis.title.x=element_blank(),
  #           axis.text=element_text(size=12),
  #           axis.title=element_text(size=12,face="bold"),
  #           #axis.title.y=element_blank(),
  #           legend.position="none",
  #           #anel.background=element_blank(),
  #           #panel.grid.major=element_blank(),
  #           #panel.grid.minor=element_blank(),
  #           # plot.background=element_blank())
  #           #plot.margin = unit(c(1,1,1,1), "cm")
  #           plot.title = element_text(size = 16),
  #           strip.text.x = element_text(size = 16, colour = "black", angle = 0),
  #           strip.background = element_rect(fill="ivory"),
  #           panel.border = element_blank(),
  #           panel.grid.major = element_blank(), 
  #           panel.grid.minor = element_blank(),
  #           panel.background = element_blank(), 
  #           axis.line = element_line(colour = "black")
  #         )
  #         
  #         g0 <- g0 + labs(
  #           caption = "The lines indicate the median, 2.5 and 97.5 percentiles" 
  #         ) 
  #         print(g0)
  
  
  #  })
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # print Efron's data
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  output$tablex <- DT::renderDataTable({
    # 
    # d <- ruben()$data.set
    # 
    # x <- d
    # 
    # foo <- x
    #  datatable(x,
    #            
    #            rownames = TRUE,
    #            
    #            options = list(
    #                searching = TRUE,
    #                #pageLength = 20,
    #                paging=FALSE,
    #                lengthMenu = FALSE ,
    #                lengthChange = FALSE,
    #                autoWidth = TRUE
    #                # colReorder = TRUE,
    #                # deferRender = TRUE,
    #                # scrollY = 200,
    #                # scroller = T
    #            ))  %>% 
    #      formatRound(
    #          columns= c("dye","efp"), digits=c(2,2)  ) 
    
  })
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
})

# Run the application 
shinyApp(ui = ui, server = server)