#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Rshiny ideas from on https://gallery.shinyapps.io/multi_regression/
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list=ls()) 
set.seed(333) # reproducible
library(directlabels)
library(shiny) 
library(shinyWidgets)
library(shinythemes)  # more funky looking apps
library(DT)
library(shinyalert)
library(Hmisc)
library(reshape)
library(rms)
library(ormPlot)
library(ordinal)
library(ggplot2)
library(tidyverse)
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
fig.widthx <- 593
fig.heightx <- 268
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
is.even <- function(x){ x %% 2 == 0 } # function to id. odd maybe useful
options(width=200)

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
                                             onclick ="window.open('https://raw.githubusercontent.com/eamonn2014/proportional-odds-model2/master/app.R', '_blank')"),    
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
                                            div(h5(tags$span(style="color:blue", "Total sample size"))), "1000"),
                                  
                                  tags$hr(),
                                  textInput('dist', 
                                            div(h5(tags$span(style="color:blue", "Approximate the distribution of the baseline version of the response by specifying
                                                             Beta shape parameters"))), "22,21"),
                                  
                                  textInput('levels', 
                                            div(h5(tags$span(style="color:blue", "Number of ordinal categories in response"))), "10"),
                                  tags$hr(), 
                                  textInput('or1', 
                                            div(h5(tags$span(style="color:blue", "Treatment odds ratio"))), "2"),
                                  
                                  
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
                                  tags$a(href = "https://rdrr.io/cran/rms/man/predict.lrm.html", tags$span(style="color:blue", "prediction of nodel mean"),),  
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
                              
                              
                              tabPanel("1 Baseline", value=7, 
                                       h4("The distribution of the baseline version of the response variable is specified here.
                                          By selecting a beta distribution using the shape parameters the
                                          expected baseline counts in categories can be approximated. The default is Beta(22,21)."),
                              
                                       
                                       ###############
                                    
                                       fluidRow(
                                         column(width = 6, offset = 0, style='padding:1px;',
                                                
                                                div(plotOutput("beta",  width=fig.width7, height=fig.height7)),
                                                 
                                         ) ,
                                        
                                         
                                         fluidRow(
                                           column(width = 5, offset = 0, style='padding:1px;',
                                            
                                                  div(plotOutput("reg.plotx",  width=fig.width7, height=fig.height7)) 
                                                  
                                           ))),
                                       h4(paste("Figures 1 & 2. Baseline distribution of outcome")), 
                                    
                              ) ,
                              
                              tabPanel("2 Outcome", value=3, 
                                            
                                       div(plotOutput("reg.plot99", width=fig.width1, height=fig.height1)),
                                       
                                       fluidRow(
                                         column(width = 7, offset = 0, style='padding:1px;',
                                                h4(paste("Figure 3. Observed responses, dictated by the user inputs on the left.")), 
                                           
                                         )),
                                       
                                       
                              ),
                              
                              #####
                              tabPanel("3 PO model", value=7, 
                                     
                                       ###############
                                       
                                       
                                       fluidRow(
                                         column(width = 6, offset = 0, style='padding:1px;',
                                                h4("Table 1 Proportional odds model"), 
                                                div( verbatimTextOutput("reg.summary1") )
                                         ) ,
                                         
                                 
                                         
                                         h4("Table 2 Proportional odds ratio summaries. Do we recover the input odds ratios...?"),
                                         fluidRow(
                                           column(width = 6, offset = 0, style='padding:1px;',

                                                 splitLayout(
                                                   textInput("bas1", div(h5("Enter a baseline low effect")), value="1", width=100),
                                                   textInput("bas2", div(h5("Enter a baseline high effect")),value="2", width=100)
                                                 ),

                                                 
                                                  div( verbatimTextOutput("reg.summary3")),

                                                  h4(htmlOutput("textWithNumber",) ),
                                           ))),
                                        
                              ) ,
                              
                             
                              #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                         
                              #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             
                              tabPanel("4 Predicted probability plot 1", value=3, 
                                       
                                       h5(paste("Enter 999 in the box below to see all the levels or enter level(s) of interest separated by a comma")), 
                                       textInput('rcat2', 
                                                 div(h5(tags$span(style="color:blue",
                                                                  ))), "999"),
                                       
                                       
                                       div(plotOutput("preds2", width=fig.width1, height=fig.height3)),
                                       
                                       
                                       
                                       fluidRow(
                                         column(width = 7, offset = 0, style='padding:1px;',
                                                h4(paste("Figure 4. Plot of the predicted probabilities")), 
                                                 
                                         )),
                              ),
                              #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                              tabPanel("5 Predicted prob. plot 2", 
                                       h4(paste("Figure 5 & 6. Plot of the predicted probabilities (reprise)")),
                                        
                                       h4("On the left you are looking at the point of view of what happens to a patient considering their baseline category. We can see the probability of response and how it depends on treatment. With the default inputs we can see a shift in the distribution to the higher categories if treated.  

On the right we can look at ALL baseline categories and see the predicted probability curves. 
Vertically all the curves will sum to 1 for a treatment group. 
For example, if a patient is in baseline group category 1 we can see the probability of the patient being in each category if they were treated (or alternatively if they were in placebo).
With the default inputs we can see horizontal lines in the treated responses (only for the default input values), telling us a patient's baseline value is not important to know.
"),
                                       fluidRow(
                                         column(width = 6, offset = 0, style='padding:1px;',
                                                
                                                div(plotOutput("preds", width=fig.width7, height=fig.height3)),
                                                
                                                fluidRow(
                                                  
                                                  textInput('base', 
                                                            div(h5(tags$span(style="color:blue", 
                                                                             "Enter a patient's baseline category and see their predicted probabilities for response outcomes"))), "1")
                                                  
                                                  
                                                ),
                                         ) ,
                                         
                                         fluidRow(
                                           
                                          
                                           column(width = 5, offset = 0, style='padding:1px;',
                                                 
                                                  div(plotOutput("predicts", width=fig.width7, height=fig.height3)),
                                                  
                                                  
                                                  
                                                  fluidRow(
                                                    
                                                    textInput('group', 
                                                              div(h5(tags$span(style="color:blue", 
                                                                               "select treatment group: 0 for placebo, 1 for treatment, 2 for both"))), "1"),
                                                    
                                                    textInput('rcat', 
                                                              div(h5(tags$span(style="color:blue", 
                                                                               "Response category, enter 999 to see all levels or enter level(s) of interest"))), "999"),
                                                    #""
                                                    
                                                    
                                                    
                                                  ),
                                               
                                           ))),
                                       
                                        
                                       width = 30 )     ,
                               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                              tabPanel("6 Tables of probabilities",
                                       h4(paste("Table 3 Predicted probabilities, the estimated mean Y (meanY) is calculated by summing values of Y multiplied by the estimated Prob(Y=j)")),
                                       fluidRow(
                                         column(width = 12, offset = 0, style='padding:1px;',
                                       
                                             div( verbatimTextOutput("reg.summaryp") ),
                                             h4(paste("Table 4 Predicted cummulative probabilities ")),
                                             div( verbatimTextOutput("reg.summaryc") ),
                                         ) ,
                                    
                                         
                                         ),

                              ),
                              
                             tabPanel("7 Linear model", value=3, 
                                       h4(" ANCOVA model Tables 5 & 6 and Figure 7"),
                                      
                                       
                                       fluidRow(
                                         column(width = 6, offset = 0, style='padding:1px;',
                                                
                                                div( verbatimTextOutput("reg.summary4") )
                                         ) ,
                                         
                                         fluidRow(
                                           column(width = 5, offset = 0, style='padding:1px;',
                                                  
                                                  div( verbatimTextOutput("reg.summary5")),
                                                  div(plotOutput("predictl", width=fig.widthx, height=fig.heightx)),
                                              
                                           ))),

                             ),
                                      
                                      tabPanel("8 Predicted Mean", value=3, 
                                            
                                               fluidRow(
                                                 column(width = 6, offset = 0, style='padding:1px;',
                                                        
                                                        textInput('kints',
                                                                  div(h5(tags$span(style="color:blue",
                                                                                   "Enter an intercept for the ordinal model"))), ""), 
                                                        
                                                        div(plotOutput("PP.plot", width=fig.width7, height=fig.height6)),
                                                        h4("Figure 8 Predictions for each trial arm by model"),
                                                        br() , 
                                                        
                                                        
                                                        #h5("For the ordinal model the estimated mean Y is calculated by summing values of 
                                                        #Y multiplied by the estimated Prob(Y=j)."),
                                                        
                                                        
                                                        h4("Table 7 Model predictions"),
                                                        div( verbatimTextOutput("predz"), width = 2), # 
                                                        
                                                        
                                                        
                                                 ),
                                                 
                                                 fluidRow(
                                                   
                                                   br(), br(), br() , br() , br() ,
                                                   
                                                   
                                                   column(width = 5, offset = 0, style='padding:0px;',
                                                          
                                                          #div(h5(tags$span(style="color:blue","test"))),  
                                                          
                                                          
                                                          div(plotOutput("PP.plot2", width=fig.width7, height=fig.height6)),
                                                          h4("Figure 9 Predictions for each model arm by trial arm to assess similarity in the two model predictions"),
                                                           
                                                   )))
                                        
                              ) ,
   
                                tabPanel("9 Data", 
                                         
                                         fluidRow(
                                           column(width = 3, offset = 0, style='padding:1px;',
                                       h4("Table 8 Data listing"),
                                       div( verbatimTextOutput("dat")),
                                           ),
                                       
                                       column(width = 9, offset = 0, style='padding:1px;',
                                              h4("Notes"),
                                              h6("We fit the baseline response as a continuous variable in the model\n"),
                                              h6("   To do check assumptions\n"),
                                              h6("    add references here to allow use of space on landing page"),
                                              h6("    feedback on ormfit not working as expected"),
                                              h6("    calculate means when intercept changes"),
                                       )
                                       
                                       
                                       )
                              )
                              
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
             "Use at your own risk", 
             type = "info")
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This is where a new sample is instigated 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  random.sample <- reactive({
    
    foo <- input$resample
    
    dis <- as.numeric(unlist(strsplit(input$dist,",")))
    
    trt <- as.numeric(unlist(strsplit(input$n,",")))
    
    ctr <- as.numeric(unlist(strsplit(input$levels,",")))
    
    #sample size for correlation
    n1y1 <- log(as.numeric(unlist(strsplit(input$or1,","))))   # user enter odds , need log for the maths
    # R
    n2y2 <- log(as.numeric(unlist(strsplit(input$or2,","))))    # user enter odds , need log for the maths
    
    
    base<- as.numeric(unlist(strsplit(input$base,",")))
     
    
    return(list(  
      n=trt[1],  
      lev=ctr[1],
      or1=n1y1[1], 
      or2=n2y2[1],
      shape1=dis[1], 
      shape2=dis[2],
      base=base[1]
   
      
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
    group  <- sample$group
    rcat  <- sample$rcat
    bas1  <- sample$bas1
    bas2  <- sample$bas2
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Parameters 
    
    # treatment assignment 1:1 
    treatment <- 1*(runif(n)<0.5)  
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
    #      sf1 <- summary(f1, antilog=TRUE, verbose=FALSE)
    
    return(list(  dat=dat )) 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  })
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~DO THE ANALYSIS~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  analysis <- reactive({
   
    bas1  <-  as.numeric(unlist(input$bas1))
    bas2  <-  as.numeric(unlist(input$bas2))
 
    dat <- mcmc()$dat
    kk <-   ( as.numeric(unlist(strsplit(input$kints,","))))
    dat$y <- as.numeric(as.character(dat$y)) 
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    options(datadist=NULL)  
    d <<- datadist(dat)              # 
    d$limits["Low:effect","baseline"]  <- bas1  #reset adjustment level to 
    d$limits["High:effect","baseline"] <- bas2  #reset adjustment level to 
    d <<- datadist(dat)  
    options(datadist="d")
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    sf1=NULL
    f2 <- orm(y ~baseline + treatment, data=dat )
    f3 <- ols(y ~baseline + treatment, data=dat )
    sf1 <- summary(f2, antilog=TRUE, verbose=FALSE)
    sf1 <- summary(f2, baseline=c(bas1,bas2),antilog=TRUE, verbose=FALSE)
 
    k <- NULL
    m <- Mean(f2, codes=FALSE)
    ml <- as.list(m)
    k <- ml$interceptRef
    
    ols.<- Predict(f3,  conf.type="mean",  baseline, treatment)
    
    # if input is empty do this 
    if(!isTruthy(kk)){
      
      K <- k
      orm. <-  Predict(f2, baseline, treatment, fun=m, kint=K )
      
    } else {
      
      # if there is a value use it.
      K <- kk
      orm. <-  Predict(f2,   baseline, treatment,fun=m, kint=K)
      
    }
    
    P  <- rbind( "Proportional odds model"=orm., "Ordinary least squares"=ols.)
    P2 <- rbind( "model"=orm., "model"=ols.)
 
    return(list( ols.=ols., orm.=orm. , kk=kk , P2=P2, k=k, K=K,dat=dat, m=m, f2=f2, f3=f3, P=P, sf1=sf1,d=d )) 
  })
   
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # beta dist plot 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
  
  output$beta <- renderPlot({        
    
    sample <- random.sample()
    
    shape1. <- sample$shape1
    shape2. <- sample$shape2
    
    x_values <- seq(0,1, length.out = 1000)
    
 
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
  #  end ggplot barplot of beta distribution
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
      
      p1 <- p1 + ggtitle( paste0("Observed dist. of baseline version of response, N=",pN), ) +
        theme(plot.title = element_text(size = 20, face = "bold")) #+
      
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
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # non cummulative predicted probabilities plot run the analysis again
  # not efficient I know
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  output$preds2 <- renderPlot({
    
    sample <- random.sample()
    # n    <- sample$n
    levz <- sample$lev
    base <- sample$base
    
    datx <- mcmc()$dat
    
    rcat <-  (as.numeric(unlist(strsplit(input$rcat2,","))))   
    
    # I can get non cummulative probabilites using clm
    # Don't know how to do it with rms? yet
    
    Res.clm <- clm(y ~treatment + baseline, data=datx)
   
    
    newdat <- data.frame(
      baseline =   (rep(1:levz)),
      treatment = rep(0:1, each = levz)
    )
    
    
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
    
    
    l$response <- as.numeric(l$response)
    l$treatment <- factor( l$treatment)
    l$treatment <- ifelse( l$treatment %in% 0,"Placebo","Treatment")
    l$response <- as.factor(l$response)
    
    br1 <- length(unique(l$baseline))
    
    
    if (rcat %in% 999) {r = 1:levz} else   {r = rcat} 
    
    l <-  l[l$response %in% r,]
   
    
  gp <-  ggplot(l, aes(x = baseline,  y = estimate, color = response)) + 
      geom_line(size = 1) + 
      geom_ribbon(aes(ymin = lower,   ymax = upper,
                      fill = response, color = response),  alpha = 0.4, linetype = 0) +
    theme_bw() + 
    scale_x_continuous( breaks=1:br1, labels=1:br1)+
      facet_grid(~treatment) +
    
    
  theme(panel.background=element_blank(),
        plot.title=element_text(), plot.margin = unit(c(5.5,12,5.5,5.5), "pt"), 
        legend.text=element_text(size=12),
        legend.title=element_text(size=14),
       # legend.title=element_blank(),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        axis.title.y=element_text(size=16),  
        axis.title.x=element_text(size=16),  
        axis.title = element_text(size = 20) , 
        plot.caption=element_text(hjust = 0, size = 7),
       strip.text.x = element_text(size = 14, colour = "black", angle = 0)
       
       ) +
 
    
    labs(title=paste0(c("Predicted probabilities of response categories"), collapse=" "), 
         x = "Baseline category",
         y = "Predicted probability",
        # subtitle =c("xxxxxxxxxxxxxx"),
         caption = "")  
    
    
    print(gp)
    
  })
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # non cummulative predicted probabilities plot run the analysis again
  # not efficient I know
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  output$preds <- renderPlot({
    
    sample <- random.sample()
   # n    <- sample$n
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
            plot.caption=element_text(hjust = 0, size = 7) ,
            
            legend.justification = c(0, 1), 
            legend.position = c(0.05, .99))  +
      
      
      
      labs(title=paste0(c("Predicted probabilities of response categories"), collapse=" "), 
           x = "Response category",
           y = "Predicted probability",
           #subtitle =c("xxxxxxxxxxxxxx"),
           caption = "")  
    # guides(fill=guide_legend(title="Treatment"))
    # 
    
    print(gp)
    
  })

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~tables of predictions
  
  predictz <- reactive({  
    
    sample <- random.sample()
    f    <- analysis()$f2
    levz <- sample$lev
    
    group <- (as.numeric(unlist(strsplit(input$group,","))))    
    rcat <-  (as.numeric(unlist(strsplit(input$rcat,","))))     
    
    require(reshape)
    
    newdat <- data.frame(
      baseline = rep(1:levz),
      treatment = rep(0:1, each = levz))
    
    xx <- predict(f, newdat, type="fitted.ind")    
    
    probs <- cbind(newdat,xx )   
    
    # adding mean to right of dataframe
    d1 <- probs[-1:-2] 
    v1 <- 1:length(names(probs)[-1:-2])
    d2 <- t(t(d1)*v1)
    meanY <- apply(d2,1,sum)
    probs <- cbind(probs, meanY)
    # end
    
    xx <- predict(f, newdat, type="fitted") 
    
    cprobs <- cbind(newdat,xx )
    
    p.with.ci <- predict_with_ci(f, np = 100, fun = stats::plogis)
    
    plotci <- plot(f, baseline, treatment, fun = stats::plogis)
 
    return(list(probs=probs, cprobs=cprobs, p.with.ci=p.with.ci , plotci=plotci$data))
    
  })
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plots of predictions
  
  output$predicts <- renderPlot({   
    
    sample <- random.sample()
    # 
    f    <- analysis()$f2
    # 
    levz <- sample$lev
    
    
    group <- (as.numeric(unlist(strsplit(input$group,","))))    
    # R
    rcat <- (as.numeric(unlist(strsplit(input$rcat,","))))     
    
    require(reshape)
    
    newdat <- data.frame(
      baseline = rep(1:levz),
      treatment = rep(0:1, each = levz))
    
    xx <- predict(f, newdat, type="fitted.ind")    #
    
    mm <- melt(data.frame(xx))
    
    mm <- cbind(newdat,mm )
    
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    mm$variable <-  gsub(".*\\.","", mm$variable)
    
    mm <- plyr::arrange(mm,treatment,variable,baseline )
    
    mm$flag <- rep(seq_along( rle(mm$variable)$values ), times = rle(mm$variable)$lengths )
    
    
    if (group %in% 0) {g = c(0)} else if (group %in% 1) {g = c(1)} else if (group %in% 2) {g = c(0,1)} 
    
    mm <-  mm[mm$treatment %in% g,]
    
    
    if (rcat %in% 999) {r = 1:levz} else   {r = rcat} 
    
    mm <-  mm[mm$variable %in% r,]
    
    A <- ifelse(mm$treatment %in% 0, "Placebo","Treatment")
    mm$grp <- paste(A, mm$variable)
    
    A <- ifelse(mm$treatment %in% 0, "Pl.","Trt.")
    mm$var2 <- paste(A, mm$variable)
    
     
    gpp <- ggplot(mm, aes(baseline, value, group=factor(grp))) +
      geom_line(aes(color=factor(A))) +
      
      scale_x_continuous( breaks=1:levz, labels=1:levz) +  
      
      theme(panel.background=element_blank(),
            plot.title=element_text(), plot.margin = unit(c(5.5,12,5.5,5.5), "pt"), 
            legend.text=element_text(size=12),
            legend.title=element_text(size=0),
            #legend.title=element_blank(),
            axis.text.x = element_text(size=10),
            axis.text.y = element_text(size=10),
            axis.line.x = element_line(color="black"),
            axis.line.y = element_line(color="black"),
            axis.title.y=element_text(size=16),  
            axis.title.x=element_text(size=16),  
            axis.title = element_text(size = 20) , 
            plot.caption=element_text(hjust = 0, size = 7) ,
            legend.position="none") +
      
      labs(title=paste0(c("Predicted probabilities of response categories"), collapse=" "), 
           x = "Baseline category",
           y = "Predicted probability",
           #subtitle =c("xxxxxxxxxxxxxx"),
           caption = "")  +
      geom_dl(aes(label = var2), method = list(dl.combine("first.points", "last.points"),
                                               cex = 0.9)) 
    # guides(fill=guide_legend(title="Treatment"))  
    # 
    
    print(gpp)
    
    
  })
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # text 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
  
  output$PP.plot <- renderPlot({   
    
    K <- analysis()$K
    
    txt <- paste0("Ordinal intercept ", K)
    
    levz <- input$levels
    
    P <- analysis()$P
    
    P$treatment <- ifelse(P$treatment %in% 0, "Placebo", "Treatment")
     
    
    ggplot(P ,  ylab='' ) +  
      scale_x_continuous( breaks=1:levz, labels=1:levz) +  
         theme(panel.background=element_blank(),
            plot.title=element_text(), plot.margin = unit(c(5.5,12,5.5,5.5), "pt"), 
            legend.text=element_text(size=12),
            legend.title=element_text(size=0),
            #legend.title=element_blank(),
            axis.text.x = element_text(size=10),
            axis.text.y = element_text(size=10),
            axis.line.x = element_line(color="black"),
            axis.line.y = element_line(color="black"),
            axis.title.y=element_text(size=16),  
            axis.title.x=element_text(size=16),  
            axis.title = element_text(size = 20) , 
            plot.caption=element_text(hjust = 0, size = 7),
            
            legend.justification = c(0, 1), 
            legend.position = c(0.05, .99))  +
      
      # legend.position="none") +
      
      
      labs(title="Comparing estimates between models", 
           x = "Baseline category",
           y = "Predicted Mean",
           subtitle =txt,
           caption = "- For the ordinal model the estimated mean Y is calculated by summing \nvalues of Y multiplied by the estimated Prob(Y=j).")
    
    
  }) 
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~baseline plots~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  output$PP.plot2 <- renderPlot({   
    
    K <- analysis()$K
    
    txt <- paste0("Ordinal intercept ", K)
    
    levz <- input$levels
    
    P <- analysis()$P
    
    P$treatment <- ifelse(P$treatment %in% 0, "Placebo", "Treatment")
    
     
    ggplot(P,  groups=".set.", ylab='' ) +  
      
      scale_x_continuous( breaks=1:levz, labels=1:levz) +  
      
       
      theme(panel.background=element_blank(),
            plot.title=element_text(), plot.margin = unit(c(5.5,12,5.5,5.5), "pt"), 
            legend.text=element_text(size=12),
            legend.title=element_text(size=0),
            #legend.title=element_blank(),
            axis.text.x = element_text(size=10),
            axis.text.y = element_text(size=10),
            axis.line.x = element_line(color="black"),
            axis.line.y = element_line(color="black"),
            axis.title.y=element_text(size=16),  
            axis.title.x=element_text(size=16),  
            axis.title = element_text(size = 20) , 
            plot.caption=element_text(hjust = 0, size = 7),
            
            legend.justification = c(0, 1), 
            legend.position = c(0.05, .99))  +
      
      # legend.position="none") +
      
      
      labs(title="Same data, comparing estimates", 
           x = "Baseline category",
           y = "Predicted Mean",
           subtitle =txt,
           caption = "")
    
    
  })
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~baseline predictions~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  predz <- reactive({
    
    dat <- mcmc()$dat
     
    
    ols. <- analysis()$ols.
    orm. <- analysis()$orm.
   
    p  <- rbind(ols=ols., orm=orm.)
    
    return(list( p=p )) 
  })  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # text 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
  
  output$textWithNumber <- renderText({ 
    
    A <- analysis()$f2     
    
    d <- analysis()$d
    
    dat <- mcmc()$dat
    
    dat$y <- as.numeric(as.character(dat$y))
    
    sample <- random.sample()
    levz <- sample$lev
    
    f <- A$coefficients
    x <-length(f) -2
    
    diff <- d$limits["High:effect","baseline"]-d$limits["Low:effect","baseline"]
    
    c1  <-expit(f[1][[1]]  + f['baseline'][[1]])
    
    c10 <- expit(f[x][[1]]  + f['baseline'][[1]])
    
    HTML(paste0( "Let's interpret the output on the left. The coefficient alongside y>=",min(dat$y)+1," is "
                 , tags$span(style="color:red", p2( f     [1][[1]]) ) ,
                 " this is the log odds of having a response in categories ",min(dat$y)+1," and above, when treatment and baseline are 0. 
                So convert this to a probability after adding the contribution on the log scale of being in baseline category 1, for example gives "
                 , tags$span(style="color:red", p3(c1)) , 
                 " and subtract from one to give the probability of being in the lowest category "
                 , tags$span(style="color:red", p3(1-  c1 )) ," for a patient in reference treatment group (placebo).",
                 br(), br(),  
                 
                 " The coefficient alongside y>=",max(dat$y)," is "
                 , tags$span(style="color:red", p2( f     [x][[1]]) ) ,
                 ", this is the log odds of having a response in the top category only, converting this to a probability 
                 after adding on the log scale the contribution of being in baseline category 1, for example gives "
                 , tags$span(style="color:red", p3(c10)) , 
                 " for a patient in reference treatment group (placebo)."
                 , tags$span(style="color:red"),  ". Check these probabilities to the top left cell and top right cell of Table 3, tab 6 !" ,
                 br(), br(),  
                 
                 " The coefficient left, 'baseline' is a log odds ratio
                 comparing an individual one baseline category higher whilst being identical in all other predictors and is "
                 , tags$span(style="color:red", p3(f['baseline'])[[1]])  , 
                  ", we can exponentiate this to give "
                     , tags$span(style="color:red", p3(exp(f['baseline']))[[1]])  ,
                 " . We can see the effect of changing multiple categories on the outcome by multiplying " 
                  , tags$span(style="color:red", p3(f['baseline'])[[1]]),  
                " unit change by the 'Diff.'= ",
                 tags$span(style="color:red", diff ) 
                 ," which is selectable above and exponentiating...to give..." 
                 , tags$span(style="color:red", p3(exp(f['baseline']*diff) ) [[1]] )  , 
                 " and so on.",
                br(), br(),  
                " The coefficient left, 'treatment' is a log odds ratio
                 comparing an individual in the treated group to the placebo group whilst 
                being identical in all other predictors and is "
                , tags$span(style="color:red", p3(f['treatment'])[[1]])  , 
                ", we can exponentiate this to give the odds ratio: "
                , tags$span(style="color:red", p3(exp(f['treatment']))[[1]])     , 
                " "
                
                ))    
                
  })
  
  
  output$textWithNumber1 <- renderText({ 
    
    A <- analysis()$f2     
    
    
  })
  
  output$dat <- renderPrint({
    
    d <- mcmc()$dat
    
    d <- plyr::arrange(d, baseline, treatment)
    
    return(print(d, digits=4))
  })
  
  
  output$predz <- renderPrint({
    
    return(print(predz()$p, digits=4))
  })
  
  output$predt <- renderPrint({
    
    return(print(predt()$pt, digits=4))
  })
  
  
  output$reg.summary1 <- renderPrint({
    
    return( (analysis()$f2 ))
    
  })
  
  output$reg.summary3 <- renderPrint({
    
    return(print(analysis()$sf1, digits=4))
    
  })
  
  output$reg.summary4 <- renderPrint({
    
    return(print(lmx()$linear, digits=4))
    
  })

  output$reg.summary5 <- renderPrint({
    
    return(print(lmx()$an, digits=4))
    
  })
  
  output$reg.summaryp <- renderPrint({
    
    return(print(predictz()$prob, digits=4))
    
  })
  
  output$reg.summaryc <- renderPrint({
    
    return(print(predictz()$cprob, digits=4))
    
  })
  
  output$reg.summaryci <- renderPrint({
    
    return(print(predictz()$plotci, digits=4))
    
  })
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  lmx <- reactive({
    
     
    dat <- mcmc()$dat
    
    dat$y <- as.numeric(as.character(dat$y))
     
    d <<- datadist(dat)
    options(datadist="d")
    
    f <- ols(y ~treatment + (baseline), data=dat)
    an <- anova(f)
    
    return(list(linear=f , an=an )) 
     
    
  })
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  output$predictl <- renderPlot({   
    
    dat <- mcmc()$dat
    
    dat$y <- as.numeric(as.character(dat$y))
    
    linear <- ols(y ~treatment + (baseline), data=dat)
    
    an1 <- anova(linear)
     
    ggplot(Predict(linear, treatment), anova=an1, pval=TRUE)  
    
   
  })
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
})

# Run the application 
shinyApp(ui = ui, server = server)