#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Rshiny ideas from on https://gallery.shinyapps.io/multi_regression/
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list=ls()) 
set.seed(333) # reproducible
library(directlabels)
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

options(width=200)


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
                                  
                                   
                                  # tabPanel("4 Predicted Mean",
                                  # 
                                  #          # div(plotOutput("predictm", width=fig.width9, height=fig.height9)),
                                  # 
                                  #          h4("Enter an intercept, default is the middle intercept, seems to correspond to ['median Y' value-1] on tab2 model output, and then find which category this refers to in the printed 'Frequencies of Responses'.
                                  #         This will closely resemble the linear regression prediction.")  ,
                                  # 
                                  #          textInput('kints',
                                  #                    div(h5(tags$span(style="color:blue",
                                  #                                     "test"))), ""),
                                  #          div(plotOutput("PP.plot", width=fig.width6, height=fig.height6)),
                                  # 
                                  # 
                                  # 
                                  # ),
                                  
                                  
                                  
                                  
                                #####
                                
                                
                                tabPanel("4 Predicted Mean b", value=3, 
                                         
                                        
                                         
                                         fluidRow(
                                           column(width = 6, offset = 0, style='padding:1px;',
                                                  
                                                  textInput('kints',
                                                            div(h5(tags$span(style="color:blue",
                                                                             "test"))), ""), 
                                                 
                                                  div(plotOutput("PP.plot", width=fig.width7, height=fig.height6))),
                                           
                                           
                                           fluidRow(
                                             column(width = 5, offset = 0, style='padding:1px;',
                                                    div( verbatimTextOutput("predt") ), # 
                                                    div( verbatimTextOutput("preds") ), # 
                                                  
                                                 
                                             )))
                                         
                                         
                                         
                                         
                                ),
                                
                                #####
                                  
                                tabPanel("4a Predicted treat", value=3, 
                                         
                                         
                                         
                                         fluidRow(
                                           column(width = 6, offset = 0, style='padding:1px;',
                                                  
                                              #    textInput('kints',
                                               #             div(h5(tags$span(style="color:blue",
                                                #                             "test"))), ""), 
                                                  
                                               div(plotOutput("r.plot", width=fig.width7, height=fig.height6))),
                                           
                                           
                                           fluidRow(
                                             column(width = 5, offset = 0, style='padding:1px;',
                                                  #  div( verbatimTextOutput("predt") ), # 
                                                   # div( verbatimTextOutput("preds") ), # 
                                                    
                                                    
                                             )))
                                         
                                         
                                         
                                         
                                ),
                                  
                                  
                                  
                                  
                                  
                                  
                                  tabPanel("5 data", 
                                           
                                        
                                          
                                           div( verbatimTextOutput("dat")),
                                           
                                           
                                           
                                  ),
                               
                                  
                                  tabPanel("5 predictions", 
                                           
                                           
                                           
                                         # div( verbatimTextOutput("preds")),
                                           
                                           
                                           
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
               "It's in progress! https://rdrr.io/cran/rms/man/residuals.lrm.html baseline categorical continuous", 
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
      
        
        return(list(  
            n=trt[1],  
            lev=ctr[1],
            or1=n1y1[1], 
            or2=n2y2[1],
            shape1=dis[1], 
            shape2=dis[2]
           
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
      
      return(list(  sf1=sf1 , dat=dat )) 
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    })
  
    
    
    analysis <- reactive({
      
      dat <- mcmc()$dat
      kk <-   ( as.numeric(unlist(strsplit(input$kints,","))))
      dat$y <- as.numeric(as.character(dat$y)) 
      
      # dat$baseline <- factor(dat$baseline)
      # use harrell's po function analyse the data
      d <<- datadist(dat)
      options(datadist="d") 
      
      f2 <- orm(y ~baseline + treatment, data=dat )
      f3 <- ols(y ~baseline + treatment, data=dat )
      
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
      
      return(list( ols.=ols., orm.=orm. , kk=kk  , k=k, K=K,dat=dat, m=m, f2=f2, f2=f3)) 
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    })
    
     
      output$PP.plot <- renderPlot({   

 
        ols. <- analysis()$ols.
        orm. <- analysis()$orm.
        
        K <- analysis()$K
        
        txt <- paste0("Ordinal intercept ", K)
        
        
        P  <- rbind( "Proportional odds model"=orm., "Ordinary least squares"=ols.)
        P$treatment <- ifelse(P$treatment %in% 0, "Placebo", "Treatment")
        
        ggplot(P , ylab='' ) +  
          
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
          
      
          labs(title=txt, 
               x = "Baseline category",
               y = "Predicted Mean",
               subtitle =c("xxxxxxxxxxxxxx"),
               caption = "")
        
        
      })
    
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   

      preds <- reactive({
        
        dat <- mcmc()$dat
    
        # d <<- datadist(dat)
        # options(datadist="d") 
        
        ols. <- analysis()$ols.
        orm. <- analysis()$orm.
        #K <- analysis()$K
        #m <- analysis()$m
        
        p  <- rbind(ols=ols., orm=orm.)

        return(list( p=p)) 
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      })  
      
    
      
      predt <- reactive({
        
        dat <- mcmc()$dat
    
        
        f2 <- orm(y ~baseline + treatment, data=dat )
        f3 <- ols(y ~baseline + treatment, data=dat )
        
        K <- analysis()$K
        m <- analysis()$m
    
        
        ols1. <- Predict(f3,  conf.type="mean",  treatment)
        orm1. <- Predict(f2,  treatment, fun=m, kint=K)
       
       pt  <- rbind(ols=ols1., orm=orm1.)
        
        return(list( pt=pt)) 
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      })  
      
      
      
      output$r.plot <- renderPlot({   
        
      #   dat <- mcmc()$dat
      #   kk <-   ( as.numeric(unlist(strsplit(input$kints,","))))
      #   dat$y <- as.numeric(as.character(dat$y))
      #   
      # 
      #   f2 <- orm(y ~baseline + treatment, data=dat )
      #   f3 <- ols(y ~baseline + treatment, data=dat )
      #   
      #   k <- NULL
      #   m <- Mean(f2, codes=FALSE)
      #   ml <- as.list(m)
      #   k <- ml$interceptRef
      #   
      #   ols.<- Predict(f3, conf.int=FALSE)
      #   
      #   # if input is empty do this 
      #   if(!isTruthy(kk)){
      #     
      #     K <- k
      #     orm. <-  Predict(f2, conf.int=FALSE, fun=m , kint=K) 
      #     
      #   } else {
      #     
      #     # if there is a value use it.
      #     K <- kk
      #     orm. <-  Predict(f2, conf.int=FALSE, fun=m , kint=K) 
      #     
      #   }
      #   
      # 
      # r <- rbind(ols      =  ols.,
      #            ordinal   = orm.
      # )
      
        pp <- predt()$pt
   
        
      plot(pp, groups='.set.')
      
})  

      
      
      
      
      
      
        
      output$dat <- renderPrint({
 
        return(print(mcmc()$dat, digits=4))
      })
      
      
      output$preds <- renderPrint({
        
        return(print(preds()$p, digits=4))
      })
      
      output$predt <- renderPrint({
        
        return(print(predt()$pt, digits=4))
      })
      
      
      
      
      
      
      # 
      # output$compare <- renderPlot({
      #   
      # #  print((mcmc()$r))
      #   mcmc()$r
      # })
      # 
      #  output$p <- renderPrint({
      #    
      #    return(print(mcmc()$pp, digits=4))
      #    
      # })
      #  
      #  output$orm.tm <- renderPrint({
      #    
      #    return(print(mcmc()$orm.tm, digits=4))
      #    
      #  })
      #  
      #  output$ols.tm <- renderPrint({
      #    
      #    return(print(mcmc()$ols.tm, digits=4))
      #    
     #  })
      #~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   output$check <- renderPlot({   
  #     
  #      
  #     Data <- mcmc()
  #     
  #     dat <- Data$dat
  #     f3 <- Data$f3  #ols
  #     f2 <- Data$f2  #orm
  #     levz <- Data$levz
  #     
  #   #  (f2 <- orm(y ~treatment + baseline, data=dat))
  #     #(f3 <- ols(y ~treatment + baseline, data=dat ))
  #    # 
  #     k <- NULL
  #     m <- Mean(f2, codes=FALSE)
  #     ml <- as.list(m)
  #     k <- ml$interceptRef
  #     i <- names(ml$intercepts)[ml$interceptRef]
  #     
  # #    ols.<- Predict(f3,  conf.type="mean",  baseline, treatment)
  #     
  #     
  #     if(!isTruthy( input$kints )) {   #if empty
  #       
  #       
  #       r <- rbind(ols      =  Predict(f3, conf.int=FALSE),
  #                  ordinal   = Predict(f2, conf.int=FALSE, fun=m, kint=k) 
  #       )
  #       
  #     } else {
  #       
  #       kk <-   ( as.numeric(unlist(strsplit(input$kints,","))))
  #       
  #       
  #       r <- rbind(ols      =  Predict(f3, conf.int=FALSE),
  #                  ordinal   = Predict(f2, conf.int=FALSE, fun=m, kint=kk) 
  #       )
  #       
  #     }
  #     
  #     plot(r, groups='.set.')
  #     
  #        
  #     
  #   })
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
      
    
})

# Run the application 
shinyApp(ui = ui, server = server)