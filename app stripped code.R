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

## convenience functions
p0 <- function(x) {formatC(x, format="f", digits=1)}
p1 <- function(x) {formatC(x, format="f", digits=1)}
p2 <- function(x) {formatC(x, format="f", digits=2)}
p3 <- function(x) {formatC(x, format="f", digits=3)}
p5 <- function(x) {formatC(x, format="f", digits=5)}
logit <- function(p) log(1/(1/p-1))
expit <- function(x) 1/(1/exp(x) + 1)
inv_logit <- function(logit) exp(logit) / (1 + exp(logit))
is.even <- function(x){ x %% 2 == 0 } # function to id. odd maybe useful
options(width=200)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

n <- 1000                           
dist <- c(22,21)
levels <- 10
b1<-or1<-2
b2<-or2 <-1
                              
bas1=1
bas2=2
rcat2<- 999 
 group <-1
 rcat <-999
 kints<- "" 
      base=1                                          
                                                 
                              
 
    
    dis <- dist
    
    trt <- n 
    
    ctr <- levels 
    
    n1y1 <-   or1 
    
    n2y2 <-  or2 
    
    
    base <- base 
    
    
   
      n=trt  
      lev=ctr 
      or1=n1y1 
      or2=n2y2 
      shape1=dis[1] 
      shape2=dis[2] 
      base=base 
      levz=ctr 
      
 
    
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
    dat99 <- dat <- data.frame(treatment, baseline, y = factor(y))
    #      sf1 <- summary(f1, antilog=TRUE, verbose=FALSE)
    
 
 ###########################
    
    bas1  <-  as.numeric( bas1)
    bas2  <-  as.numeric(bas2)

    kk <-   kints 
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
    
 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # beta dist plot 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
  
 
    shape1. <- shape1
    shape2. <- shape2
    
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
    

  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  end ggplot barplot of beta distribution
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    

    levz <- lev
    n   <- n
    
    dat <- dat
    
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
    
 
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    levz <-lev
    n   <-n
    

    
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
    
 
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # non cummulative predicted probabilities plot run the analysis again
  # not efficient I know
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
    
    # n    <- sample$n
    levz <-  lev
    base <-  base
    
    datx <-  dat99
    
    rcat <-  rcat2    
    
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
    
 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # non cummulative predicted probabilities plot run the analysis again
  # not efficient I know
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
    # n    <- sample$n
    levz <-  lev
    base <-  base
    
    datx <- dat99
    
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
    
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~tables of predictions
 
    
   
    f    <-  f2
    levz <-  lev
    
    group <- group   
    rcat <-   rcat  
    
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
    
   #plotci <- plot(f, baseline, treatment, fun = stats::plogis)
    
    
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plots of predictions
  
   
    # 
    f    <- f2
    # 
    levz <-  lev
    
    
    group <- group   
    rcat <-   rcat  
    
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
    
    print(gpp)
    
 
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # text 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
  
 
    
    K <-  K
    
    txt <- paste0("Ordinal intercept ", K)
    
    levz <- levz
    
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
    
    
 
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~baseline plots~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
    
    
    K <-  K
    
    txt <- paste0("Ordinal intercept ", K)
    
    levz <- levz
    
    P <-  P
    
   # P$treatment <- ifelse(P$treatment %in% 0, "Placebo", "Treatment")
    
    
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
    
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # a plot of coef from series of logistic regression models. checking assumptions
  
 
    #levz <- sample$lev
    dat <-  dat99
    
    f1    <-  f2
    
  #  or1 #<- log(as.numeric(unlist(strsplit(input$or1,","))))   # user enter odds , need log for the maths
    # R
 #   or2 <- log(as.numeric(unlist(strsplit(input$or2,","))))    # user enter odds , need log for the maths
    
    # x <- plyr::arrange(x, variable, value)
    # I recommend partial residual plots using the rms package's 
    # lrm and residuals.lrm functions. You can also fit a series of binary models using different 
    # #cutoffs for Y and plot the log odds ratios vs. cutoff.
    #https://stats.stackexchange.com/questions/25988/proportional-odds-assumption-in-ordinal-logistic-regression-in-r-with-the-packag?noredirect=1&lq=1
    
    ps <- sort(unique(dat$y))  # get cut offs
    
    d <<- datadist(dat)
    options(datadist="d")
    A <- d$limits["Low:effect",]$baseline -0
    B <- d$limits["High:effect",]$baseline +0
    ps <- ps[1:(length(ps)-1)] # not using highest one
    ps <- ps[A:B] # n
    
    x<- matrix(NA, ncol = length(ps), nrow = 3)  # matrix to collect estimates
    
    Y <- as.numeric(dat$y) 
    # function to run logistic regression and capture coefs
    for ( i in 1:length(ps)) {
      
      dat$y1 <- ifelse(Y <=  as.numeric(ps[i]) ,0,1)
      x[1,i] <- coef(orm(y1 ~treatment + baseline, data=dat))["baseline"][[1]]
      x[2,i] <- coef(orm(y1 ~treatment + baseline, data=dat))["treatment"][[1]]
      x[3,i] <- ps[i]
    }
    
    # manage 
    x <- t(as.data.frame(x))
    x <-  as.data.frame(x)
    
    names(x) <- c("baseline","treatment","cutoff")
    x<- melt(x, id.vars= c("cutoff"))
    x$true <- ifelse(x$variable %in% "baseline",  or2,  or1)  #baslein
    
    # bring in the PO estimates
    co <- confint(f1)[c("baseline","treatment"),]
    ceb <- coef(f1)[c("baseline" )]
    cet <- coef(f1)[c( "treatment")]
    
    xx <- rbind(ceb, cet) 
    xx <- cbind(xx,co)
    xx <- as.data.frame(xx)
    xx$variable <- c("baseline","treatment")
    xx <- merge(x,xx)
    
    names(xx) <- c("variable","cutoff", "value", "true","est","lower", "upper")
    require(ggrepel)
    p <- ggplot(xx, aes(x=cutoff , y= value  , colour=variable, label=cutoff))+
      geom_point(size = 5)+ ylab("Logit")+xlab("A series of binary models using different cutoffs for Y") +
      facet_wrap(.~variable, scales="free") +
      geom_text_repel(aes(label = cutoff),
                      size = 5) +
      geom_hline(data = xx, aes(yintercept = true), linetype="dashed", color="black", size=.7) +
      geom_hline(data = xx, aes(yintercept = lower), linetype="solid", color="blue", size=.7) +
      geom_hline(data = xx, aes(yintercept = upper), linetype="solid", color="blue", size=.7) +
      geom_hline(data = xx, aes(yintercept = est), linetype="solid", color="blue", size=.7)   
    
    
    p <- p + theme(panel.background=element_blank(),
                   plot.title=element_text(), plot.margin = unit(c(5.5,12,5.5,5.5), "pt"), 
                   legend.text=element_text(size=12),
                   legend.title=element_text(size=14),
                   axis.text.x = element_text(size=12),
                   axis.text.y = element_text(size=10),
                   axis.line.x = element_line(color="black"),
                   axis.line.y = element_line(color="black"),
                   axis.title = element_text(size = 20) , 
                   plot.caption=element_text(hjust = 0, size = 13),  #left align
                   strip.text = element_text(size=15),
                   axis.title.x = element_text(color = "grey20",size = 15) ,
                   axis.title.y = element_text(color = "grey20",size = 15)
    )
    
    g <- p + theme(legend.position="none") 
    g
    
    
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~assumption plot~~~~~~~~~~~~~~~~~~~~~~~~    
  # on the fly plot harrell's PO assumption plot...
  
   
    
    dat <-  dat99
    levz <- levz
    l2 <- as.numeric(levz)-1
    y <- as.integer(dat$y)  
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # here turn text into a function to allow flexibilty in changing levels
    start <- "function(y)
   c("
    
    res <-  vector()
    
    x <- for(i in 2:levz) {
      res[i] <- paste0(" 'y>=" ,i, "'=qlogis(mean(y>=",i,")),") 
    }
    
    txt <- paste(res[-1], collapse = '')
    z <- gsub(".{1}$", "", txt)
    v <- z[]
    
    this <- paste0(start, v, ")",collapse='')
    
    #~~~~~~~~~~~~~~~~~~
    Assumption <-  eval(parse(text= (this)))
    
    s<- summary(y ~ treatment + baseline, fun=Assumption, data=dat)
     
    
    is.na(s) <- do.call(cbind,lapply(s, is.infinite))  # need to remove infinity as plot will return xlim error
    
    plot(s, which=1:l2, pch=1:l2, xlab='Logit',
         main ="Checking the proportional odds assumption")
    
  
 
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~baseline predictions~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
    
    ols. <-  ols.
    orm. <-  orm.
    
    p  <- rbind(ols=ols., orm=orm.)
    plot(P)
   
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # text 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
  
 
    A <-  f2     
    
    d <-  d
    
    dat <-  dat99
    
    dat$y <- as.numeric(as.character(dat$y))
    
   
    levz <- leva
    
    f <- A$coefficients
    x <-length(f) -2
    
    diff <- d$limits["High:effect","baseline"]-d$limits["Low:effect","baseline"]
    
    c1  <-expit(f[1][[1]]  + f['baseline'][[1]])
    
    c10 <- expit(f[x][[1]]  + f['baseline'][[1]])
    
   paste0( "Let's interpret the output on the left. The coefficient alongside y>=",min(dat$y)+1," is "
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
                 
    )     
    
 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
    
    dat <-  dat99
    
    dat$y <- as.numeric(as.character(dat$y))
    
    d <<- datadist(dat)
    options(datadist="d")
    
    f <- ols(y ~treatment + (baseline), data=dat)
    an <- anova(f)
    
 
    
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
    
    linear <- ols(y ~treatment + (baseline), data=dat)
    
    an1 <- anova(linear)
    
    ggplot(Predict(linear, treatment), anova=an1, pval=TRUE)  
    
    
 
 