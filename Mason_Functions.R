########## Functions Mason Garrison Commonly Uses. I do not claim I wrote these functions. I typically have modified existing code to suite my purposes.

# cor_1star
## Prints out clean correlation matrix. If you want to replace the upper triangle with the sample size, toggle include.n, whose default is FALSE.
cor_1star <- function(x,digit=3,sig=.05,include.n=FALSE){ 
  require(Hmisc) 
  x <- as.matrix(x) 
  R <- rcorr(x)$r 
  p <- rcorr(x)$P 
  
  ## define notions for significance levels; spacing is important.
  mystars <- ifelse(p < sig, "* ", " ")
  
  ## trunctuate the matrix that holds the correlations to two decimal
  R <- format(round(cbind(rep(-1.11, ncol(x)), R), digit))[,-1] 
  
  ## build a new matrix that includes the correlations with their apropriate stars 
  Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x)) 
  diag(Rnew) <- paste(diag(R), " ", sep="") 
  rownames(Rnew) <- colnames(x) 
  colnames(Rnew) <- paste(colnames(x), "", sep="") 
  
  Rnew <- as.matrix(Rnew)
  if(!include.n){
    ## remove upper triangle
    Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
    Rnew <- as.data.frame(Rnew) 
    Rnew <- cbind(Rnew[1:length(Rnew)-1])
  }
  if(include.n){
    #repalce upper triange with sample size
    require(psych)
    ct<-corr.test(x)
    ct<-as.matrix(ct$n)
    Rnew[upper.tri(Rnew, diag = TRUE)] <- ct[upper.tri(ct, diag = TRUE)]
    Rnew <- as.data.frame(Rnew) 
    Rnew <- cbind(Rnew[1:length(Rnew)])#-1])
  }
  
  return(Rnew) 
}



# detachAllData
## Modified detachAllData function from defunct epicalc package. Allows user to detach all datasets at once.
detachAllData<-function () {
  pos.to.detach <- (1:length(search()))[
    substring(search(),first = 1, last = 8)
    !="package:"&search() 
    !=".GlobalEnv"&search()
    !="Autoloads"&search() 
    !="CheckExEnv"&search() 
    !="tools:rstudio"&search() 
    !="TempEnv"]
  for (i in 1:length(pos.to.detach)) {
    if (length(pos.to.detach) > 0) {
      detach(pos = pos.to.detach[1])
      pos.to.detach <- (1:length(search()))[
        substring(search(),first = 1, last = 8) 
        !="package:"&search() 
        !=".GlobalEnv"&search() 
        !="Autoloads"&search() 
        !="CheckExEnv"&search() 
        !="tools:rstudio"&search() 
        !="TempEnv"]
    }
  }
}

# Reverse Substring
substrRight <- function(x, n, fromend=0, end=nchar(x)-fromend){
  substr(x, end-n+1,end)}

# ROUNDING FUNCTIONS
## Function rounds to three digits and defaults to na.rm = TRUE
## Cor
Cor <- function(x,digit=3, use="pairwise.complete.obs") {
  round(cor(x, use = use),digit)}

## Max
Max <- function(x,digit=3,na.rm = TRUE) {
  round(max(x, na.rm=na.rm),digit)}

## Mean
Mean <- function(x,digit=3,na.rm = TRUE) {
  round(mean(x, na.rm=na.rm),digit)}

## Median
Median <- function(x,digit=3,na.rm = TRUE) {
  round(median(x, na.rm=na.rm),digit)}

## Min
Min <- function(x,digit=3,na.rm = TRUE) {
  round(min(x, na.rm=na.rm),digit)} 

## RowMedians
RowMedians <-function(x,digit=3,na.rm = TRUE) {
  round(rowMedians(as.matrix(x), na.rm=na.rm),digit)}

## Sd
Sd <- function(x,digit=3,na.rm = TRUE) {
  round(sd(x, na.rm=na.rm),digit)}

# seed.alpha
## Function converts string into hex. Set 'set.seed' to TRUE if you want the seed set. Set 'keep.seed' to FALSE if you don't want the seed value returned.
seed.alpha <- function(x,set.seed=FALSE,keep.seed=TRUE) {
  require("digest")
  hexval <- paste0("0x",digest(x,"crc32"))
  intval <- type.convert(hexval) %% .Machine$integer.max
  if(set.seed){
    set.seed(intval)
  }
  if(keep.seed){
    return(intval)
  }
}


df_AceEstimate=function(ACE){
  data<-data.frame(ACE@ASquared,ACE@CSquared,ACE@ESquared,ACE@CaseCount)
  names(data)<-c("ASquared","CSquared","ESquared","CaseCount")
  data}



cor_3stars <- function(x,digit=3){ 
  require(Hmisc) 
  x <- as.matrix(x) 
  R <- rcorr(x)$r 
  p <- rcorr(x)$P 
  
  ## define notions for significance levels; spacing is important.
  mystars <- ifelse(p < .001, "***", ifelse(p < .01, "** ", ifelse(p < .05, "* ", " ")))
  
  ## trunctuate the matrix that holds the correlations to two decimal
  R <- format(round(cbind(rep(-1.11, ncol(x)), R), digit))[,-1] 
  
  ## build a new matrix that includes the correlations with their apropriate stars 
  Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x)) 
  diag(Rnew) <- paste(diag(R), " ", sep="") 
  rownames(Rnew) <- colnames(x) 
  colnames(Rnew) <- paste(colnames(x), "", sep="") 
  
  ## remove upper triangle
  Rnew <- as.matrix(Rnew)
  Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
  Rnew <- as.data.frame(Rnew) 
  
  ## remove last column and return the matrix (which is now a data frame)
  Rnew <- cbind(Rnew[1:length(Rnew)-1])
  return(Rnew) 
}




#### Remove row if missing value in specific collumn

completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

#### Legacy Functions
Mean_0 <- function(x) base::round(mean(x, na.rm=TRUE),0)
Mean_1 <- function(x) base::round(mean(x, na.rm=TRUE),1)
Mean_2 <- function(x) base::round(mean(x, na.rm=TRUE),2)
Mean_3 <- function(x) base::round(mean(x, na.rm=TRUE),3)
Mean_4 <- function(x) base::round(mean(x, na.rm=TRUE),4)

Sd_0 <- function(x) round(stats::sd(x, na.rm=TRUE),0)
Sd_1 <- function(x) round(stats::sd(x, na.rm=TRUE),1)
Sd_2 <- function(x) round(stats::sd(x, na.rm=TRUE),2)
Sd_3 <- function(x) round(stats::sd(x, na.rm=TRUE),3)
Sd_4 <- function(x) round(stats::sd(x, na.rm=TRUE),4)

##
lm_eqn = function(df){
  m = lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, list(a = format(coef(m)[1], digits = 2), 
                                                                                  b = format(coef(m)[2], digits = 2), 
                                                                                  r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}

## Extract Values from aalen
aalen_simple_summary <- function(model=NULL){
  iv<-names(model$obs.testBeq0)
if(length(dimnames(model$gamma)[[1]]) >0){
  iv = c(iv, dimnames(model$gamma)[[1]])
}
  diag<-as.data.frame(diag(model$robvar.gamma))
  rownames(diag)<-colnames(model$robvar.gamma)
  data<- data.frame(matrix(as.numeric(NA),nrow=length(iv),ncol = 10))
  names(data)<-c("Var","obs.testBeq0","pval.testBeq0","obs.testBeqC","pval.testBeqC","obs.testBeqC.is","pval.testBeqC.is","gamma","zobs.gamma","p.gamma")
  data$Var<-as.character(data$Var)
  for(i in 1:length(iv)){
    if(substr(iv[i],1,6)!="const("){
  data[i,]<-I(c(iv[i],model$obs.testBeq0[paste0(iv[i])],model$pval.testBeq0[paste0(iv[i])],model$obs.testBeqC[paste0(iv[i])],model$pval.testBeqC[paste0(iv[i])],model$obs.testBeqC.is[paste0(iv[i])],model$pval.testBeqC.is[paste0(iv[i])],NA,NA,NA))
  }
    if(substr(iv[i],1,6)=="const("){
      data[i,]<-I(c(iv[i], NA,NA,NA,NA,NA,NA,model$gamma[paste0(iv[i]),],model$gamma[paste0(iv[i]),]/sqrt(diag[paste0(iv[i]),1]),2*pnorm(-abs(model$gamma[paste0(iv[i]),]/sqrt(diag[paste0(iv[i]),1])))))
                    }
  }
  for(i in c(2:ncol(data))) {
    data[,i] <- as.numeric(as.character(data[,i]))
  }
  names(data)<-c("Var","t obs","p","Kolmogorov-Smirnov","p(K-S)","Cramer von Mises ","p(CvM)", "Parametric Gamma","z(gamma)","p(z(gamma))")
  data}
### Transpose Data Frame
transpose_dataframe<- function(dataframe=NULL)
{tdataframe <- as.data.frame(t(dataframe[,-1]))
colnames(tdataframe) <-  dataframe[,1]
tdataframe}



#### Find True Price

True.Price<-function(price=NULL,stock=NULL,base=100) {
  return(price*base/stock)
}


# Functions
# Functions to make ggplot KM survivor curves made with survfit() in library(survival)
#
# code written by Ramon Saccilotto and modified by Mason Garrison
# and included in his ggplot2 tutorial
# 2010-12-08

# define custom function to create a survival data.frame
createSurvivalFrame <- function(f.survfit){
  # initialise frame variable
  f.frame <- NULL
  
  # check if more then one strata
  if(length(names(f.survfit$strata)) == 0){
    # create data.frame with data from survfit
    f.frame <- data.frame(time=f.survfit$time, n.risk=f.survfit$n.risk, n.event=f.survfit$n.event, 
                          n.censor = f.survfit$n.censor, surv=f.survfit$surv, upper=f.survfit$upper, 
                          lower=f.survfit$lower)
    # create first two rows (start at 1)
    f.start <- data.frame(time=c(0, f.frame$time[1]), n.risk=c(f.survfit$n, f.survfit$n), n.event=c(0,0), 
                          n.censor=c(0,0), surv=c(1,1), upper=c(1,1), lower=c(1,1))
    # add first row to dataset
    f.frame <- rbind(f.start, f.frame)
    # remove temporary data
    rm(f.start)
  } 
  else {
    # create vector for strata identification
    f.strata <- NULL
    for(f.i in 1:length(f.survfit$strata)){
      # add vector for one strata according to number of rows of strata
      f.strata <- c(f.strata, rep(names(f.survfit$strata)[f.i], f.survfit$strata[f.i]))
    }
    # create data.frame with data from survfit (create column for strata)
    f.frame <- data.frame(time=f.survfit$time, n.risk=f.survfit$n.risk, n.event=f.survfit$n.event, 
                          n.censor = f.survfit$n.censor, surv=f.survfit$surv, upper=f.survfit$upper, 
                          lower=f.survfit$lower, strata=factor(f.strata))
    # remove temporary data
    rm(f.strata)
    # create first two rows (start at 1) for each strata
    for(f.i in 1:length(f.survfit$strata)){
      # take only subset for this strata from data
      f.subset <- subset(f.frame, strata==names(f.survfit$strata)[f.i])
      # create first two rows (time: 0, time of first event)
      f.start <- data.frame(time=c(0, f.subset$time[1]), n.risk=rep(f.survfit[f.i]$n, 2), 
                            n.event=c(0,0), n.censor=c(0,0), surv=c(1,1), upper=c(1,1), lower=c(1,1), 
                            strata=rep(names(f.survfit$strata)[f.i],2))
      # add first two rows to dataset
      f.frame <- rbind(f.start, f.frame)
      # remove temporary data
      rm(f.start, f.subset)
    }
    # reorder data
    f.frame <- f.frame[order(f.frame$strata, f.frame$time), ]
    # rename row.names
    rownames(f.frame) <- NULL
  }
  # return frame
  return(f.frame)
}

# define custom function to draw kaplan-meier curve with ggplot
qplot_survival <- function(f.frame, f.CI="default", f.shape=3){
  # use different plotting commands dependig whether or not strata's are given
  if("strata" %in% names(f.frame) == FALSE){
    # confidence intervals are drawn if not specified otherwise
    if(f.CI=="default" | f.CI==TRUE ){
      # create plot with 4 layers (first 3 layers only events, last layer only censored)
      # hint: censoring data for multiple censoring events at timepoint are overplotted
      # (unlike in plot.survfit in survival package)
      ggplot(data=f.frame) + geom_step(aes(x=time, y=surv), direction="hv") + geom_step(aes(x=time, 
                                                                                            y=upper), directions="hv", linetype=2) + geom_step(aes(x=time,y=lower), 
                                                                                                                                               direction="hv", linetype=2)  + 
        geom_point(data=subset(f.frame, n.censor==1), aes(x=time, y=surv), shape=1, position=position_jitter(width=0, height=0.05))  + 
        geom_point(data=subset(f.frame, n.censor==0), aes(x=time, y=surv), shape=16,size = 2, position=position_jitter(width=0, height=0.025))
    }
    else {
      # create plot without confidence intervalls
      ggplot(data=f.frame,) + geom_step(aes(x=time, y=surv), direction="hv") + 
        geom_point(data=subset(f.frame, n.censor==1), aes(x=time, y=surv), shape=1, position=position_jitter(width=0, height=0.05))  + 
        geom_point(data=subset(f.frame, n.censor==0), aes(x=time, y=surv), size = 2, shape=16, position=position_jitter(width=0, height=0.025))
    }
  }
  else {
    if(f.CI=="default" | f.CI==FALSE){
      # without CI 
      ggplot(data=f.frame, aes(group=strata, colour=strata)) + geom_step(aes(x=time, y=surv), 
                                                                         direction="hv") + 
        geom_point(data=subset(f.frame, n.censor==1), aes(x=time, y=surv), shape=1, position=position_jitter(width=0, height=0.05))  + 
        geom_point(data=subset(f.frame, n.censor==0), aes(x=time, y=surv), shape=16,size = 500, position=position_jitter(width=0, height=0.025))
    }
    else {
      # with CI (hint: use alpha for CI)
      
      ggplot(data=f.frame, aes(colour=strata, group=strata)) + geom_step(aes(x=time, y=surv), 
                                                                         direction="hv") + geom_step(aes(x=time, y=upper), directions="hv", linetype=2, 
                                                                                                     alpha=0.5) + geom_step(aes(x=time,y=lower), direction="hv", linetype=2, alpha=0.5) + 
        geom_point(data=subset(f.frame, n.censor==1), aes(x=time, y=surv), shape=f.shape)
    }
  }
}