# Missing Data Removed Default, Rounding Added
Mean <- function(x,digit=3) {
  round(mean(x, na.rm=TRUE),digit)}
Sd <- function(x,digit=3) {
  round(sd(x, na.rm=TRUE),digit)}
Median <- function(x,digit=3) {
  round(median(x, na.rm=TRUE),digit)}
Min <- function(x,digit=3) {
  round(min(x, na.rm=TRUE),digit)} 
Max <- function(x,digit=3) {
  round(max(x, na.rm=TRUE),digit)}

RowMedians <-function(x,digit=3,na.rm = TRUE) {
  round(rowMedians(as.matrix(x), na.rm=na.rm),digit)
  }

Cor <- function(x,digit=3) {
  round(cor(x, use = "pairwise.complete.obs"),digit)}
# Reverse Substrinf
substrRight <- function(x, n, fromend=0, end=nchar(x)-fromend){
  substr(x, end-n+1,end)}

df_AceEstimate=function(ACE){
  data<-data.frame(ACE@ASquared,ACE@CSquared,ACE@ESquared,ACE@CaseCount)
  names(data)<-c("ASquared","CSquared","ESquared","CaseCount")
  data}

# Correlation Matrix Nice
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
  if(include.n==FALSE){
    ## remove upper triangle
    Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
  }
  if(include.n==TRUE){
    #repalce upper triange with sample size
    require(psych)
    ct<-corr.test(x)
    ct<-as.matrix(ct$n)
    Rnew[upper.tri(Rnew, diag = TRUE)] <- ct[upper.tri(ct, diag = TRUE)]
  }
  Rnew <- as.data.frame(Rnew) 
  ## remove last column and return the matrix (which is now a data frame)

  Rnew <- cbind(Rnew[1:length(Rnew)-1])

  return(Rnew) 
}

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


## Knitr/Sweave Stuff


listing<- function(x, options) {
 paste0("\\begin{lstlisting}[breaklines=true]\n",
        x, "\\end{lstlisting}\n")
}

### Modified detachAllData from defunct epicalc package
detachAllData<-function () {
  pos.to.detach <- (1:length(search()))[substring(search(), 
                                                  first = 1, last = 8) != "package:" & search() != ".GlobalEnv" & 
                                          search() != "Autoloads" & search() != "CheckExEnv" & 
                                          search() != "tools:rstudio" & search() != "TempEnv"]
  for (i in 1:length(pos.to.detach)) {
    if (length(pos.to.detach) > 0) {
      detach(pos = pos.to.detach[1])
      pos.to.detach <- (1:length(search()))[substring(search(), 
                                                      first = 1, last = 8) != "package:" & search() != 
                                              ".GlobalEnv" & search() != "Autoloads" & search() != 
                                              "CheckExEnv" & search() != "tools:rstudio" & 
                                              search() != "TempEnv"]
    }
  }
}

#### Remove row if missing value in specific collumn

completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}


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