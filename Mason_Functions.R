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
cor_1star <- function(x,digit=3,sig=.05){ 
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
  
  ## remove upper triangle
  Rnew <- as.matrix(Rnew)
  Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
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
