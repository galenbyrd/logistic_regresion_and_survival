library(dplyr)
library(readxl)
library(survival)
library(survMisc)
library(ggplot2)
library(magrittr)
library(mfp)

fixp <- function(x, dig=3){
  x <- as.data.frame(x)
  if(substr(names(x)[ncol(x)],1,2) != "Pr")
    warning("The name of the last column didn't start with Pr. This may indicate that p-values weren't in the last row, and thus, that this function is inappropriate.")
  x[,ncol(x)] <- round(x[,ncol(x)], dig)
  for(i in 1:nrow(x)){
    if(x[i,ncol(x)] == 0)
      x[i,ncol(x)] <- paste0("< 0.", paste0(rep(0,dig-1), collapse=""), "1")
  }
  x
}
## for use with models formed from coxph
summary_table.coxph <- function(model) {
  require(xtable)
  summ <- summary(model)$coef[,-2] # remove exp(coef column)
  summ <- fixp(summ)
  intervals <- confint(model)
  combined <- merge(summ, intervals, by="row.names", all.x=TRUE)
  xt <- xtable(combined, digits=c(1,1,3,4,2,3,3,3), align=c("l","l", "S", "S", "S", "S[table-format = <2.3]", "S", "S"))
  ## xt <- fixp(xt)
  names(xt) <- c("", "\\mc{Estimate}", "\\mc{Std. Error}", "\\mc{$z\\text{-value}$}", "\\mc{Pr$\\parens{> |z|}$}", "\\mc{2.5 \\%}", "\\mc{97.5 \\%}")
  ## align(xt) <- "lSSSS[table-format = <2.3]SS"
  print(xt, include.rownames=FALSE, booktabs=TRUE, sanitize.colnames.function = identity, sanitize.text.function = identity, sanitize.numbers=FALSE)
}
#summary_table.coxph(fit2)

test_assumptions <- function(model) {
  testI <- cox.zph(model, transform=I, global=TRUE)
  testLog <- cox.zph(model, transform=log, global=TRUE)
  testKM <- cox.zph(model, transform="km", global=TRUE)
  testRank <- cox.zph(model, transform="rank", global=TRUE)
  data.frame(chisqI=testI$table[,"chisq"], pI=testI$table[,"p"],
             chisqLog=testLog$table[,"chisq"], pLog=testLog$table[,"p"],
             chisqKM=testKM$table[,"chisq"], pKM=testKM$table[,"p"],
             chisqRank=testRank$table[,"chisq"], pRank=testRank$table[,"p"])
}
fracpoly.table <- function(time, event, x, data, others=c()) {
  extra <- ""
  for (cov in others) {
    extra <- paste(extra, " + ", "data[,'", cov, "']", sep="")
  }
  ## print(extra)
  
  getFormula <- function(form) {
    as.formula(paste(form, extra))
  }
  
  model.null <- coxph(getFormula("Surv(data[, time], data[,event]) ~ 1"))
  model.main <- coxph(getFormula("Surv(data[, time], data[,event]) ~ data[,x]"))
  model.mfp2 <- mfp(getFormula("Surv(data[, time], data[,event]) ~ fp(data[,x], alpha=1, df=2)"), family="cox")
  model.mfp4 <- mfp(getFormula("Surv(data[, time], data[,event]) ~ fp(data[,x], alpha=1, df=4)"), family="cox")
  
  if (length(others) == 0) {
    null.deviance <- -2 * model.null$loglik
  } else {
    null.deviance <- -2 * model.null$loglik[2]
  }
  main.deviance <- -2 * model.main$loglik[2]
  mfp2.deviance <- model.mfp2$dev
  mfp4.deviance <- model.mfp4$dev
  
  diff0 <- null.deviance - mfp4.deviance
  diff1 <- main.deviance - mfp4.deviance
  diff2 <- mfp2.deviance - mfp4.deviance
  
  ## for p-values of these Likelihood ratio tests, use df = 4 - current
  pval0 <- pchisq(diff0, 4, lower.tail=FALSE)
  pval1 <- pchisq(diff1, 3, lower.tail=FALSE)
  pval2 <- pchisq(diff2, 2, lower.tail=FALSE)
  
  ## make a table
  mfptest <- data.frame(deviance=c(null.deviance, main.deviance, mfp2.deviance, mfp4.deviance),
                        dev_diff=c(diff0,diff1,diff2,NA),
                        p_val=c(pval0,pval1,pval2,NA),
                        powers=c("","1",
                                 # get the power of each model
                                 # for mfp2, just get the first term since its a 1 term model
                                 model.mfp2$powers["data[, x]",][1],
                                 # for mfp4, collapse the vector into a single string so you get both
                                 paste(model.mfp4$powers["data[, x]",], collapse=" ")))
  
  row.names(mfptest) <- c("Not in Model", "Linear", "m=1", "m=2")
  
  mfptest
}
test_assumptions <- function(model) {
  testI <- cox.zph(model, transform=I, global=TRUE)
  testLog <- cox.zph(model, transform=log, global=TRUE)
  testKM <- cox.zph(model, transform="km", global=TRUE)
  testRank <- cox.zph(model, transform="rank", global=TRUE)
  data.frame(chisqI=testI$table[,"chisq"], pI=testI$table[,"p"],
             chisqLog=testLog$table[,"chisq"], pLog=testLog$table[,"p"],
             chisqKM=testKM$table[,"chisq"], pKM=testKM$table[,"p"],
             chisqRank=testRank$table[,"chisq"], pRank=testRank$table[,"p"])
}
ratioConfint <- function(model, x, y, a, b, level=0.95) {
  ## finds point estimate and confidence interval for a linear combination of two variables: aX + bY
  ## model is some kind of model (glm or coxph)
  ## x and y are strings of the two coefficients of the model
  ## a and b are what x and y get multiplied by
  ## level is the size of the confidence interval (default = 0.95 i.e. 95% confidence interval)
  ## sample call for the interaction of bmi (continuous) with gender (dichotomous),
  ## comparing gender = 0 to gender = 1 at set levels of bmi with 90% confidence intervals:
  ## ratioConfint(model=model, x="gender", y="gender:bmi", a=1, b=c(20,25,30), level=0.9)
  X <- coef(model)[x]
  Y <- coef(model)[y]
  varX <- vcov(model)[x, x]
  varY <- vcov(model)[y, y]
  covXY <- vcov(model)[x, y]
  
  pm <- function(x,y) {
    data.frame(estimate=x, low=(x - y), high=(x + y))
  }
  stderr <- sqrt(a^2 * varX + b^2 * varY + 2 * a * b * covXY)
  lowBound <- (1-level)/2
  highBound <- 1 - lowBound
  df <- exp(pm(a*X + b*Y, qnorm(highBound) * stderr))
  colnames(df) <- c("estimate", paste(lowBound*100, "%", sep=""), paste(highBound*100, "%", sep=""))
  df
}


DATA <- read.table("~/Desktop/STAT 229/Byrd Framingham Survival Problem.out", header = TRUE)

#DATA %<>% mutate_at("CIGPDAY", factor)
DATA %<>% mutate(CIGS21 = ifelse(CIGPDAY >= 21, 0, 1))


fit <- coxph(Surv(DEATHYRS, DEATH) ~ SEX + AGE + SYSBP + CIGPDAY + BMI + ANYCHD + CVD, data = DATA)
summary(fit)
########################### 1 ######################################
summary_table.coxph(fit)

########################### 2 ######################################
fracpoly.table('DEATHYRS','DEATH','AGE',DATA,c('SEX','AGE','SYSBP','CIGPDAY','BMI','ANYCHD','CVD'))
fracpoly.table('DEATHYRS','DEATH','BMI',DATA,c('SEX','AGE','SYSBP','CIGPDAY','BMI','ANYCHD','CVD'))
fracpoly.table('DEATHYRS','DEATH','SYSBP',DATA,c('SEX','AGE','SYSBP','CIGPDAY','BMI','ANYCHD','CVD'))
#fracpoly.table('DEATHYRS','DEATH','CIGPDAY',DATA,c('SEX','AGE','SYSBP','CIGPDAY','BMI','ANYCHD','CVD'))

########################### 3 ######################################
test_assumptions(fit)

sresids <- residuals( fit, type="scaledsch" )
colnames( sresids ) <- names( fit$coef )
time <- as.numeric( rownames( sresids ) )

plot( log(time), sresids[,'SEX'], xlab=" Log Time",ylab="Scaled Schoenfeld Residual",main="Residuals for SEX")
lines( lowess( log(time), sresids[,'SEX'] ), col="red", lwd=2 )

########################### 4 ######################################
resids <- residuals( fit, type="score" )
plot(DATA[,"AGE"], resids[,"AGE"],  xlab="AGE", ylab="Score Residual")
plot(DATA[,"SYSBP"], resids[,"SYSBP"],  xlab="SYSBP", ylab="Score Residual")
plot(DATA[,"BMI"], resids[,"BMI"],  xlab="BMI", ylab="Score Residual")
plot(DATA[,"CIGPDAY"], resids[,"CIGPDAY"],  xlab="CIGPDAY", ylab="Score Residual")

boxplot(resids[,"SEX"]~SEX,data=DATA, xlab="SEX", ylab="Score Residuals")
boxplot(resids[,"ANYCHD"]~ANYCHD,data=DATA, xlab="ANYCHD", ylab="Score Residuals")
boxplot(resids[,"CVD"]~CVD,data=DATA, xlab="CVD", ylab="Score Residuals")


########################### 5 ######################################
gof(fit)

########################### 6 ######################################
ratioConfint(fit, 'CIGPDAY','CIGPDAY',a=20,b=0)
ratioConfint(fit, 'CIGPDAY','CIGPDAY',a=10,b=0)

########################### 7 ######################################
summary(fit)
ratioConfint(fit, 'AGE','AGE',a=5,b=0)
IQR(DATA$AGE)
ratioConfint(fit, 'BMI','BMI',a=2,b=0)
IQR(DATA$BMI)
ratioConfint(fit, 'SYSBP','SYSBP',a=10,b=0)
IQR(DATA$SYSBP)

ratioConfint(fit, 'CIGPDAY','CIGPDAY',a=10,b=0)
IQR(DATA$CIGPDAY)



