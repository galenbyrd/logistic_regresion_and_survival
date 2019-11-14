library(dplyr)
library(readxl)
library(survival)
library(survMisc)
library(ggplot2)
#library(tidyverse)
library(magrittr)

DATA <- read_excel("~/Desktop/STAT 229/Stat 229a WHAS 500.xls")
DATA %<>% mutate(age65 = age - 65, hr85 = hr-85,lenfolyr=lenfol/365)


fita <- coxph(Surv(lenfol, fstat) ~ gender, data = DATA)
summary(fita)
summary_table.coxph(fita)
exp(confint(fita, level=0.9))

fitb <- coxph(Surv(lenfol, fstat) ~ gender+bmi, data = DATA)
summary(fitb)
exp(confint(fitb, level=0.9))
exp(5*fitb$coefficients['bmi'])
exp(5*fitb$coefficients['bmi']+qnorm(.95)*5*0.01491)
exp(5*fitb$coefficients['bmi']-qnorm(.95)*5*0.01491)
summary_table.coxph(fitb)


fitd <- coxph(Surv(lenfol, fstat) ~ gender+bmi+gender*bmi, data = DATA)
summary(fitd)
summary_table.coxph(fitd)


fit2 <- coxph(Surv(lenfol, fstat) ~ gender + age65 + gender*age65 + mitype + miord + mitype*miord + hr85 + chf, data = DATA)
summary(fit2)
summary_table.coxph(fit2)
baseSurv<-survfit(fit2,newdata = data.frame(gender=0,age65=0,mitype=0,miord=0,hr85=0,chf=0))
plot(baseSurv, xlab='Time (days)',ylab='Survival Probability',main='Survival Function')

rscore <-predict(fit2, type='risk')
modRisk<- rscore-(0.784220*DATA$chf)

chf0 <- exp(median(modRisk))
chf1 <- exp(median(modRisk)-0.784220)

df1 <- rbind(data.frame(time = baseSurv$time, surv = exp(-baseSurv$cumhaz)^chf0,chf=0),data.frame(time = baseSurv$time,surv = exp(-baseSurv$cumhaz)^chf1,chf=1))
ggplot(df1, aes(x = time, y=surv, color = factor(chf))) +
  geom_step() +
  labs(title = 'Survival Function at Two levels of congestive heart failure', x='Time', y = 'Survival Probability', color = 'chf')




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

## gender 0 vs 1 at various bmis
ratioConfint(model=fitd, x="gender", y="gender:bmi", a=1, b=c(20,25,30), level=0.9)

## bmi x vs  x + 5 at gender 0 //// y = bmi and b=0 since we're not using it
ratioConfint(model=fitd, x="bmi", y="bmi", a=5, b=0, level=0.9)

## bmi x vs x + 5 at gender 1
ratioConfint(model=fitd, "bmi", "gender:bmi", a=5, b=5, level=0.9)








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
summary_table.coxph(fit2)










ggplot(mapping = aes(baseSurv$time))+
  geom_step(aes(y=exp(-baseSurv$cumhaz)^chf0)) +
  geom_step(aes(y=exp(-baseSurv$cumhaz)^chf1), color='red') +
  labs(title = 'Survival Function at Two levels of chf', x='Time', y = 'Survival Probability') +
  scale_colour_discrete(name  ="Payer",breaks=c("Congestive Heart Failure", "No chf"), labels=c(1,0))

         


