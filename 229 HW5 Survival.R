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

DATA <- read.table("~/Desktop/STAT 229/actg320.dat", col.names = c('id','time','censor','time_d','censor_d','tx','txgrp','strat2','sex','race','ivdrug','hemophil','karnof','cd4','priorzdv','age'))

factor_cols <- c("txgrp", "race", "ivdrug", "karnof","sex")
DATA %<>% mutate_at(factor_cols, factor)
DATA %<>% mutate(karnofmid = ifelse(karnof == 80|karnof == 70, 0, 1))

fit2 <- coxph(Surv(time, censor) ~ tx + karnofmid + cd4 + age, data = DATA)
summary(fit2)
########################### 1.1 ######################################
test_assumptions(fit2)

sresids <- residuals( fit2, type="scaledsch" )
colnames( sresids ) <- names( fit2$coef )
time <- as.numeric( rownames( sresids ) )

plot( log(time), sresids[,'cd4'], xlab=" Log Time",ylab="Scaled Schoenfeld Residual",main="Residuals for CD4")
lines( lowess( log(time), sresids[,'cd4'] ), col="red", lwd=2 )

plot( log(time), sresids[,'age'], xlab=" Log Time",ylab="Scaled Schoenfeld Residual",main="Residuals for age")
lines( lowess( log(time), sresids[,'age'] ), col="red", lwd=2 )

plot( log(time), sresids[,'tx'], xlab=" Log Time",ylab="Scaled Schoenfeld Residual",main="Residuals for tx")
lines( lowess( log(time), sresids[,'tx'] ), col="red", lwd=2 )

plot( log(time), sresids[,'karnofmid'], xlab=" Log Time",ylab="Scaled Schoenfeld Residual",main="Residuals for karnofmid")
lines( lowess( log(time), sresids[,'karnofmid'] ), col="red", lwd=2 )

########################### 1.2 ######################################
resids <- residuals( fit2, type="score" )
plot(DATA[,"cd4"], resids[,"cd4"],  xlab="cd4", ylab="Score Residual")
plot(DATA[,"age"], resids[,"age"],  xlab="age", ylab="Score Residual")
boxplot(resids[,"karnofmid"]~karnofmid,data=DATA, xlab="Karnofmid", ylab="Score Residuals")
boxplot(resids[,"tx"]~tx,data=DATA, xlab="tx", ylab="Score Residuals")


########################### 1.3 ######################################
gof(fit2)

########################### 1.4 ######################################
summary(fit2)
ratioConfint(fit2, 'age','age',a=5,b=0)
ratioConfint(fit2, 'cd4','cd4',a=50,b=0)

summary(DATA$age)
IQR(DATA$age)

########################### 1.6 ######################################
rscore <-predict(fit2, type='risk')
modRisk<- rscore-(-0.643987*DATA$tx)
tx0 <- exp(median(modRisk))
tx1 <- exp(median(modRisk)-0.643987)


baseSurv<-survfit(fit2,newdata = data.frame(tx=0,age=median(DATA$age),cd4=median(DATA$cd4),karnofmid=1))
df1 <- rbind(data.frame(time = baseSurv$time, surv = exp(-baseSurv$cumhaz)^tx0,tx=0),data.frame(time = baseSurv$time,surv = exp(-baseSurv$cumhaz)^tx1,tx=1))
ggplot(df1, aes(x = time, y=surv, color = factor(tx))) +
  geom_step() +
  ylim(0,1)
  labs(title = 'Survival Function for differing Treatments', x='Time', y = 'Survival Probability', color = 'chf')


