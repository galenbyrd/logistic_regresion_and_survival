library(dplyr)
library(readxl)
library(survival)
library(survMisc)
library(ggplot2)
#library(tidyverse)
library(magrittr)
library(mfp)

DATA <- read.table("~/Desktop/STAT 229/whas100.dat", col.names = c('id','admitDate','folDate','los','lenfol','fstat','age','gender','bmi'))
#DATA %<>% mutate(lenfolyr=lenfol/365, age2 = age^2, bmi2 = bmi^2)

fit1 <- mfp(Surv(lenfol, fstat) ~ fp(age), data = DATA, family = 'cox', verbose = TRUE)

fit1 <- mfp(Surv(lenfol, fstat) ~ fp(bmi), data = DATA, family = 'cox', verbose = TRUE)

fracpoly.table('lenfol','fstat','age',DATA)
fracpoly.table('lenfol','fstat','bmi',DATA)



###################################### QUESTION 2 #########################################
DATA <- read.table("~/Desktop/STAT 229/actg320.dat", col.names = c('id','time','censor','time_d','censor_d','tx','txgrp','strat2','sex','race','ivdrug','hemophil','karnof','cd4','priorzdv','age'))

factor_cols <- c("txgrp", "race", "ivdrug", "karnof","sex")
DATA %<>% mutate_at(factor_cols, factor)
DATA %<>% mutate(karnofmid = factor(ifelse(karnof == 80|karnof == 70, 0, 1)))

# PART 1 Purposeful Selection
fit2 <- coxph(Surv(time, censor) ~ tx + sex + ivdrug + karnof + cd4 + priorzdv + age, data = DATA)
summary(fit2)

# PART 2/3 Purposeful Selection
# remove priorzdv and we see no change in coefs, so not a confounder and we can leave it out.
fit2 <- coxph(Surv(time, censor) ~ tx + sex + ivdrug + karnof + cd4 + age, data = DATA)
summary(fit2)
summary_table.coxph(fit2)

# remove sex and we see no change in coefs, so not a confounder and we can leave it out.
fit2 <- coxph(Surv(time, censor) ~ tx + ivdrug + karnof + cd4 + age, data = DATA)
summary(fit2)

# remove ivdrug and we see no change in coefs, so not a confounder and we can leave it out.
fit2 <- coxph(Surv(time, censor) ~ tx + karnof + cd4 + age, data = DATA)
summary(fit2)

# remove karnof and we see no change in coefs, so not a confounder and we can leave it out.
fit2 <- coxph(Surv(time, censor) ~ tx + cd4 + age, data = DATA)
summary(fit2)
summary_table.coxph(fit2)


fit2 <- coxph(Surv(time, censor) ~ tx + karnofmid + cd4 + age, data = DATA)
summary(fit2)

# PART 4 Purposeful Selection
covariates <- c('priorzdv','sex','ivdrug',"txgrp", "strat2",  "race", "hemophil")
for (cov in covariates) {
  form <- paste("Surv(time, censor) ~ tx + karnofmid + cd4 + age+", cov)
  model <- coxph(as.formula(form), data = DATA)
  print(summary(model))
  cat("\n")
}
# Dont add in any initially removed variables. This is our prelim main effects model.

# PART 5 Purposeful Selection
# LOOK At LINEARITY OF CONTINUOUS VARIABLES
fit2 <- mfp(Surv(time, censor) ~ tx + karnofmid + fp(cd4) + age, data = DATA, family = 'cox', verbose = TRUE)
fit2 <- mfp(Surv(time, censor) ~ tx + karnofmid + cd4 + fp(age), data = DATA, family = 'cox', verbose = TRUE)


# PART 6 Purposeful Selection
# consider interactions
covariates <- c("tx*karnofmid",  "tx*cd4", "tx*age",'karnofmid*cd4','karnofmid*age','cd4*age')
for (cov in covariates) {
  form <- paste("Surv(time, censor) ~ tx + karnofmid + cd4 + age+", cov)
  model <- coxph(as.formula(form), data = DATA)
  print(summary(model))
  cat("\n")
}

fit2 <- coxph(Surv(time, censor) ~ tx + karnofmid + cd4 + age, data = DATA)
summary(fit2)




fracpoly.table <- function(time, event, x, data) {
  model.null <- coxph(Surv(data[, time], data[,event]) ~ 1, data=data)
  model.main <- coxph(Surv(data[, time], data[,event]) ~ data[,x], data=data)
  model.mfp2 <- mfp(Surv(data[, time], data[,event]) ~ fp(data[,x], alpha=1, df=2), family="cox")
  model.mfp4 <- mfp(Surv(data[, time], data[,event]) ~ fp(data[,x], alpha=1, df=4), family="cox")
  
  null.deviance <- -2 * model.null$loglik
  main.deviance <- -2 * model.main$loglik[2]
  
  diff0 <- null.deviance - model.mfp4$dev
  diff1 <- main.deviance - model.mfp4$dev
  diff2 <- model.mfp2$dev - model.mfp4$dev
  
  ## for p-values of these Likelihood ratio tests, use df = 4 - current
  pval0 <- pchisq(diff0, 4, lower.tail=FALSE)
  pval1 <- pchisq(diff1, 3, lower.tail=FALSE)
  pval2 <- pchisq(diff2, 2, lower.tail=FALSE)
  
  ## make a table
  mfptest <- data.frame(deviance=c(null.deviance, main.deviance, model.mfp2$dev, model.mfp4$dev),
                        dev_diff=c(diff0,diff1,diff2,NA),
                        p_val=c(pval0,pval1,pval2,NA),
                        powers=c("","1",
                                 levels(model.mfp2$fptable$power1),
                                 paste(levels(model.mfp4$fptable$power1),
                                       levels(model.mfp4$fptable$power2))))
  
  row.names(mfptest) <- c("Not in Model", "Linear", "m=1", "m=2")
  
  mfptest
}


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

