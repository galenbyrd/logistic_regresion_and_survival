library(mfp)
library(readr)
library(splines)
library(dplyr)
library(tidyverse)
library(magrittr)
library(rms)
library(mfp)
library(MASS)
library(ResourceSelection)
library(oddsratio)

DATA <- read_delim("~/Desktop/STAT 229/Byrd4.out","\t", escape_double = FALSE, trim_ws = TRUE) 
summary(DATA)
str(DATA)

## columns to be converted to factors
factor_cols <- c("female", "insured", "hosp_region", "hosp_trauma", "hosp_teach", "shock", "zipinc_qrtl", "gsw_int", "gun_typ", "cannabis", "gsw_ass", "gsw_suicide")

## convert columns
DATA %<>% mutate(id = as.character(id), hosp_id = as.character(hosp_id))
DATA %<>% mutate_at(factor_cols, factor)

covariates <- c("age", "female", "insured", "hosp_region", "hosp_trauma", "shock", "year", "zipinc_qrtl", "gsw_int", "gun_typ", "PrDeathICD9", "TMPM", "cannabis", "gsw_ass", "gsw_suicide")
for (cov in covariates) {
  form <- paste("died ~", cov)
  model <- glm(as.formula(form), data=DATA, family=binomial)
  print(summary(model)$coef)
  #print(confint(model))
  #summary_table(model)
  cat("\n")
}

table(DATA$hosp_teach,DATA$died)

# DO NOT INCLUDE: hosp_region, hosp_teach, year, gsw_int, gun_typ, gsw_ass AND TMPM as it is same as PrDeath
f1 = glm(died ~ age+female+insured+hosp_trauma+shock+zipinc_qrtl+PrDeathICD9+cannabis+gsw_suicide,family = 'binomial',data=DATA)
summary(f1)
# Remove zipinc_qrtl

f2 = glm(died ~ age+female+insured+hosp_trauma+shock+PrDeathICD9+cannabis+gsw_suicide,family = 'binomial',data=DATA)
summary(f2)
# Remove hosp_trauma

f3 = glm(died ~ age+female+insured+shock+PrDeathICD9+cannabis+gsw_suicide,family = 'binomial',data=DATA)
summary(f3)
# ALL SIGNIFICANT SO -> try adding in covariates deleted in univariate case

covariates <- c("hosp_region", "year",  "gsw_int", "gun_typ", "gsw_ass")
for (cov in covariates) {
  form <- paste("died ~ age+female+insured+shock+PrDeathICD9+cannabis+gsw_suicide+", cov)
  model <- glm(as.formula(form), data=DATA, family=binomial)
  print(summary(model)$coef)
  cat("\n")
}
# ADD gsw_ass
f3 = glm(died ~ age+female+insured+shock+PrDeathICD9+cannabis+gsw_suicide+gsw_ass,family = 'binomial',data=DATA)
summary(f3)
summary_table(f3)


DATA %<>% mutate(hosp_west = factor(ifelse(hosp_region == 4, 0, 1)),
                 hosp_trauma_d = factor(ifelse(hosp_trauma == 0 | hosp_trauma == 1, 0, 1)),
                 zipinc_qrtl_high = factor(ifelse(zipinc_qrtl == 4, 0, 1)),
                 PrDeathICD9log = log((PrDeathICD9/0.1)),
                 PrDeathICD9logSquared = I(log((PrDeathICD9/0.1))^2),
                 year_d = factor(ifelse(year <2013 , 0, 1)))

f4 = glm(died ~ age+female+insured+shock+PrDeathICD9+cannabis+gsw_suicide+gsw_ass+hosp_west,family = 'binomial',data=DATA)
summary(f4)
# Add hosp_west

f5 = glm(died ~ age+female+insured+shock+PrDeathICD9+cannabis+gsw_suicide+gsw_ass+hosp_trauma_d,family = 'binomial',data=DATA)
summary(f5)
# Add hosp_trauma_d

f6 = glm(died ~ age+female+insured+shock+PrDeathICD9+cannabis+gsw_suicide+gsw_ass+zipinc_qrtl_high,family = 'binomial',data=DATA)
summary(f6)
# Add zipinc_qrtl_high

f7 = glm(died ~ age+female+insured+shock+PrDeathICD9+cannabis+gsw_suicide+gsw_ass+year_d,family = 'binomial',data=DATA)
summary(f7)
# Add year_d

f8 = glm(died ~ age+female+insured+shock+PrDeathICD9+cannabis+gsw_suicide+gsw_ass+hosp_west+hosp_trauma_d+zipinc_qrtl_high+year_d,family = 'binomial',data=DATA)
summary(f8)
# Remove hosp_west

f9 = glm(died ~ age+female+insured+shock+PrDeathICD9+cannabis+gsw_suicide+gsw_ass+hosp_trauma_d+zipinc_qrtl_high+year_d,family = 'binomial',data=DATA)
summary(f9)
summary_table(f9)
# Good, now try frac poly

f10 = mfp(died ~ fp(age)+female+insured+shock+PrDeathICD9+cannabis+gsw_suicide+gsw_ass+hosp_trauma_d+zipinc_qrtl_high+year_d,family = 'binomial',data=DATA,verbose = TRUE)
summary(f10)
# Not needed for age

f10 = mfp(died ~ age+female+insured+shock+fp(PrDeathICD9)+cannabis+gsw_suicide+gsw_ass+hosp_trauma_d+zipinc_qrtl_high+year_d,family = 'binomial',data=DATA,verbose = TRUE)
summary(f10)
# Use frac poly for PrDeathICD9

f11 = glm(died ~ age+female+insured+shock+PrDeathICD9log+PrDeathICD9logSquared+cannabis+gsw_suicide+gsw_ass+hosp_trauma_d+zipinc_qrtl_high+year_d,family = 'binomial',data=DATA)
summary(f11)

# NOW TEST INTERACTIONS
covariates <- c("age*female", "age*insured",  "age*shock", "age*PrDeathICD9log", "age*PrDeathICD9logSquared", "age*cannabis", "age*gsw_suicide", "age*gsw_ass", "age*hosp_trauma_d", "age*zipinc_qrtl_high", "age*year_d",
                "female*insured", "female*shock", "female*PrDeathICD9log", "female*PrDeathICD9logSquared", "female*cannabis", "female*gsw_suicide", "female*gsw_ass", "female*hosp_trauma_d", "female*zipinc_qrtl_high", "female*year_d")
for (cov in covariates) {
  form <- paste("died ~ age+female+insured+shock+PrDeathICD9log+PrDeathICD9logSquared+cannabis+gsw_suicide+gsw_ass+hosp_trauma_d+zipinc_qrtl_high+year_d+", cov)
  model <- glm(as.formula(form), data=DATA, family=binomial)
  print(summary(model)$coef)
  cat("\n")
}
# Significant interactions: age*insured, age*gsw_suicide, age*gsw_ass

f12 = glm(died ~ age+female+insured+shock+PrDeathICD9log+PrDeathICD9logSquared+cannabis+gsw_suicide+gsw_ass+hosp_trauma_d+zipinc_qrtl_high+year_d+age*insured+age*gsw_suicide+age*gsw_ass,family = 'binomial',data=DATA)
summary(f12)

# FULL MODEL:
f150 = glm(died ~ age+female+insured+shock+PrDeathICD9log+PrDeathICD9logSquared+cannabis+gsw_suicide+gsw_ass+hosp_trauma_d+zipinc_qrtl_high+year_d+age*insured+age*gsw_suicide+age*gsw_ass,family = 'binomial',data=DATA)
summary(f150)
summary_table(f150)

# Odds ratio tables:
or_glm(data = DATA, model = f150, incr = list(age = 2, PrDeathICD9log=1 ,PrDeathICD9logSquared=4 ))

# For contunious variable PrDeathICD9
or_table <- function(fit, values) {
  b1<- 1.17145
  b2<- -0.13276
  a <- log((values + 0.1)/values)
  b <- ((log(values + 0.1))^2 - (log(values))^2)
  or <- b1 * a + b2 * b
  v1 <- vcov(fit)[6,6]
  v3 <- vcov(fit)[7,7]
  covb1b3 <- vcov(fit)[6,7]
  se <- sqrt( v1*a^2 + v3*b^2 + 2*a*b*covb1b3)
  ## print(se)
  pm <- function(x,y) {
    data.frame(high=(x - y),low=( x + y))
  }
  ##    exp(pm(log(or), 1.96*se))
  data.frame(OR=or, CONF=exp(pm(log(or), 1.96*se)))
}
or_table(f150,c(0.1,0.3,0.5,0.7,0.9))



# AGE AND INSURED
or_table <- function(fit, values) {
  or <- exp(coef(fit)[4] + values*coef(fit)[14])
  insured <- vcov(fit)[4,4]
  interact <- vcov(fit)[14,14]
  covb1b3 <- vcov(fit)[[4,14]]
  se <- sqrt(insured + values^2 * interact + 2*values*covb1b3)
  ## print(se)
  pm <- function(x,y) {
    data.frame(high=(x - y),low=( x + y))
  }
  ##    exp(pm(log(or), 1.96*se))
  data.frame(OR=or, CONF=exp(pm(log(or), 1.96*se)))
}
or_table(f150,c(8,10,12,14,16,18))
or_table <- function(fit) {
  or <- exp(2*(coef(fit)[2] + coef(fit)[14]))
  age <- vcov(fit)[2,2]
  interact <- vcov(fit)[14,14]
  covb1b3 <- vcov(fit)[[2,14]]
  se <- sqrt(4 * (age+interact + 2*covb1b3))
  ## print(se)
  pm <- function(x,y) {
    data.frame(high=(x - y),low=( x + y))
  }
  ##    exp(pm(log(or), 1.96*se))
  data.frame(OR=or, CONF=exp(pm(log(or), 1.96*se)))
}
or_table(f150)

vcov(f150)

# AGE AND gsw_suicide
or_table <- function(fit, values) {
  or <- exp(coef(fit)[9] + values*coef(fit)[15])
  suicide <- vcov(fit)[9,9]
  interact <- vcov(fit)[15,15]
  covb1b3 <- vcov(fit)[[9,15]]
  se <- sqrt(suicide + values^2 * interact + 2*values*covb1b3)
  ## print(se)
  pm <- function(x,y) {
    data.frame(high=(x - y),low=( x + y))
  }
  ##    exp(pm(log(or), 1.96*se))
  data.frame(OR=or, CONF=exp(pm(log(or), 1.96*se)))
}
or_table(f150,c(8,10,12,14,16,18))
or_table <- function(fit) {
  or <- exp(2*(coef(fit)[2] + coef(fit)[15]))
  age <- vcov(fit)[2,2]
  interact <- vcov(fit)[15,15]
  covb1b3 <- vcov(fit)[[2,15]]
  se <- sqrt(4 * (age+interact + 2*covb1b3))
  ## print(se)
  pm <- function(x,y) {
    data.frame(high=(x - y),low=( x + y))
  }
  ##    exp(pm(log(or), 1.96*se))
  data.frame(OR=or, CONF=exp(pm(log(or), 1.96*se)))
}
or_table(f150)


# AGE AND gsw_ass
or_table <- function(fit, values) {
  or <- exp(coef(fit)[10] + values*coef(fit)[16])
  assault <- vcov(fit)[10,10]
  interact <- vcov(fit)[16,16]
  covb1b3 <- vcov(fit)[[10,16]]
  se <- sqrt(assault + values^2 * interact + 2*values*covb1b3)
  ## print(se)
  pm <- function(x,y) {
    data.frame(high=(x - y),low=( x + y))
  }
  ##    exp(pm(log(or), 1.96*se))
  data.frame(OR=or, CONF=exp(pm(log(or), 1.96*se)))
}
or_table(f150,c(8,10,12,14,16,18))
or_table <- function(fit) {
  or <- exp(2*(coef(fit)[2] + coef(fit)[16]))
  age <- vcov(fit)[2,2]
  interact <- vcov(fit)[16,16]
  covb1b3 <- vcov(fit)[[2,16]]
  se <- sqrt(4 * (age+interact + 2*covb1b3))
  ## print(se)
  pm <- function(x,y) {
    data.frame(high=(x - y),low=( x + y))
  }
  ##    exp(pm(log(or), 1.96*se))
  data.frame(OR=or, CONF=exp(pm(log(or), 1.96*se)))
}
or_table(f150)




# EVALUATE
library(ResourceSelection)
## bartley (p. 140)
## not looking too hot, should be less than 10, but usually around 50
rejected <- 0
for(i in 1:100){
  newD <- sample_n(DATA, 5000, replace=TRUE)
  predicted <- predict(f150, newdata=newD, type="response")
  chi <- hoslem.test(newD$died, predicted)$statistic
  p <- pchisq(chi, 10, lower.tail=FALSE)
  #print(paste(chi, p))
  if (p < 0.05) {
    rejected <- rejected + 1
  }
}
print(rejected)


prob=predict(f150,type=c("response"))
DATA$prob=prob
library(pROC)
g <- roc(died ~ prob, data = DATA)
plot(g)
g$auc
# Area under the curve: 0.9231
# excelent discrimination


## need to have \mc macro defined as well as siunitx package: \newcommand{\mc}[1]{\multicolumn{1}{c}{#1}}
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
summary_table <- function(model) {
  require(xtable)
  summ <- summary(model)$coef
  summ <- fixp(summ)
  intervals <- confint(model)
  combined <- merge(summ, intervals, by="row.names", all.x=TRUE)
  xt <- xtable(combined, digits=c(1,1,3,4,2,3,3,3), align=c("l","l", "S", "S", "S", "S[table-format = <2.3]", "S", "S"))
  ## xt <- fixp(xt)
  names(xt) <- c("", "\\mc{Estimate}", "\\mc{Std. Error}", "\\mc{$z\\text{-value}$}", "\\mc{Pr$\\parens{> |z|}$}", "\\mc{2.5 \\%}", "\\mc{97.5 \\%}")
  ## align(xt) <- "lSSSS[table-format = <2.3]SS"
  print(xt, include.rownames=FALSE, booktabs=TRUE, sanitize.colnames.function = identity, sanitize.text.function = identity, sanitize.numbers=FALSE)
}
summary_table(f150)
