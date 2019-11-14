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

APS <- read_delim("~/Desktop/STAT 229/APS/APS.txt","\t", escape_double = FALSE, trim_ws = TRUE)
## columns to be converted to factors
factor_cols <- c("PLACE", "PLACE3", "RACE", "GENDER", "NEURO", "EMOT", "DANGER", "ELOPE", "CUSTD", "VIOL")
## convert columns
APS %<>% mutate(ID = as.character(ID),
                PD = ifelse(PLACE < 2, 0, 1),
                BEHAV_D = factor(ifelse(BEHAV < 7, 0, 1)),
                DANGER_D = factor(ifelse(DANGER < 1, 0, 1)),
                NEURO_3 = factor(ifelse(NEURO == 3 , 1,0)))
APS %<>% mutate_at(factor_cols, factor)

covariates <- c("AGE", "RACE", "GENDER", "NEURO_3", "EMOT", "DANGER_D", "ELOPE", "LOS", "BEHAV_D", "CUSTD", "VIOL")
for (cov in covariates) {
  form <- paste("PD ~", cov)
  model <- glm(as.formula(form), data=APS, family=binomial)
  print(summary(model)$coef)
  cat("\n")
}

f1 = glm(PD ~ AGE+RACE+GENDER+EMOT+DANGER_D+ELOPE+BEHAV_D+CUSTD+VIOL,family = 'binomial',data=APS)
summary(f1)
# Remove ELOPE
f2 = glm(PD ~ AGE+RACE+GENDER+EMOT+DANGER_D+BEHAV_D+CUSTD+VIOL,family = 'binomial',data=APS)
summary(f2)
# Remove VIOL
f3 = glm(PD ~ AGE+RACE+GENDER+EMOT+DANGER_D+BEHAV_D+CUSTD,family = 'binomial',data=APS)
summary(f3)
# Remove EMOT
f3 = glm(PD ~ AGE+RACE+GENDER+DANGER_D+BEHAV_D+CUSTD,family = 'binomial',data=APS)
summary(f3)

# ADD NEURO_3
f4 = glm(PD ~ AGE+RACE+GENDER+DANGER_D+BEHAV_D+CUSTD+NEURO_3,family = 'binomial',data=APS,verbose = TRUE)
summary(f4)

# ADD NEURO_3
f4 = mfp(PD ~ fp(AGE)+RACE+GENDER+DANGER_D+BEHAV_D+CUSTD+NEURO_3,family = 'binomial',data=APS,verbose = TRUE)
summary(f4)
# keep age as linear, as the deviance is not significantly different from m=2

f4 = glm(PD ~ AGE+RACE+GENDER+DANGER_D+BEHAV_D+CUSTD+NEURO_3,family = 'binomial',data=APS)
summary(f4)

# test all interactions with Age/race/gender
# add in one at a time to all main effects at at least .05
# those significant for univariate interacrtions are added to the main model, then remove one by one.
# hopefully only 1 or 2 interactions in final model



covariates <- c("AGE*RACE", "AGE*GENDER", "AGE*DANGER_D", "AGE*BEHAV_D", "AGE*CUSTD", "AGE*NEURO_3", "RACE*GENDER", "RACE*DANGER_D", "RACE*BEHAV_D","RACE*CUSTD", "RACE*NEURO_3", "GENDER*DANGER_D", "GENDER*BEHAV_D","GENDER*CUSTD","GENDER*NEURO_3")
for (cov in covariates) {
  form <- paste("PD ~ AGE+RACE+GENDER+DANGER_D+BEHAV_D+CUSTD+NEURO_3+", cov)
  model <- glm(as.formula(form), data=APS, family=binomial)
  print(summary(model)$coef)
  cat("\n")
}

# Only add interaction of age*race as it is the only interaction that comes back significant
f4 = glm(PD ~ AGE+RACE+GENDER+DANGER_D+BEHAV_D+CUSTD+NEURO_3+AGE*RACE,family = 'binomial',data=APS)
summary(f4)


# Odds ratio for age and interaction as age is the only continuous
or_glm(data = APS, model = f4, incr = list(AGE = 2))

or_table <- function(fit, values) {
  or <- exp(coef(fit)[3] + values*coef(fit)[9])
  varb1 <- vcov(fit)[3,3]
  varb3 <- vcov(fit)[9,9]
  covb1b3 <- vcov(fit)[[3,9]]
  
  se <- sqrt(varb1 + values^2 * varb3 + 2*values*covb1b3)
  ## print(se)
  pm <- function(x,y) {
    data.frame(high=(x - y),low=( x + y))
  }
  
  ##    exp(pm(log(or), 1.96*se))
  data.frame(OR=or, CONF=exp(pm(log(or), 1.96*se)))
}

or_table(f4,c(12,14,16,18))


vcov(f4)

or_table <- function(fit) {
  or <- exp(2*(coef(fit)[2] + coef(fit)[9]))

  varb1 <- vcov(fit)[2,2]
  varb3 <- vcov(fit)[9,9]
  covb1b3 <- vcov(fit)[[2,9]]
  se <- sqrt(4 * (varb1+varb3 + 2*covb1b3))
  ## print(se)
  pm <- function(x,y) {
    data.frame(high=(x - y),low=( x + y))
  }
  
  ##    exp(pm(log(or), 1.96*se))
  data.frame(OR=or, CONF=exp(pm(log(or), 1.96*se)))
}

or_table(f4)



set.seed(123)
n <- 508
x <- rnorm(n)
y <- rbinom(n,1,plogis(0.1+.5*x))
hos <- hoslem.test(f4$y, fitted(f4))
hos
hos$observed
hos$expected
# X-squared = 7.947, df = 8, p-value = 0.4387
# High p-value means accept our null hypothesis that the model is good. There is no evidence to suggest it is bad



prob=predict(f4,type=c("response"))
APS$prob=prob
library(pROC)
g <- roc(PD ~ prob, data = APS)
plot(g)
g$auc
# Area under the curve: 0.8717
# excelent discrimination




















