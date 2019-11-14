library(dplyr)
library(readxl)
library(survival)
library(survMisc)
library(ggplot2)

DATA <- read_excel("~/Desktop/STAT 229/Stat 229a WHAS 500.xls")
DATA <- mutate(DATA, lenfolyr=lenfol/365)


fit <- coxph(Surv(lenfol, fstat) ~ age + gender + hr + sysbp + diasbp + bmi + cvd + miord + mitype, data = DATA)
summary(fit)


fit2 <- coxph(Surv(lenfol, fstat) ~ age + hr + diasbp + bmi, data = DATA)
summary(fit2)

plot(survfit(fit2, newdata = data.frame(age=0,hr=0,diasbp=0,bmi=0)), data = DATA, xlab='Time (days)',ylab='Survival Probability',main='Survival Function')
plot(survfit(fit2), xlab='Time (days)',ylab='Survival Probability',main='Survival Function', centered= FALSE)

DATA <- mutate(DATA, agemed = median(age), hrmed = median(hr), diasbpmed = median(diasbp), bmimed = median(bmi))
fit3 <- coxph(Surv(lenfol, fstat) ~ agemed + gender + hr + diasbp + bmi, data = DATA)
summary(fit3)
ggsurvplot(survfit(fit3), data = DATA, pval = TRUE)
