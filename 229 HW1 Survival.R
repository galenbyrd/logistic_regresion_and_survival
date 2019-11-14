library(dplyr)
library(survminer)
library(readxl)
library(survival)
library(survMisc)

DATA <- read_excel("~/Desktop/STAT 229/Stat 229a WHAS 500.xls")

fit <- survfit(Surv(lenfol, fstat) ~ gender, data = DATA)
summary(fit)
quantile(fit)


by(DATA$lenfol,DATA$gender,quantile)

quantile(DATA$lenfol)
ggsurvplot(fit, data = DATA, pval = TRUE)
plot(fit, lty = 1:3)

survdiff(Surv(lenfol, fstat) ~ gender, data = DATA)
logrank_test(Surv(lenfol, fstat) ~ gender, data = DATA)
coxph <- coxph(Surv(lenfol, fstat) ~ gender, data = DATA)
summary(coxph)

t1 <-ten(fit)
capture.output(comp(t1))
d <- head(data.frame(chisq=attr(t1, "lrt")$chiSq,p=attr(t1,"lrt")$pChisq),4)
row.names(d) <-c('Log-Rank','Wilcoxon','Tarone-Ware','Peto-Prentice')
d










shitSample <- data.frame(ID =c(1:12),
                         time = c(1.2,3.4,5.0,5.1,6.1,7.1,0.4,1.2,4.3,4.9,5.0,5.1),
                         gender = c(0,0,0,0,0,0,1,1,1,1,1,1),
                         died = c(1,1,0,1,1,1,1,1,1,1,1,0),
                         censor = c(0,0,1,0,0,0,0,0,0,0,0,1))
fit2 <- survfit(Surv(time, died) ~ gender, data = shitSample)
summary(fit2)
quantile(fit2)
ggsurvplot(fit2, data = shitSample)






men <- c(1.2,3.4,5.0,5.1,6.1,7.1)
menp <- c(1.0,0.833,0.666,0.5,0.333,0.167,0)
quantile(men)
women <- c(0.4,1.2,4.3,4.9,5.0,5.1)
womenp <- c(1.0,0.833,0.666,0.5,0.333,0.167,0.167)
quantile(women)



f <- stepfun(men,menp)
plot(f)
f <- stepfun(women,womenp)
plot(f)