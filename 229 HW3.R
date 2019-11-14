library(mfp)
library(readr)
library(splines)
library(dplyr)

BURN <- read_delim("~/Desktop/STAT 229/BURN/BURN1000.txt","\t", escape_double = FALSE, trim_ws = TRUE)

factor_cols <- c("DEATH","RACEC","INH_INJ")
BURN <- BURN %>% mutate_at(factor_cols,factor)

f = glm(DEATH ~ AGE+TBSA+RACEC+INH_INJ+AGE*TBSA,family = 'binomial',data=BURN)
summary(f)
f = glm(DEATH ~ AGE+TBSA+RACEC+INH_INJ+AGE*RACEC,family = 'binomial',data=BURN)
summary(f)
f = glm(DEATH ~ AGE+TBSA+RACEC+INH_INJ+AGE*INH_INJ,family = 'binomial',data=BURN)
summary(f)
f = glm(DEATH ~ AGE+TBSA+RACEC+INH_INJ+TBSA*RACEC,family = 'binomial',data=BURN)
summary(f)
f = glm(DEATH ~ AGE+TBSA+RACEC+INH_INJ+TBSA*INH_INJ,family = 'binomial',data=BURN)
summary(f)
f = glm(DEATH ~ AGE+TBSA+RACEC+INH_INJ+RACEC*INH_INJ,family = 'binomial',data=BURN)
summary(f)

f = glm(DEATH ~ AGE+TBSA+RACEC+INH_INJ+AGE*INH_INJ,family = 'binomial',data=BURN)
summary(f)


knots <- quantile(BURN$AGE, probs = c(.05,.35,.65,.95))
splineFunc <- glm(BURN$DEATH ~ ns(BURN$AGE, knots = knots)+BURN$TBSA+ BURN$RACEC+ BURN$INH_INJ)
summary(splineFunc)
