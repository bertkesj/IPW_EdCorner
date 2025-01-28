# example: risk under a single hypothetical intervention, and getting hazard ratio compared with no intervention

library(dplyr)
library(Hmisc)
library(survival)
load("data/sim_cohort.RData")



# data process
source('sim/1-data_process.R')

# checking exposure at work
quantile(filter(sim_cohort, atwork==1)$x, 0:20/20)

#1-0) set limit (limits) to generalize to a list of limits
limits <- c(2)

# Create Clones
source('sim/2-create_clones.R')

#STEVE edit: rename 
cens_data <- combined_wtd_data

#6) estimate risk under an exposure limit

#  6a) Use Aalen-Johansen estimator to get risk for each outcome at limit
tt = survfit(Surv(agein, age, event)~1, data=cens_data, id=cloneid, weight=ipw)

riskdf1 = data.frame(tt$time, tt$pstate)
names(riskdf1) = c("age", "surv", "risk_other", "risk_lc")
tail(riskdf1)

#  6a) Use Aalen-Johansen estimator to get risk for each outcome at the "natural course"
# note that this could require IPCW for censoring due to loss to follow-up, but these data do not have that issue

sim_cohort$event <- factor(sim_cohort$d2 + sim_cohort$d1*2, 0:2, labels=c("censor", "d_other", "d_lc"))
tt0 = survfit(Surv(agein, age, event)~1, data=sim_cohort, id=id)

riskdf0 = data.frame(tt0$time, tt0$pstate)
names(riskdf0) = c("age0", "surv0", "risk_other0", "risk_lc0")
tail(riskdf0)

#7) estimate effect of introducing an exposure limit (impact of hypothetical intervention vs. doing nothing)
# 7a) Hazard ratio from weighted Cox model (note use of robust variance)
sim_cohort$exposed = 0
cens_data$exposed = 1
sim_cohort$ipw = 1
sim_cohort$cloneid = paste0("cloneobs_", sim_cohort$id)
combined_data <- bind_rows(sim_cohort, cens_data)

coxph(Surv(agein, age, d1)~exposed, 
      data=filter(combined_data, ipw != 0), weight=ipw, 
      id=cloneid, cluster=id) %>% 
  summary %>% print


# 7b) Risk difference at a few illustrative ages from weighted Aalen-Johansen estimator
# Note: this will be negative if the exposure is harmful
# risk difference at age 80 
riskdf1[riskdf1$age==80, "risk_lc"] - riskdf0[riskdf0$age==80, "risk_lc0"]
# risk difference at age 90
riskdf1[riskdf1$age==90, "risk_lc"] - riskdf0[riskdf0$age==90, "risk_lc0"]


######### BOOTSTRAPPING #########
# We can optionally do bootstrapping for confidence intervals on the HR which will generally be better
# than robust confidence intervals except in very large sample sizes. Bootstrapping
# is needed for the Aalen-Johansen estimator
# Code to do this is not explicitly done here but it is a straightforward set of steps:
# 1. sample, with replacement, 10,000 individuals from the population (original size = 10,000)
# 2. Carry out clone-censor-weight + estimate effects on the sample from #1
# 3. Repeat 1 and 2 200+ times, recording effect estimate for each sample
# 4. Bootstrap standard error is the standard deviation of the effect estimates across the 200+ bootstrap samples
