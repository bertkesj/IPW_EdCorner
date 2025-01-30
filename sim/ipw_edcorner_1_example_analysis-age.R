# example: risk under a single hypothetical intervention, and getting hazard ratio compared with no intervention

library(dplyr)
library(Hmisc)
library(survival)
load("data/sim_cohort.RData")

# Read in data
source('sim/1-data_process.R')

# Create Clones
source('sim/2-create_clones.R')
combined_wtd_data <- clones(data = sim_cohort,
                             limit_var = x,
                             limits = c(Inf, 2),
                             baseline_formula = ~ age + rcspline.eval(age, nk=4) + 
                               wagestatus + male + race,
                             fu_formula = ~ age + rcspline.eval(age, nk=4) + 
                               mxl + cumatworkl + 
                               wagestatus + male + race,
                            pass_thru_vars=vars(event,wtcontr,cumx,atwork,cens,wagestatus, male, race))

#6) estimate risk under an exposure limit
#  6a) Use Aalen-Johansen estimator to get risk for each outcome at limit
tt = survfit(Surv(agein, age, event) ~ 1, 
             data=filter(combined_wtd_data,
                         limit == 2),
             id=cloneid, 
             weight=ipw)

riskdf1 = data.frame(tt$time, tt$pstate)
names(riskdf1) = c("age", "surv", "risk_other", "risk_lc")
tail(riskdf1)

#  6a) Use Aalen-Johansen estimator to get risk for each outcome at the "natural course"
#      note that this could require IPCW for censoring due to loss to follow-up, 
#      but these data do not have that issue
tt0 = survfit(Surv(agein, age, event) ~ 1, 
              data=filter(combined_wtd_data,
                          limit == Inf),
              id=id)

riskdf0 = data.frame(tt0$time, tt0$pstate)
names(riskdf0) = c("age0", "surv0", "risk_other0", "risk_lc0")
tail(riskdf0)

# 7) estimate effect of introducing an exposure limit (impact of hypothetical intervention vs. doing nothing)
# 7a) Hazard ratio from weighted Cox model (note use of robust variance)
coxph(Surv(agein, age, event == 'd_lc') ~ factor(limit,
                                                 c(Inf, 2)), 
      data=filter(combined_wtd_data, ipw != 0), 
      weight=ipw, 
      id=cloneid, 
      cluster=id) %>% 
  summary %>% print


# 7b) Risk difference at a few illustrative ages from weighted Aalen-Johansen estimator
#     Note: this will be negative if the exposure is harmful

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
