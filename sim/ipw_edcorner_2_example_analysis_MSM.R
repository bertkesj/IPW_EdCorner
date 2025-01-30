# example: a series of hypothetical interventions, getting hazard ratio and estimating risk at specified intervention from the MSM
library(dplyr)
library(Hmisc)
library(survival)
load("data/sim_cohort.RData")

#1) data process
source('sim/1-data_process.R')

#1b) define reasonable limits
sourc('sim/1b-check_define_reasonable_limits')

# do all analysis: note this is wrapped in a function to facilitate bootstrapping
doipw = function(sim_cohort, limits){
  
  # Create Clones
  source('sim/2-create_clones.R')
  combined_wtd_data <- clones(data = sim_cohort,
                              limit_var = x,
                              limits = c(Inf, limits),
                              baseline_formula = ~ age + rcspline.eval(age, nk=4) + 
                                wagestatus + male + race,
                              fu_formula = ~ age + rcspline.eval(age, nk=4) + 
                                mxl + cumatworkl + 
                                wagestatus + male + race,
                              pass_thru_vars=vars(event,wtcontr,cumx,atwork,cens,wagestatus, male, race))
  
  # marginal structural policy Cox model: effect of incremental change in the limit
  # robust variance is done by ID level, but observations are done by CLONEID level
  ##########################################################################################
  # NOTE that the linear ln-HR may not be appropriate here we may chose to model with a spline: rcs(limit)
  # this should be a user option
  ##########################################################################################
  ttr1 = coxph(Surv(agein, age, event == 'd_lc') ~ limit, 
               data=filter(combined_wtd_data, 
                           limit != Inf, 
                           ipw>0), 
               id=cloneid, 
               weight=ipw, 
               cluster=id, 
               x=FALSE, y=FALSE)
  ttr2 = coxph(Surv(agein, age, event == 'd_other') ~ limit, 
               data=filter(combined_wtd_data, 
                           limit != Inf, 
                           ipw>0), 
               id=cloneid, 
               weight=ipw, 
               cluster=id, 
               x=FALSE, y=FALSE)
  
  source('sim/3-clone_summary.R')
  # multi-state model: here this just fits both cox models at once: the output is a little different, 
  # and it's only needed for risk difference calculation
  # we could likely just include this and pull out individual level model estimates out of it, 
  # which gives some flexibility
  ttrms = coxph(Surv(agein, age, event) ~ limit,
                data=filter(combined_wtd_data, 
                            limit != Inf, 
                            ipw>0), 
                id=cloneid, 
                weight=ipw, 
                cluster=id, 
                x=FALSE, y=FALSE)
  # Output: two marginal structural models (one for each outcome) and some basic diagnostics
  list(msmr1=ttr1, msmr2=ttr2, msmms=ttrms, dx=wtdx)
}



# 4) fit the MSM
# Note: this is a big model: on my computer I have to increase the "vector memory limit" for the robust variance estimator
mem.maxVSize(vsize = 32768)
msm_estimates = doipw(sim_cohort, limits)
# check diagnostics
msm_estimates$dx

# ln-HR: effect of one unit change in the exposure limit (Marginal structural policy model)
ttr = msm_estimates$msmr1
summary(ttr)

# 5)  survival curve at specified limits based on predicting from MSM (this is more efficient than the Aalen-Johansen estimator in the other analysis at a single intervention)
#  5a) Use Aalen-Johansen estimator to get risk for each outcome at the "natural course"
# note that this could require IPCW for censoring due to loss to follow-up, but these data do not have that issue
tt0 = survfit(Surv(agein, age, event)~1, data=sim_cohort, id=id)

riskdf0 = data.frame(tt0$time, tt0$pstate)
names(riskdf0) = c("age0", "surv0", "risk_other0", "risk_lc0")
tail(riskdf0)

#5b) Risk from the MSM (an MS Cox model, or MS multistate model to account for competing risks)
# risk at limit of 2.0
ccr = survfit(msm_estimates$msmms, newdata = data.frame(limit=2.0))
(riskdf1 = data.frame(
  age=ccr$time, 
  survint = ccr$pstate[,,1],
  risk_otherint = ccr$pstate[,,2],
  risk_lcint = ccr$pstate[,,3]
))
tail(riskdf1)

# risk difference across age range (may need to match the two datasets on age first)
rdlc = riskdf1$risk_lcint-riskdf0$risk_lc0
rddf = cbind((riskdf0), (riskdf1), risk_difference_lc = (rdlc))
tail(rddf)



# risk curves for lung cancer (intervention in red)
plot(rddf$age, rddf$risk_lc0, type="s", ylab="risk", xlab="age", xlim=c(16,90))
lines(rddf$age, rddf$risk_lcint, col=2, type="s")

# plotting risk difference over time
plot(rddf$age, rddf$risk_difference_lc, ylab="risk difference", xlab="age", type="s", xlim=c(16,90))