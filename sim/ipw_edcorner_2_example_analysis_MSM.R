# example: a series of hypothetical interventions, getting hazard ratio and estimating risk at specified intervention from the MSM
library(dplyr)
library(Hmisc)
library(survival)
load("data/sim_cohort.RData")

#1) data process
source('sim/1-data_process.R')


#1) set grid of limits based on regular values of empirical distribution
# the lowest limit must include some uncensored deaths
# find the max exposure over a career: 
#  we can only explore interventions in which some individuals are observed to follow limits
#  so here we find a regular grid of values for the maximum annual exposure
maxexposures = as.numeric(tapply(sim_cohort$x, sim_cohort$id, max))

# Finding 30 potential occupational limits for a marginal structural model
# strategy: pick limits that are well within the bounds of the data
(limits = quantile(maxexposures[maxexposures>0], seq(.20,0.95,length.out=30))) # 30 limits all within the range of observed lifetime max exposures
# want to have at least 20% remain uncensored,

# check number of deaths by limit
deaths = data.frame(
  maxexposures = maxexposures,
  numd1 = as.numeric(tapply(sim_cohort$d1, sim_cohort$id, max)),
  numd2 = as.numeric(tapply(sim_cohort$d2, sim_cohort$id, max))
)
tot_events = do.call(rbind, lapply(limits, function(x){
  data.frame(
    obs=sum(deaths$maxexposures < x), 
    limit = x,
    tot_d1=sum(deaths$numd1[deaths$maxexposures < x]), 
    tot_d2=sum(deaths$numd2[deaths$maxexposures < x]),
    r_d1=sum(deaths$numd1[deaths$maxexposures < x])/sum(deaths$maxexposures < x), 
    r_d2=sum(deaths$numd2[deaths$maxexposures < x])/sum(deaths$maxexposures < x)
  )
}))

tot_events

# feasible limits
limits = filter(tot_events, tot_d1 > 10 & tot_d2 > 10)$limit

# do all analysis: note this is wrapped in a function to facilitate bootstrapping
doipw = function(sim_cohort, limits){
  
  # Create Clones
  source('sim/2-create_clones.R')
  
  # marginal structural policy Cox model: effect of incremental change in the limit
  # robust variance is done by ID level, but observations are done by CLONEID level
  ##########################################################################################
  # NOTE that the linear ln-HR may not be appropriate here we may chose to model with a spline: rcs(limit)
  # this should be a user option
  ##########################################################################################
  ttr1 = coxph(Surv(agein, age, d1)~limit, data=filter(combined_wtd_data, ipw>0), 
               id=cloneid, weight=ipw, cluster=id, x=FALSE, y=FALSE)
  ttr2 = coxph(Surv(agein, age, d2)~limit, data=filter(combined_wtd_data, ipw>0), 
               id=cloneid, weight=ipw, cluster=id, x=FALSE, y=FALSE)
  # multi-state model: here this just fits both cox models at once: the output is a little different, and it's only needed for risk difference calculation
  # we could likely just include this and pull out individual level model estimates out of it, which gives some flexibility
  ttrms = coxph(Surv(agein, age, event)~limit, data=filter(combined_wtd_data, ipw>0), 
                id=cloneid, weight=ipw, cluster=id, x=FALSE, y=FALSE)
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

sim_cohort$event <- factor(sim_cohort$d2 + sim_cohort$d1*2, 0:2, labels=c("censor", "d_other", "d_lc"))
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