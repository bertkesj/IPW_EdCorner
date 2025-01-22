# example: a series of hypothetical interventions, getting hazard ratio and estimating risk at specified intervention from the MSM
library(dplyr)
library(Hmisc)
library(survival)
load("data/sim_cohort.RData")

#1) data process
sim_cohort <- sim_cohort %>%
  group_by(id) %>%
  mutate(
    one = 1,
    time = cumsum(one),
    timein = time-1,
    agein = age-1,
    xl = lag(x, default=0),
    cxl = lag(cumx, default=0),
    cumatwork = cumsum(atwork),
    cumatworkl = lag(cumatwork, default=0),
    atworkl = lag(atwork, default=0),  
    leftwork = as.numeric(atworkl==1 & atwork == 0),
    cxl2 = lag(cxl, default=0),  
  ) %>%
  select(-one) %>%
  ungroup()
sim_cohort <- sim_cohort %>%
  group_by(id) %>%
  mutate(
    male = as.numeric(gender == 'M')
  ) %>%
  ungroup()


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
  
  #1b) create copy of the dataset for each limit
  clones <- do.call(rbind,
                    lapply(1:length(limits), function(x) {  
                      df = sim_cohort                  
                      df$limit = as.numeric(limits[x])
                      df$cloneid = paste0(x, "_", df$id)
                      df                    
                    }
                    ))
  #2) artificially censor data (keeping first observation that is censored)
  cens_data <- clones %>%
    group_by(cloneid) %>%
    mutate(
      cens = pmin(1,cumsum(x > limit)),
      drop = pmin(2, cumsum(cens))
    ) %>%
    group_by() %>%
    filter(
      drop < 2
    )
  
  # 3) create 1/0 weights to use for confounding/censoring during follow-up (alternative is to create new datasets)
  cens_data <- group_by(cens_data, cloneid) %>%
    mutate(
      one = 1,
      nobs = cumsum(one),
      fobs____ = as.numeric(nobs == 1)
    ) %>%
    mutate(
      conf_weight = as.numeric(nobs  == 1), # this duplicates fobs____ but is kept for clarity
      fu_weight = as.numeric((nobs  > 1) & (atwork == 1)),    
      dconf = 0,
      dcens = 0,
      intervened = 0,
    ) %>%
    select(-c(one, nobs)) %>%
    ungroup()
  # check: which(tapply(sim_cohort$x, sim_cohort$id, max) < limit & tapply(sim_cohort$atwork, sim_cohort$id, sum)>3)[2] # index 427
  # check: names(tapply(sim_cohort$x, sim_cohort$id, max))[427]
  # check: print(select(cens_data, c(id, x, limit, conf_weight, fu_weight, cens, fobs____)) %>% filter(id==427), n=50)
  # check: print(select(filter(sim_cohort, id==427), c(id, x, atwork)))
  
  
  #4) fit censoring models (can include pre-baseline exposure in practice)
  
  # fit models if there is more than a small amount of censoring (seems to work either way, but this avoids convergence problems)
  for (l in limits){
    limidx = which(cens_data$limit == l)
    tempdat = cens_data[limidx,]
    tempdat$mxl = tempdat$cxl / (tempdat$cumatworkl + as.numeric(tempdat$time==1))
    agekn0 = attr(rcspline.eval(filter(tempdat, conf_weight==1)$age, nk = 4), "knots")
    yearkn0 = attr(rcspline.eval(filter(tempdat, conf_weight==1)$year, nk = 4), "knots")
    agekn = attr(rcspline.eval(filter(tempdat, fu_weight==1)$age, nk = 4), "knots")
    yearkn = attr(rcspline.eval(filter(tempdat, fu_weight==1)$year, nk = 4), "knots")
    #
    if(sum(tempdat[tempdat$conf_weight==1,"cens"]) > 10){
      summary(confnmod <- glm(cens ~ age + rcspline.eval(age, knots=agekn0) , data = tempdat, weight=conf_weight, family=binomial()))
      summary(confdmod <- glm(cens ~ age + rcspline.eval(age, knots=agekn0) + wagestatus + male + race, data = tempdat, weight=conf_weight, family=binomial()))
      # need to use base R operations here
      # only get non-zero predictions if "conf_weight" == 1 (eligible to be "censored" at baseline)
      cens_data[limidx,"nconf"] = tempdat$conf_weight*as.numeric(predict(confnmod, type="response"))
      cens_data[limidx,"dconf"] = tempdat$conf_weight*as.numeric(predict(confdmod, type="response"))
    } else{
      cens_data[limidx,"nconf"] = 0
      cens_data[limidx,"dconf"] = 0
    }
    if(sum(tempdat[tempdat$fu_weight==1,"cens"]) > 1){
      summary(censnmod <- glm(cens ~ age + rcspline.eval(age, knots=agekn), data = tempdat, weight=fu_weight, family=binomial()))
      summary(censdmod <- glm(cens ~  mxl + cumatworkl + age + rcspline.eval(age, knots=agekn) + wagestatus + male + race, data = tempdat, weight=fu_weight, family=binomial()))
      # only get non-zero predictions if "fu_weight" == 1 (eligible to be censored during follow-up)
      cens_data[limidx,"ncens"] = tempdat$fu_weight*as.numeric(predict(censnmod, type="response")) 
      cens_data[limidx,"dcens"] = tempdat$fu_weight*as.numeric(predict(censdmod, type="response")) 
    } else{
      cens_data[limidx,"ncens"] = 0
      cens_data[limidx,"dcens"] = 0
    }
  }
  
  
  # create final weights
  combined_wtd_data <- cens_data %>% 
    mutate(
      wtcontr = case_when(
        # weight: (stabilized) inverse probability of NOT being censored
        # note regarding unstabilized weights: Cain et al 2010 show that weight stabilization is not guaranteed to reduce variance
        ((fobs____ == 1) & (atwork==1)) ~ (1-cens)*(1-nconf)/(1-dconf),
        ((fobs____ == 0) & (atwork==1)) ~ (1-cens)*(1-ncens)/(1-dcens),
        .default=1
      )
    ) %>%
    select(id,cloneid,agein,age,d1,d2,wtcontr,x,cumx,atwork,cens,limit, wagestatus, male, race, fobs____) %>% # optional step here to reduce comp. burden
    group_by(cloneid) %>%
    mutate(
      ipw = cumprod(wtcontr)
    ) %>%
    group_by() %>%
    mutate(
      # weight truncation at 10.0 following Cain et al 2010 (should be a user option)
      ipw_trunc = pmin(ipw, 10.0)
    )
  
  # check: summary(select(combined_wtd_data, c(ipw, ipw_trunc)))
  # check: print(select(combined_wtd_data, c(id, age, atwork, x, cens, wtcontr, ipw, ipw_trunc)) %>% filter(id==427), n=10)
  
  # check mean weights by intervention
  # note mean is taken across all possible observations, even those with weights = 0
  N = nrow(sim_cohort)
  Nid = length(unique(sim_cohort$id))
  wtdx = data.frame(
    limit = limits,
    mean_cumx = do.call(c,lapply(limits, function(x) mean(combined_wtd_data[combined_wtd_data$limit == x  & combined_wtd_data$cens == 0,]$cumx))),
    n_conf = do.call(c,lapply(limits, function(x) sum(combined_wtd_data[combined_wtd_data$limit == x & combined_wtd_data$fobs____==1,]$cens))),
    n_cens = do.call(c,lapply(limits, function(x) sum(combined_wtd_data[combined_wtd_data$limit == x & combined_wtd_data$fobs____==0,]$cens))),
    n_d1w = do.call(c,lapply(limits, function(x) sum(combined_wtd_data[combined_wtd_data$limit == x & combined_wtd_data$cens == 0,]$ipw*combined_wtd_data[combined_wtd_data$limit == x & combined_wtd_data$cens == 0,]$d1))),
    n_d2w = do.call(c,lapply(limits, function(x) sum(combined_wtd_data[combined_wtd_data$limit == x & combined_wtd_data$cens == 0,]$ipw*combined_wtd_data[combined_wtd_data$limit == x & combined_wtd_data$cens == 0,]$d2))),
    mean_ipw = do.call(c,lapply(limits, function(x) mean(combined_wtd_data[combined_wtd_data$limit == x & combined_wtd_data$cens == 0,]$ipw)))
  )
  
  combined_wtd_data$event <- factor(combined_wtd_data$d2 + combined_wtd_data$d1*2, 0:2, labels=c("censor", "d_other", "d_lc"))
  
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