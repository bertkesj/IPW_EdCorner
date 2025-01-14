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
    cxl = lag(cumx, default=0)
  ) %>%
  select(-one) %>%
  group_by()

#1) set grid of limits based on regular values of empirical distribution
# the lowest limit must include some uncensored deaths
# find the max exposure over a career: 
#  we can only explore interventions in which some individuals are observed to follow limits
#  so here we find a regular grid of values for the maximum annual exposure
maxexposures = as.numeric(tapply(sim_cohort$x, sim_cohort$id, max))


limits = quantile(maxexposures, seq(0.033,1,by=0.033)) # 30 limits all within the range of observed lifetime max exposures

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
    group_by()
  # check: which(tapply(sim_cohort$x, sim_cohort$id, max) < limit & tapply(sim_cohort$atwork, sim_cohort$id, sum)>3)[2] # index 427
  # check: names(tapply(sim_cohort$x, sim_cohort$id, max))[427]
  # check: print(select(cens_data, c(id, x, limit, conf_weight, fu_weight, cens, fobs____)) %>% filter(id==427), n=50)
  # check: print(select(filter(sim_cohort, id==427), c(id, x, atwork)))
  
  
  #3) fit censoring models (can include pre-baseline exposure in practice)
  agekn0 = attr(rcspline.eval(filter(cens_data, conf_weight==1)$age, nk = 4), "knots")
  yearkn0 = attr(rcspline.eval(filter(cens_data, conf_weight==1)$year, nk = 4), "knots")
  agekn = attr(rcspline.eval(filter(cens_data, fu_weight==1)$age, nk = 4), "knots")
  yearkn = attr(rcspline.eval(filter(cens_data, fu_weight==1)$year, nk = 4), "knots")
  
  # fit models if there is more than a small amount of censoring (seems to work either way, but this avoids convergence problems)
  for (l in limits){
    limidx = which(cens_data$limit == l)
    tempdat = cens_data[limidx,]
    if(sum(tempdat[tempdat$conf_weight==1,"cens"]) > 10){
      summary(confdmod <- glm(cens ~rcspline.eval(tempdat$year, knots=yearkn0) + wagestatus + race + gender , data = tempdat, weight=conf_weight, family=binomial()))
      # need to use base R operations here
      # only get non-zero predictions if "conf_weight" == 1 (eligible to be "censored" at baseline)
      cens_data[limidx,"dconf"] = tempdat$conf_weight*as.numeric(predict(confdmod, type="response"))
    } else{
      cens_data[limidx,"dconf"] = 0
    }
    if(sum(tempdat[tempdat$fu_weight==1,"cens"]) > 1){
      summary(censdmod <- glm(cens ~ rcspline.eval(year, knots=yearkn) + wagestatus + xl + race + gender, data = tempdat, weight=fu_weight, family=binomial()))
      # only get non-zero predictions if "fu_weight" == 1 (eligible to be censored during follow-up)
      cens_data[limidx,"dcens"] = tempdat$fu_weight*as.numeric(predict(censdmod, type="response")) 
    } else{
      cens_data[limidx,"dcens"] = 0
    }
  }
  
  # create final weights
  cens_data_simple <- cens_data %>% 
    mutate(
      wtcontr = case_when(
        # weight: (unstabilized) inverse probability of NOT being censored
        # note regarding unstabilized weights: Cain et al 2010 show that weight stabilization doesn't work in IPCW as well as it does for static regime MSMs
        # notably the numerator model should not be fit to the 
        ((fobs____ == 1) & (atwork==1)) ~ (1-cens)*(1)/(1-dconf),
        ((fobs____ == 0) & (atwork==1)) ~ (1-cens)*(1)/(1-dcens),
        .default=1
      )
    ) %>%
    select(id,cloneid,agein,age,d1,d2,wtcontr,x,cumx,atwork,cens,limit, wagestatus, gender, race, fobs____) %>% # optional step here to reduce comp. burden
    group_by(cloneid) %>%
    mutate(
      ipw = cumprod(wtcontr)
    ) %>%
    group_by() %>%
    mutate(
      # weight truncation at 10.0 following Cain et al 2010 (should be a user option)
      ipw_trunc = pmin(ipw, 10.0)
    )
  # check: summary(select(cens_data_simple, c(ipw, ipw_trunc)))
  # check: print(select(cens_data_simple, c(id, age, atwork, x, cens, wtcontr, ipw, ipw_trunc)) %>% filter(id==427), n=10)
  
  
  
  # check mean weights by intervention
  # note mean is taken across all possible observations, even those with weights = 0
  N = nrow(sim_cohort)
  wtdx = data.frame(
    limit = limits,
    mean_cumx = do.call(c,lapply(limits, function(x) mean(cens_data_simple[cens_data_simple$limit == x  & cens_data_simple$cens == 0,]$cumx))),
    n_conf = do.call(c,lapply(limits, function(x) sum(cens_data_simple[cens_data_simple$limit == x & cens_data_simple$fobs____==1,]$cens))),
    n_cens = do.call(c,lapply(limits, function(x) sum(cens_data_simple[cens_data_simple$limit == x & cens_data_simple$fobs____==0,]$cens))),
    n_d1 = do.call(c,lapply(limits, function(x) sum(cens_data_simple[cens_data_simple$limit == x & cens_data_simple$cens == 0,]$d1))),
    n_d2 = do.call(c,lapply(limits, function(x) sum(cens_data_simple[cens_data_simple$limit == x & cens_data_simple$cens == 0,]$d2))),
    n_d1w = do.call(c,lapply(limits, function(x) sum(cens_data_simple[cens_data_simple$limit == x & cens_data_simple$cens == 0,]$ipw*cens_data_simple[cens_data_simple$limit == x & cens_data_simple$cens == 0,]$d1))),
    n_d2w = do.call(c,lapply(limits, function(x) sum(cens_data_simple[cens_data_simple$limit == x & cens_data_simple$cens == 0,]$ipw*cens_data_simple[cens_data_simple$limit == x & cens_data_simple$cens == 0,]$d2))),
    mean_ipw = do.call(c,lapply(limits, function(x) sum(cens_data_simple[cens_data_simple$limit == x & cens_data_simple$cens == 0,]$ipw)))/N,
    mean_ipwt = do.call(c,lapply(limits, function(x) sum(cens_data_simple[cens_data_simple$limit == x & cens_data_simple$cens == 0,]$ipw_trunc)))/N
  )
  
  sim_cohort$event <- factor(sim_cohort$d2 + sim_cohort$d1*2, 0:2, labels=c("censor", "d_other", "d_lc"))
  sim_cohort$intercept <- 1
  cens_data_simple$event <- factor(cens_data_simple$d2 + cens_data_simple$d1*2, 0:2, labels=c("censor", "d_other", "d_lc"))

  # marginal structural policy Cox model: effect of incremental change in the limit
  (nc = coxph(Surv(agein, age, event)~intercept, data=sim_cohort, id=id)) # trick to get survival curve from cohort
  (tt = coxph(Surv(agein, age, event)~limit, data=filter(cens_data_simple, ipw>0), id=cloneid, weight=ipw)) # original variance (invalid) need bootstrap to get correct SEs
  (ttr = coxph(Surv(agein, age, event)~limit, data=filter(cens_data_simple, ipw>0), id=cloneid, weight=ipw, cluster=id)) # Robust variance (valid in large samples)
  list(nc=nc, msm=tt, msmr=ttr, dx=wtdx)
}


# traditional estimates:
(std_hr1 <- coxph(Surv(agein, age, d1)~cxl + gender + race + wagestatus, data=sim_cohort))
(std_hr2 <- coxph(Surv(agein, age, d2)~cxl + gender + race + wagestatus, data=sim_cohort))
# this model is not needed: just showing that the multi-state model is equivalent to two Cox models
sim_cohort$event <- factor(sim_cohort$d2 + sim_cohort$d1*2, 0:2, labels=c("censor", "d_other", "d_lc"))
(std_hr <- coxph(Surv(agein, age, event)~cxl + gender + race + wagestatus, data=sim_cohort, id=id))
# if exposure is harmful, then coefficient should be positive (can be biased from HWSB)

# IPW estimates, addressing bias from HWSB
mem.maxVSize(vsize = 32768)
point_estimates = doipw(sim_cohort, limits)
# check diagnostics
point_estimates$dx

# ln-HR: effect of one unit change in the exposure limit (Marginal structural policy model)
(ttr = point_estimates$msmr)


# survival curve at specified limits based on predicting from MSM
# NOTE that the linear ln-HR may not be appropriate here and a model with a spline
# on the "limit" variable may be warranted
cc_nc = survfit(point_estimates$nc, newdata = data.frame(intercept=1))  # can alternatively get predictions for a very high limit that approximates the natural course
ccr = survfit(point_estimates$msmr, newdata = data.frame(limit=3.0))

commontimes = intersect(cc_nc$time, ccr$time)
ccncindex = which(cc_nc$time %in% commontimes)
ccindex = which(ccr$time %in% commontimes)

# risk difference, natural course vs. limit of 0.1
rddf = data.frame(
  age=ccr$time[ccncindex], 
  risklc_nc=cc_nc$pstate[,,3][ccindex], 
  risklc_int = ccr$pstate[,,3][ccindex]
  )
rddf$rdlc = with(rddf, risklc_int-risklc_nc)


# risk difference at age 90 (older ages may be of interest or too sparse)
tail(filter(rddf, age<=90))


# risk curves for lung cancer (intervention in red)
plot(rddf$age, rddf$risklc_nc, type="s", ylab="risk", xlab="age", xlim=c(16,90))
lines(rddf$age, rddf$risklc_int, col=2, type="s")

# plotting risk difference over time
plot(rddf$age, rddf$rdlc, ylab="risk difference", xlab="age", type="s", xlim=c(16,90))
