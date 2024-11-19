# example: risk under a single hypothetical intervention, and getting hazard ratio compared with no intervention

library(dplyr)
library(Hmisc)
library(survival)
load("data/sim_cohort.RData")



# data process
sim_cohort <- sim_cohort %>%
  group_by(id) %>%
  mutate(
    one = 1,
    time = cumsum(1),
    timein = time-1,
    agein = age-1,
    xl = lag(x, default=0),
    cxl = lag(cumx, default=0)
  ) %>%
  select(-one) %>%
  group_by()

# checking exposure at work
quantile(filter(sim_cohort, atwork==1)$x, 0:20/20)

#1) set limit (can use limit > max exposure, which should return a "natural course" estimate)

limit = 3

#2) artificially censor data (keeping first observation that is censored)

cens_data <- sim_cohort %>%
  group_by(id) %>%
  mutate(
    cs = cumsum(x > limit),
    cens = pmin(1,cumsum(x > limit)),
    drop = pmin(2, cumsum(cens))
  ) %>%
  group_by() %>%
  filter(
    drop < 2
  )


# 3) create 1/0 weights to use for confounding/censoring during follow-up (alternative is to create new datasets)
cens_data = cens_data %>%
  group_by(id) %>%
  mutate(
    one = 1,
    nobs = cumsum(one),
    fobs____ = as.numeric(nobs==1)
  ) %>%
  mutate(
   conf_weight = as.numeric(nobs  == 1),
   fu_weight = as.numeric((nobs  > 1) & (atwork == 1))
  ) %>%
  select(
    -c(one, nobs)
  ) %>%
  group_by()

length(unique(cens_data$id))
sum(filter(cens_data, conf_weight == 1)$cens)
sum(filter(cens_data, fu_weight == 1)$cens)

sum(filter(cens_data, drop == 0)$d1)
sum(filter(cens_data, drop == 0)$d2)

  
#3) fit weight models (can include pre-baseline exposure in practice)
agekn0 = attr(rcspline.eval(filter(cens_data, conf_weight==1)$age, nk = 4), "knots")
yearkn0 = attr(rcspline.eval(filter(cens_data, conf_weight==1)$year, nk = 4), "knots")
agekn = attr(rcspline.eval(filter(cens_data, fu_weight==1)$age, nk = 4), "knots")
yearkn = attr(rcspline.eval(filter(cens_data, fu_weight==1)$year, nk = 4), "knots")

# fit models if there is more than a small amount of censoring (seems to work either way, but this avoids convergence problems)
if(sum(cens_data[cens_data$conf_weight==1,"cens"]) > 10){
#  confdmod <- glm(cens ~ rcspline.eval(age, knots=agekn0) + rcspline.eval(year, knots=yearkn0) + wagestatus + gender + race, data = cens_data, weight=conf_weight, family=binomial())
  confdmod <- glm(cens ~ rcspline.eval(year, knots=yearkn0) + wagestatus, data = cens_data, weight=conf_weight, family=binomial())
  #confnmod <- glm(cens ~ rcspline.eval(age, knots=agekn), data = cens_data, weight=conf_weight, family=binomial())
  cens_data = cens_data %>% 
    mutate(
      dconf = as.numeric(predict(confdmod, type="response")),
    )
} else{
  cens_data = cens_data %>% 
    mutate(
      dconf = 0,
    )
}

if(sum(cens_data[cens_data$fu_weight==1,"cens"]) > 10){
  #censdmod <- glm(cens ~ rcspline.eval(age, knots=agekn) + rcspline.eval(year, knots=yearkn) + wagestatus + gender + race + cxl, data = cens_data, weight=fu_weight, family=binomial())
  censdmod <- glm(cens ~ year + wagestatus, data = cens_data, weight=fu_weight, family=binomial())
  #censnmod <- glm(cens ~ rcspline.eval(age, knots=agekn), data = cens_data, weight=fu_weight, family=binomial())
  cens_data = cens_data %>% 
    mutate(
      dcens = as.numeric(predict(censdmod, type="response")),
    )
} else{
  cens_data = cens_data %>% 
    mutate(
      dcens = 0,
    )
}

# calculating weights
cens_data = cens_data %>% 
  mutate(
    wtcontr = case_when(
      # note using unstabilized weights: Cain et al 2010 show that weight stabilization doesn't work in IPCW as well as it does for static regime MSMs
      # in particular if marginal censoring is used in the numerator, it becomes a function of the time-varying exposure (which it should not be) because it is
      # restricted to those who are consistent with the exposure regime in the past
      ((fobs____ == 1) & (atwork==1)) ~ (1-cens)*(1)/(1-dconf),
      ((fobs____ == 0) & (atwork==1)) ~ (1-cens)*(1)/(1-dcens),
      .default = 1
    )
  ) %>%
  group_by(id) %>%
  mutate(
    ipw = cumprod(wtcontr)
  ) %>%
  group_by()
  
summary(filter(cens_data)$ipw)
mean(filter(cens_data, cens==0)$ipw)
sum(filter(cens_data, cens==0)$ipw)
nrow(sim_cohort)
sum(filter(cens_data, cens==0)$ipw)/nrow(sim_cohort)
nrow(cens_data)
nrow(filter(cens_data, cens==0))


sum(filter(cens_data, cens==0)$d1)
sum(filter(cens_data, cens==0)$d2)
length(unique(cens_data$id))

# Use Aalen-Johansen estimator to get risk for each outcome
cens_data$event <- factor(cens_data$d2 + cens_data$d1*2, 0:2, labels=c("censor", "d_other", "d_lc"))
tt = survfit(Surv(agein, age, event)~1, data=cens_data, id=id, weight=ipw)

riskdf = data.frame(tt$time, tt$pstate)
names(riskdf) = c("age", "surv", "risk_other", "risk_lc")



# note the calculated risk will not equal the observed risk due to late entry
# this could be ameliorated by using time-on-study as the time-scale
sum(sim_cohort$d1)/length(unique(sim_cohort$id))
sum(sim_cohort$d2)/length(unique(sim_cohort$id))


# now limited to  age 90
tail(filter(riskdf, age<=90))


# getting a Hazard ratio to estimate impact of hypothetical intervention (vs. doing nothing)
ncdf = sim_cohort
ncdf$ipw = 1
ncdf$intervention = 0
ncdf$cloneid = paste0("nc_",ncdf$id)
intdf = cens_data
intdf$intervention = 1
intdf$cloneid = paste0("int_",intdf$id)

compdat = merge(ncdf,intdf, on="cloneid", all=TRUE) %>%
  filter(ipw > 0) # concatenate data (not really a merge)
compdat$event <- factor(compdat$d2 + compdat$d1*2, 0:2, labels=c("censor", "d_other", "d_lc"))
# use Robust variance to get confidence intervals (bootstrap would be preferred)
# this multi-state model is equivalent to fitting two Cox models
ft <- coxph(Surv(agein, age, event)~intervention, data=compdat, id=cloneid, weight=ipw, cluster=id)
## interpretation: intervention to reduce exposure should reduce mortality, so coefficient(s) should be negative if exposure is harmful

# Use Aalen-Johansen-like estimator to get risk for each outcome from multi-state model
cens_data$event <- factor(cens_data$d2 + cens_data$d1*2, 0:2, labels=c("censor", "d_other", "d_lc"))
cc = survfit(ft, newdata = data.frame(intervention=c(0,1)))
rddf = data.frame(tt$time, do.call(cbind, lapply(1:3, function(x) tt$pstate[,,x])))
names(rddf) <- c("age", "surv_nc", "surv_int", "riskother_nc", "riskother_int", "risklc_nc", "risklc_int")
rddf$rdlc = with(rddf, risklc_int-risklc_nc)

# risk difference at age 90 (older ages may be of interest or too sparse)
tail(filter(rddf, age<=90))

# risk curves for lung cancer (intervention in red)
plot(rddf$age, rddf$risklc_nc, type="s", ylab="risk", xlab="age", xlim=c(16,90))
lines(rddf$age, rddf$risklc_int, col=2, type="s")

# plotting risk difference over time
plot(rddf$age, rddf$rdlc, ylab="risk difference", xlab="age", type="s", xlim=c(16,90))



