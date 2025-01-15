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

# checking exposure at work
quantile(filter(sim_cohort, atwork==1)$x, 0:20/20)

#1) clone the data for an intervention of an annual limit at 2
# create a "clone id" that is distinct to each clone
clones <- sim_cohort %>%
  mutate(
    limit = as.numeric(2.0),
    cloneid = paste0("clone2.0_", id)
  )

#2) artificially censor data (keeping first observation that is censored)

cens_data <- clones %>%
  group_by(cloneid) %>%
  mutate(
    cens = pmin(1,cumsum(x > limit)),
    drop = pmin(2, cumsum(cens))
  ) %>%
  ungroup() %>%
  filter(
    drop < 2
  )

# 2b) Checking for potential positivity issues by age
plot(density(filter(cens_data, cens==0 & atwork==1)$age), col="black", main="Age of active work", xlab="Age", ylab="Kernel Density", lwd=2)
lines(density(filter(sim_cohort, atwork==1)$age), col="red", lwd=2)
legend("topright", legend=c("Uncensored (limit=2)", "Cohort"), col=c("black", "red"), lty=1, lwd=2)

# also check by year
plot(density(filter(cens_data, cens==0 & atwork==1)$year), col="black", main="Year of active work", xlab="Year", ylab="Kernel Density", lwd=2)
lines(density(filter(sim_cohort, atwork==1)$year), col="red", lwd=2)
legend("topright", legend=c("Uncensored (limit=2)", "Cohort"), col=c("black", "red"), lty=1, lwd=2)


# 3) create 1/0 weights to use for confounding/censoring during follow-up (alternative is to create new datasets with some observations dropped)
#  note: these are not the inverse probability of censoring weights that will be used later
cens_data <- group_by(cens_data, cloneid) %>%
  mutate(
    one = 1,
    nobs = cumsum(one),
    fobs____ = as.numeric(nobs == 1)
  ) %>%
  mutate(
    conf_weight = as.numeric(nobs  == 1), # all of these are at work
    fu_weight = as.numeric((nobs  > 1) & (atwork == 1)),    
    dconf = 0,
    dcens = 0,
    intervened = 0,
  ) %>%
  select(-c(one, nobs)) %>%
  ungroup()


#4) fit censoring models (can include pre-baseline exposure in practice)
# use flexible functions of age and time based on restricted cubic splines
agekn0 = attr(rcspline.eval(filter(cens_data, time==1 & atwork == 1)$age, nk = 4), "knots")
yearkn0 = attr(rcspline.eval(filter(cens_data, time==1 & atwork == 1)$year, nk = 4), "knots")
agekn = attr(rcspline.eval(filter(cens_data, time>1 & atwork == 1)$age, nk = 4), "knots")
yearkn = attr(rcspline.eval(filter(cens_data, time>1 & atwork == 1)$year, nk = 4), "knots")
cxkn = attr(rcspline.eval(filter(cens_data, time>1 & atwork == 1)$cxl, nk = 4), "knots")
clkn = attr(rcspline.eval(filter(cens_data, time>1 & atwork == 1)$cumatworkl, nk = 4), "knots")
cens_data$mxl = cens_data$cxl / (cens_data$cumatworkl + as.numeric(cens_data$time==1))

# only fit models if there is more than a small amount of censoring (seems to work either way, but this avoids convergence problems)
# 4a) "censored at baseline" model 
if(sum(cens_data[cens_data$conf_weight==1,"cens"]) > 10){
  confdmod <- glm(cens~ 
                    age + rcspline.eval(age, knots=agekn0) + 
                    wagestatus + male + race, 
                  data = cens_data, weight=conf_weight, family=binomial())
  
  # numerator for stabilizing weight
  confnmod <- glm(cens~ age + rcspline.eval(age, knots=agekn0) , data = cens_data, weight=conf_weight, family=binomial())
} 

# 4b) "censored during follow-up" model 
if(sum(cens_data[cens_data$fu_weight==1,"cens"]) > 10){
  censdmod <- glm(cens~ mxl + cumatworkl + 
                    age + rcspline.eval(age, knots=agekn) + 
                    wagestatus + male + race , data = cens_data, weight=fu_weight, 
                  family=binomial())
  # numerator for stabilizing weight
  censnmod <- glm(cens ~ age + rcspline.eval(age, knots=agekn), data = cens_data, weight=fu_weight, family=binomial())
  } 

# 5) create weights
# 5a) "censored at baseline" model 
if(sum(cens_data[cens_data$conf_weight==1,"cens"]) > 10){
  # weight components
  cens_data$nconf = cens_data$conf_weight*as.numeric(predict(confnmod, type="response"))
  cens_data$dconf = cens_data$conf_weight*as.numeric(predict(confdmod, type="response"))
} else{
  cens_data = cens_data %>% 
    mutate(
      dconf = 0,
      nconf = 0,
    )
}

# 5b) "censored during follow-up" model 
if(sum(cens_data[cens_data$fu_weight==1,"cens"]) > 10){
  # weight components
  cens_data$ncens = cens_data$fu_weight*as.numeric(predict(censnmod, type="response"))
  cens_data$dcens = cens_data$fu_weight*as.numeric(predict(censdmod, type="response")) 
} else{
  cens_data = cens_data %>% 
    mutate(
      dcens = 0,
      ncens = 0,
    )
}
# 5c) calculate weights combining predictions from both censoring models and accrue over time
cens_data <- cens_data %>% 
  mutate(
    wtcontr = case_when(
      ((fobs____ == 1) & (atwork==1)) ~ (1-cens)*(1-nconf)/(1-dconf),
      ((fobs____ == 0) & (atwork==1)) ~ (1-cens)*(1-ncens)/(1-dcens),
      .default=1
    )
  ) %>%
  group_by(cloneid) %>%
  mutate(
    ipw = cumprod(wtcontr),
    ipwt = pmin(10, cumprod(wtcontr)) # sometimes truncated weights are advocated
  ) %>%
  ungroup()

  
#6) estimate risk under an exposure limit

#  6a) Use Aalen-Johansen estimator to get risk for each outcome at limit
cens_data$event <- factor(cens_data$d2 + cens_data$d1*2, 0:2, labels=c("censor", "d_other", "d_lc"))
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

