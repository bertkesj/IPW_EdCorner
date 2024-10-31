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

# checking exposure at work
quantile(filter(sim_cohort, atwork==1)$x, 0:20/20)

#1) set limit (example uses limit > max exposure, which should return a "natural course" estimate)

limit = 69

#2) artificially censor data

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
  confdmod <- glm(cens ~ rcspline.eval(age, knots=agekn0) + rcspline.eval(year, knots=yearkn0) + wagestatus + gender + race, data = cens_data, weight=conf_weight, family=binomial())
  confnmod <- glm(cens ~ 1, data = cens_data, weight=conf_weight, family=binomial())
  cens_data = cens_data %>% 
    mutate(
      dconf = as.numeric(predict(confdmod, type="response")),
      nconf = as.numeric(predict(confnmod, type="response")),
    )
} else{
  cens_data = cens_data %>% 
    mutate(
      dconf = 0,
      nconf = 0,
    )
}

if(sum(cens_data[cens_data$fu_weight==1,"cens"]) > 10){
  censdmod <- glm(cens ~ rcspline.eval(age, knots=agekn) + rcspline.eval(year, knots=yearkn) + wagestatus + gender + race + cxl, data = cens_data, weight=fu_weight, family=binomial())
  censnmod <- glm(cens ~ 1, data = cens_data, weight=fu_weight, family=binomial())
  cens_data = cens_data %>% 
    mutate(
      dcens = as.numeric(predict(censdmod, type="response")),
      ncens = as.numeric(predict(censnmod, type="response"))
    )
} else{
  cens_data = cens_data %>% 
    mutate(
      dcens = 0,
      ncens = 0,
    )
}




cens_data = cens_data %>% 
  mutate(
    wtcontr = case_when(
      ((fobs____ == 1) & (atwork==1)) ~ (1-nconf)/(1-dconf),
      ((fobs____ == 0) & (atwork==1)) ~ (1-ncens)/(1-dcens),
      TRUE ~ 1
    )
  ) %>%
  group_by(id) %>%
  mutate(
    ipw = cumprod(wtcontr)
  ) %>%
  group_by()
  
mean(filter(cens_data, cens==0)$ipw)
sum(filter(cens_data, cens==0)$d1)
sum(filter(cens_data, cens==0)$d2)

#
cens_data$event <- factor(cens_data$d2 + cens_data$d1*2, 0:2, labels=c("censor", "d_other", "d_lc"))

#tt = survfit(Surv(ageout_new, d)~1, data=dat)
#tail(cbind(tt$time, tt$surv))
#sum(dat$ageout_new >= tt$time[62])

tt = survfit(Surv(timein, time, event)~1, data=cens_data, id=id, weight=ipw)

riskdf = data.frame(time = tt$time, cumhaz = as.matrix(tt$cumhaz)) %>%
  mutate(
    lim = limit,
    haz1 = c(0,diff(cumhaz.1.2)),
    haz2 = c(0,diff(cumhaz.1.3)),
    s = exp(-(cumhaz.1.2 + cumhaz.1.3)),
    sprev = lag(s, default=1),
    risk1 = cumsum(sprev*haz1),
    risk2 = cumsum(sprev*haz2),
  )
tail(riskdf)

# observed
# note the calculated risk will not equal the observed risk due to differential total possible years of follow-up by hire date
sum(sim_cohort$d1)/length(unique(sim_cohort$id))
sum(sim_cohort$d2)/length(unique(sim_cohort$id))

# now limited to first 50 years of follow-up (min total possible time)
tail(filter(riskdf, time<=50))

# note the calculated risk should approximately equal the observed risk with truncation at 50 years of potential follow-up
sim_cohort50 = filter(sim_cohort, time<=50)
sum(sim_cohort50$d1)/length(unique(sim_cohort50$id))
sum(sim_cohort50$d2)/length(unique(sim_cohort50$id))
