#1) clone the data for an intervention of an annual limit at 2
# create a "clone id" that is distinct to each clone
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
  ungroup() %>%
  filter(
    drop < 2
  )

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
