library(dplyr) # TODO: change to specific tidyverse packages (dplyr?)
library(lubridate)
library(LTASR)
library(rms)

##function
invlogit <- function(x) {
  mx <- pmin(500, x)
    1/(1 + exp(-mx)) %>%
    return()
} 

#Simulate Full Cohort
source('sim/ipw_edcorner_0_datagen.R')

#Delete person time after first death (d1c, d2c)
fc_del <- full_cohort %>%
  group_by(id) %>%
  mutate(d1c = cumsum(d1),
         d2c = cumsum(d2)) %>%
  group_by(id, d1c) %>%
  filter(d1c == 0 |
           (d1c == 1 & row_number() == 1)) %>%
  group_by(id, d2c) %>%
  filter(d2c == 0 |
           (d2c == 1 & row_number() == 1)) 

# Save csv for package
#sim_cohort <- fc_del %>% group_by()
#save(sim_cohort, file="data/sim_cohort.RData")

###Function to create cohort copy
cohort <- function(limit){
  ################################################################################
  # Censored
  per_his <- fc_del %>%
    group_by(id) %>%
    mutate(censored = 1*(x > limit) %>%
             cumsum() %>%
             `>=`(1))
  
  ################################################################################
  # Employed Model
  emp_c <- per_his %>%
    filter(atwork == 1) %>%
    group_by(id)  %>%
    filter(censored == 0 |
             (censored == 1 & (lag(censored) == 0 |
                                 row_number() == 1))) %>%
    mutate(cumx_prev = cumx - 1)
  
  emp_c$Pc <- glm((censored==1) ~ rcs(age,4) + rcs(year,4) + 
                    factor(race) + factor(gender) +
                    factor(wagestatus) + rcs(cumx_prev,4) , 
                  data=emp_c,
                  family=binomial(link='logit')) %>%
    predict(newdata=emp_c) %>%
    invlogit()
  
  emp_c$Pca <- glm((censored==1) ~ rcs(age,4), 
                   data=emp_c,
                   family=binomial(link='logit')) %>%
    predict(newdata=emp_c) %>%
    invlogit()
  
  Mn <- emp_c %>%
    select(id, year, age, Pc, Pca)
  
  ################################################################################
  # Weighted data
  per_weight <- per_his %>%
    group_by(id)  %>%
    filter(censored == 0) %>%
    
    left_join(Mn, by=c('id', 'year', 'age')) %>%
    mutate(Pca = if_else(is.na(Pca), 0, Pca),
           Pc = if_else(is.na(Pc), 0, Pc)) %>%
    group_by(id) %>%
    mutate(IPW = cumprod((1-Pca)/(1-Pc)),

           IPWunadj = cumprod(1/(1-Pc)),
           cohort = limit) %>%
    ungroup() %>%
    mutate(IPW = pmin(quantile(IPW, .99), IPW),
           IPWunadj = pmin(quantile(IPWunadj, .99), IPWunadj))
  return(per_weight)
}

library(survival)
dat <- map(c(5:10 / 10, 2:10, 1000),
           cohort) %>%
  reduce(bind_rows)  %>%
  mutate(age_beg = age,
         age_end = age + .99) %>%
  select(id, age_beg, age_end, d1, IPW, IPWunadj, cohort)
cox <- coxph(Surv(age_beg, age_end, d1) ~ factor(cohort,
                                                 c(1000, 5:10 / 10, 2:10)), 
             weight=IPW, 
             data=dat)
summary(cox)$coefficients %>%
  as_tibble(rownames = 'var') %>%
  mutate(var = str_remove(var, 
                          fixed('factor(cohort, c(1000, 5:10/10, 2:10))')) %>%
           as.numeric()) %>%
  #View()
  ggplot() +
  geom_point(aes(x=var, 
                 y=exp(coef))) +
  scale_x_continuous(trans='log10') +
  ylab('HR (referent: no threshold)') +
  xlab('Threshold') +
  geom_hline(yintercept=1) +
  theme_bw()


dat %>%
  group_by(cohort) %>%
  dplyr::summarize(Ncases = sum(d1),
            Nw = sum(IPW*d1),
            Pw = sum(IPW), 
            Nwunadj = sum(IPWunadj*d1),
            Pwunadj = sum(IPWunadj)) %>%
  View()


