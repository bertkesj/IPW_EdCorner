library(tidyverse)
library(lubridate)
library(LTASR)
library(rms)

##function
invlogit <- function(x) {
  mx <- pmin(500, x)
  exp(mx) / (1 + exp(mx)) %>%
    return()
} 

full_cohort <- read_csv('sim/cohort.csv')
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
dat <- map(c(5:10 / 10, 2:10, 75),
           cohort) %>%
  reduce(bind_rows)  %>%
  mutate(age_beg = age,
         age_end = age + .99) %>%
  select(id, age_beg, age_end, d1, IPW, IPWunadj, cohort)
cox <- coxph(Surv(age_beg, age_end, d1) ~ factor(cohort,
                                                 c(75, 5:10 / 10, 2:10)), 
             weight=IPW, 
             data=dat)
summary(cox)$coefficients %>%
  as_tibble(rownames = 'var') %>%
  mutate(var = str_remove(var, 
                          fixed('factor(cohort, c(75, 5:10/10, 2:10))')) %>%
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




#######################################
# Alternative Cox Model Function
#######################################
coxphreg <- function(lin = ~ 1, loglin = ~ 1, data){
  ll <- function(b) {
    
    hloglin <- rep(0, nrow(data))
    hlin <- rep(0, nrow(data))
    if (ncol(Xloglin) != 0) hloglin <- Xloglin %*% b[1:ncol(Xloglin)]
    if (ncol(Xlin) != 0) hlin <-  Xlin %*% b[(ncol(Xloglin) + 1):tot]
    data$h <- exp(hloglin)*(1 + hlin)
    
    ll <- data %>%
      group_by(case_id) %>%
      dplyr::summarize(den = sum(IPW*h),
                       num = sum(IPW*c*h)) %>%
      mutate(L = num / den) %>%
      dplyr::summarize(x = sum(log(L)))
    
    return(ll$x)
  }
  g <- function(b) {
    if (ncol(Xloglin) != 0) hloglin <- exp(as.numeric(Xloglin %*% b[1:ncol(Xloglin)])) else hloglin <- rep(1, nrow(data))
    if (ncol(Xlin) != 0) hlin <-  (1 + as.numeric(Xlin %*% b[(ncol(Xloglin) + 1):tot])) else hlin <- rep(1, nrow(data))
    data$hloglin <- hloglin
    data$hlin <- hlin
    
    bbloglin <- map_dbl(asplit(Xloglin, 2),
                        ~ {
                          data$x <- as.numeric(.x)
                          data %>%
                            group_by(case_id) %>%
                            mutate(h = hloglin*hlin) %>%
                            dplyr::summarize(den = sum(IPW*h),
                                             num = sum(IPW*h*x),
                                             cx = sum(IPW*c*x),
                                             .groups = 'drop') %>%
                            mutate(L = cx - num/den) %>%
                            dplyr::summarize(x = sum(L)) %>%
                            `$`(x)
                        })
    bblin <- map_dbl(asplit(Xlin, 2),
                     ~ {
                       data$x <- as.numeric(.x)
                       data %>%
                         group_by(case_id) %>%
                         dplyr::summarize(den = sum(hloglin*hlin),
                                          num = sum(hloglin*x),
                                          cxn = sum(c*x),
                                          cxd = sum(c*hlin),
                                          .groups = 'drop') %>%
                         mutate(L = cxn/cxd - num/den) %>%
                         dplyr::summarize(x = sum(L)) %>%
                         `$`(x)
                     })
    return(c(bbloglin, bblin))
  }
  
  
  Xloglin <- model.matrix(loglin, data = data) %>%
    `[`(,-1, drop = F)
  Xlin <- model.matrix(lin, data = data) %>%
    `[`(,-1, drop = F)
  
  
  tot <- ncol(Xlin) + ncol(Xloglin)
  b <- rep(0, tot)
  
  o <- optim(b, ll, hessian = T, method = "BFGS", control=list(fnscale=-1))#
  #se <- sqrt(diag(solve(-o$hessian)))
  #lower <- o$par - 1.96*se
  #upper <- o$par + 1.96*se
  out <- tibble(var = c(colnames(Xloglin), colnames(Xlin)),
                est = o$par,
                #se = se,
                #lower = lower,
                #upper = upper
  )
  return(list(output=out, 
              neg2LL_diff = -2*ll(b) - - 2*o$value, 
              AIC = - 2*o$value + 2*tot,
              hessian = o$hessian))
}



RRs <- c()
cohs <- c(5:10 / 10, 2:10)
for (coh in cohs){
  print(coh)
  dat75 <- dat %>%
    filter(cohort %in% c(75, coh)) %>%
    mutate(id = paste0(id, '_', cohort))
  
  cases <- dat75 %>%
    filter(d1 == 1) %>%
    select(case_id = id,
           age = age_beg)
  
  risk_sets <- left_join(cases, 
                         dat75, 
                         by=c('age' = 'age_beg'),
                         relationship = "many-to-many") %>%
    mutate(c = case_id == id) %>%
    select(case_id, id, c, cohort, IPW) %>%
    distinct()
  rr <- coxphreg(loglin = ~ factor(cohort,
                                   c(75, coh)),
                 data=risk_sets) %>%
    `$`(output) %>%
    `$`(est)
  RRs <- c(RRs, rr)
  print(rr)
}
     