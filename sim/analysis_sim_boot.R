
library(survival)

ids <- unique(full_cohort$id)


boot <- map(1:10,
    ~{
      samp_id <- sample(ids,
                        replace = T)
      fc_del <- tibble(id = samp_id) %>%
        group_by(id) %>%
        mutate(rep = row_number()) %>%
        arrange(id, rep) %>%
        left_join(full_cohort,
                  by = 'id',
                  relationship = 'many-to-many') %>%
        mutate(id = paste0(id, rep))
    
      dat <- map(c(5:10 / 10, 2:10, 1000),
                 ~cohort(., fc_del)) %>%
        reduce(bind_rows)  %>%
        mutate(age_beg = age,
               age_end = age + .99) %>%
        select(id, age_beg, age_end, d1, IPW, IPWunadj, cohort)
      cox <- coxph(Surv(age_beg, age_end, d1) ~ factor(cohort,
                                                       c(1000, 5:10 / 10, 2:10)), 
                   #id = paste(id, cohort),
                   weight=IPW, 
                   data=dat,
                   ties='breslow')
      summary(cox)$coefficients %>%
        as_tibble(rownames = 'var') %>%
        mutate(var = str_remove(var, 
                                fixed('factor(cohort, c(1000, 5:10/10, 2:10))')) %>%
                 as.numeric(),
               rep = .x) %>%
        return()
}) %>%
  reduce(bind_rows)






fc_del <- full_cohort
dat <- map(c(5:10 / 10, 2:10, 1000),
           ~cohort(., fc_del)) %>%
  reduce(bind_rows)  %>%
  mutate(age_beg = age,
         age_end = age + .99) %>%
  select(id, age_beg, age_end, d1, IPW, IPWunadj, cohort)
cox <- coxph(Surv(age_beg, age_end, d1) ~ factor(cohort,
                                                 c(1000, 5:10 / 10, 2:10)), 
             #id = paste(id, cohort),
             weight=IPW, 
             data=dat,
             ties='breslow')
orig <- summary(cox)$coefficients %>%
  as_tibble(rownames = 'var') %>%
  mutate(var = str_remove(var, 
                          fixed('factor(cohort, c(1000, 5:10/10, 2:10))')) %>%
           as.numeric()) 



boot %>%
  group_by(var) %>%
  dplyr::summarize(med = median(`exp(coef)`),
                   lower = quantile(`exp(coef)`, .025),
                   upper = quantile(`exp(coef)`, .975)) %>%
  ggplot() +
  geom_point(aes(x=var, 
                 y=med)) +
  geom_point(data=orig,
             aes(x=var, 
                 y=`exp(coef)`),
             color='red') +
  geom_errorbar(aes(x=var, 
                    ymin =lower,
                    ymax = upper)) +
  scale_x_continuous(trans='log10') +
  ylab('HR (referent: no threshold)') +
  xlab('Threshold') +
  geom_hline(yintercept=1) +
  theme_bw()

orig %>%
  mutate(lower = exp(coef - 1.96*`robust se`),
         upper = exp(coef + 1.96*`robust se`)) %>%
  ggplot() +
  geom_point(aes(x=var, 
                 y=exp(coef))) +
  geom_errorbar(aes(x=var, 
                    ymin =lower,
                    ymax = upper)) +
  scale_x_continuous(trans='log10') +
  ylab('HR (referent: no threshold)') +
  xlab('Threshold') +
  geom_hline(yintercept=1) +
  theme_bw()

