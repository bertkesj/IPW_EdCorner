# input: person-year data
#        limit variable and limits to define censoring scenarios
#           I include Inf as a limit which outputs original data
#        one sided baseline and follow-up formulas to use when modeling censoring with a glm
#        pass through variables to include in output data
#        truncation value for final IPWs
# output: single data frame which includes all cloned data with calculated IPWs

clones <- function(data,
                   limit_var,
                   limits,
                   baseline_formula,
                   fu_formula,
                   pass_thru_vars,
                   ipw_trunc = 10){
  limit_var <- enquo(limit_var)
  #1) clone the data for an intervention of an annual limit at 2
  # create a "clone id" that is distinct to each clone
  #1b) create copy of the dataset for each limit
  clones <- do.call(rbind,
                    lapply(1:length(limits), function(x) {  
                      df = data                  
                      df$limit = as.numeric(limits[x])
                      df$cloneid = paste0(x, "_", df$id)
                      df                    
                    }
                    ))
  #2) artificially censor data (keeping first observation that is censored)
  cens_data <- clones %>%
    group_by(cloneid) %>%
    mutate(
      cens = pmin(1,cumsum(!!limit_var > limit)),
      drop = pmin(2, cumsum(cens))
    ) %>%
    group_by() %>%
    filter(
      drop < 2
    )
  
  # 3) create 1/0 weights to use for confounding/censoring during follow-up 
  #    (alternative is to create new datasets)
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
  
  # Censoring models 
  for (l in limits){
    
    ### Baseline model. Check if there is sufficient (> 10) amount of censoring
    
    #filter data by censor limits scenario (limidx)
    limidx = which(cens_data$limit == l)
    tempdat = cens_data[limidx,]
    
    if(sum(tempdat[tempdat$conf_weight==1,"cens"]) > 10){
      confnmod <- glm(cens ~ age + rcspline.eval(age, nk=4) , 
                      data = tempdat, 
                      weight=conf_weight, 
                      family=binomial())
      confdmod <- glm(paste('cens', 
                            paste(as.character(baseline_formula),
                                  collapse = ' ')) %>%
                        formula(), 
                      data = tempdat, 
                      weight=conf_weight, 
                      family=binomial())
      cens_data[limidx,"nconf"] = tempdat$conf_weight*as.numeric(predict(confnmod, type="response"))
      cens_data[limidx,"dconf"] = tempdat$conf_weight*as.numeric(predict(confdmod, type="response"))
    } else{
      cens_data[limidx,"nconf"] = 0
      cens_data[limidx,"dconf"] = 0
    }
    
    ## Follow-up model
    if(sum(tempdat[tempdat$fu_weight==1,"cens"]) > 1){
      censnmod <- glm(cens ~ age + rcspline.eval(age, nk=4), 
                      data = tempdat, 
                      weight=fu_weight, 
                      family=binomial())
      censdmod <- glm(paste('cens', 
                            paste(as.character(fu_formula),
                                  collapse = ' ')) %>%
                        formula(), 
                      data = tempdat, 
                      weight=fu_weight, 
                      family=binomial())
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
        # note regarding unstabilized weights: Cain et al 2010 show that weight 
        # stabilization is not guaranteed to reduce variance
        ((fobs____ == 1) & (atwork==1)) ~ (1-cens)*(1-nconf)/(1-dconf),
        ((fobs____ == 0) & (atwork==1)) ~ (1-cens)*(1-ncens)/(1-dcens),
        .default=1
      )
    ) %>%
    select(id,cloneid,agein,age,
           !!limit_var,!!!pass_thru_vars,
           wtcontr,limit, fobs____) %>% # optional step here to reduce comp. burden
    group_by(cloneid) %>%
    mutate(
      ipw = cumprod(wtcontr)
    ) %>%
    group_by() %>%
    mutate(
      ipw_trunc = pmin(ipw, ipw_trunc)
    )
  
  return(combined_wtd_data)
}



