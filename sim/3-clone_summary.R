
# check mean weights by intervention
# note mean is taken across all possible observations, even those with weights = 0
N = nrow(combined_wtd_data)
Nid = length(unique(combined_wtd_data$id))
wtdx = data.frame(
  limit = limits,
  mean_cumx = do.call(c,lapply(limits, function(x) mean(combined_wtd_data[combined_wtd_data$limit == x  & combined_wtd_data$cens == 0,]$cumx))),
  n_conf = do.call(c,lapply(limits, function(x) sum(combined_wtd_data[combined_wtd_data$limit == x & combined_wtd_data$fobs____==1,]$cens))),
  n_cens = do.call(c,lapply(limits, function(x) sum(combined_wtd_data[combined_wtd_data$limit == x & combined_wtd_data$fobs____==0,]$cens))),
  n_d1w = do.call(c,lapply(limits, function(x) sum(combined_wtd_data[combined_wtd_data$limit == x & combined_wtd_data$cens == 0,]$ipw*combined_wtd_data[combined_wtd_data$limit == x & combined_wtd_data$cens == 0,]$d1))),
  n_d2w = do.call(c,lapply(limits, function(x) sum(combined_wtd_data[combined_wtd_data$limit == x & combined_wtd_data$cens == 0,]$ipw*combined_wtd_data[combined_wtd_data$limit == x & combined_wtd_data$cens == 0,]$d2))),
  mean_ipw = do.call(c,lapply(limits, function(x) mean(combined_wtd_data[combined_wtd_data$limit == x & combined_wtd_data$cens == 0,]$ipw)))
)
