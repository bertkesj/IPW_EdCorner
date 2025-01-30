
#1) set grid of limits based on regular values of empirical distribution
# the lowest limit must include some uncensored deaths
# find the max exposure over a career: 
#  we can only explore interventions in which some individuals are observed to follow limits
#  so here we find a regular grid of values for the maximum annual exposure
maxexposures = as.numeric(tapply(sim_cohort$x, sim_cohort$id, max))

# Finding 30 potential occupational limits for a marginal structural model
# strategy: pick limits that are well within the bounds of the data
(limits = quantile(maxexposures[maxexposures>0], seq(.20,0.95,length.out=30))) # 30 limits all within the range of observed lifetime max exposures
# want to have at least 20% remain uncensored,

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

tot_events

# feasible limits
limits = filter(tot_events, tot_d1 > 10 & tot_d2 > 10)$limit
