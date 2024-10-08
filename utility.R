##################################################################
# weight_aggregator.R
#  Internal functions for miscellaneous tasks
#  Functions: 
##################################################################

safelog <- function(x,eps=1e-10){
  ifelse(x<eps,eps,x)
}
expit <- function(mu) 1/(1+exp(-mu))
logit <- function(p) safelog(p) - safelog(1-p)