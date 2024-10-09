##################################################################
# weight_aggregator.R
#  Internal functions for miscellaneous tasks
#  Functions: 
##################################################################

safelog <- function(x,eps=.Machine$double.eps){
  ifelse(x<eps,log(eps),log(x))
}
safeexp <- function(x,absmax=min(log(.Machine$double.xmax), abs(log(.Machine$double.xmin)))){
  ifelse(abs(x)<absmax,exp(x),exp(sign(x)*absmax))
}
safe_expit <- function(mu)1/(1+safeexp(-mu))
safe_logit <- function(p) safelog(p) - safelog(1-p)