##################################################################
# api.R
#  User facing functions for inverse probability weighting and
#  analysis of weighted models
#  Functions: 
#   cox_msm_ipcw
#   cumrisk_ipcw
##################################################################


cox_msm_ipcw = function(
 formula,
 ipw_conf_formulas,
 ipw_cens_formulas,
 data,
 exposure="",
 limits=NULL,
 id=NULL,
 bootstrap=FALSE
  ...
){
  #' @title Fit a Cox marginal structural model for a dynamic treatment regime
  #' @description 
  #' @details 
  #' @param formula an r formula representing the conditional model for the outcome (survival syntax)
  #' @param ipw_conf_formulas a list of r glm-style formulas representing the conditional model for censoring indicator: list(denominator_formula, numerator_formula)
  #' @param ipw_cens_formulas a list of r glm-style formulas representing the conditional model for censoring indicator: list(denominator_formula, numerator_formula)
  #' @param data a data frame
  #' @param interventions
  #' @param id NULL, or a variable name indexing individual units
  #' @param time NULL, or a variable name indexing time. If NULL, then data are assumed to be sorted by id, time (ascending)
  #' @param bootstrap 
  #' @param ... arguments to \code{\link[survival]{coxph}}
  #' @seealso \code{\link[survival]{coxph}}
  #' @concept occupational interventions
  #' @import survival stats
  #' @export
  #' @examples
  #' set.seed(50)
  #' # more here
  
  
  #0. bookkeeping (background work to ensure scoping works) ----

    origcall <- thecall <- match.call(expand.dots = FALSE)
    names(thecall) <- gsub("f", "formula", names(thecall))
    m <- match(c("formula", "data", "weights", "offset"), names(thecall), 0L)
    thecall <- thecall[c(1L, m)]
    thecall$drop.unused.levels <- TRUE
    
    # TODO: cut this?
    thecall[[1L]] <- quote(stats::model.frame)
    thecalle <- eval(thecall, parent.frame()) # a model frame pulled in from the environment in which the function was called

  #1. error catching
  if(is.null(id)){
    # note: we could allow 1 obs per id and choose to ignore this, which
    # would turn this into a special case of standard IPW
    stop("id must be specified")
  }
  #1. data handling
    if(!is.null(time)){
      # TODO: check whether sorting in place would affect data in global scope
      # check sorting by time, id
      datasorted = sort_by_time()  # TODO: create this function sort dataset by id, time
    } else {
      datasorted = data
    }

  #3. create copies of data with censoring indicators ----
    datalist = copy_and_censor(data, id, limits) # TODO: create this function, create a dataset with only the first observations
    # each dataset in the list has a new created variable INTERVENTION 

  #4. estimate weight contribution: censoring at baseline/confounding ----
    bldatalist = grab_first_obs(datalist, id, time) # TODO: create this function, create a dataset with only the first observations
    fudatalist = grab_fu_obs(datalist, id, time) # TODO: create this function, create a dataset with only the follow-up observations
  
  
  #5. estimate weights: censoring during follow-up ----
  lapply(bldatalist, function(x){
    x$deltaIPConfW = estimate_confounding_weights(x) # TODO: create function that takes in a dataframe, fits weight models, and adds the weight to the dataframe
  })
  lapply(fudatalist, function(x){
    x$deltaIPCensW = estimate_censoring_weights(x) # TODO: create function that takes in a dataframe, fits weight models, and adds the weight to the dataframe
  })
  # aggregate weights
  analysisdatalist = NULL # TODO: check radon example for how this was done

  #6. fit weighted Cox model ----
  fulldata = do.call(rbind, analysisdatalist)
  survival::coxph(formula, data=fulldata, weights=IPCW, id=id, robust=!bootstrap, ...)

}