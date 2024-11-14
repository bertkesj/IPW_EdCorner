##################################################################
# weight_api.R
#  api functions for estimating weights
#  Functions: 
##################################################################


# TODO: create function that takes in a dataframe, fits weight models, and adds the weight to the dataframe
# TODO: make user visible with documentation
# assumption: data sorted by time
estimate_confounding_weights <- function(formulas, df, id, censvar, method="glm"){
  bldat <- df %>% 
    group_by(id) %>% 
    slice_head(n = 1) %>%
    ungroup()
  
  confounding_model_glm(f, bldat) # TODO: create this function
}
dd = estimate_confounding_weights(df, "id", "cens")



# TODO: create function that takes in a dataframe, fits weight models, and adds the weight to the dataframe
# TODO: make user visible with documentation
estimate_censoring_weights <- function(formulas, df, id, censvar, method="glm"){
  
}