source("sim/ipw_edcorner_0_datagen_fns.R")
# data generation
data(us_119ucod_19602021, package="LTASR")
rateset = getallcauserates(us_119ucod_19602021$rates)


word2number <- function(word){
  wordvec = tolower(do.call(c,strsplit(gsub(" |_|-", "", word),"")))
  #as.numeric(paste0(match(wordvec, letters), collapse=""))
  sum(match(wordvec, letters))
}

set.seed(word2number("bertke keil kelly-reif"))

# the data
sim_cohort = cohort_gen(n=10000,            # number of participants
                         exposure_limit=10, # set the annual exposure limit in the observed data (can be Inf, too) - useful for establishing true risk estimates if needed
                                            # preferably set to a large finite value, rather than Inf (generally there will be some operational limits)
                         cause=16,          # Cause of death of interest (following LTAS numbering of minor causes)
                         # these should make sense
                         age_firsthire = 16, 
                         age_lasthire = 40, 
                         year_firsthire = 1950, 
                         year_lasthire = 1980, 
                         year_eof=2030, 
                         rateset=rateset
                         )

# save(sim_cohort, file="data/sim_cohort.RData")
# save(sim_cohort, file="ipcw.limits/data/sim_cohort.RData")

# basic data dictionary:
# 'data.frame':	5243 obs. of  18 variables:          Each record comprises 1 person-year
# $ id              : int  1 1 1 1 1 1 1 1 1 1 ...   Unique participant identifier
# $ gender          : chr  "F" "F" "F" "F" ...       Gender
# $ race            : chr  "W" "W" "W" "W" ...       Race (W or N)
# $ wagestatus      : int  0 0 0 0 0 0 0 0 0 0 ...   Wage status (1=wage, 0 = salary)
# $ agecat          : chr  "[30,35)" "[30,35)" "[30,35)" "[30,35)" ...      Observation specific age category
# $ yearcat         : chr  "[1960,1965)" "[1960,1965)" "[1960,1965)" "[1960,1965)" ...  Observation specific year category
# $ age             : num  31 32 33 34 35 36 37 38 39 40 ...                         Observation specific age    
# $ year            : int  1957 1958 1959 1960 1961 1962 1963 1964 1965 1966 ...   Observation specific year
# $ minor           : num  16 16 16 16 16 16 16 16 16 16 ...                       Cause of death of interest (following LTAS numbering of minor causes)
# $ rate            : num  9.09e-06 9.11e-06 9.14e-06 9.06e-06 2.70e-05 ...        Annual rate of cause of death of interest for participant (Based on LTAS rates, but updated with exposure/confounder)
# $ rate_allcause   : num  0.00114 0.00114 0.00115 0.00114 0.00181 ...            Annual rate of cause of death other than cause of interest  for participant (Based on LTAS rates, but updated with exposure/confounder)
# $ baseline_frailty: num  0.171 0.171 0.171 0.171 0.171 ...                      Baseline frailty variable: predicts employment status and death
# $ rate_term       : num  0.0218 0.0217 0.0219 0.0222 0.0218 ...                 Annual rate of leaving work (not used in analysis)
# $ atwork          : num  1 1 1 0 0 0 0 0 0 0 ...                                Participant worked during the person year (1=year, 0=no)
# $ d1              : num  0 0 0 0 0 0 0 0 0 0 ...                                Indicator: death from cause of interest
# $ d2              : num  0 0 0 0 0 0 0 0 0 0 ...                                Indicator: death from any other cause
# $ x               : num  0.232 0.31 0.333 0 0 ...                               Exposure
# $ cumx            : num  0.232 0.542 0.875 0 0 ...                              Cumulative exposure
