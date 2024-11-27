
cohort_gen <- function(n=100, exposure_limit=1, age_firsthire = 16, age_lasthire = 40, year_firsthire = 1950, year_lasthire = 1980, year_eof=2030, rateset=NULL,cause=16){
  if(is.null(rateset)) stop("Please specify rateset")
  minorrateset = getonecauserates(rateset,cause)
  
  
  dat = data.frame(id = 1:n)
  dat = within(dat, {
    age_hire = round(runif(n, age_firsthire, age_lasthire))
    year_hire = round(runif(n, year_firsthire, year_lasthire))
    termrate = runif(n, .01, .1)  # randomly chosen baseline exponential termination rate (will account for time varying factors below)
    year_eof = rep(year_eof, n)
    wage = rbinom(n, 1, 0.7) # wage or salary worker 1 = wage, 0 = salary
    gender = sample(c("M", "F"), n, replace=TRUE)
    race = sample(c("N", "W", "W", "W"), n, replace=TRUE)
    #baseline_frailty = runif(n, 0.1, 0.9)
    baseline_frailty = runif(n, -1.0, 1.0)
  })
  
  # full worker data up to the point of end of follow up, excluding employment, exposure, mortality data
  workers_full =lapply(1:n, function(x){
    tdat = worker_gen0row(dat[x,], minorrateset)
    tdat$rate_term = dat[x,]$termrate
    tdat
  }
  )
  # employ, expose, and permit mortality to create observed person-time
  workers_alive = do.call(rbind, lapply(workers_full, worker_genx, limit=exposure_limit))
  
  # partially re-order some columns to make it look more typical
  workers_alive = workers_alive[,unique(c("id", "gender", "race", "wagestatus", names(workers_alive)))]
  # drop some names
  dropnames = c("minor", "rate", "rate_allcause", "baseline_frailty", "rate_term")
  currnames = unique(c("id", "x", "cumx", names(workers_alive), "d1", "d2"))
  keepnames = setdiff(currnames, dropnames)
  workers_alive[,keepnames,drop=FALSE]
}


sim_exposure <- function(mn,sd,limit, maxit=30){
  # mn = 1.0
  # sd = 1.0
  #limit = 1.0
  if(limit==0) return(0)
  (x = rnorm(1, mn, sd))
  it = 0
  while(x > log(limit) && it <= maxit){
    x = rnorm(1, mn, sd)
    if(x > log(limit) && it == maxit){
      x = log(limit)
    }
  }
  exp(x)
}

logit <- function(x){
  log(x/(1-x))
}

worker_genx <- function(worker, limit=Inf){
  # worker = workers_full[[2]]
  # worker$rate_term = worker$rate_term[1]
  # worker$rate = worker$rate[1]
  # worker$rate_allcause = worker$rate[1]
  nobs = nrow(worker)
  iswaged = worker$wagestatus[1]
  U = worker$baseline_frailty[1]
  worker$atwork = c(1,rep(0,nobs-1))
  worker$d1 = 0
  worker$d2 = 0
  worker$x = 0
  worker$cumx = 0
  i = 0
  leftwork = 0
  alive =1 
  cumx = 0;
  xl = 0;
  while((i < nobs) && alive){
    i = i+1
    # employment
    if(leftwork == 0 && i > 1) {
      worker$rate_term[i] = plogis(logit(worker$rate_term[i]) + 0.01*worker$cumx[i-1] + 0.1*(worker$age[i]-50)/10 + 2*U+ escalate(worker$age[i], 65, 0.5)) 
      worker$atwork[i] = rbinom(1, 1,  1-worker$rate_term[i])
      leftwork = 1 - worker$atwork[i]
    }
    # exposure
    if(leftwork == 0){
      # testx = runif(10000)*10
      # mu = 0; sd = 0.6
      # sum(sapply(testx, function(x) x*dlnorm(x,mu,sd)))/sum(sapply(testx, function(x) dlnorm(x,mu,sd)))
      mn = 0.0 + 0.1*(xl-1.2) + 0.3*(iswaged-0.7)  + 0.3*((worker$gender[i]=="M")-0.5) + 0.02*(worker$year[i]-1980)/10
      sd = 0.6;
      worker$x[i] = sim_exposure(mn, sd, limit)  
      cumx = cumx + worker$x[i]
    }
    worker$cumx[i] = cumx;
    # mortality
    # cause of interest
    worker$rate[i] =              plogis(logit(worker$rate[i])           + 0.02*(worker$cumx[i]-12) + U + 2*(iswaged-0.7) + escalate(worker$age[i], 95, 0.5))
    #worker$rate[i] = plogis(logit(worker$rate[i]) + 0.01*worker$cumx[i] + U + -iswaged)
    worker$d1[i] = rbinom(1, 1,  worker$rate[i])
    # all other causes
    if(worker$d1[i] == 0){
      worker$rate_allcause[i] = plogis(logit(worker$rate_allcause[i]) + 0.01*(worker$cumx[i]-12) + U + 2*(iswaged-0.7) + escalate(worker$age[i], 95, 0.5))
      #worker$rate_allcause[i] = plogis(logit(worker$rate_allcause[i]) + 0.01*worker$cumx[i] + U + -iswaged)
      worker$d2[i] = rbinom(1, 1,  worker$rate_allcause[i])
    }
    if(worker$d1[i]==1 || worker$d2[i]==1) alive = 0
    xl = worker$x[i]
  }
  worker[1:i,]
}


worker_gen0row <- function(row, minorrateset){
  worker_gen0(row$id, row$age_hire, row$age_term, row$year_hire, row$year_eof, row$wage, row$gender, row$race, row$baseline_frailty, minorrateset)
}

#' worker_gen0
#' @description Generating longitudinal worker data under no exposure
#' 
#' @param age_hire 
#'
#' @param age_term 
#' @param year_hire 
#' @param year_eof 
#' @param wage 
#'
#' @examples 
#' age_hire = 25
#' age_term = 26
#' year_hire = 1932
#' year_eof = 1939
#' wage = 1
#' id = 10
#' gender = "M"
#' race = "W"
worker_gen0 <- function(id, age_hire, age_term, year_hire, year_eof, wage, gender, race, baseline_frailty, minorrateset){
  ageendpts = c("[15,20)","[20,25)","[25,30)","[30,35)","[35,40)","[40,45)",
                "[45,50)","[50,55)","[55,60)","[60,65)","[65,70)","[70,75)",
                "[75,80)","[80,85)","[85, Inf)")
  
  yearendpts = c("[1960,1965)","[1965,1970)","[1970,1975)","[1975,1980)",
                 "[1980,1985)","[1985,1990)","[1990,1995)","[1995,2000)",
                 "[2000,2005)","[2005,2010)","[2010,2015)","[2015,2020)",
                 "[2020,2025)","[2025,2030)")
  # worker history, if unexposed
  # need to fix: risk function for termination will depend on u, and exposure will affect employment
  (year = seq(year_hire, year_eof))
  
  (gender = rep(gender, length(year)))
  (race = rep(race, length(year)))
  (age = age_hire + year - year_hire)
  (agecat = sapply(age, function(x) catify(x, ageendpts)))
  (yearcat = sapply(year, function(x) catify(x, yearendpts)))
  worker = data.frame(gender,race,age,year,agecat,yearcat)
  worker = merge(worker, minorrateset, 
                 by.x=c("agecat", "yearcat", "gender", "race"), 
                 by.y=c("ageCat", "CPCat", "gender", "race"), 
                 all.x=TRUE, all.y=FALSE, sort=FALSE
  )
  #
  worker$id = (id = rep(id, length(year)))
  worker$wagestatus = rep(wage, length(year))
  #worker$atwork = as.numeric(age <= age_term) 
  worker$baseline_frailty = rep(baseline_frailty, length(year)) 
  worker 
}




catify <- function(val, cats){
  # note: this will give a category for any valid value (ignores start of first interval and end of last interval)
  (endpts = as.numeric(gsub("\\[([0-9Inf]{2,4}),([ 0-9Inf]{2,4})\\)","\\2" ,cats)))
  (idx = which(val < endpts))
  if(length(idx)==0) idx = c(length(endpts))
  cats[idx[1]]
}

getonecauserates <- function(rateset, cause=16){
  rateset[rateset$minor==cause,]
}


getallcauserates <- function(rateset){
  nms = names(rateset)
  rateset$allcats = with(rateset, paste(gender, race, CPCat, ageCat))
  sumset = data.frame(tapply(rateset$rate, rateset$allcats, sum))
  names(sumset) <- c("rate_allcause")
  sumset$allcats = rownames(sumset)
  rownames(sumset) = NULL
  nmset = rateset[rateset$minor==119,c("allcats", "gender", "race", "CPCat", "ageCat")]
  #sumset = merge(nmset, sumset, by=c("allcats"), all=TRUE)
  fullrates = merge(rateset, sumset, by="allcats")
  fullrates[,c(nms, "rate_allcause")]
}


escalate <- function(age, threshold, rate){
  # accelerate baseline risk of death at much older ages: outputs a log odds ratio # years over over threshold
  (age > threshold)*(age-threshold)*rate
}
