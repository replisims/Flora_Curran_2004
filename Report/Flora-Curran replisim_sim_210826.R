# Load packages
if(!require(SimDesign)){install.packages('SimDesign')}
library(SimDesign)


#____________________DESIGN______________________________________

Design <- createDesign(
  N = c(100, 200, 500, 1000),
  categories = c(2, 5),
  skewness_kurtosis = list(c(0, 0), c(.75, 1.75), c(.75, 3.75),
                           c(1.25, 1.75), c(1.25, 3.75)),
  factors = c(1, 2),
  indicators = c(5, 10),
  estimator = c("WLSMV", "WLS"),
  
  # remove known problematic conditions
  subset = !(estimator == 'WLS' & N %in% c(100, 200) &
               factors == 2 & indicators == 10))
Design

# include syntax generation function
source('FloraCurran2004-functions.R')



#________________________GENERATE_______________________________

Generate <- function(condition, fixed_objects = NULL) {
  
  Attach(condition)
  
  syntax <- genLavaanSyntax(factors = factors, indicators = indicators)
  
  cdat <- simulateData(syntax, model.type = 'cfa', sample.nobs = N,
                       skewness = skewness_kurtosis[1L],
                       kurtosis = skewness_kurtosis[2L])
  
  tau <- if(categories == 5)
    c(-1.645, -0.643, 0.643, 1.645) else 0
  # data generation fix described in Flora's (2002) unpublished dissertation
  if (categories == 5 && all(skewness_kurtosis == c(1.25, 1.75)))
    tau[1] <- -1.125
  
  dat <- apply(cdat, 2, function(x, tau){
    dat <- numeric(length(x))
    for (i in 1:length(tau))
      dat[x > tau[i]] <- i
    dat
  }, tau = tau)
  
  # throw error if number of categories not correct
  # Not included here because it is not defined as improper solution by F&C
  #  if (!all( apply(dat, 2, function(x) length(unique(x))) == categories))
  #    stop('Number of categories generated is incorrect')
  
  dat
  
}



#________________________ANALYSE_______________________________

Analyse <- function(condition, dat, fixed_objects = NULL) {
  
  Attach(condition)
  
  syntax <- genLavaanSyntax(factors = factors, indicators = indicators, 
                            analyse = TRUE)
  
  mod <- cfa(syntax, dat, ordered = colnames(dat), estimator = estimator)
  
  # Check that model and coefficients are reasonable
  if (!lavInspect(mod, 'converged')) stop('Model did not converge')
  
  # Identify factor loading estimates
  pick_lambdas <- matrix(TRUE, indicators * factors, factors)
  if(factors == 2)
    pick_lambdas[(indicators + 1):(indicators * 3)] <- F
  cfs <- lavInspect(mod, what = "std")$lambda[pick_lambdas]
  
  # Check Heywood cases
  if(any(cfs > 1 | cfs < -1))
    stop('Model contains Heywood cases')
  if(factors > 2 && abs(lavInspect(mod, what="std")$psi[2,1]) >= 1)
    stop('Latent variable psi matrix not positive definite')
  
  # Extract desired results
  fit <- fitMeasures(mod)
  ses <- lavInspect(mod, what="se")$lambda[pick_lambdas]
  stat_names <- extract_stats <- c("chisq", "df", "pvalue")
  if(estimator == "WLSMV")
    extract_stats <- paste0(extract_stats, ".scaled")
  fitstats <- fit[extract_stats]
  names(fitstats) <- stat_names
  
  # Factor correlation
  phi21 <- if(factors == 2)
    lavInspect(mod, what = "std")$psi[1,2] else NULL
  
  # Compile
  ret <- c(fitstats, mean_ses = mean(ses), lambda = cfs, phi21 = phi21)
  ret
}

#________________________SUMMARISE_______________________________

Summarise <- function(condition, results, fixed_objects = NULL) {
  
  # Factor loadings
  lambdas <- results[ , grepl('lambda', colnames(results))]
  pool_mean_lambdas <- mean(apply(lambdas, 2, mean)) # Equation 10
  pool_SD_lambdas <- sqrt(mean(apply(lambdas, 2, var))) # Equation 11
  pool_RB_lambdas <- bias(pool_mean_lambdas, parameter = .7, 
                          type = 'relative', percent = T)
  
  # Factor correlations
  mean_phi21 <- if (condition$factors == 2) mean(results$phi21) else NULL
  SD_phi21 <- if (condition$factors == 2) sd(results$phi21) else NULL
  RB_phi21 <- if (condition$factors == 2)
    bias(results$phi21, parameter = .3, type = 'relative', percent = T) else NULL
  
  mean_se <- mean(results$mean_ses)
  
  # Goodness-of-fit
  edr_05 <- EDR(results$pvalue, alpha = .05)
  mean_X2 <- mean(results$chisq)
  sd_X2 <- sd(results$chisq)
  RB_X2 <- bias(results$chisq, parameter = results$df, type = 'relative',
                percent = T, unname = T)
  
  # Compile
  ret <- c(edr_05 = edr_05, M_X2 = mean_X2, SD_X2 = sd_X2, RB_X2 = RB_X2, 
           M_lambda = pool_mean_lambdas, SD_lambda = pool_SD_lambdas, 
           RB_lambda = pool_RB_lambdas, M_phi21 = mean_phi21,
           SD_phi21 = SD_phi21, RB_phi21 = RB_phi21, mean_se = mean_se)
  ret
}


#________________________RUNSIMULATION_______________________________

res <- runSimulation(design = Design, replications = 500, generate = Generate,
                     analyse = Analyse, summarise = Summarise, max_errors = 100,
                     packages = 'lavaan', parallel = T,
                     filename = 'FloraCurran2004_210826', save_results = T)
res
