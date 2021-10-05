# repliSims replication of Flora & Curran ("F&C", 2004)
# Additional simulation study to reproduce Table 2: Data generation
# Written by Y. Andre Wang, based on Chalmer & Adkins (2020)
# Last update: October 5, 2021


# Prepare Packages --------------------------------------------------------

# Name the packages needed
packages <- c("psych", "SimDesign")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == F)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = T))

# Source function used to generate data with skewness and kurtosis
source('ValeMaurelli1983.R')



# DESIGN ------------------------------------------------------------------

Design <- createDesign(
  sample.size = c(100, 200, 500, 1000),
  categories = c(2, 5),
  skewness_kurtosis = list(c(0, 0), c(.75, 1.75), c(.75, 3.75), c(1.25, 1.75), 
                           c(1.25, 3.75)),
  rho = c(.147, .49))
Design



# GENERATE ----------------------------------------------------------------

Generate <- function(condition, fixed_objects = NULL) {
  
  Attach(condition)
  
  cdat <- ValeMaurelli1983(
    n = sample.size,
    COR = matrix(c(1, rho, rho, 1), nrow = 2, ncol = 2, byrow = T),
    skewness = skewness_kurtosis[1L],
    kurtosis = skewness_kurtosis[2L])
  
  tau <- if(categories == 5) c(-1.645, -0.643, 0.643, 1.645) else 0
  
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
  if (!all(apply(dat, 2, function(x) length(unique(x))) == categories))
    stop('Number of categories generated is incorrect')
  
  dat
}



# ANALYSE -----------------------------------------------------------------

Analyse <- function(condition, dat, fixed_objects = NULL) {
  
  rho_est <- polychoric(dat)$rho[1, 2]
  
  ret <- c(rho_est = rho_est)
}
  


# SUMMARISE ---------------------------------------------------------------

Summarise <- function(condition, results, fixed_objects = NULL) {

  Attach(condition)
  
  rho_mean <- mean(results$rho_est)
  rho_sd <- sd(results$rho_est)
  rho_RB <- bias(results$rho_est, parameter = rho, type = 'relative',
                percent = TRUE, unname = TRUE)
  
  ret <- c(rho_mean = rho_mean, rho_sd = rho_sd, rho_RB = rho_RB)
  ret
}



# RUN SIMULATION ----------------------------------------------------------

res <- runSimulation(design = Design, replications = 500, generate = Generate,
                     analyse = Analyse, summarise = Summarise, max_errors = 100,
                     packages = 'psych', parallel = F,
                     filename = 'FloraCurran2004-results-table2', 
                     save_results = T)
res
