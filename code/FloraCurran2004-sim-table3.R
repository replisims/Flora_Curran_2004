# repliSims replication of Flora & Curran ("F&C", 2004)
# Additional simulation study to reproduce Table 3: Data generation
# Written by Y. Andre Wang
# Last update: October 5, 2021


# Prepare Packages --------------------------------------------------------

# Name the packages needed
packages <- c("lavaan", "SimDesign", "tidyverse")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == F)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = T))

# include syntax generation function
source('FloraCurran2004-functions.R')



# Create a data frame of simulation conditions ----------------------------

Design <- createDesign(
  N = c(100, 200),
  categories = c(2, 5),
  skewness_kurtosis = list(c(0, 0), c(.75, 1.75), c(.75, 3.75),
                           c(1.25, 1.75), c(1.25, 3.75)),
  factors = c(1, 2),
  indicators = c(5, 10),
  estimator = "WLS")

Design <- mutate(Design, tau = ifelse(
  Design$categories == 2, 0, list(c(-1.645, -0.643, 0.643, 1.645))
))

# data generation fix described in Flora's (2002) unpublished dissertation
Design$tau[c(15:16, 35:36, 55:56, 75:76)] <- list(c(-1.125, -0.643, 0.643, 1.645))

# View design
Design



# Simulation study function -----------------------------------------------

fc3 <- function(ksim = 500, cond, seed = 42) {
  
  # Data generation model syntax
  syntaxGen <- genLavaanSyntax(factors = Design$factors[cond], 
                               indicators = Design$indicators[cond])
  
  # Data analysis model syntax
  syntaxAnz <- genLavaanSyntax(factors = Design$factors[cond], 
                               indicators = Design$indicators[cond], analyse = T)
  
  # Define tau (thresholds for transforming continuous data into ordinal data)
  Tau <- Design$tau[[cond]]
  
  # Set simulation seed
  set.seed(seed)
  
  # Define results objects
  nonconv <- 0 # nonconverging model count
  heywood <- 0 # Heywood case-containing model count
  
  
  # Loop by iteration
  for (i in 1:ksim) {
    
    # Simulate and store data
    cdat <- simulateData(syntaxGen, model.type = 'cfa', 
                         sample.nobs = Design$N[cond],
                         skewness = Design$skewness_kurtosis[[cond]][1L],
                         kurtosis = Design$skewness_kurtosis[[cond]][2L])
    
    
    # Transform continuous data into ordinal data
    dat <- apply(cdat, 2, function(x, Tau){
      dat <- numeric(length(x))
      for (i in 1:length(Tau))
        dat[x > Tau[i]] <- i
      dat
    }, Tau = Tau)
    
    # Test-fitting CFA to data
    test_model <- try(cfa(syntaxAnz, dat, ordered = colnames(dat), 
                          estimator = Design$estimator[cond]), silent = T)
    
    # If CFA was not successfully fitted (i.e., test_model did not run) 
    if (inherits(test_model, "try-error")) {
      
      nonconv <- nonconv + 1 # nonconvergence count + 1
      
    } else {
      
      mod <- cfa(syntaxAnz, dat, ordered = colnames(dat), 
                 estimator = Design$estimator[cond])
      
      # If CFA was fitted but mod did not converge
      if (!lavInspect(mod, 'converged')) {
        
        nonconv <- nonconv + 1 # nonconvergence count + 1
        
      } else {
        
        # Heywood cases
        pick_lambdas <- matrix(T, Design$indicators[cond]*Design$factors[cond], 
                               Design$factors[cond])
        
        if(Design$factors[cond] == 2)
          pick_lambdas[(
            Design$indicators[cond] + 1):(Design$indicators[cond]*3)] <- F
        
        cfs <- lavInspect(mod, what = "std")$lambda[pick_lambdas]
        
        if(any(cfs > 1 | cfs < -1)){
          
          heywood <- heywood + 1
          
        } else if(
          
          Design$factors[cond] > 2 && abs(
            lavInspect(mod, what = "std")$psi[2,1]) >= 1){
          
          heywood <- heywood + 1
          
        }
        
      }
      
    }
    
  }
  
  # record output
  return(data.frame(cond = cond, nonconv = nonconv, heywood = heywood))
  
}



# Run simulations and store data ------------------------------------------


# Create shell for storing outcomes
results_compiled <- NULL

# Loop through conditions and record outcomes in .csv files with condition #
for(i in 1:nrow(Design)) {
  results <- fc3(cond = i)
  results_compiled <- rbind(results_compiled, results)
  write.table(results_compiled, sep = ",", row.names = F, 
              paste("results_compiled_", i, ".csv", sep = ""))
}



# Compile data into table -------------------------------------------------

# Read compiled data file with all conditions
tab3 <- read.csv("results_compiled_80.csv")

# Combine data with design conditions
tab3 <- cbind(Design, tab3)

# Split skewness_kurtosis
tab3$skewness <- as.character(lapply(tab3$skewness_kurtosis, `[[`, 1))
tab3$kurtosis <- as.character(lapply(tab3$skewness_kurtosis, `[[`, 2))

# Calculate number of improper solutions
tab3$ip <- tab3$heywood + tab3$nonconv

# Rename nonconvergence
tab3 <- rename(tab3, nc = nonconv)

# Identify model number
tab3 <- mutate(tab3, mod = paste(factors, indicators))
tab3$mod <- recode(tab3$mod, 
                   "1 5" = "mod1", "1 10" = "mod2", 
                   "2 5" = "mod3", "2 10" = "mod4")

# Keep only variables needed
tab3 <- tab3[, c("N", "skewness", "kurtosis", "categories", "mod", "ip", "nc")]

# Convert numbers to percentages (k = 500)
tab3$ip <- tab3$ip / 500 * 100
tab3$nc <- tab3$nc / 500 * 100

# Two-category model
tab3_cat2 <- tab3 %>%
  filter(categories == 2) %>%
  select(N, skewness, kurtosis, mod, ip, nc) %>% 
  pivot_wider(names_from = mod, names_sep = ".", values_from = c(ip, nc)) %>%
  arrange(N)

# Reorder columns
tab3_cat2 <- tab3_cat2[, c(
  "N", "skewness", "kurtosis", "ip.mod1", "nc.mod1", 
  "ip.mod2", "nc.mod2", "ip.mod3", "nc.mod3",
  "ip.mod4", "nc.mod4")]

# Five-category model
tab3_cat5 <- tab3 %>%
  filter(categories == 5) %>%
  select(N, skewness, kurtosis, mod, ip, nc) %>% 
  pivot_wider(names_from = mod, names_sep = ".", values_from = c(ip, nc)) %>%
  arrange(N)

# Reorder columns
tab3_cat5 <- tab3_cat5[, c(
  "N", "skewness", "kurtosis", "ip.mod1", "nc.mod1", 
  "ip.mod2", "nc.mod2", "ip.mod3", "nc.mod3",
  "ip.mod4", "nc.mod4")]

# Combine into one R object
tab3_list <- list(tab3_cat2, tab3_cat5)

# saveRDS(tab3_list, "FloraCurran2004-results-table3.rds")
