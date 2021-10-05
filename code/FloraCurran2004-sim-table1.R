# repliSims replication of Flora & Curran ("F&C", 2004)
# Additional simulation study to reproduce Table 1: Data generation
# Written by Y. Andre Wang
# Last update: October 5, 2021


# Prepare Packages --------------------------------------------------------

# Name the packages needed
packages <- c("moments", "lavaan", "SimDesign")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == F)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = T))



# Cretae a data frme of simulation conditions -----------------------------

Design <- data.frame(cond = 1:5,
                     skew = c(0, .75, .75, 1.25, 1.25),
                     kurt = c(0, 1.75, 3.75, 1.75, 3.75))
Design



# Create function implementing Fleishman's (1978) method ------------------

# For simulating univariate data with nonzero skewness and kurtosis
# Adapted from lavaan::simulateData()

fleishman1978_abcd <- function(skewness, kurtosis) {
  system.function <- function(x, skewness, kurtosis) {
    b.=x[1L]; c.=x[2L]; d.=x[3L]
    eq1 <- b.^2 + 6*b.*d. + 2*c.^2 + 15*d.^2 - 1
    eq2 <- 2*c.*(b.^2 + 24*b.*d. + 105*d.^2 + 2) - skewness
    eq3 <- 24*(b.*d. + c.^2*(1 + b.^2 + 28*b.*d.) +
                 d.^2*(12 + 48*b.*d. + 141*c.^2 + 225*d.^2)) - kurtosis
    eq <- c(eq1,eq2,eq3)
    sum(eq^2)
  }
  
  out <- nlminb(start = c(1, 0, 0), objective = system.function,
                scale = 10, control = list(trace = 0),
                skewness = skewness, kurtosis = kurtosis)
  if(out$convergence != 0) warning("no convergence")
  b. <- out$par[1L]; c. <- out$par[2L]; d. <- out$par[3L]; a. <- -c.
  c(a.,b.,c.,d.)
}



# Create simulation function ----------------------------------------------

fc1 <- function(ksim = 50000, cond, seed = 3012) {
  
  # Set seed
  set.seed(seed)
  
  # Generate normal data of sample size 50000
  dat_x <- rnorm(n = ksim, mean = 0, sd = 1)
  
  # Skewness and kurtosis
  coef_f <- fleishman1978_abcd(skewness = Design$skew[cond], 
                               kurtosis = Design$kurt[cond])
  dat_y <- coef_f[1] + coef_f[2]*dat_x + coef_f[3]*dat_x^2 + coef_f[4]*dat_x^3
  
  # Create orginal data
  dat <- numeric(50000)
  
  if (cond != 4) {
    tau <- c(-1.645, -0.643, 0.643, 1.645)
    for (i in 1:5) {
      dat[dat_y > tau[i]] <- i
    }
  } 
  
  else {
    tau <- c(-1.125, -0.643, 0.643, 1.645)
    for (i in 1:5) {
      dat[dat_y > tau[i]] <- i
    }
  }
  
  # Record skewness and kurtosis of ordinal data
  return(data.frame(skewness_con = Design$skew[cond],
                    kurtosis_con = Design$kurt[cond],
                    skewness_ord = moments::skewness(dat),
                    kurtosis_ord = moments::kurtosis(dat) - 3))
}



# Run simulations and store data ------------------------------------------

results_compiled <- NULL
for (i in 1:nrow(Design)){
  results <- fc1(cond = i)
  results_compiled <- rbind(results_compiled, results)
}
results_compiled
# saveRDS(results_compiled, "FloraCurran2004-results-table1.rds")
