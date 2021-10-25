# repliSims results on Flora & Curran ("F&C", 2004)
# Written by Y. Andre Wang & Udi Alter
# Last update: September 1, 2021

# Load packages
if(!require(dplyr)){install.packages('dplyr')}
if(!require(ggplot2)){install.packages('ggplot2')}
if(!require(SimDesign)){install.packages('SimDesign')}
library(dplyr); library(ggplot2); library(SimDesign)

# Read primary simulation replication data (Tables 4-11 and Figure 6)
res <- readRDS("FloraCurran2004_210826.rds")
head(res)

# Format data
res %>% 
  
  # Label skewness (s), kurtosis (k), Type I error rate; calculate RB for lambda
  mutate(s = sapply(res$skewness_kurtosis, function(x) x[1]),
         k = sapply(res$skewness_kurtosis, function(x) x[2]),
         "% reject" = edr_05 * 100) %>%
  
  # Keep only 5-category ordinal data for constructing tables, per F&C
  filter(categories == 5) -> 
  
  res_dat

# Read additional simulation replication data (Tables 1-3)
dat_tab1 <- readRDS("FloraCurran2004_Table1.rds")
dat_tab2 <- readRDS("FloraCurran2004_Table2.rds")
dat_tab3 <- readRDS("FloraCurran2004_Table3.rds")



# Table 1 -----------------------------------------------------------------

# Skewness and Kurtosis of Univariate Latent Response (y*) Distributions and
# Five-Category Ordinal (y) Distributions

dat_tab1$skewness_ord <- round(dat_tab1$skewness_ord, 2)
dat_tab1$kurtosis_ord <- round(dat_tab1$kurtosis_ord, 2)
dat_tab1


# Table 2 -----------------------------------------------------------------

# Means, Standard Deviations, and Relative Bias of 
# Polychoric Correlation Estimates

tab2_rho1 <- dat_tab2 %>% 
  mutate(s = sapply(skewness_kurtosis, function(x) x[1]),
         k = sapply(skewness_kurtosis, function(x) x[2])) %>%
  arrange(categories, sample.size) %>%
  filter(rho == .147) %>%
  select(sample.size, s, k, categories, rho, rho_mean, rho_sd, rho_RB)

tab2_rho2 <- dat_tab2 %>% 
  mutate(s = sapply(skewness_kurtosis, function(x) x[1]),
         k = sapply(skewness_kurtosis, function(x) x[2])) %>%
  arrange(categories, sample.size) %>%
  filter(rho == .49) %>%
  select(sample.size, s, k, categories, rho, rho_mean, rho_sd, rho_RB)

# rho = .147 for 2 categories and 5 categories
tab2_rho1

# rho = .49 for 2 categories and 5 categories
tab2_rho2



# Table 3 -----------------------------------------------------------------

# Rates of Improper Solutions Obtained With Full WLS Estimation

# ip = rate of improper solutions
# nc = rate of nonconverged solutions

# 2 categories model
dat_tab3[[1]]

# 5 categories model
dat_tab3[[2]]



# *- Table 4 --------------------------------------------------------------

# Chi-Square Test Statistics for Model 1 (Five Indicators, One Factor)
res_dat %>%
  filter(indicators == 5, factors == 1, estimator == 'WLS') %>%
  select(N, s, k, M_X2, SD_X2, RB_X2, "% reject") -> WLS_4
res_dat %>%
  filter(indicators == 5, factors == 1, estimator == 'WLSMV') %>%
  select(N, s, k, RB_X2, "% reject") -> WLSMV_4
full_join(WLS_4, WLSMV_4, by = c("N", "s", "k"), suffix= c('', '.DWLS')) %>%
  arrange(N, s, k) -> tab4
tab4



# *- Table 5 --------------------------------------------------------------

# Chi-Square Test Statistics for Model 2 (10 Indicators, One Factor)
res_dat %>%
  filter(indicators == 10, factors == 1, estimator == 'WLS') %>%
  select(N, s, k, M_X2, SD_X2, RB_X2, "% reject") -> WLS_5
res_dat %>%
  filter(indicators == 10, factors == 1, estimator == 'WLSMV') %>%
  select(N, s, k, RB_X2, "% reject") -> WLSMV_5
full_join(WLS_5, WLSMV_5, by = c("N", "s", "k"), suffix= c('', '.DWLS')) %>%
  arrange(N, s, k) -> tab5
tab5



# *- Table 6 --------------------------------------------------------------

# Chi-Square Test Statistics for Model 3 (10 Indicators, Two Correlated Factors)
res_dat %>%
  filter(indicators == 5, factors == 2, estimator == 'WLS') %>%
  select(N, s, k, M_X2, SD_X2, RB_X2, "% reject") -> WLS_6
res_dat %>%
  filter(indicators == 5, factors == 2, estimator == 'WLSMV') %>%
  select(N, s, k, RB_X2, "% reject") -> WLSMV_6
full_join(WLS_6, WLSMV_6, by = c("N", "s", "k"), suffix= c('', '.DWLS')) %>%
  arrange(N, s, k) -> tab6
tab6



# *- Table 7 --------------------------------------------------------------

# Chi-Square Test Statistics for Model 4 (20 Indicators, Two Correlated Factors)
res_dat %>%
  filter(indicators == 10, factors == 2, estimator == 'WLS') %>%
  select(N, s, k, M_X2, SD_X2, RB_X2, "% reject") -> WLS_7
res_dat %>%
  filter(indicators == 10, factors == 2, estimator == 'WLSMV') %>%
  select(N, s, k, RB_X2, "% reject") -> WLSMV_7
full_join(WLS_7, WLSMV_7, by = c("N", "s", "k"), suffix= c('', '.DWLS')) %>%
  arrange(N, s, k) -> tab7
tab7

# F&B reported results for N = 200 but they are not part of the replication here
# due to nonconvergence (i.e., no model converged in our replications)
# See replication of Table 3 for details on improper solutions



# *- Figure 6 -------------------------------------------------------------

# Generate summary data for Figure 6 in F&C
res_dat %>%
  filter(indicators == 5) %>%
  mutate(TypeI = edr_05 * 100,
         est_fact = factor(estimator):factor(factors)) %>%
  group_by(N, est_fact) %>%
  summarise(mean_TypeI = mean(TypeI)) -> fig6dat
labels <- c("Model 1 (full WLS estimation)", "Model 3 (full WLS estimation)",
            "Model 1 (robust WLS estimation)", "Model 3 (robust WLS estimation)")

# Reproduce Figure 6 using the ggplot2 package
ggplot(fig6dat,
       aes(x = factor(N), y = mean_TypeI, linetype = est_fact,
            shape = est_fact, group = est_fact)) +
  geom_line() + geom_point(size = 4) +
  geom_hline(yintercept = 5, colour = 'red', linetype = 'dashed') +
  xlab("Sample Size") + ylab("Type I Error Rate") +
  scale_linetype_discrete(labels = labels) +
  scale_shape_discrete(labels = labels) +
  ylim(0, 80) + theme_bw() +
  theme(legend.title= element_blank(),
        legend.position = c(.8, .7),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()
  )



# Table 8 -----------------------------------------------------------------

# Mean Factor Loadings for Model 1 (Five Indicators, One Factor)
res_dat %>%
  filter(indicators == 5, factors == 1, estimator == 'WLS') %>%
  select(N, s, k, M_lambda, SD_lambda, RB_lambda) -> WLS_8
res_dat %>%
  filter(indicators == 5, factors == 1, estimator == 'WLSMV') %>%
  select(N, s, k, M_lambda, SD_lambda, RB_lambda) -> WLSMV_8
full_join(WLS_8, WLSMV_8, by = c("N", "s", "k"), suffix= c('', '.DWLS')) %>%
  arrange(N, s, k) -> tab8
tab8



# Table 9 -----------------------------------------------------------------

# Mean Factor Loadings for Model 2 (10 Indicators, One Factor)
res_dat %>%
  filter(indicators == 10, factors == 1, estimator == 'WLS') %>%
  select(N, s, k, M_lambda, SD_lambda, RB_lambda) -> WLS_9
res_dat %>%
  filter(indicators == 10, factors == 1, estimator == 'WLSMV') %>%
  select(N, s, k, M_lambda, SD_lambda, RB_lambda) -> WLSMV_9
full_join(WLS_9, WLSMV_9, by = c("N", "s", "k"), suffix= c('', '.DWLS')) %>%
  arrange(N, s, k) -> tab9
tab9



# Table 10 ----------------------------------------------------------------

# Mean Parameter Estimates for Model 3 (10 Indicators, Two Factors)
res_dat %>%
  filter(indicators == 5, factors == 2, estimator == 'WLS') %>%
  select(N, s, k, M_lambda, SD_lambda, RB_lambda, 
         M_phi21, SD_phi21, RB_phi21) -> WLS_10
res_dat %>%
  filter(indicators == 5, factors == 2, estimator == 'WLSMV') %>%
  select(N, s, k, M_lambda, SD_lambda, RB_lambda,
         M_phi21, SD_phi21, RB_phi21) -> WLSMV_10
full_join(WLS_10, WLSMV_10, by = c("N", "s", "k"), suffix= c('', '.DWLS')) %>%
  arrange(N, s, k) -> tab10
tab10



# Table 11 ----------------------------------------------------------------

# Mean Parameter Estimates for Model 4 (20 Indicators, Two Factors)
res_dat %>%
  filter(indicators == 10, factors == 2, estimator == 'WLS') %>%
  select(N, s, k, M_lambda, SD_lambda, RB_lambda,
         M_phi21, SD_phi21, RB_phi21) -> WLS_11
res_dat %>%
  filter(indicators == 10, factors == 2, estimator == 'WLSMV') %>%
  select(N, s, k, M_lambda, SD_lambda, RB_lambda,
         M_phi21, SD_phi21, RB_phi21) -> WLSMV_11
full_join(WLS_11, WLSMV_11, by = c("N", "s", "k"), suffix= c('', '.DWLS')) %>%
  arrange(N, s, k) -> tab11
tab11

# F&B reported results for N = 200 but they are not part of the replication here
# due to nonconvergence (i.e., no model converged in our replications)
# See replication of Table 3 for details on improper solutions
