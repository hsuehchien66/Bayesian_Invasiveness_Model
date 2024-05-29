## Generating simulation data 
setwd("/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/simulation")
library(tibble)
library(progressionEstimation)
#library(rstan)
library(cmdstanr)
library(dplyr)
require(tidyverse)
require(magrittr)
library(posterior)
library(bayesplot)

process_2feature_input_data <- function(input_df, main_feature = "feature1", feature2 = "feature2"){
  
  if (!(main_feature %in% colnames(input_df))) {
    stop("Type column not in input data")
  }
  
  # Convert to factors
  input_df %<>%
    dplyr::mutate(study = factor(study)) %>%
    dplyr::mutate((!!main_feature) := factor(!!! dplyr::syms(main_feature))) %>%
    dplyr::mutate((!!feature2) := factor(!!! dplyr::syms(feature2))) 

  # Calculate input
  input_df <- as.data.frame(input_df)
  i_values <- as.integer(input_df$study)
  j_values <- as.integer(input_df[, which(colnames(input_df) == main_feature)])
  g1_values <- as.integer(input_df[, which(colnames(input_df) == feature2)])
  c_ijg1 <- input_df$carriage
  d_ijg1 <- input_df$disease
  n_i <- input_df$carriage_samples
  N_i <- input_df$surveillance_population
  t_i <- input_df$time_interval
  progression_rate_data <- list(
    i_max = max(i_values),
    j_max = max(j_values),
    g1_max = max(g1_values),
    n_obs = length(c_ijg1),
    i_values = i_values,
    j_values = j_values,
    g1_values = g1_values,
    c_ijg1 = c_ijg1,
    d_ijg1 = d_ijg1,
    n_i = n_i,
    N_i = N_i,
    t_i = t_i
  )
  return(progression_rate_data)
}

get_mean<-function(parameter,model) {
  #return(as.numeric(rstan::summary(model,pars=c(parameter))$summary[,"mean"]))
  return(model$summary(variables = parameter, "mean")$mean)
}

get_lower<-function(parameter,model) {
  #return(as.numeric(rstan::summary(model,pars=c(parameter))$summary[,"97.5%"]))
  return(model$summary(variables = parameter, 
                       extra_quantiles = ~posterior::quantile2(., probs = .0275))$q2.75)
}

get_upper<-function(parameter,model) {
  #return(as.numeric(rstan::summary(model,pars=c(parameter))$summary[,"2.5%"]))
  return(model$summary(variables = parameter, 
                       extra_quantiles = ~posterior::quantile2(., probs = .975))$q97.5)
  
}

get_median<-function(parameter,model) {
  #return(as.numeric(rstan::summary(model,pars=c(parameter))$summary[,"50%"]))
  return(model$summary(variables = parameter, "median")$median)
}

'''
fit$summary(
  variables = NULL,
  posterior::default_summary_measures(),
  extra_quantiles = ~posterior::quantile2(., probs = c(.0275, .975))
)
'''

setwd("/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/simulation")
# ------------------------- Carriage Rate -------------------------------

'''
Prevalance / Carriage rate
All simulations use the same randomly-generated set of carriage prevalance (rho_ij) for different types different locations
'''

n_types <- 20 
n_studies <- 20
n_gene1 <- 2

n_types <- 5 
n_studies <- 5
n_gene1 <- 2

# Simulate rho_ijg1g2
set.seed(284)
random_prevalance <- runif(n_types*n_studies*n_gene1, min = 0, max = 0.05) # 0.05 because there are 20 types (n_types * prevalance cannot exceed 1)

# generate table 
carriage_prevalance <- 
  tibble(
    "serotype" = factor(rep(paste0("Serotype_", 1:n_types), n_studies*n_gene1), levels = paste0("Serotype_", 1:n_types)), 
    "study" = factor(rep(paste0("Population_", LETTERS[1:n_studies]), each = n_types*n_gene1)),
    "gene1" = factor(rep(paste0("Gene1_", LETTERS[1:n_gene1]), each = n_types*n_studies)), 
    "true_rho_ijg1" = random_prevalance
  )

carriage_prevalance %>%
  kableExtra::kable()  %>%
  kableExtra::kable_styling(latex_options = "scale_down")



# ------------------------- Serotype Progression Rate -------------------------------

'''
The progresssion rates for the types (nu_j) are set to be vary by type
'''
set.seed(284)
random_progression_rates <- 10**(-1*runif(n_types, min = 1, max = 5))

serotype_progression_rates <-
  tibble(
    "serotype" = factor(paste0("Serotype_", 1:n_types), levels = paste0("Serotype_", 1:n_types)), 
    "true_nu_j" = random_progression_rates
  )

serotype_progression_rates %>%
  kableExtra::kable()  %>%
  kableExtra::kable_styling(latex_options = "scale_down")


# ------------------------- Population Factor -------------------------------

'''
The population factors are randomly generated. To simulate different levels of disease surveillance between locations, 
the (gamma_i) parameters were set such that the first location was used as the standard, reltative to which others varied
'''

# Sample size: n_i
set.seed(284)
n_i <- round(rnorm(n_studies, 1000, 250))
# Population size: N_i
set.seed(284)
N_i <- round(10**(runif(n_studies, min = 4, max = 6)))
# Population factor: r_i
set.seed(284)
random_r_i <- 10**rnorm(n_studies-1, mean = 0, sd=2) 
random_r_i <- c(1, random_r_i) # first population factor is 1

set.seed(284)
random_t <- round(runif(n_studies, min = 1, max = 4))

population_factor <-
  tibble(
    "study" = factor(paste0("Population_", LETTERS[1:n_studies])), 
    "true_r_i" = random_r_i, 
    "carriage_samples" = n_i,
    "surveillance_population" = N_i, 
    "time_interval" = random_t
  )

population_factor %>%
  kableExtra::kable()  %>%
  kableExtra::kable_styling(latex_options = "scale_down")


# ------------------------- Gene1 relative progression rate -------------------------------
set.seed(284)
random_g1 <- c(1, 1.5)

gene1_progression_rate <-
  tibble(
    "gene1" = factor(paste0("Gene1_", LETTERS[1:n_gene1])), 
    "true_nu_g1" = random_g1
  )

gene1_progression_rate %>%
  kableExtra::kable()  %>%
  kableExtra::kable_styling(latex_options = "scale_down")


# ------------------------- Simulated dataframe -------------------------------
simulation_df <- 
  carriage_prevalance %>%
  dplyr::left_join(serotype_progression_rates, by=c("serotype")) %>%
  dplyr::left_join(population_factor, by=c("study")) %>%
  dplyr::left_join(gene1_progression_rate, by=c("gene1"))


# ------------------------- Serotype-based_gene1_adjusted_model -------------------------------
## generate simulated observed data
serotype_based_gene1_adjusted_poisson_simulation_df <-
  simulation_df %>%
  dplyr::rowwise() %>%
  dplyr::mutate(carriage = rbinom(1, carriage_samples, prob = true_rho_ijg1)) %>%
  dplyr::mutate(disease = rpois(1, surveillance_population*true_rho_ijg1*time_interval*true_nu_j*true_nu_g1*true_r_i)) %>%
  dplyr::ungroup()

## making into a list for BIM 
serotype_based_gene1_adjusted_poisson_simulation_BIM_list <-
  process_2feature_input_data(input_df = serotype_based_gene1_adjusted_poisson_simulation_df, 
                              main_feature = "serotype", feature2 = "gene1")


# ------------------------- Model Fitting -------------------------------
## general hyperparameters
n_chains <- 2
n_iter <- 1e4
n_core <- 2

#options(mc.cores = parallel::detectCores())

## serotype_determined_gene_adjusted_poisson
#serotype_based_gene1_adjusted_poisson_stan <- stan_model("/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM/stan/serotype_determined_gene_adjusted_poisson.stan")
serotype_based_gene1_adjusted_poisson_stan <- cmdstanr::cmdstan_model("/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM/stan/serotype_determined_gene_adjusted_poisson.stan")

## if using MCMC
#serotype_based_gene1_adjusted_poisson_stan$sample(data = serotype_based_gene1_adjusted_poisson_simulation_BIM_list, 
                                                  iter_sampling = iter, chains = n_chains, parallel_chains = n_core)

## if using variational inference
serotype_based_gene1_adjusted_poisson_simulation_fit <- serotype_based_gene1_adjusted_poisson_stan$variational(data = serotype_based_gene1_adjusted_poisson_simulation_BIM_list)
#variational_fit$summary(variables = c("rho_ijg1"))

variational_fit_df <- variational_fit$draws(format = "df")
variational_fit_summary <- variational_fit$summary(variables = "nu_j", "mean", "sd")


serotype_based_gene1_adjusted_poisson_simulation_fit <- 
  rstan::sampling(serotype_jgbased_gene1_adjusted_poisson_stan, 
                  serotype_based_gene1_adjusted_poisson_simulation_BIM_list, 
                  iter=n_iter, 
                  chains=n_chains, 
                  cores=n_core)

summary_fit <- summary(serotype_based_gene1_adjusted_poisson_simulation_fit)

rhat_values <- summary_fit$summary[, "Rhat"]
ess_values <- summary_fit$summary[, "n_eff"]
mcmc_rhat(rhat_values)
mcmc_neff(ess_values)


output_df <- serotype_based_gene1_adjusted_poisson_simulation_df
i_levels = levels(output_df %>% dplyr::pull(study))
j_levels = levels(output_df %>% dplyr::pull(serotype))
g1_levels = levels(output_df %>% dplyr::pull(gene1))


output_carriage_df <- data.frame(
  #"study" = i_levels
  #"type" = j_levels
  "rho" = get_median("rho_ijg1",serotype_based_gene1_adjusted_poisson_simulation_fit),
  "rho_lower" = get_lower("rho_ijg1",serotype_based_gene1_adjusted_poisson_simulation_fit),
  "rho_upper" = get_upper("rho_ijg1",serotype_based_gene1_adjusted_poisson_simulation_fit)
)


output_carriage_df <- data.frame(
  #"study" = i_levels
  #"type" = j_levels
  "rho" = get_median("rho_ijg1",serotype_based_gene1_adjusted_poisson_simulation_fit),
  "rho_lower" = get_lower("rho_ijg1",serotype_based_gene1_adjusted_poisson_simulation_fit),
  "rho_upper" = get_upper("rho_ijg1",serotype_based_gene1_adjusted_poisson_simulation_fit)
)


output_df %<>% dplyr::bind_cols(output_carriage_df)


output_location_parameters <- data.frame(
  "study" = i_levels,
  "gamma" = get_median("gamma_i",serotype_based_gene1_adjusted_poisson_simulation_fit),
  "gamma_lower" = get_lower("gamma_i",serotype_based_gene1_adjusted_poisson_simulation_fit),
  "gamma_upper" = get_upper("gamma_i",serotype_based_gene1_adjusted_poisson_simulation_fit)
)

output_df %<>% dplyr::left_join(output_location_parameters, by = c("study"="study"))

output_serotype_progression_rates_df <- data.frame(
  "serotype" = j_levels,
  "nu_j" = get_median("nu_j",serotype_based_gene1_adjusted_poisson_simulation_fit),
  "nu_j_lower" = get_lower("nu_j",serotype_based_gene1_adjusted_poisson_simulation_fit),
  "nu_j_upper" = get_upper("nu_j",serotype_based_gene1_adjusted_poisson_simulation_fit)
)

output_df %<>% dplyr::left_join(output_serotype_progression_rates_df, by = setNames("serotype","serotype"))


output_gene1_progression_rates_df <- data.frame(
  "gene1" = g1_levels,
  "nu_g1" = get_median("nu_g1",serotype_based_gene1_adjusted_poisson_simulation_fit),
  "nu_g1_lower" = get_lower("nu_g1",serotype_based_gene1_adjusted_poisson_simulation_fit),
  "nu_g1_upper" = get_upper("nu_g1",serotype_based_gene1_adjusted_poisson_simulation_fit)
)

output_df %<>% dplyr::left_join(output_gene1_progression_rates_df, by = setNames("gene1","gene1"))


output_diff_df <- output_df %>%
  dplyr::mutate(rho_diff = abs(true_rho_ijg1 - rho)/true_rho_ijg1) %>%
  dplyr::mutate(nu_j_diff = abs(true_nu_j - nu_j)/true_nu_j) %>%
  dplyr::mutate(nu_g1_diff = abs(true_nu_g1 - nu_g1)/true_nu_g1) %>%
  dplyr::mutate(r_diff = abs(true_r_i - gamma)/true_r_i) 




write.table(output_diff_df, "Serotype_specific_popadjusted_onegene_poisson_simulation_df_20types_5pops_2g1.txt", sep="\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

serotype_prate_plot <-
  ggplot(output_diff_df,
         aes(x = true_nu_j,
             y = nu_j,
             ymin = nu_j_lower,
             ymax = nu_j_upper)) +
  geom_abline(slope = 1, intercept = 0, lty = 2, colour = "coral") +
  ylab("Predicted serotype progression rate") +
  xlab("Observed serotype progression rate") +
  theme_bw()

serotype_prate_plot <- serotype_prate_plot +
  geom_point(aes(color = serotype)) +
  geom_errorbar(aes(color = serotype), alpha = 0.75)+
  theme(
    axis.title = element_text(size = 25), 
    axis.text = element_text(size = 20), 
    legend.text = element_text(size = 20), 
    legend.title = element_text(size = 25)
  )

png("Serotype_specific_popadjusted_onegene_poisson_simulation_serotype_prate.png", width = 1000, height = 1000)
serotype_prate_plot
dev.off()


gene1_prate_plot <-
  ggplot(output_diff_df,
         aes(x = true_nu_g1,
             y = nu_g1,
             ymin = nu_g1_lower,
             ymax = nu_g1_upper)) +
  geom_abline(slope = 1, intercept = 0, lty = 2, colour = "coral") +
  ylab("Predicted gene1 progression rate") +
  xlab("Observed gene1 progression rate") +
  theme_bw()

gene1_prate_plot <- gene1_prate_plot +
  geom_point(aes(color = gene1)) +
  geom_errorbar(aes(color = gene1), alpha = 0.75)+
  theme(
    axis.title = element_text(size = 25), 
    axis.text = element_text(size = 20), 
    legend.text = element_text(size = 20), 
    legend.title = element_text(size = 25)
  )

png("Serotype_specific_popadjusted_onegene_poisson_simulation_gene1_prate.png", width = 1000, height = 1000)
gene1_prate_plot
dev.off()


rho_plot <-
  ggplot(output_diff_df,
         aes(x = true_rho_ijg1,
             y = rho,
             ymin = rho_lower,
             ymax = rho_upper)) +
  geom_abline(slope = 1, intercept = 0, lty = 2, colour = "coral") +
  ylab("Predicted carriage rate") +
  xlab("Observed carriage rate") +
  theme_bw()

rho_plot <- rho_plot +
  geom_point(aes(color = study)) +
  geom_errorbar(aes(color = study), alpha = 0.75)+
  theme(
    axis.title = element_text(size = 25), 
    axis.text = element_text(size = 20), 
    legend.text = element_text(size = 20), 
    legend.title = element_text(size = 25)
  )

png("Serotype_specific_popadjusted_onegene_poisson_simulation_rho.png", width = 1000, height = 1000)
rho_plot
dev.off()

gamma_plot <-
  ggplot(output_diff_df,
         aes(x = true_r_i,
             y = gamma,
             ymin = gamma_lower,
             ymax = gamma_upper)) +
  geom_abline(slope = 1, intercept = 0, lty = 2, colour = "coral") +
  ylab("Predicted gamma rate") +
  xlab("Observed gamma rate") +
  theme_bw()

gamma_plot <- gamma_plot +
  geom_point(aes(color = study)) +
  geom_errorbar(aes(color = study), alpha = 0.75)+
  theme(
    axis.title = element_text(size = 25), 
    axis.text = element_text(size = 20), 
    legend.text = element_text(size = 20), 
    legend.title = element_text(size = 25)
  )

png("Serotype_specific_popadjusted_onegene_poisson_simulation_gamma.png", width = 1000, height = 1000)
gamma_plot
dev.off()
