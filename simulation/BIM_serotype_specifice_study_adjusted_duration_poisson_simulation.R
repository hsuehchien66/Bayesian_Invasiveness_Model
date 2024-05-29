## Generating simulation data 
setwd("/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/simulation")
library(tibble)
library(progressionEstimation)
library(rstan)
library(dplyr)
require(tidyverse)
require(magrittr)

get_mean<-function(parameter,model) {
  return(as.numeric(rstan::summary(model,pars=c(parameter))$summary[,"mean"]))
}

get_upper<-function(parameter,model) {
  return(as.numeric(rstan::summary(model,pars=c(parameter))$summary[,"97.5%"]))
}

get_lower<-function(parameter,model) {
  return(as.numeric(rstan::summary(model,pars=c(parameter))$summary[,"2.5%"]))
}

get_median<-function(parameter,model) {
  return(as.numeric(rstan::summary(model,pars=c(parameter))$summary[,"50%"]))
}
# ------------------------- Carriage Rate -------------------------------

'''
Prevalance / Carriage rate
All simulations use the same randomly-generated set of carriage prevalance (rho_ij) for different types different locations
'''

n_types <- 20 
n_studies <- 5
# Simulate rho_ij
set.seed(284)
random_prevalance <- runif(n_types*n_studies, min = 0, max = 0.05) # 0.05 because there are 20 types (n_types * prevalance cannot exceed 1)

# generate table 
carriage_prevalance <- 
  tibble(
    "serotype" = factor(rep(paste0("Serotype_", 1:n_types), n_studies), levels = paste0("Serotype_", 1:n_types)), 
    "study" = factor(rep(paste0("Population_", LETTERS[1:n_studies]), each = n_types)), 
    "true_rho_ij" = random_prevalance
  )

carriage_prevalance %>%
  kableExtra::kable()  %>%
  kableExtra::kable_styling(latex_options = "scale_down")



# ------------------------- Progression Rate -------------------------------

'''
The progresssion rates for the types (nu_j) are set to be vary by type

set.seed(284)
random_progression_rates <- 10**(-1*runif(n_types, min = 1, max = 5))

progression_rates <-
  tibble(
    "serotype" = factor(paste0("Serotype_", 1:n_types), levels = paste0("Serotype_", 1:n_types)), 
    "true_nu_j" = random_progression_rates
  )

progression_rates %>%
  kableExtra::kable()  %>%
  kableExtra::kable_styling(latex_options = "scale_down")
'''
# ------------------------- Probability of progressing to disease -------------------------------
set.seed(284)
random_disease_prob <- 10**(-1*runif(n_types, min = 0, max = 5))

disease_prob <-
  tibble(
    "serotype" = factor(paste0("Serotype_", 1:n_types), levels = paste0("Serotype_", 1:n_types)), 
    "true_psi_j" = random_disease_prob
  )

disease_prob %>%
  kableExtra::kable()  %>%
  kableExtra::kable_styling(latex_options = "scale_down")

# ------------------------- Carriage duration -------------------------------
set.seed(284)
library(truncnorm)
random_carriage_duration <- rtruncnorm(n_types, a=0, b=1, mean = 0.08, sd = 0.1)

carriage_duration <-
  tibble(
    "serotype" = factor(paste0("Serotype_", 1:n_types), levels = paste0("Serotype_", 1:n_types)), 
    "true_dur_j" = random_carriage_duration
  )

carriage_duration %>%
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

population_factor <-
  tibble(
    "study" = factor(paste0("Population_", LETTERS[1:n_studies])), 
    "true_r_i" = random_r_i, 
    "carriage_samples" = n_i,
    "surveillance_population" = N_i, 
    "time_interval" = 1
  )

population_factor %>%
  kableExtra::kable()  %>%
  kableExtra::kable_styling(latex_options = "scale_down")


# ------------------------- Simulated dataframe -------------------------------
simulation_df <- 
  carriage_prevalance %>%
  dplyr::left_join(disease_prob, by=c("serotype")) %>%
  dplyr::left_join(carriage_duration, by=c("serotype")) %>%
  dplyr::left_join(population_factor, by=c("study"))

# ------------------------- Serotype-specific_model -------------------------------
## generate simulated observed data
Serotype_specific_duration_poisson_simulation_df <-
  simulation_df %>%
  dplyr::rowwise() %>%
  dplyr::mutate(carriage = rbinom(1, carriage_samples, prob = true_rho_ij)) %>%
  dplyr::mutate(disease = rpois(1, surveillance_population*true_rho_ij*time_interval*(true_psi_j/true_dur_j)*true_r_i)) %>%
  dplyr::ungroup()



## making into a list for BIM 
Serotype_specific_duration_poisson_simulation_BIM_df <-
  progressionEstimation::process_input_data(Serotype_specific_duration_poisson_simulation_df, type = "serotype")


# ------------------------- Model Fitting -------------------------------
## general hyperparameters
n_chains <- 3
n_iter <- 2e4
n_core <- 2

#options(mc.cores = parallel::detectCores())

## study-adjusted serotype-specific Poisson simulated data

serotype_specific_poisson_pop_adjusted_duration_stan <- stan_model("/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/simulation/serotype_specific_study_adjusted_with_duration.stan")
study_adjusted_serotype_specific_poission_simulations_fit <- 
  rstan::sampling(serotype_specific_poisson_pop_adjusted_duration_stan, 
                  Serotype_specific_duration_poisson_simulation_BIM_df, 
                  iter=n_iter, 
                  chains=n_chains, 
                  cores=n_core)

launch_shinystan(study_adjusted_serotype_specific_poission_simulations_fit)

saveRDS(study_adjusted_serotype_specific_poission_simulations_fit, "study_adjusted_serotype_specific_duration_poission_simulations_fit.rds")
sim_duration_res<-readRDS("study_adjusted_serotype_specific_duration_poission_simulations_fit.rds")




input_df <- Serotype_specific_duration_poisson_simulation_df
i_levels = levels(input_df %>% dplyr::pull(study))
j_levels = levels(input_df %>% dplyr::pull(serotype))
model_output <- sim_duration_res

output_df <- input_df

carriage_df <- data.frame(
  #"study" = i_levels
  #"type" = j_levels
  "rho" = get_median("rho_ij",model_output),
  "rho_lower" = get_lower("rho_ij",model_output),
  "rho_upper" = get_upper("rho_ij",model_output)
)

output_df %<>% dplyr::bind_cols(carriage_df)


location_parameters <- data.frame(
  "study" = i_levels,
  "gamma" = get_median("gamma_i",model_output),
  "gamma_lower" = get_lower("gamma_i",model_output),
  "gamma_upper" = get_upper("gamma_i",model_output)
)

output_df %<>% dplyr::left_join(location_parameters, by = c("study"="study"))

disease_prob_df <- data.frame(
  "serotype" = j_levels,
  "psi" = as.numeric(get_median("psi_j",model_output)),
  "psi_lower" = as.numeric(get_lower("psi_j",model_output)),
  "psi_upper" = as.numeric(get_upper("psi_j",model_output))
)

output_df %<>% dplyr::left_join(disease_prob_df, by = setNames("serotype","serotype"))


duration_df <- data.frame(
  "serotype" = j_levels,
  "dur" = as.numeric(get_median("dur_j",model_output)),
  "dur_lower" = as.numeric(get_lower("dur_j",model_output)),
  "dur_upper" = as.numeric(get_upper("dur_j",model_output))
)

output_df %<>% dplyr::left_join(duration_df, by = setNames("serotype","serotype"))


output_diff_df <- output_df %>%
  dplyr::mutate(rho_diff = abs(true_rho_ij - rho)/true_rho_ij) %>%
  dplyr::mutate(psi_diff = abs(true_psi_j - psi)/true_psi_j) %>%
  dplyr::mutate(dur_diff = abs(true_dur_j - dur)/true_dur_j) %>%
  dplyr::mutate(r_diff = abs(true_r_i - gamma)/true_r_i) %>%
  dplyr::mutate(true_nu_j = true_psi_j/true_dur_j) %>%
  dplyr::mutate(nu = psi/dur) %>%
  dplyr::mutate(nu_diff = abs(true_nu_j - nu)/true_nu_j) 
  
  



write.table(output_diff_df, "Serotype_specific_duration_poisson_simulation_df.txt", sep="\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

prate_plot <-
  ggplot(output_diff_df,
         aes(x = true_nu_j,
             y = nu)) +
  geom_abline(slope = 1, intercept = 0, lty = 2, colour = "coral") +
  ylab("Predicted progression rate") +
  xlab("Observed progression rate") +
  theme_bw()

prate_plot <- prate_plot +
  geom_point(aes(color = study)) +
  #geom_errorbar(aes(color = study), alpha = 0.75)+
  theme(
    axis.title = element_text(size = 25), 
    axis.text = element_text(size = 20), 
    legend.text = element_text(size = 20), 
    legend.title = element_text(size = 25)
  )



png("Serotype_specific_duration_poisson_simulation_prate.png", width = 1000, height = 1000)
prate_plot
dev.off()

rho_plot <-
  ggplot(output_diff_df,
         aes(x = true_rho_ij,
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

png("Serotype_specific_duration_poisson_simulation_rho.png", width = 1000, height = 1000)
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

png("Serotype_specific_duration_poisson_simulation_gamma.png", width = 1000, height = 1000)
gamma_plot
dev.off()

psi_plot <-
  ggplot(output_diff_df,
         aes(x = true_psi_j,
             y = psi,
             ymin = psi_lower,
             ymax = psi_upper)) +
  geom_abline(slope = 1, intercept = 0, lty = 2, colour = "coral") +
  ylab("Predicted disease probability") +
  xlab("Observed disease probability") +
  theme_bw()

psi_plot <- psi_plot +
  geom_point(aes(color = serotype)) +
  geom_errorbar(aes(color = serotype), alpha = 0.75)+
  theme(
    axis.title = element_text(size = 25), 
    axis.text = element_text(size = 20), 
    legend.text = element_text(size = 20), 
    legend.title = element_text(size = 25)
  )

png("Serotype_specific_duration_poisson_simulation_psi.png", width = 1000, height = 1000)
psi_plot
dev.off()

dur_plot <-
  ggplot(output_diff_df,
         aes(x = true_dur_j,
             y = dur,
             ymin = dur_lower,
             ymax = dur_upper)) +
  geom_abline(slope = 1, intercept = 0, lty = 2, colour = "coral") +
  ylab("Predicted carriage duration") +
  xlab("Observed carriage duration") +
  theme_bw()

dur_plot <- dur_plot +
  geom_point(aes(color = serotype)) +
  geom_errorbar(aes(color = serotype), alpha = 0.75)+
  theme(
    axis.title = element_text(size = 25), 
    axis.text = element_text(size = 20), 
    legend.text = element_text(size = 20), 
    legend.title = element_text(size = 25)
  )

png("Serotype_specific_duration_poisson_simulation_dur.png", width = 1000, height = 1000)
dur_plot
dev.off()



