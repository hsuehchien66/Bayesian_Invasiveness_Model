## Generating simulation data 
setwd("/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/simulation")
library(tibble)
library(progressionEstimation)
library(rstan)
library(dplyr)
require(tidyverse)
require(magrittr)

process_3feature_input_data <- function(input_df, main_feature = "feature1", feature2 = "feature2", feature3 = "feature3"){
  
  if (!(main_feature %in% colnames(input_df))) {
    stop("Type column not in input data")
  }
  
  # Convert to factors
  input_df %<>%
    dplyr::mutate(study = factor(study)) %>%
    dplyr::mutate((!!main_feature) := factor(!!! dplyr::syms(main_feature))) %>%
    dplyr::mutate((!!feature2) := factor(!!! dplyr::syms(feature2))) %>%
    dplyr::mutate((!!feature3) := factor(!!! dplyr::syms(feature3)))
  
  # Calculate input
  input_df <- as.data.frame(input_df)
  i_values <- as.integer(input_df$study)
  j_values <- as.integer(input_df[, which(colnames(input_df) == main_feature)])
  g1_values <- as.integer(input_df[, which(colnames(input_df) == feature2)])
  g2_values <- as.integer(input_df[, which(colnames(input_df) == feature3)])
  c_ij <- input_df$carriage
  d_ij <- input_df$disease
  n_i <- input_df$carriage_samples
  N_i <- input_df$surveillance_population
  t_i <- input_df$time_interval
  progression_rate_data <- list(
    i_max = max(i_values),
    j_max = max(j_values),
    g1_max = max(g1_values),
    g2_max = max(g2_values),
    n_obs = length(c_ij),
    i_values = i_values,
    j_values = j_values,
    g1_values = g1_values,
    g2_values = g2_values,
    c_ijg1g2 = c_ij,
    d_ijg1g2 = d_ij,
    n_i = n_i,
    N_i = N_i,
    t_i = t_i
  )
  return(progression_rate_data)
}

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
setwd("/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/simulation")
# ------------------------- Carriage Rate -------------------------------

'''
Prevalance / Carriage rate
All simulations use the same randomly-generated set of carriage prevalance (rho_ij) for different types different locations
'''

n_types <- 5 
n_studies <- 5
n_gene1 <- 2
n_gene2 <- 2

# Simulate rho_ijg1g2
set.seed(284)
random_prevalance <- runif(n_types*n_studies*n_gene1*n_gene2, min = 0, max = 0.05) # 0.05 because there are 20 types (n_types * prevalance cannot exceed 1)

# generate table 
carriage_prevalance <- 
  tibble(
    "serotype" = factor(rep(paste0("Serotype_", 1:n_types), n_studies*n_gene1*n_gene2), levels = paste0("Serotype_", 1:n_types)), 
    "study" = factor(rep(paste0("Population_", LETTERS[1:n_studies]), each = n_types*n_gene1*n_gene2)),
    "gene1" = factor(rep(paste0("Gene1_", LETTERS[1:n_gene1]), each = n_types*n_studies*n_gene2)), 
    "gene2" = factor(rep(paste0("Gene2_", LETTERS[1:n_gene2]), each = n_types*n_studies*n_gene1)),
    "true_rho_ijg1g2" = random_prevalance
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


# ------------------------- Gene1 relative progression rate -------------------------------
set.seed(284)
random_g1 <- 10**rnorm(n_gene1, mean = 0, sd=3) 

gene1_progression_rate <-
  tibble(
    "gene1" = factor(paste0("Gene1_", LETTERS[1:n_gene1])), 
    "true_nu_g1" = random_g1
  )

gene1_progression_rate %>%
  kableExtra::kable()  %>%
  kableExtra::kable_styling(latex_options = "scale_down")

# ------------------------- Gene2 relative progression rate -------------------------------
set.seed(285)
random_g2 <- 10**rnorm(n_gene2, mean = 0, sd=2) 

gene2_progression_rate <-
  tibble(
    "gene2" = factor(paste0("Gene2_", LETTERS[1:n_gene2])), 
    "true_nu_g2" = random_g2
  )

gene2_progression_rate %>%
  kableExtra::kable()  %>%
  kableExtra::kable_styling(latex_options = "scale_down")

# ------------------------- Simulated dataframe -------------------------------
simulation_df <- 
  carriage_prevalance %>%
  dplyr::left_join(serotype_progression_rates, by=c("serotype")) %>%
  dplyr::left_join(population_factor, by=c("study")) %>%
  dplyr::left_join(gene1_progression_rate, by=c("gene1")) %>%
  dplyr::left_join(gene2_progression_rate, by=c("gene2"))
  

# ------------------------- Serotype-based_gene1_gene2_adjusted_model -------------------------------
## generate simulated observed data
serotype_based_gene1_gene2_adjusted_poisson_simulation_df <-
  simulation_df %>%
  dplyr::rowwise() %>%
  dplyr::mutate(carriage = rbinom(1, carriage_samples, prob = true_rho_ijg1g2)) %>%
  dplyr::mutate(disease = rpois(1, surveillance_population*true_rho_ijg1g2*time_interval*true_nu_j*true_nu_g1*true_nu_g2*true_r_i)) %>%
  dplyr::ungroup()

## making into a list for BIM 
serotype_based_gene1_gene2_adjusted_poisson_simulation_BIM_list <-
  process_3feature_input_data(input_df = serotype_based_gene1_gene2_adjusted_poisson_simulation_df, 
                                main_feature = "serotype", feature2 = "gene1", feature3 = "gene2")


# ------------------------- Model Fitting -------------------------------
## general hyperparameters
n_chains <- 2
n_iter <- 1e4
n_core <- 2

#options(mc.cores = parallel::detectCores())

## Serotype-based_gene1_gene2_adjusted_model
serotype_based_gene1_gene2_adjusted_poisson_stan <- stan_model("/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM/stan/serotype_based_g1_g2_adjusted.stan")
serotype_based_gene1_gene2_adjusted_poisson_stan <- stan_model("serotype_based_g1_g2_adjusted.stan")

serotype_based_gene1_gene2_adjusted_poisson_simulation_fit <- 
    rstan::sampling(serotype_based_gene1_gene2_adjusted_poisson_stan, 
                    serotype_based_gene1_gene2_adjusted_poisson_simulation_BIM_list, 
                  iter=n_iter, 
                  chains=n_chains, 
                  cores=n_core)

parameter_list <- as.data.frame(serotype_based_gene1_gene2_adjusted_poisson_simulation_fit)

#launch_shinystan(study_adjusted_serotype_specific_poission_simulations_fit)

saveRDS(serotype_based_gene1_gene2_adjusted_poisson_simulation_fit, "serotype_based_gene1_gene2_adjusted_poisson_simulation_fit.rds")


gene_sim <- readRDS("serotype_based_gene1_gene2_adjusted_poisson_simulation_fit.rds")





input_df <- serotype_based_gene1_gene2_adjusted_poisson_simulation_df
i_levels = levels(input_df %>% dplyr::pull(study))
j_levels = levels(input_df %>% dplyr::pull(serotype))
g1_levels = levels(input_df %>% dplyr::pull(gene1))
g2_levels = levels(input_df %>% dplyr::pull(gene2))
model_output <- gene_sim

output_df <- input_df

carriage_df <- data.frame(
  #"study" = i_levels
  #"type" = j_levels
  "rho" = get_median("rho_ijg1g2",model_output),
  "rho_lower" = get_lower("rho_ijg1g2",model_output),
  "rho_upper" = get_upper("rho_ijg1g2",model_output)
)

output_df %<>% dplyr::bind_cols(carriage_df)


location_parameters <- data.frame(
  "study" = i_levels,
  "gamma" = get_median("gamma_i",model_output),
  "gamma_lower" = get_lower("gamma_i",model_output),
  "gamma_upper" = get_upper("gamma_i",model_output)
)

output_df %<>% dplyr::left_join(location_parameters, by = c("study"="study"))

serotype_progression_rates_df <- data.frame(
  "serotype" = j_levels,
  "nu_j" = as.numeric(get_median("nu_j",model_output)),
  "nu_j_lower" = as.numeric(get_lower("nu_j",model_output)),
  "nu_j_upper" = as.numeric(get_upper("nu_j",model_output))
)

output_df %<>% dplyr::left_join(serotype_progression_rates_df, by = setNames("serotype","serotype"))


gene1_progression_rates_df <- data.frame(
  "gene1" = g1_levels,
  "nu_g1" = as.numeric(get_median("nu_g1",model_output)),
  "nu_g1_lower" = as.numeric(get_lower("nu_g1",model_output)),
  "nu_g1_upper" = as.numeric(get_upper("nu_g1",model_output))
)

output_df %<>% dplyr::left_join(gene1_progression_rates_df, by = setNames("gene1","gene1"))

get_median<-function(parameter,model) {
  return(as.numeric(rstan::summary(model,pars=c(parameter))$summary[,"50%"]))
}

gene2_progression_rates_df <- data.frame(
  "gene2" = g2_levels,
  "nu_g2" = as.numeric(get_median("nu_g2",model_output)),
  "nu_g2_lower" = as.numeric(get_lower("nu_g2",model_output)),
  "nu_g2_upper" = as.numeric(get_upper("nu_g2",model_output))
)

output_df %<>% dplyr::left_join(gene2_progression_rates_df, by = setNames("gene2","gene2"))

output_diff_df <- output_df %>%
  dplyr::mutate(rho_diff = abs(true_rho_ijg1g2 - rho)/true_rho_ijg1g2) %>%
  dplyr::mutate(nu_j_diff = abs(true_nu_j - nu_j)/true_nu_j) %>%
  dplyr::mutate(nu_g1_diff = abs(true_nu_g1 - nu_g1)/true_nu_g1) %>%
  dplyr::mutate(nu_g2_diff = abs(true_nu_g2 - nu_g2)/true_nu_g2) %>%
  dplyr::mutate(r_diff = abs(true_r_i - gamma)/true_r_i) 




write.table(output_diff_df, "Serotype_specific_popadjusted_twogenes_poisson_simulation_df.txt", sep="\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

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



png("Serotype_specific_poisson_simulation_prate.png", width = 1000, height = 1000)
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

png("Serotype_specific_poisson_simulation_rho.png", width = 1000, height = 1000)
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

png("Serotype_specific_poisson_simulation_gamma.png", width = 1000, height = 1000)
gamma_plot
dev.off()
