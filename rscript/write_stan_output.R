---
title: "write stan model output"
author: "Raymond Cheng"
date: "11/03/2024"
output: model output table
---
library(readxl)
library(progressionEstimation)
library(imputeTS)
require(tidyverse)
require(magrittr)
require(rstan)
require(bridgesampling)
require(loo)
require(cowplot)
require(ggrepel)
require(xlsx)
require(roxygen2)
library(ggpubr) # for multiple ggplot
library(gridExtra) # for multiple ggplot
library(dplyr)  

get_mean<-function(parameter,model) {
    return(as.numeric(rstan::summary(model,pars=c(parameter))$summary[,1]))
  }

get_upper<-function(parameter,model) {
  return(as.numeric(rstan::summary(model,pars=c(parameter))$summary[,8]))
}

get_lower<-function(parameter,model) {
  return(as.numeric(rstan::summary(model,pars=c(parameter))$summary[,4]))
}

get_median<-function(parameter,model) {
  return(as.numeric(rstan::summary(model,pars=c(parameter))$summary[,6]))
}

### write output function
two_features_output <- function(input_df, model_output){
  input_df[, "feature1"] <- as.factor(input_df[, "feature1"])
  input_df[, "feature2"] <- as.factor(input_df[, "feature2"])
  
  i_levels = levels(input_df %>% dplyr::pull("study"))
  j_levels = levels(input_df %>% dplyr::pull("feature1"))
  g1_levels = levels(input_df %>% dplyr::pull("feature2"))
  
  output_df <- input_df
  carriage_df <- data.frame(
    #"study" = i_levels
    #"type" = j_levels
    "rho" = get_median("rho_ijg1",model_output),
    "rho_lower" = get_lower("rho_ijg1",model_output),
    "rho_upper" = get_upper("rho_ijg1",model_output)
  )
  
  output_df %<>% dplyr::bind_cols(carriage_df)
  
  
  location_parameters <- data.frame(
    "study" = i_levels,
    "gamma" = get_median("gamma_i",model_output),
    "gamma_lower" = get_lower("gamma_i",model_output),
    "gamma_upper" = get_upper("gamma_i",model_output)
  )
  
  output_df %<>% dplyr::left_join(location_parameters, by = c("study"="study"))
  
  feature1_progression_rate_df <- data.frame(
    "feature1" = j_levels,
    "nu_j" = as.numeric(get_median("nu_j",model_output)),
    "nu_j_lower" = as.numeric(get_lower("nu_j",model_output)),
    "nu_j_upper" = as.numeric(get_upper("nu_j",model_output))
  )
  
  output_df %<>% dplyr::left_join(feature1_progression_rate_df, by = setNames("feature1","feature1"))
  
  feature2_progression_rate_df <- data.frame(
    "feature2" = g1_levels,
    "nu_g1" = as.numeric(get_median("nu_g1",model_output)),
    "nu_g1_lower" = as.numeric(get_lower("nu_g1",model_output)),
    "nu_g1_upper" = as.numeric(get_upper("nu_g1",model_output))
  )
  
  output_df %<>% dplyr::left_join(feature2_progression_rate_df, by = setNames("feature2","feature2"))
  
  output_df %<>%
    dplyr::mutate(carriage_prediction = get_median("c_ijg1_pred", model_output)) %>%
    dplyr::mutate(carriage_prediction_lower = get_lower("c_ijg1_pred", model_output)) %>%
    dplyr::mutate(carriage_prediction_upper =  get_upper("c_ijg1_pred", model_output)) %>%
    dplyr::mutate(disease_prediction = get_median("d_ijg1_pred", model_output)) %>%
    dplyr::mutate(disease_prediction_lower = get_lower("d_ijg1_pred", model_output)) %>%
    dplyr::mutate(disease_prediction_upper =  get_upper("d_ijg1_pred", model_output))
  
  # Add in absolute deviation
  output_df %<>%
    dplyr::mutate(carriage_abs_dev = abs(carriage - carriage_prediction)) %>%
    dplyr::mutate(disease_abs_dev = abs(disease - disease_prediction))
  
  return(output_df)
}
  
'''
This is the script to transform BIM stan output to tab delimited output table, so that it is easy to check

Input: 
1. input table for stan model
2. stan model fit output
'''

### read input table for stan model
# feature1: serotype; feature2: iga variant
BIM_sero_iga_input <- read.table("/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_input/BIM_sero_iga_input.txt", sep = "\t", header = TRUE, stringsAsFactors = TRUE)
BIM_iga_sero_input <- BIM_sero_iga_input
colnames(BIM_iga_sero_input) <- c("study", "feature2", "carriage", "disease", "carriage_samples", "surveillance_population", "time_interval", "feature1")
# feature1: iga variant; feature2: serotype

### iga_based_serotype_adjusted
iga_based_serotype_adjusted_model_fit <- readRDS("/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_output_rdata/small_pop/iga_based_serotype_adjusted_model_fit.rds")
iga_based_serotype_adjusted_model_fit_output <- two_features_output(input_df = BIM_iga_sero_input,
                                                                    model_output = iga_based_serotype_adjusted_model_fit)
write.table(iga_based_serotype_adjusted_model_fit_output, "/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_output_results/iga_based_serotype_adjusted_model_fit_output.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

### serotype_based_iga_adjusted
serotype_based_iga_adjusted_model_fit <- readRDS("/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_output_rdata/small_pop/serotype_based_iga_adjusted_model_fit.rds")
serotype_based_iga_adjusted_model_fit_output <- two_features_output(input_df = BIM_sero_iga_input, 
                                                                    model_output = serotype_based_iga_adjusted_model_fit)

write.table(serotype_based_iga_adjusted_model_fit_output, "/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_output_results/serotype_based_iga_adjusted_model_fit_output.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


### iga_based
serotype_based_iga_adjusted_model_fit <- readRDS("/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_output_rdata/small_pop/serotype_based_iga_adjusted_model_fit.rds")
serotype_based_iga_adjusted_model_fit_output <- two_features_output(input_df = BIM_sero_iga_input, 
                                                                    model_output = serotype_based_iga_adjusted_model_fit)

write.table(serotype_based_iga_adjusted_model_fit_output, "/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_output_results/serotype_based_iga_adjusted_model_fit_output.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)












input_df$feature1 <- as.factor(input_df$feature1)
i_levels = levels(input_df %>% dplyr::pull("study"))
j_levels = levels(input_df %>% dplyr::pull("feature1"))
g1_levels = levels(input_df %>% dplyr::pull("feature2"))
model_output <- iga_based_serotype_adjusted_model_fit

output_df <- input_df

# Extract predictions and intervals
output_df %<>%
  dplyr::mutate(carriage_prediction = get_median("c_ijg1_pred", model_output)) %>%
  dplyr::mutate(carriage_prediction_lower = get_lower("c_ijg1_pred", model_output)) %>%
  dplyr::mutate(carriage_prediction_upper =  get_upper("c_ijg1_pred", model_output)) %>%
  dplyr::mutate(disease_prediction = get_median("d_ijg1_pred", model_output)) %>%
  dplyr::mutate(disease_prediction_lower = get_lower("d_ijg1_pred", model_output)) %>%
  dplyr::mutate(disease_prediction_upper =  get_upper("d_ijg1_pred", model_output))

# Add in absolute deviation
output_df %<>%
  dplyr::mutate(carriage_abs_dev = abs(carriage - carriage_prediction)) %>%
  dplyr::mutate(disease_abs_dev = abs(disease - disease_prediction))

write.table(output_df, "/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_output_results/iga_based_serotype_adjusted_model_output.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)



### ---------------
BIM_sero_gpsc_input <- read.table("/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_input/BIM_sero_gpsc_input.txt", sep = "\t", header = TRUE, stringsAsFactors = TRUE)
BIM_gpsc_sero_input <- BIM_sero_gpsc_input
colnames(BIM_gpsc_sero_input) <- c("study", "feature2", "carriage", "disease", "carriage_samples", "surveillance_population", "time_interval", "feature1")

gpsc_based_serotype_adjusted_model_fit <- readRDS("/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_output_rdata/small_pop/s_pneumoniae_gpsc_only_model_fit.rds")

input_df <- BIM_gpsc_sero_input
input_df$feature1 <- as.factor(input_df$feature1)

i_levels = levels(input_df %>% dplyr::pull("study"))
j_levels = levels(input_df %>% dplyr::pull("feature1"))
g1_levels = levels(input_df %>% dplyr::pull("feature2"))
model_output <- gpsc_based_serotype_adjusted_model_fit

output_df <- input_df

# Extract predictions and intervals
output_df %<>%
  dplyr::mutate(carriage_prediction = get_median("c_ij_pred", model_output)) %>%
  dplyr::mutate(carriage_prediction_lower = get_lower("c_ij_pred", model_output)) %>%
  dplyr::mutate(carriage_prediction_upper =  get_upper("c_ij_pred", model_output)) %>%
  dplyr::mutate(disease_prediction = get_median("d_ij_pred", model_output)) %>%
  dplyr::mutate(disease_prediction_lower = get_lower("d_ij_pred", model_output)) %>%
  dplyr::mutate(disease_prediction_upper =  get_upper("d_ij_pred", model_output))

# Add in absolute deviation
output_df %<>%
  dplyr::mutate(carriage_abs_dev = abs(carriage - carriage_prediction)) %>%
  dplyr::mutate(disease_abs_dev = abs(disease - disease_prediction))

write.table(output_df, "/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_output_results/gpsc_based_serotype_adjusted_model_output.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

### ---------------
BIM_sero_gpsc_input <- read.table("/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_input/BIM_sero_gpsc_input.txt", sep = "\t", header = TRUE, stringsAsFactors = TRUE)

gpsc_based_model_fit <- readRDS("/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_output_rdata/small_pop/s_pneumoniae_gpsc_only_model_fit.rds")
model_output <- gpsc_based_model_fit

input_df <- BIM_sero_gpsc_input
input_df$feature1 <- as.factor(input_df$feature1)

output_df <- input_df

# Extract predictions and intervals
output_df %<>%
  dplyr::mutate(carriage_prediction = get_median("c_ij_pred", model_output)) %>%
  dplyr::mutate(carriage_prediction_lower = get_lower("c_ij_pred", model_output)) %>%
  dplyr::mutate(carriage_prediction_upper =  get_upper("c_ij_pred", model_output)) %>%
  dplyr::mutate(disease_prediction = get_median("d_ij_pred", model_output)) %>%
  dplyr::mutate(disease_prediction_lower = get_lower("d_ij_pred", model_output)) %>%
  dplyr::mutate(disease_prediction_upper =  get_upper("d_ij_pred", model_output))

# Add in absolute deviation
output_df %<>%
  dplyr::mutate(carriage_abs_dev = abs(carriage - carriage_prediction)) %>%
  dplyr::mutate(disease_abs_dev = abs(disease - disease_prediction))

write.table(output_df, "/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_output_results/serotype_based_model_output.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


### ---------------
BIM_sero_iga_input <- read.table("/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_input/BIM_sero_iga_input.txt", sep = "\t", header = TRUE, stringsAsFactors = TRUE)

sero_based_iga_adjusted_model_fit <- readRDS("/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_output_rdata/small_pop/serotype_based_iga_adjusted_model_fit.rds")

input_df <- BIM_sero_iga_input
input_df$feature1 <- as.factor(input_df$feature2)
model_output <- sero_based_iga_adjusted_model_fit

output_df <- input_df

# Extract predictions and intervals
output_df %<>%
  dplyr::mutate(carriage_prediction = get_median("c_ijg1_pred", model_output)) %>%
  dplyr::mutate(carriage_prediction_lower = get_lower("c_ijg1_pred", model_output)) %>%
  dplyr::mutate(carriage_prediction_upper =  get_upper("c_ijg1_pred", model_output)) %>%
  dplyr::mutate(disease_prediction = get_median("d_ijg1_pred", model_output)) %>%
  dplyr::mutate(disease_prediction_lower = get_lower("d_ijg1_pred", model_output)) %>%
  dplyr::mutate(disease_prediction_upper =  get_upper("d_ijg1_pred", model_output))

# Add in absolute deviation
output_df %<>%
  dplyr::mutate(carriage_abs_dev = abs(carriage - carriage_prediction)) %>%
  dplyr::mutate(disease_abs_dev = abs(disease - disease_prediction))

write.table(output_df, "/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_output_results/serotype_based_iga_adjusted_model_output.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
