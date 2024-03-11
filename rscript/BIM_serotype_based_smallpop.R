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
# ------------------------------------------------------
### FUNCTIONs
## read data from multiple sheets excel file
multiplesheets = function(fname){
  ## getting info about all excel sheets
  # excel_sheet(path): fetch all the worksheet names 
  sheets = readxl::excel_sheets(fname)
  ## use lapply method applies the read_excel method over every sheet of the workbook
  # read_excel(path, sheet)
  tibble = lapply(sheets, function(x) readxl::read_excel(fname, sheet=x))
  data_frame = lapply(tibble, as.data.frame)
  
  ## assigning names to dataframes
  names(data_frame) = sheets
  
  return(data_frame)
} 
## process metadata to BIM input format
metadata_to_BIM_imput <- function(metadata_df, carriage_rate_df,
                                  pop_colname="Country.x", pair_colname1="In_Silico_serotype", pair_colname2="GPSC_PoPUNK2", 
                                  BIM_colnames = c("study", "feature1", "carriage", "disease","carriage_samples", "surveillance_population", "time_interval", "feature2", "pair_name")){
  for(isolate in 1:nrow(metadata_df)){
    metadata_df[isolate, "pair"] = paste(metadata_df[isolate, pop_colname], 
                                         metadata_df[isolate, pair_colname1],
                                         metadata_df[isolate, pair_colname2], sep = "_")
  }
  pair_freq <- as.data.frame(table(metadata_df$pair))
  colnames(pair_freq) <- c("name", "freq")
  
  dsize <- matrix(0, nrow = nrow(pair_freq), ncol = 9)
  pop_pair_df = as.data.frame(dsize)
  colnames(pop_pair_df) <- BIM_colnames
  
  for (i in 1:nrow(pair_freq)) {
    pop = strsplit(as.character(pair_freq[i,1]), "_")[[1]][1]
    pair_freq$pop[i] = pop
    feature1 = strsplit(as.character(pair_freq[i,1]),"_")[[1]][2]
    feature2 = strsplit(as.character(pair_freq[i,1]),"_")[[1]][length(strsplit(as.character(pair_freq[i,1]),"_")[[1]])]
    pop_pair_df[i, "study"] = pop
    pop_pair_df[i, "feature1"] = feature1
    pop_pair_df[i, "feature2"] = feature2
    pop_pair_df[i, "pair_name"] = as.character(pair_freq[i, "name"])
    
  }
  uniq_pair = unique(metadata_df$pair)
  ### each row is a population-serotype-gpsc pair
  for (p in uniq_pair){
    q_pop = metadata_df[which(metadata_df$pair == p),]
    t_count = nrow(q_pop)
    q_disease = q_pop[which(q_pop$Manifest_type == "IPD"),]
    disease_count = nrow(q_disease)
    q_carriage = q_pop[which(q_pop$Manifest_type != "IPD"),]
    # carriage_count = nrow(q_carriage)
    carriage_count = t_count - disease_count
    pop_pair_df[which(pop_pair_df$pair_name == p), "carriage"] = carriage_count
    pop_pair_df[which(pop_pair_df$pair_name == p), "disease"] = disease_count
  }
  ### we do not know the number of samples we collected in the GPSC table
  ### carriage rate is collected from literature reviews
  ### surveillance population is retrieved from UN stats
  for (study in unique(pop_pair_df$study)){
    GPS_subset = metadata_df[which(metadata_df$Country.x == study), ]
    GPS_time = max(GPS_subset$Year_collection, na.rm = TRUE) - min(GPS_subset$Year_collection, na.rm = TRUE) +1
    study_subset = pop_pair_df[which(pop_pair_df$study == study),]
    study_sample = sum(study_subset$carriage) + sum(study_subset$disease)
    study_ref = carriage_rate_df[which(carriage_rate_df$Country == study),]
    study_carriage_rate = study_ref$Carriage_rate
    study_pop_size = study_ref$Pop_size
    pop_pair_df[which(pop_pair_df$study == study), "carriage_samples"] = round(study_sample/study_carriage_rate)
    pop_pair_df[which(pop_pair_df$study == study), "surveillance_population"] = round(study_pop_size * 1000000)
    pop_pair_df[which(pop_pair_df$study == study), "time_interval"] = GPS_time
  }
  ### make South Africa as reference population (make it as the first level)
  pop_pair_df$study = as.factor(pop_pair_df$study)
  pop_pair_df_study_level = levels(pop_pair_df$study)
  pop_pair_df_study_level = append(pop_pair_df_study_level, "SOUTH AFRICA", after = 0)
  pop_pair_df_study_level = pop_pair_df_study_level[-which(pop_pair_df_study_level == "SOUTH AFRICA")[2]] ## South Africa is duplicated, remove the second South Africa 
  pop_pair_df$study = factor(pop_pair_df$study, levels = pop_pair_df_study_level)
  # factorise character type column
  pop_pair_df$feature1 = as.factor(pop_pair_df$feature1)
  pop_pair_df$feature2 = as.factor(pop_pair_df$feature2)
  pair_freq$pop <- as.factor(pair_freq$pop)
  BIM_sero_variant_input = pop_pair_df[, c("study", "feature1", "carriage", "disease", "carriage_samples", "surveillance_population", "time_interval", "feature2")]
  
  res_list = list("Pair_Freq_table" = pair_freq,
                  "BIM_Input_table" = BIM_sero_variant_input)
  return(res_list)
}

## process input table to the input for BIM model
process_input_data <- function(input_df, main_feature = "feature1", use_feature2 = TRUE, combine_feature2 = FALSE, condense = FALSE) {
  if (!(main_feature %in% colnames(input_df))) {
    stop("Type column not in input data")
  }
  # Process input data
  if (combine_feature2 | use_feature2) { # make a dataframe with one more column of feature-pair
    input_df %<>%
      #input_df_feature2
      tidyr::unite("combined",!!main_feature,feature2, sep='_', remove = FALSE) %>%
      dplyr::mutate(combined = factor(combined))   ## add a column to combine the feature-pair
    if (condense) {
      input_df <- combine_rows(input_df, col_name = "combined")
      #input_df_feature2_condense <- combine_rows(input_df_feature2, col_name = "combined")
    }
    if (combine_feature2) { # change main feature to the feature-pair
      main_feature = "combined"
    } else if (use_feature2) { # the input dataframe remain the same as the original one
      input_df %<>% dplyr::select(-combined)
      #input_df_feature2  %<>% dplyr::select(-combined)
    }
  } else if (condense) { # aggregate the input dataframe to feature1 level
    input_df <- combine_rows(input_df %>% dplyr::select(study,
                                                        !!main_feature,
                                                        carriage,
                                                        disease,
                                                        carriage_samples,
                                                        surveillance_population,
                                                        time_interval),
                             col_name = main_feature)
  }
  # Convert to factors
  input_df %<>%
    dplyr::mutate(study = factor(study)) %>%
    dplyr::mutate((!!main_feature) := factor(!!! dplyr::syms(main_feature)))
  if ("feature2" %in% colnames(input_df)) {
    input_df %<>%
      dplyr::mutate(feature2 = factor(feature2))
  }
  # Calculate input
  input_df <- as.data.frame(input_df)
  i_values <- as.integer(input_df$study)
  j_values <- as.integer(input_df[, which(colnames(input_df) == main_feature)])
  c_ij <- input_df$carriage
  d_ij <- input_df$disease
  n_i <- input_df$carriage_samples
  N_i <- input_df$surveillance_population
  t_i <- input_df$time_interval
  if (use_feature2) {
    g1_values <- as.integer(input_df$feature2)
    progression_rate_data <- list(
      i_max = max(i_values),
      j_max = max(j_values),
      g1_max = max(g1_values),
      n_obs = length(c_ij),
      i_values = i_values,
      j_values = j_values,
      g1_values = g1_values,
      c_ijg1 = c_ij,
      d_ijg1 = d_ij,
      n_i = n_i,
      N_i = N_i,
      t_i = t_i
    )
  } else {
    progression_rate_data <- list(
      i_max = max(i_values),
      j_max = max(j_values),
      n_obs = length(c_ij),
      i_values = i_values,
      j_values = j_values,
      c_ij = c_ij,
      d_ij = d_ij,
      n_i = n_i,
      N_i = N_i,
      t_i = t_i
    )
  }
  return(progression_rate_data)
}


## Bayesian Invasiveness Model: Poisson (SeroAdjusted & PopAdjusted to be set for serotype-dependent and population/study dependent)
Bim_Poisson = function(SeroAdjusted=FALSE, PopAdjusted=FALSE, 
                       NumChains=2, NumIter=1e4, 
                       Model="poisson",
                       GPSCMajor=FALSE,
                       GPSCMinor=FALSE,
                       Bframe){
  ## Bframe is a data.frame with 8 features
  s_pneumoniae_data <- progressionEstimation::process_input_data(Bframe) ## convert into a list
  ## Fit the Baysian progression model
  
  # We can fit a model that assumes disease occurs as a Poisson process with a fixed progression rate (nu) per unit time, for each serotype j.
  # The unit time will be number of year
  # Hence, the expected number of disease isolates will be dependent only on 
  # 1. the carriage prevalence of serotype j (rho)
  # 2. the size of the population under surveillance for disease in the study (N)
  # 3. the time interval over which disease surveillance was conducted
  # 4. scaling factor accounts for variation in different populations (gamma). No scaling factor here since there is only one studied population. location_adjustment=FALSE
  s_pneumoniae_poisson_fit <- progressionEstimation::fit_progression_rate_model(input_data = s_pneumoniae_data,
                                                                                type_specific = SeroAdjusted, ## nu varies among serotypes
                                                                                location_adjustment = PopAdjusted,
                                                                                stat_model = Model,
                                                                                strain_as_primary_type = GPSCMajor,
                                                                                strain_as_secondary_type = GPSCMinor,
                                                                                num_chains = NumChains,
                                                                                num_iter = NumIter)  
  ## The parameter nu_j refers to the progression rate estimates for each serotype j.
  ## We can assess the convergence between the two MCMCs run above for these progression rates
  if(SeroAdjusted==TRUE){
    poisson_nu_trace = rstan::traceplot(s_pneumoniae_poisson_fit,
                                        pars= "nu_j")
  }else{
    poisson_nu_trace = rstan::traceplot(s_pneumoniae_poisson_fit,
                                        pars= "nu")
  }
  
  ## Both chains of MCMC have converged on similar values, suggesting we have identified genuine variation between serotypes.
  ## This can be formally tested by estimating rhat values across MCMCs
  poisson_rhat = plot(s_pneumoniae_poisson_fit, plotfun= "rhat", binwidth= 0.00005)
  
  ## We can now combine the results of the model fit with the original data to enable the MCMCs outputs to be interpreted.
  s_pneumoniae_poisson_output_df <- progressionEstimation::process_progression_rate_model_output(s_pneumoniae_poisson_fit, Bframe)
  case_carrier_pred = progressionEstimation::plot_case_carrier_predictions(s_pneumoniae_poisson_output_df, n_label = 3)
  serotype_time_distribution = progressionEstimation::plot_progression_rates(s_pneumoniae_poisson_output_df,
                                                                             unit_time= "year",
                                                                             type_name= "Serotype")
  
  resList = list("Model_Fit"=s_pneumoniae_poisson_fit, 
                 "Nu_Trace"= poisson_nu_trace, 
                 "Rhat"=poisson_rhat, 
                 "Model_Output"=s_pneumoniae_poisson_output_df,
                 "Case_Carrier_Prediction"=case_carrier_pred,
                 "Serotype_Distribution"=serotype_time_distribution)
}

plot_ordered_population_scale_factors <- function(model_output_df) {
  if (!("carriage_prediction" %in% colnames(model_output_df))) {
    stop("Need to include model output in data frame for plotting")
  }
  
  model_output_df %<>%
    dplyr::group_by(study) %>%
    dplyr::slice_head(n=1) %>%
    dplyr::ungroup()
  
  ggplot(model_output_df,
         aes(x = reorder(study, -gamma), y = gamma, ymin = gamma_lower, ymax = gamma_upper)) +
    geom_point() +
    geom_errorbar() +
    ylab(paste0("Population scale factor")) +
    xlab("Population") +
    scale_y_continuous(trans = "log10") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
}

plot_ordered_progression_rates <- function(model_output_df, type = "type", unit_time = "unit time", type_name = "type",
                                           colour_col = NULL, colour_palette = NULL, use_sample_size = FALSE) {
  if (!(any(grepl("nu",colnames(model_output_df))))) {
    stop("Need to include model output in data frame for plotting")
  }
  progression_rate_values <- c("nu", "nu_lower", "nu_upper")
  if (type == "strain" & "secondary_nu" %in% colnames(model_output_df)) {
    progression_rate_values <- paste0("secondary_", progression_rate_values)
  }
  y_label_text = paste0("Progression rate (disease per carrier per ",unit_time,")")
  if ("secondary_nu" %in% colnames(model_output_df)) {
    y_label_text = paste0("Progression rate contribution (disease per carrier per ",unit_time,")")
  }
  
  if (use_sample_size) {
    model_output_df %<>%
      dplyr::group_by(!!! dplyr::syms(type)) %>%
      dplyr::mutate(num_observations = sum(carriage+disease)) %>%
      dplyr::ungroup()
  }
  
  model_output_df %<>%
    dplyr::group_by(!!! dplyr::syms(type)) %>%
    dplyr::slice_head(n=1) %>%
    dplyr::ungroup()
  
  if (is.null(colour_col)) {
    if (use_sample_size) {
      base_graph <-
        ggplot(model_output_df,
               aes(x = reorder(get(!!type), -get(progression_rate_values[1])),
                   y = get(progression_rate_values[1]),
                   ymin = get(progression_rate_values[2]),
                   ymax = get(progression_rate_values[3]),
                   shape = num_observations))
    } else {
      base_graph <-
        ggplot(model_output_df,
               aes(x = reorder(get(!!type), -get(progression_rate_values[1])),
                   y = get(progression_rate_values[1]),
                   ymin = get(progression_rate_values[2]),
                   ymax = get(progression_rate_values[3])))
    }
  } else {
    if (use_sample_size) {
      base_graph <-
        ggplot(model_output_df,
               aes(x = reorder(get(!!type), -get(progression_rate_values[1])),
                   y = get(progression_rate_values[1]),
                   ymin = get(progression_rate_values[2]),
                   ymax = get(progression_rate_values[3]),
                   colour = get(colour_col),
                   fill = get(colour_col),
                   shape = num_observations))
    } else {
      base_graph <-
        ggplot(model_output_df,
               aes(x = reorder(get(!!type), -get(progression_rate_values[1])),
                   y = get(progression_rate_values[1]),
                   ymin = get(progression_rate_values[2]),
                   ymax = get(progression_rate_values[3]),
                   colour = get(colour_col),
                   fill = get(colour_col)))
    }
  }
  
  point_graph <-
    base_graph +
    geom_point() +
    geom_errorbar() +
    ylab(y_label_text) +
    xlab(type_name) +
    scale_y_continuous(trans ="log10") +
    theme_bw()
  
  if (!is.null(colour_col) & !is.null(colour_palette)) {
    point_graph <- point_graph +
      scale_colour_manual(values = colour_palette,
                          name = colour_col) +
      scale_fill_manual(values = colour_palette,
                        name = colour_col) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            legend.position = "bottom")
  } else if (!is.null(colour_col)) {
    point_graph <- point_graph +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            legend.position = "bottom")
  } else {
    point_graph <- point_graph +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  }
  
  if (use_sample_size) {
    point_graph <- point_graph +
      scale_shape_binned(name = "Number of\nisolates",
                         breaks = c(5,10,25,50,100))
  }
  
  return(point_graph)
}

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


# only include South Africa, Malawi, The Gambia, Nepal, Peru, India
# 12414 isolates, 6 countries, 5417 carriages, 6997 IPD

BIM_sero_gpsc_input <- read.table("BIM_sero_gpsc_input.txt", sep = "\t", header = TRUE)


## run BIM -------------------------------------------------
# feature1: serotype 
s_pneumoniae_sero_only_data <- process_input_data(BIM_sero_gpsc_input, main_feature = "feature1", 
                                                  use_feature2 = FALSE, combine_feature2 = FALSE, condense = FALSE)

one_feature_model_name = "adjusted_type_specific_poisson.stan"
one_feature_model_stan = stan_model(one_feature_model_name)
num_chains=3
num_iter=1e4
num_cores=parallel::detectCores()

s_pneumoniae_sero_only_model_fit <-rstan::sampling(one_feature_model_stan,
                                                   data = s_pneumoniae_sero_only_data,
                                                   iter = num_iter,
                                                   cores = num_cores,
                                                   chains = num_chains)
s_pneumoniae_sero_only_model_fit@model_name <- "serotype_based_model"
saveRDS(s_pneumoniae_sero_only_model_fit, "s_pneumoniae_sero_only_model_fit.rds")







