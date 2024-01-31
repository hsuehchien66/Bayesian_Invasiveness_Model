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
                                  BIM_colnames = c("study", "type", "carriage", "disease","carriage_samples", "surveillance_population", "time_interval", "strain", "pair_name")){
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
    sero = strsplit(as.character(pair_freq[i,1]),"_")[[1]][2]
    gpsc = strsplit(as.character(pair_freq[i,1]),"_")[[1]][length(strsplit(as.character(pair_freq[i,1]),"_")[[1]])]
    pop_pair_df[i, "study"] = pop
    pop_pair_df[i, "type"] = sero
    pop_pair_df[i, "strain"] = gpsc
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
  pop_pair_df$type = as.factor(pop_pair_df$type)
  pop_pair_df$strain = as.factor(pop_pair_df$strain)
  pair_freq$pop <- as.factor(pair_freq$pop)
  BIM_sero_gpsc_input = pop_pair_df[, c("study", "type", "carriage", "disease", "carriage_samples", "surveillance_population", "time_interval", "strain")]
  
  res_list = list("Pair_Freq_table" = pair_freq,
                  "BIM_Input_table" = BIM_sero_gpsc_input)
  return(res_list)
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
# ------------------------------------------------------
### Analysis
## Read data
{
  gene_path = "/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_input/gwas_res_lrtpvalue_sorted_unitigs_e04_pres_matched.rtab"
  variant_pres = read.table(gene_path, header = T, row.names = 1, sep = "\t", comment.char = "$", check.names = F)
  variant_pres = as.data.frame(t(variant_pres))
  variant_pres <- variant_pres %>%
    rownames_to_column(var = "RowName")
  colnames(variant_pres)[1] = "Lane_id"
  
  # top 1 unitig: GATTATAATGTTACACCGAATTTTGTAGACC
  top1unitig = variant_pres[,c("Lane_id","GATTATAATGTTACACCGAATTTTGTAGACC")]
  
  # load GPS dataset
  GPS_path = "/Users/hc14/Documents/PhD_project/PhD_datasets/GPS/GPS1_database_v3.3_selected_row.xlsx"
  GPS_dataset = multiplesheets(GPS_path)
  
  ## Merge Tables
  GPS_table1_Metadata_v3 = GPS_dataset$table1_Metadata_v3
  GPS_table2_QC_v3 = GPS_dataset$table2_QC_v3
  GPS_table3_analysis_v3 = GPS_dataset$table3_analysis_v3
  ## Use the "Public_name" to link 
  GPS_table3_Uni = GPS_table3_analysis_v3[which(GPS_table3_analysis_v3$Duplicate == "UNIQUE"), ]
  GPS_merge = merge(GPS_table3_Uni, GPS_table1_Metadata_v3, by="Public_name")
  GPS_Bayes_selec = GPS_merge[,c("Lane_id","In_Silico_serotype", "GPSC_PoPUNK2", "Country.x", "Month_collection", "Year_collection", "Manifest_type")]
}

GPS_Bayes_selec_candidate_variant_top1 <- merge(GPS_Bayes_selec, top1unitig, by="Lane_id")
serotype_profile <- as.data.frame(table(GPS_Bayes_selec_candidate_variant_top1$In_Silico_serotype))
gpsc_profile <- as.data.frame(table(GPS_Bayes_selec_candidate_variant_top1$GPSC_PoPUNK2))
country_profile <- as.data.frame(table(GPS_Bayes_selec_candidate_variant_top1$Country.x))

carriage_rate_path = "/Users/hc14/Documents/PhD_project/PhD_datasets/GPS/GPS_carriage_rate.xlsx"
carriage_rate_file = read_excel(carriage_rate_path)
carriage_rate_df = as.data.frame(carriage_rate_file)
### Generating input data ---------------------------------------------------------------
### Generating the dataframe with each row as a Population-Serotype-GPSC pair
#---------------------------------------------------------------
input_data_list <- metadata_to_BIM_imput(GPS_Bayes_selec_candidate_variant_top1, carriage_rate_df)
BIM_sero_gpsc_input <- input_data_list$BIM_Input_table

sero_variant_input_data_list <- metadata_to_BIM_imput(GPS_Bayes_selec_candidate_variant_top1, carriage_rate_df, 
                                                      pair_colname1 = "In_Silico_serotype", pair_colname2 = "GATTATAATGTTACACCGAATTTTGTAGACC", 
                                                      BIM_colnames = c("study", "type", "carriage", "disease","carriage_samples", "surveillance_population", "time_interval", "variant", "pair_name"))

### BIM columns function doesn't work
gpsc_variant_input_data_list <- metadata_to_BIM_imput(GPS_Bayes_selec_candidate_variant_top1, carriage_rate_df, 
                                                      pair_colname1 = "GPSC_PoPUNK2", pair_colname2 = "GATTATAATGTTACACCGAATTTTGTAGACC", 
                                                      BIM_colnames = c("study", "strain", "carriage", "disease","carriage_samples", "surveillance_population", "time_interval", "variant", "pair_name"))

BIM_gpsc_variant_data_list <- gpsc_variant_input_data_list$BIM_Input_table
colnames(BIM_gpsc_variant_data_list) <- c("study", "strain", "carriage", "disease","carriage_samples", "surveillance_population", "time_interval", "variant")
write.table(BIM_gpsc_variant_data_list, "/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_input/BIM_gpsc_variant_data_list.txt",
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

BIM_sero_variant_input <- sero_variant_input_data_list$BIM_Input_table  # type: serotype ; strain: variant
write.table(BIM_sero_gpsc_input, "/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_input/BIM_sero_gpsc_input.txt",
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
write.table(BIM_sero_variant_input, "/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_input/BIM_sero_variant_input.txt",
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

## run BIM -------------------------------------------------

### run serotype-based gpsc-adjustedPoisson Bayesian model
s_pneumoniae_sero_gpsc_data <- progressionEstimation::process_input_data(BIM_sero_gpsc_input, type = "type", use_strain = TRUE, combine_strain = FALSE, condense = FALSE)
s_pneumoniae_poisson_serobased_gpsc_adjust_fit <- progressionEstimation::fit_progression_rate_model(input_data = s_pneumoniae_sero_gpsc_data,
                                                                              type_specific = TRUE,
                                                                              location_adjustment = TRUE,
                                                                              stat_model = "poisson",
                                                                              strain_as_primary_type = FALSE,
                                                                              strain_as_secondary_type = TRUE,
                                                                              num_chains = 2,
                                                                              num_iter = 1e4)
### run gpsc-based serotype-adjusted Poisson Bayesian model
s_pneumoniae_gpsc_sero_data <- progressionEstimation::process_input_data(BIM_sero_gpsc_input, type = "type", use_strain = TRUE, combine_strain = FALSE, condense = FALSE)
s_pneumoniae_poisson_gpscbased_seroadjust_fit <- progressionEstimation::fit_progression_rate_model(input_data = s_pneumoniae_gpsc_sero_data,
                                                                                                    type_specific = TRUE,
                                                                                                    location_adjustment = TRUE,
                                                                                                    stat_model = "poisson",
                                                                                                    strain_as_primary_type = TRUE,
                                                                                                    strain_as_secondary_type = FALSE,
                                                                                                    num_chains = 2,
                                                                                                    num_iter = 1e4)
### run gpsc based Poisson Bayesian model
s_pneumoniae_gpsc_data <- progressionEstimation::process_input_data(BIM_sero_gpsc_input, type = "strain", use_strain = FALSE, combine_strain = FALSE, condense = FALSE)
s_pneumoniae_poisson_gpsc_fit <- progressionEstimation::fit_progression_rate_model(input_data = s_pneumoniae_gpsc_data,
                                                                                                    type_specific = TRUE,
                                                                                                    location_adjustment = TRUE,
                                                                                                    stat_model = "poisson",
                                                                                                    strain_as_primary_type = FALSE,
                                                                                                    strain_as_secondary_type = FALSE,
                                                                                                    num_chains = 2,
                                                                                                    num_iter = 1e4)

### run serotype based Poisson Bayesian model
s_pneumoniae_sero_data <- progressionEstimation::process_input_data(BIM_sero_gpsc_input, type = "type", use_strain = FALSE, combine_strain = FALSE, condense = FALSE)
s_pneumoniae_poisson_sero_fit <- progressionEstimation::fit_progression_rate_model(input_data = s_pneumoniae_sero_data,
                                                                                   type_specific = TRUE,
                                                                                   location_adjustment = TRUE,
                                                                                   stat_model = "poisson",
                                                                                   strain_as_primary_type = FALSE,
                                                                                   strain_as_secondary_type = FALSE,
                                                                                   num_chains = 2,
                                                                                   num_iter = 1e4)

### run serotype based variant adjusted Poisson Bayesian model
s_pneumoniae_serobased_variantadjusted_data <- progressionEstimation::process_input_data(BIM_sero_variant_input, type = "type", use_strain = TRUE, combine_strain = FALSE, condense = FALSE)
s_pneumoniae_poisson_serobased_variantadjusted_fit <- progressionEstimation::fit_progression_rate_model(input_data = s_pneumoniae_serobased_variantadjusted_data,
                                                                                   type_specific = TRUE,
                                                                                   location_adjustment = TRUE,
                                                                                   stat_model = "poisson",
                                                                                   strain_as_primary_type = FALSE,
                                                                                   strain_as_secondary_type = TRUE,
                                                                                   num_chains = 2,
                                                                                   num_iter = 1e4)

### run variant based serotype adjusted Poisson Bayesian model
s_pneumoniae_variantbased_serotypeadjusted_data <- progressionEstimation::process_input_data(BIM_sero_variant_input, type = "type", use_strain = TRUE, combine_strain = FALSE, condense = FALSE)
s_pneumoniae_poisson_variantbased_serotypeadjusted_fit <- progressionEstimation::fit_progression_rate_model(input_data = s_pneumoniae_variantbased_serotypeadjusted_data,
                                                                                                        type_specific = TRUE,
                                                                                                        location_adjustment = TRUE,
                                                                                                        stat_model = "poisson",
                                                                                                        strain_as_primary_type = TRUE,
                                                                                                        strain_as_secondary_type = FALSE,
                                                                                                        num_chains = 2,
                                                                                                        num_iter = 1e4)

### run sero based v dataframe Poisson Bayesian model
colnames(BIM_sero_variant_input) <- c("study", "serotype", "carriage", "disease", "carriage_samples", "surveillance_population", "time_interval", "variant")
s_pneumoniae_serobased_v_data <- progressionEstimation::process_input_data(BIM_sero_variant_input, type = "serotype", use_strain = FALSE, combine_strain = FALSE, condense = FALSE)
s_pneumoniae_poisson_serobased_v_fit <- progressionEstimation::fit_progression_rate_model(input_data = s_pneumoniae_serobased_v_data,
                                                                                           type_specific = TRUE,
                                                                                           location_adjustment = TRUE,
                                                                                           stat_model = "poisson",
                                                                                           strain_as_primary_type = FALSE,
                                                                                           strain_as_secondary_type = FALSE,
                                                                                           num_chains = 2,
                                                                                           num_iter = 1e4)


### run variant based Poisson Bayesian model
colnames(BIM_sero_variant_input) <- c("study", "serotype", "carriage", "disease", "carriage_samples", "surveillance_population", "time_interval", "type")
s_pneumoniae_variantbased_data <- progressionEstimation::process_input_data(BIM_sero_variant_input, type = "type", use_strain = FALSE, combine_strain = FALSE, condense = FALSE)
s_pneumoniae_poisson_variantbased_fit <- progressionEstimation::fit_progression_rate_model(input_data = s_pneumoniae_variantbased_data,
                                                                                                            type_specific = TRUE,
                                                                                                            location_adjustment = TRUE,
                                                                                                            stat_model = "poisson",
                                                                                                            strain_as_primary_type = FALSE,
                                                                                                            strain_as_secondary_type = FALSE,
                                                                                                            num_chains = 2,
                                                                                                            num_iter = 1e4)

### run variant based no-pop-adjusted Poisson Bayesian model
colnames(BIM_sero_variant_input) <- c("study", "serotype", "carriage", "disease", "carriage_samples", "surveillance_population", "time_interval", "type")
s_pneumoniae_variantbased_data <- progressionEstimation::process_input_data(BIM_sero_variant_input, type = "type", use_strain = FALSE, combine_strain = FALSE, condense = TRUE)
s_pneumoniae_poisson_variantbased_nopopadjusted_fit <- progressionEstimation::fit_progression_rate_model(input_data = s_pneumoniae_variantbased_data,
                                                                                           type_specific = TRUE,
                                                                                           location_adjustment = FALSE,
                                                                                           stat_model = "poisson",
                                                                                           strain_as_primary_type = FALSE,
                                                                                           strain_as_secondary_type = FALSE,
                                                                                           num_chains = 2,
                                                                                           num_iter = 1e4)
### run gpsc_based variant_adjusted Poisson Bayesian model
# strain ==> variant
colnames(BIM_gpsc_variant_data_list) <- c("study", "strain_type", "carriage", "disease", "carriage_samples", "surveillance_population", "time_interval", "strain")
s_pneumoniae_gpscbased_variantadjusted_data <- progressionEstimation::process_input_data(BIM_gpsc_variant_data_list, type = "strain_type", use_strain = TRUE, combine_strain = FALSE, condense = FALSE)
s_pneumoniae_poisson_gpscbased_variantadjusted_fit <- progressionEstimation::fit_progression_rate_model(input_data = s_pneumoniae_gpscbased_variantadjusted_data,
                                                                                                         type_specific = TRUE,
                                                                                                         location_adjustment = TRUE,
                                                                                                         stat_model = "poisson",
                                                                                                         strain_as_primary_type = FALSE,
                                                                                                         strain_as_secondary_type = TRUE,
                                                                                                         num_chains = 2,
                                                                                                         num_iter = 1e4)

### run gpsc_based v Poisson Bayesian model
# strain ==> variant
colnames(BIM_gpsc_variant_data_list) <- c("study", "strain_type", "carriage", "disease", "carriage_samples", "surveillance_population", "time_interval", "strain")
s_pneumoniae_gpscbased_v_data <- progressionEstimation::process_input_data(BIM_gpsc_variant_data_list, type = "strain_type", use_strain = FALSE, combine_strain = FALSE, condense = FALSE)
s_pneumoniae_poisson_gpscbased_v_fit <- progressionEstimation::fit_progression_rate_model(input_data = s_pneumoniae_gpscbased_v_data,
                                                                                                        type_specific = TRUE,
                                                                                                        location_adjustment = TRUE,
                                                                                                        stat_model = "poisson",
                                                                                                        strain_as_primary_type = FALSE,
                                                                                                        strain_as_secondary_type = FALSE,
                                                                                                        num_chains = 2,
                                                                                                        num_iter = 1e4)

### run gpsc_based v Poisson Bayesian model
# strain ==> variant
colnames(BIM_gpsc_variant_data_list) <- c("study", "strain_type", "carriage", "disease", "carriage_samples", "surveillance_population", "time_interval", "variant")
s_pneumoniae_variantbased_g_data <- progressionEstimation::process_input_data(BIM_gpsc_variant_data_list, type = "variant", use_strain = FALSE, combine_strain = FALSE, condense = FALSE)
s_pneumoniae_poisson_variantbased_g_fit <- progressionEstimation::fit_progression_rate_model(input_data = s_pneumoniae_variantbased_g_data,
                                                                                          type_specific = TRUE,
                                                                                          location_adjustment = TRUE,
                                                                                          stat_model = "poisson",
                                                                                          strain_as_primary_type = FALSE,
                                                                                          strain_as_secondary_type = FALSE,
                                                                                          num_chains = 2,
                                                                                          num_iter = 1e4)

save(s_pneumoniae_variantbased_g_data, s_pneumoniae_poisson_variantbased_g_fit, 
     file = "s_pneumoniae_variantbased_g.RData")

save(s_pneumoniae_gpscbased_v_data, s_pneumoniae_poisson_gpscbased_v_fit, 
     file = "s_pneumoniae_gpscbased_v.RData")

save(s_pneumoniae_gpscbased_variantadjusted_data, s_pneumoniae_poisson_gpscbased_variantadjusted_fit, 
     file = "s_pneumoniae_gpscbased_variantadjusted.RData")


save(s_pneumoniae_serobased_v_data, s_pneumoniae_poisson_serobased_v_fit, 
     file = "s_pneumoniae_serobased_v.RData")

save(s_pneumoniae_variantbased_data, s_pneumoniae_poisson_variantbased_nopopadjusted_fit, 
     file = "variantbased_nopopadjusted.RData")

save(s_pneumoniae_sero_gpsc_data, s_pneumoniae_poisson_serobased_gpsc_adjust_fit,
     s_pneumoniae_gpsc_sero_data, s_pneumoniae_poisson_gpscbased_seroadjust_fit,file = "serotype_gpsc_BIM.RData")

save(s_pneumoniae_gpsc_data, s_pneumoniae_poisson_gpsc_fit, file = "gpsc_BIM.RData")

save(s_pneumoniae_sero_data, s_pneumoniae_poisson_sero_fit,
     s_pneumoniae_gpsc_data, s_pneumoniae_poisson_gpsc_fit, file = "serotype_gpsc_only_BIM.RData")

save(s_pneumoniae_serobased_variantadjusted_data, s_pneumoniae_poisson_serobased_variantadjusted_fit,
     file = "serotype_gpsc_BIM.RData")

save(s_pneumoniae_variantbased_data, s_pneumoniae_poisson_variantbased_fit,
     file = "/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_output_rdata/variant_only_BIM.RData")
### model comparison ----------------------------

## Model comparison to avoid overfitting: LOO-CV & Bayes Factor

## LOO-CV uses likelihoods. It is not directly affected by prior distributions. (loo package in R)

## LOO-CV  
progressionEstimation::compare_model_fits_with_loo(list(s_pneumoniae_poisson_sero_fit, s_pneumoniae_poisson_gpsc_fit)) %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling(latex_options = "scale_down")
## The output is a table. Each row is a different mode. The rows are ordered with the best-fitting model on the top row. 
## ELPD: expected log pointwise predictive density, negative value indicates worse fits
## One model is regarded as outperforming another when "elpd_diff" is greater than 4, and larger than the corresponding "se_diff"

## LOO-CV doesn't work for different models with different number of parameters

## Bayes Factors are calculated using marginal likelihoods, which involves directly sampling from the prior. (bridgesampling package in Stan)
sero_vs_gpsc_comparison <- progressionEstimation::compare_model_fits_with_bf(list(s_pneumoniae_poisson_sero_fit,
                                                       s_pneumoniae_poisson_gpsc_fit)) %>%
  dplyr::rename("log(Bayes Factor)" = log_Bayes_factor) %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling(latex_options = "scale_down")

## post-analysis ---------------------------------

#variantbased_serotypeadjusted start 
rstan::traceplot(s_pneumoniae_poisson_variantbased_serotypeadjusted_fit,
                 pars= "nu_j")

s_pneumoniae_poisson_variantbased_serotypeadjusted_output_df <- progressionEstimation::process_progression_rate_model_output(s_pneumoniae_poisson_variantbased_serotypeadjusted_fit, 
                                                                                                                             BIM_sero_variant_input,
                                                                                                                             strain_as_secondary_type = TRUE)
write.table(s_pneumoniae_poisson_variantbased_serotypeadjusted_output_df, 
            file = "/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_output_results/variant_based_serotype_adjusted/s_pneumoniae_poisson_variantbased_serotypeadjusted_output_df.txt", 
            sep="\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

case_carrier_pred_variantbased_serotypeadjusted = progressionEstimation::plot_case_carrier_predictions(s_pneumoniae_poisson_variantbased_serotypeadjusted_output_df , n_label = 3)

variant_prate_variantbased_serotypeadjusted = progressionEstimation::plot_progression_rates(s_pneumoniae_poisson_variantbased_serotypeadjusted_output_df,
                                                                      type="strain",
                                                                      unit_time= "year",
                                                                      type_name= "Variant")
#variantbased_serotypeadjusted end

#variantbased start
s_pneumoniae_poisson_variantbased_output_df <- progressionEstimation::process_progression_rate_model_output(s_pneumoniae_poisson_variantbased_fit, 
                                                                                                                             BIM_sero_variant_input,
                                                                                                                             strain_as_secondary_type = FALSE, 
                                                                                                            condense = TRUE)
write.table(s_pneumoniae_poisson_variantbased_output_df, 
            file = "/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_output_results/variant_based/s_pneumoniae_poisson_variantbased_output_df.txt", 
            sep="\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

case_carrier_pred_variantbased = progressionEstimation::plot_case_carrier_predictions(s_pneumoniae_poisson_variantbased_output_df , n_label = 3)

variant_prate_variantbased = progressionEstimation::plot_progression_rates(s_pneumoniae_poisson_variantbased_output_df,
                                                                                            type="type",
                                                                                            unit_time= "year",
                                                                                            type_name= "Variant")
pdf(file = "/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_output_results/variant_based/case_carrier_pred_variantbased.pdf", width=16, height = 8)
case_carrier_pred_variantbased
dev.off()

pdf(file = "/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_output_results/variant_based/variant_prate_variantbased.pdf", width=16, height = 8)
variant_prate_variantbased
dev.off()
#variantbased end


# Check convergence
## All the parameters are estimated with rhat below the generally accepted threshold of 1.05.
## Most of the rhat values are close to 1.0
## It confirms convergence of the parameter estimates
poisson_nu_trace = rstan::traceplot(s_pneumoniae_poisson_fit,
                                    pars= "nu_j")

plot(s_pneumoniae_poisson_serobased_gpsc_adjust_fit, plotfun= "rhat", binwidth= 0.00005)










### supplementary ----------------------------------
### for microreact
{
  GPS_microreact = GPS_merge
  GPS_microreact$Manifest_type[is.na(GPS_microreact$Manifest_type)] = 0
  for(i in 1:nrow(GPS_merge)){
    if(GPS_microreact[i,"Manifest_type"]=="IPD"){
      GPS_microreact[i,"DiseaseCarrier"]="IPD"
    }else{
      GPS_microreact[i,"DiseaseCarrier"]="Carrier"
    }
  }
  write.csv(GPS_microreact, file="/Users/hc14/Documents/PhD_project/PhD_datasets/GPS/GPS1_microreact.csv", quote=F, row.names = FALSE)
}

## Descriptive Stat
{
  
  ### A. plot for number of sample in each population
  ggplot(GPS_merge, aes(x=reorder(Country.x, Country.x, function(x)-length(x))))+
    geom_bar(fill="red")+
    labs(x="Country")+
    theme(panel.background = element_blank(),text = element_text(family = "Helvetica"), axis.title = element_text(size=4), axis.text.x=element_text(size=4,angle=90,vjust=0.6), axis.text.y=element_text(size=4))
  ###
  
  
  ### B. plot for disease and carrier count in each population  
  GPS_country_carrier_disease_count = data.frame()
  GPS_country_manifest_ggplot_df = data.frame()
  num_country = 1
  for(ct in unique(GPS_merge$Country.x)){
    country = ct
    GPS_country = GPS_merge[which(GPS_merge$"Country.x" == ct),]
    #carriage_count = dim(GPS_country[which(GPS_country$Manifest_type != "IPD"),])[1]
    disease_count = dim(GPS_country[which(GPS_country$Manifest_type == "IPD"),])[1]
    carriage_count = dim(GPS_country)[1] - disease_count
    GPS_country_carrier_disease_count[num_country, "ID"] = paste0("$",num_country)
    GPS_country_carrier_disease_count[num_country, "Country"] = ct
    GPS_country_carrier_disease_count[num_country, "CarriageCount"] = carriage_count
    GPS_country_carrier_disease_count[num_country, "DiseaseCount"] = disease_count
    
    GPS_country_manifest_ggplot_df[c(num_country*2-1,num_country*2), "Country"] = ct
    GPS_country_manifest_ggplot_df[num_country*2-1, "Manifest"] = "Carrier"
    GPS_country_manifest_ggplot_df[num_country*2-1, "Count"] = carriage_count
    GPS_country_manifest_ggplot_df[num_country*2, "Manifest"] = "Disease"
    GPS_country_manifest_ggplot_df[num_country*2, "Count"] = disease_count
    
    num_country = num_country + 1
  }
  
  ggplot(GPS_country_manifest_ggplot_df, aes(fill=Manifest, y=Count, x=Country))+
    geom_bar(position="stack", stat="identity")+
    theme(panel.background = element_blank(),text = element_text(family = "Helvetica"), axis.title = element_text(size=8), axis.text.x=element_text(size=8,angle=90,vjust=0.6), axis.text.y=element_text(size=8))
  ###
}

### Generating input data ---------------------------------------------------------------
### Generating the dataframe with each row as a Population-Serotype pair
### Population-serotype-GPSC pair
{
  ### Population-serotype-GPSC pair
  for(isolate in 1:nrow(GPS_Bayes_selec_candidate_variant_top1)){
    GPS_Bayes_selec_candidate_variant_top1[isolate, "Pair"] = paste(GPS_Bayes_selec_candidate_variant_top1[isolate, "Country.x"],GPS_Bayes_selec_candidate_variant_top1[isolate, "In_Silico_serotype"], GPS_Bayes_selec_candidate_variant_top1[isolate, "GPSC_PoPUNK2"],sep="_")
  }
  
  pair_freq = as.data.frame(table(GPS_Bayes_selec_candidate_variant_top1$Pair))
  colnames(pair_freq) = c("name", "freq")
  for(i in 1:nrow(pair_freq)){
    pop = strsplit(as.character(pair_freq[i,1]),"_")[[1]][1]
    pair_freq$pop[i] = pop
  }
  pair_freq$pop = as.factor(pair_freq$pop)
}

### Population-serotype-variant pair
{
  ### Population-serotype-variant pair
  for(isolate in 1:nrow(GPS_Bayes_selec_candidate_variant_top1)){
    GPS_Bayes_selec_candidate_variant_top1[isolate, "sero_variant_Pair"] = paste(GPS_Bayes_selec_candidate_variant_top1[isolate, "Country.x"],GPS_Bayes_selec_candidate_variant_top1[isolate, "In_Silico_serotype"], GPS_Bayes_selec_candidate_variant_top1[isolate, "GATTATAATGTTACACCGAATTTTGTAGACC"],sep="_")
  }
  
  sero_variant_pair_freq = as.data.frame(table(GPS_Bayes_selec_candidate_variant_top1$sero_variant_Pair))
  colnames(sero_variant_pair_freq) = c("name", "freq")
  for(i in 1:nrow(sero_variant_pair_freq)){
    pop = strsplit(as.character(sero_variant_pair_freq[i,1]),"_")[[1]][1]
    sero_variant_pair_freq$pop[i] = pop
  }
  sero_variant_pair_freq$pop = as.factor(sero_variant_pair_freq$pop)
}


### generate input data for Bayesian Model 
### study, type, carriage, disease, time_interval can be directly extracted from GPSC table
### carriage and disease information are extracted from manifest in GPSC table
### carriage samples and surveillance population information are from literature reviews
### pair_name is used to create population-serotype pair, and is not designed in the input format so it should be deleted before feeding into the model
dsize = matrix(0, nrow=nrow(pair_freq), ncol=9)
pop_sero_gpsc_df = as.data.frame(dsize)
colnames(pop_sero_gpsc_df) = c("study", "type", "carriage", "disease","carriage_samples", "surveillance_population", "time_interval", "strain", "pair_name")

for(i in 1:nrow(pair_freq)){
  pop = strsplit(as.character(pair_freq[i,1]),"_")[[1]][1]
  sero = strsplit(as.character(pair_freq[i,1]),"_")[[1]][2]
  gpsc = strsplit(as.character(pair_freq[i,1]),"_")[[1]][3]
  pop_sero_gpsc_df[i, "study"] = pop
  pop_sero_gpsc_df[i, "type"] = sero
  pop_sero_gpsc_df[i, "strain"] = gpsc
  pop_sero_gpsc_df[i, "pair_name"] = as.character(pair_freq[i, "name"])
}

uniq_pair = unique(GPS_Bayes_selec_candidate_variant_top1$Pair)
## each row is a population-serotype-gpsc pair
for(p in uniq_pair){
  q_pop = GPS_Bayes_selec_candidate_variant_top1[which(GPS_Bayes_selec_candidate_variant_top1$Pair == p),]
  t_count = nrow(q_pop)
  q_disease = q_pop[which(q_pop$Manifest_type == "IPD"),]
  disease_count = nrow(q_disease)
  q_carriage = q_pop[which(q_pop$Manifest_type != "IPD"),]
  #carriage_count = nrow(q_carriage)
  carriage_count = t_count - disease_count
  pop_sero_gpsc_df[which(pop_sero_gpsc_df$pair_name == p), "carriage"] = carriage_count
  pop_sero_gpsc_df[which(pop_sero_gpsc_df$pair_name == p), "disease"] = disease_count
}

### we do not know the number of samples we collected in the GPSC table
### carriage rate is collected from literature reviews
### surveillance population is retrieved from UN stats
carriage_rate_path = "/Users/hc14/Documents/PhD_project/PhD_datasets/GPS/GPS_carriage_rate.xlsx"
carriage_rate_file = read_excel(carriage_rate_path)
carriage_rate_df = as.data.frame(carriage_rate_file)

pop_select_gps_df = pop_sero_gpsc_df[which(pop_sero_gpsc_df$study %in% carriage_rate_df$Country), ]

for(study in unique(pop_select_gps_df$study)){
  GPS_subset = GPS_Bayes_selec_candidate_variant_top1[which(GPS_Bayes_selec_candidate_variant_top1$Country.x == study), ]
  GPS_time = max(GPS_subset$Year_collection, na.rm = TRUE) - min(GPS_subset$Year_collection, na.rm = TRUE) + 1
  study_subset = pop_select_gps_df[which(pop_select_gps_df$study == study), ]
  study_sample = sum(study_subset$carriage) + sum(study_subset$disease)
  study_ref = carriage_rate_df[which(carriage_rate_df$Country == study),]
  study_carriage_rate = study_ref$Carriage_rate
  study_pop_size = study_ref$Pop_size
  pop_select_gps_df[which(pop_select_gps_df$study == study), "carriage_samples"] = round(study_sample/study_carriage_rate)
  pop_select_gps_df[which(pop_select_gps_df$study == study), "surveillance_population"] = round(study_pop_size * 1000000)
  pop_select_gps_df[which(pop_select_gps_df$study == study), "time_interval"] = GPS_time
}

### Indonesia has no time data and Netherland has less than five samples

remove_pop_list = c("INDONESIA")
bigpop_select_gps_df = pop_select_gps_df[which(!pop_select_gps_df$study %in% remove_pop_list),]

### make South Africa as reference population (make it as the first level)
bigpop_select_gps_df$study = as.factor(bigpop_select_gps_df$study)
bigpop_select_gps_df_study_level = levels(bigpop_select_gps_df$study)
bigpop_select_gps_df_study_level = append(bigpop_select_gps_df_study_level, "SOUTH AFRICA", after = 0)
bigpop_select_gps_df_study_level = bigpop_select_gps_df_study_level[-which(bigpop_select_gps_df_study_level == "SOUTH AFRICA")[2]] ## South Africa is duplicated, remove the second South Africa 
bigpop_select_gps_df$study = factor(bigpop_select_gps_df$study, levels = bigpop_select_gps_df_study_level)
### factorise character type column
bigpop_select_gps_df$type = as.factor(bigpop_select_gps_df$type)
bigpop_select_gps_df$strain = as.factor(bigpop_select_gps_df$strain)
BIM_sero_gpsc_input = bigpop_select_gps_df[,c("study", "type", "carriage", "disease", "carriage_samples", "surveillance_population", "time_interval", "strain")]



