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
metadata_to_BIM_input_3feature <- function(metadata_df, carriage_rate_df,
                                           pop_colname="Country.x", pair_colname1="In_Silico_serotype", pair_colname2="GPSC_PoPUNK2", pair_colname3 = "GATTTTCATTGCCGTTATGCCAAGCATAGCA",
                                           BIM_colnames = c("study", "feature1", "carriage", "disease","carriage_samples", "surveillance_population", "time_interval", "feature2", "feature3", "pair_name")){
  
  metadata_df <- metadata_df %>%
    tidyr::unite("combined_name", c(all_of(pop_name), all_of(pair_colname1), all_of(pair_colname2), all_of(pair_colname3)), sep = "_", remove = FALSE) %>%
    dplyr::mutate(combined_name = factor(combined_name))
  
  pair_freq <- as.data.frame(table(metadata_df$combined_name))
  colnames(pair_freq) <- c("name", "freq")
  
  dsize <- matrix(0, nrow = nrow(pair_freq), ncol = length(BIM_colnames))
  pop_pair_df = as.data.frame(dsize)
  colnames(pop_pair_df) <- BIM_colnames
  
  for (i in 1:nrow(pair_freq)) {
    pop = strsplit(as.character(pair_freq[i,"name"]), "_")[[1]][1]
    pair_freq$pop[i] = pop
    feature1 = strsplit(as.character(pair_freq[i,"name"]),"_")[[1]][2]
    feature2 = strsplit(as.character(pair_freq[i,"name"]),"_")[[1]][3]
    feature3 = strsplit(as.character(pair_freq[i,"name"]),"_")[[1]][4]
    pop_pair_df[i, "study"] = pop
    pop_pair_df[i, "feature1"] = feature1
    pop_pair_df[i, "feature2"] = feature2
    pop_pair_df[i, "feature3"] = feature3
    pop_pair_df[i, "combined_name"] = as.character(pair_freq[i, "name"])
    
  }
  
  uniq_pair = unique(metadata_df$combined_name)
  ### each row is a population-serotype-gpsc pair
  for (p in uniq_pair){
    q_pop = metadata_df[which(metadata_df$combined_name == p), ]
    t_count = nrow(q_pop)
    q_disease = q_pop[which(q_pop$Manifest_type == "IPD"), ]
    disease_count = nrow(q_disease)
    q_carriage = q_pop[which(q_pop$Manifest_type == "Carriage"), ]
    # carriage_count = nrow(q_carriage)
    carriage_count = t_count - disease_count
    pop_pair_df[which(pop_pair_df$combined_name == p), "carriage"] = carriage_count
    pop_pair_df[which(pop_pair_df$combined_name == p), "disease"] = disease_count
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
  pop_pair_df$feature3 = as.factor(pop_pair_df$feature3)
  pair_freq$pop <- as.factor(pair_freq$pop)
  BIM_sero_variant_input = pop_pair_df[, c("study", "feature1", "carriage", "disease", "carriage_samples", "surveillance_population", "time_interval", "feature2", "feature3")]
  
  res_list = list("Pair_Freq_table" = pair_freq,
                  "BIM_Input_table" = BIM_sero_variant_input)
  
  return(res_list)
}

process_3feature_input_data <- function(input_df, main_feature = "feature1", feature2 = "feature2", feature3 = "feature3"){
  
  if (!(main_feature %in% colnames(input_df))) {
    stop("Type column not in input data")
  }
  
  # Convert to factors
  input_df %<>%
    dplyr::mutate(study = factor(study)) %>%
    dplyr::mutate((!!main_feature) := factor(!!! dplyr::syms(main_feature))) %>%
    dplyr::mutate(feature2 = factor(feature2)) %>%
    dplyr::mutate(feature3 = factor(feature3))
  
  # Calculate input
  input_df <- as.data.frame(input_df)
  i_values <- as.integer(input_df$study)
  j_values <- as.integer(input_df[, which(colnames(input_df) == main_feature)])
  g1_values <- as.integer(input_df$feature2)
  g2_values <- as.integer(input_df$feature3)
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

process_3feature_output <- function(model_output, 
                                    input_df, 
                                    main_feature = "feature1", 
                                    feature2 = "feature2", 
                                    feature3 = "feature3"){
  
  input_df %<>%
    dplyr::mutate(model_name = model_output@model_name)
  
  i_levels = levels(input_df %>% dplyr::pull(study))
  j_levels = levels(input_df %>% dplyr::pull(main_feature))
  g1_levels = levels(input_df %>% dplyr::pull(feature2))
  g2_levels = levels(input_df %>% dplyr::pull(feature3))
  
  # Carriage prevalence estimates
  carriage_df <- data.frame(
    "rho" = get_median("rho_ijg1g2", model_output),
    "rho_lower" = get_lower("rho_ijg1g2", model_output),
    "rho_upper" = get_upper("rho_ijg1g2", model_output)
  )
  
  input_df %<>% dplyr::bind_cols(carriage_df)
  
  # Variation by location
  scale_parameter <- 1
  if ("gamma_i" %in%  model_output@model_pars) {
    location_parameters <- data.frame(
      "study" = i_levels,
      "gamma" = get_median("gamma_i", model_output),
      "gamma_lower" = get_lower("gamma_i", model_output),
      "gamma_upper" = get_upper("gamma_i", model_output)
    )
  } else {
    location_parameters <- data.frame(
      "study" = i_levels,
      "gamma" = 1,
      "gamma_lower" = 1,
      "gamma_upper" = 1
    )
  }
  input_df %<>% dplyr::left_join(location_parameters, by = c("study"="study"))
  
  # Calculate serotype invasiveness values
  nu_name = "nu"
  if ("nu_j" %in%  model_output@model_pars) {
    nu_name = "nu_j"
  }
  
  progression_rates_df <- data.frame(
    "serotype" = j_levels,
    "nu" = as.numeric(get_median(nu_name, model_output)),
    "nu_lower" = as.numeric(get_lower(nu_name, model_output)),
    "nu_upper" = as.numeric(get_upper(nu_name, model_output))
  )
  input_df %<>% dplyr::left_join(progression_rates_df, by = setNames("serotype", main_feature))
  
  ## gene1 relative progression rate
  if ("nu_g1" %in%  model_output@model_pars) {
    gene1_progression_rates_df <- data.frame(
      "gene1" = g1_levels,
      "gene1_nu" = get_median("nu_g1", model_output),
      "gene1_nu_lower" = get_lower("nu_g1", model_output),
      "gene1_nu_upper" = get_upper("nu_g1", model_output)
    )
    input_df %<>% dplyr::left_join(gene1_progression_rates_df, by = setNames("gene1", feature2))
  }
  
  ## gene2 relative progression rate
  if ("nu_g2" %in%  model_output@model_pars) {
    gene2_progression_rates_df <- data.frame(
      "gene2" = g2_levels,
      "gene2_nu" = get_median("nu_g2", model_output),
      "gene2_nu_lower" = get_lower("nu_g2", model_output),
      "gene2_nu_upper" = get_upper("nu_g2", model_output)
    )
    input_df %<>% dplyr::left_join(gene2_progression_rates_df, by = setNames("gene2", feature3))
  }
  
  
  # Extract predictions and intervals
  input_df %<>%
    dplyr::mutate(carriage_prediction = get_median("c_ijg1g2_pred", model_output)) %>%
    dplyr::mutate(carriage_prediction_lower = get_lower("c_ijg1g2_pred", model_output)) %>%
    dplyr::mutate(carriage_prediction_upper =  get_upper("c_ijg1g2_pred", model_output)) %>%
    dplyr::mutate(disease_prediction = get_median("d_ijg1g2_pred", model_output)) %>%
    dplyr::mutate(disease_prediction_lower = get_lower("d_ijg1g2_pred", model_output)) %>%
    dplyr::mutate(disease_prediction_upper =  get_upper("d_ijg1g2_pred", model_output))
  
  # Add in absolute deviation
  input_df %<>%
    dplyr::mutate(carriage_abs_dev = abs(carriage - carriage_prediction)) %>%
    dplyr::mutate(disease_abs_dev = abs(disease - disease_prediction))
  
  return(input_df)
}

######

gene_path = "/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_input/gwas_res_lrtpvalue_sorted_unitigs_e04_pres_matched.rtab"
variant_pres = read.table(gene_path, header = T, row.names = 1, sep = "\t", comment.char = "$", check.names = F)
variant_pres = as.data.frame(t(variant_pres))
variant_pres <- variant_pres %>%
  rownames_to_column(var = "RowName")
  colnames(variant_pres)[1] = "Lane_id"

# iga unitig: GATTATAATGTTACACCGAATTTTGTAGACC
iga_kmers <- "GATTATAATGTTACACCGAATTTTGTAGACC"
truA_kmers <- "GATTTTCATTGCCGTTATGCCAAGCATAGCA"
iga_unitig = variant_pres[,c("Lane_id", iga_kmers)]
truA_unitig <- variant_pres[, c("Lane_id", truA_kmers)]

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

GPS_Bayes_selec_candidate_variant_top1 <- merge(GPS_Bayes_selec, iga_unitig, by="Lane_id")
GPS_Bayes_selec_candidate_variants <- merge(GPS_Bayes_selec_candidate_variant_top1, truA_unitig, by="Lane_id")

serotype_profile <- as.data.frame(table(GPS_Bayes_selec_candidate_variants$In_Silico_serotype))
gpsc_profile <- as.data.frame(table(GPS_Bayes_selec_candidate_variants$GPSC_PoPUNK2))
country_profile <- as.data.frame(table(GPS_Bayes_selec_candidate_variants$Country.x))

carriage_rate_path = "/Users/hc14/Documents/PhD_project/PhD_datasets/GPS/GPS_carriage_rate.xlsx"
carriage_rate_file = read_excel(carriage_rate_path)
carriage_rate_df = as.data.frame(carriage_rate_file)

### Select suitable countries
## Country with population factor wihtin 10 scale: Slovenia, South Africa, Argentina (all disease), Brazil, USA (All disease), Malawi, The Gambia, Mozambique, PNG (all disease), Poland, Peru
## Country with realtively even counts of carrier and disease:
## South Africa, Malawi, The Gambia, Nepal, Peru, India

SA_GPS <- GPS_Bayes_selec_candidate_variants[which(GPS_Bayes_selec_candidate_variants$Country.x == "SOUTH AFRICA"), ]
SA_GPS <- SA_GPS[which(SA_GPS$Manifest_type == "Carriage" | SA_GPS$Manifest_type == "IPD"), ]

Malawi_GPS <- GPS_Bayes_selec_candidate_variants[which(GPS_Bayes_selec_candidate_variants$Country.x == "MALAWI"), ]
Malawi_GPS <- Malawi_GPS[which(Malawi_GPS$Manifest_type == "Carriage" | Malawi_GPS$Manifest_type == "IPD"), ]

Gambia_GPS <- GPS_Bayes_selec_candidate_variants[which(GPS_Bayes_selec_candidate_variants$Country.x == "THE GAMBIA"), ]
Gambia_GPS <- Gambia_GPS[which(Gambia_GPS$Manifest_type == "Carriage" | Gambia_GPS$Manifest_type == "IPD"), ]

Peru_GPS <- GPS_Bayes_selec_candidate_variants[which(GPS_Bayes_selec_candidate_variants$Country.x == "PERU"), ]
Peru_GPS <- Peru_GPS[which(Peru_GPS$Manifest_type == "Carriage" | Peru_GPS$Manifest_type == "IPD"), ]

Nepal_GPS <- GPS_Bayes_selec_candidate_variants[which(GPS_Bayes_selec_candidate_variants$Country.x == "NEPAL"), ]
Nepal_GPS <- Nepal_GPS[which(Nepal_GPS$Manifest_type == "Carriage" | Nepal_GPS$Manifest_type == "IPD"), ]

India_GPS <- GPS_Bayes_selec_candidate_variants[which(GPS_Bayes_selec_candidate_variants$Country.x == "INDIA"), ]
India_GPS <- India_GPS[which(India_GPS$Manifest_type == "Carriage" | India_GPS$Manifest_type == "IPD"), ]

curated_country <- c("SOUTH AFRICA", "MALAWI", "THE GAMBIA", "PERU", "NEPAL", "INDIA")
## 12414 isolates, 6 countries, 5417 carriages, 6997 IPD
GPS_curated_table <- rbind(SA_GPS, Malawi_GPS, Gambia_GPS, Peru_GPS, Nepal_GPS, India_GPS)
GPS_curated_table[which(GPS_curated_table$In_Silico_serotype == "SWISS_NT"), "In_Silico_serotype"] <- "SWISS"

GPS_curated_carriagerate <- carriage_rate_df[which(carriage_rate_df$Country %in% curated_country), ]

serotype_iga_truA_input_data_list <- metadata_to_BIM_input_3feature(metadata_df = GPS_curated_table, 
                                                                    carriage_rate_df = GPS_curated_carriagerate, 
                                                                    pair_colname1 = "In_Silico_serotype", 
                                                                    pair_colname2 = iga_kmers,
                                                                    pair_colname3 = truA_kmers)

BIM_sero_iga_truA_input <- serotype_iga_truA_input_data_list$BIM_Input_table

## run BIM ----------------
# feature1: serotype
# feature2: iga
# feature3: truA

BIM_sero_iga_truA_data <- process_3feature_input_data(input_df = BIM_sero_iga_truA_input, 
                                                      main_feature = "feature1", 
                                                      feature2 = "feature2", 
                                                      feature3 = "feature3")

three_feature_model_name = "/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM/stan/serotype_based_g1_g2_adjusted.stan"
three_feature_model_stan = stan_model(three_feature_model_name)
num_chains=2
num_iter=1e4
num_cores=parallel::detectCores()

serotype_based_iga_truA_adjusted_model_fit<-rstan::sampling(three_feature_model_stan,
                                                        data = BIM_sero_iga_truA_data,
                                                        iter = num_iter,
                                                        cores = num_cores,
                                                        chains = num_chains)
serotype_based_iga_truA_adjusted_model_fit@model_name <- "serotype_based_iga_truA_adjusted_model_fit"
saveRDS(serotype_based_iga_truA_adjusted_model_fit, "/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_output_rdata/serotype_based_iga_truA_adjusted_model_fit.rds")

serotype_based_iga_truA_adjusted_model_output <- process_3feature_output(model_output = serotype_based_iga_truA_adjusted_model_fit, 
                                                                         input_df = BIM_sero_iga_truA_input, 
                                                                         main_feature = "feature1", 
                                                                         feature2 = "feature2", 
                                                                         feature3 = "feature3")
write.table(serotype_based_iga_truA_adjusted_model_output, "/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_output_results/serotype_based_iga_truA_adjusted_model_output.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)




