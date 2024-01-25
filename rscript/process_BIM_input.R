library(readxl)
library(progressionEstimation)
library("imputeTS")
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
# Input: 
#  (1)initial metadata (with Serotype, GPSC, Variant, Population, Disease status, Collection time)
#     GPS_Bayes_selec
#     GPS_Bayes_selec_candidate_variant_top1
#  (2)carriage rate
# Output: 
#  (1)pair frequency  
#  (2)input for Bayesian model

### Population-serotype-variant pair
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
    gpsc = strsplit(as.character(pair_freq[i,1]),"_")[[1]][3]
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

### -------------------------------------------

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

carriage_rate_path = "/Users/hc14/Documents/PhD_project/PhD_datasets/GPS/GPS_carriage_rate.xlsx"
carriage_rate_file = read_excel(carriage_rate_path)
carriage_rate_df = as.data.frame(carriage_rate_file)

#### ---------------------------------------------

input_data_list <- metadata_to_BIM_imput(GPS_Bayes_selec_candidate_variant_top1, carriage_rate_df)
BIM_sero_gpsc_input <- input_data_list$BIM_Input_table

save(BIM_sero_gpsc_input,
     s_pneumoniae_sero_gpsc_data, s_pneumoniae_poisson_serobased_gpsc_adjust_fit,
     s_pneumoniae_gpsc_sero_data, s_pneumoniae_poisson_gpscbased_seroadjust_fit,file = "serotype_gpsc_BIM.RData")

#### post analysis --------------------------------
load("/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_output_rdata/serotype_gpsc_BIM.RData")

s_pneumoniae_poisson_gpscbased_seroadjust_output_df <- progressionEstimation::process_progression_rate_model_output(s_pneumoniae_poisson_gpscbased_seroadjust_output_df, Bframe)

case_carrier_pred = progressionEstimation::plot_case_carrier_predictions(s_pneumoniae_poisson_output_df, n_label = 3)
serotype_time_distribution = progressionEstimation::plot_progression_rates(s_pneumoniae_poisson_output_df,
                                                                           unit_time= "year",
                                                                           type_name= "Serotype")

#### model comparison
## Bayes Factors are calculated using marginal likelihoods, which involves directly sampling from the prior. (bridgesampling package in Stan)
serobased_gpscadjusted_vs_gpscbased_seroadjusted_comparison <- progressionEstimation::compare_model_fits_with_bf(list(s_pneumoniae_poisson_serobased_gpsc_adjust_fit,
                                                                                  s_pneumoniae_poisson_gpscbased_seroadjust_fit)) %>%
  dplyr::rename("log(Bayes Factor)" = log_Bayes_factor) %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling(latex_options = "scale_down")









