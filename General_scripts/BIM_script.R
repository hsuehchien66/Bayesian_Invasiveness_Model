## load required packages
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
## generate input dataframe for Bayesian Invasiveness model
BayesInvInputDframe = function(studyn="SouthAfrica", type="In_Silico_serotype", time="Year_collection", manifest="Manifest_type", 
                               pop_size, Dframe){
  min_year = min(Dframe[,time])
  max_year = max(Dframe[,time])
  m = matrix(0, nrow=length(unique(Dframe[,type])), ncol=7)
  sero_count_df = as.data.frame(m)
  colnames(sero_count_df) = c("study", "type", "carriage", "disease","carriage_samples", "surveillance_population", "time_interval")
  sero_count_df$study = studyn
  sero_count_df$carriage_samples = nrow(Dframe)
  sero_count_df$surveillance_population = pop_size
  sero_count_df$time_interval = max_year - min_year
  
  num_sero = 0
  for (sero in unique(Dframe[,type])){
    num_sero = num_sero + 1
    specific_serotype_df = Dframe[which(Dframe[,type] == sero),]
    
    sero_count = nrow(specific_serotype_df)
    disease_count = nrow(specific_serotype_df[which(specific_serotype_df[,manifest] == "IPD"),])
    carriage_count = sero_count - disease_count
    
    sero_count_df[num_sero, "type"] = sero
    sero_count_df[num_sero, "carriage"] = carriage_count
    sero_count_df[num_sero, "disease"] = disease_count
  }
  sero_count_df$study = as.factor(sero_count_df$study)
  sero_count_df$type = as.factor(sero_count_df$type)
  
  return(sero_count_df)
} 
## Bayesian Invasiveness Model: Poisson (SeroAdjusted & PopAdjusted to be set for serotype-dependent and population/study dependent)
Bim_Poisson = function(SeroAdjusted=TRUE, PopAdjusted=FALSE, NumChains=2, NumIter=1e4, Bframe){
  ## Bframe is a data.frame with 7 features
  s_pneumoniae_data <- progressionEstimation::process_input_data(Bframe) ## convert into a list
  ## Fit the serotype-specific model
  
  # We can fit a model that assumes disease occurs as a Poisson process with a fixed progression rate (nu) per unit time, for each serotype j.
  # The unit time will be number of year
  # Hence, the expected number of disease isolates will be dependent only on 
  # 1. the carriage prevalence of serotype j (rho)
  # 2. the size of the population under surveillance for disease in the study (N)
  # 3. the time interval over which disease surveillance was conducted
  # 4. scaling factor accounts for variation in different populations (gamma). No scaling factor here since there is only one studied population. location_adjustment=FALSE
  s_pneumoniae_poisson_fit <- progressionEstimation::fit_progression_rate_model(s_pneumoniae_data,
                                                                                type_specific = SeroAdjusted, ## nu varies among serotypes
                                                                                location_adjustment = PopAdjusted,
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

# ------------------------------------------------------
### Analysis
## Read data
{
  GPS_path = "/Users/hc14/Documents/PhD_project/PhD_datasets/GPS/GPS1_database_v3.3_selected_row.xlsx"
  GPS_dataset = multiplesheets(GPS_path)
  
  ## Merge Tables
  GPS_table1_Metadata_v3 = GPS_dataset$table1_Metadata_v3
  GPS_table2_QC_v3 = GPS_dataset$table2_QC_v3
  GPS_table3_analysis_v3 = GPS_dataset$table3_analysis_v3
  ## Use the "Public_name" to link 
  GPS_table3_Uni = GPS_table3_analysis_v3[which(GPS_table3_analysis_v3$Duplicate == "UNIQUE"), ]
  GPS_merge = merge(GPS_table3_Uni, GPS_table1_Metadata_v3, by="Public_name")
  GPS_Bayes_selec = GPS_merge[,c("In_Silico_serotype", "GPSC_PoPUNK2", "Country.x", "Month_collection", "Year_collection", "Manifest_type")]
}


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
## Descriptive Stat
{
  ggplot(GPS_merge, aes(x=reorder(Country.x, Country.x, function(x)-length(x))))+
    geom_bar(fill="red")+
    labs(x="Country")+
    theme(panel.background = element_blank(),text = element_text(family = "Helvetica"), axis.title = element_text(size=4), axis.text.x=element_text(size=4,angle=90,vjust=0.6), axis.text.y=element_text(size=4))
}

## Select populations
{
  ## select South Africa
  GPS_Bayes_selec_SouthAfrica = GPS_Bayes_selec[which(GPS_Bayes_selec$"Country.x" == "SOUTH AFRICA"), ]
  serotype_freq = as.data.frame(table(GPS_Bayes_selec_SouthAfrica$In_Silico_serotype))
}
country_IPD_count = dim(GPS_Bayes_selec_SouthAfrica[which(GPS_Bayes_selec_SouthAfrica$Manifest_type == "IPD"),])[1] 
country_carriage_count = dim(GPS_Bayes_selec_SouthAfrica[which(GPS_Bayes_selec_SouthAfrica$Manifest_type != "IPD"),])[1]

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

## Generate Input Dataframe for Bayesian Invasiveness model (Serotype)
s_pneumoniae_GPS_SA_df = BayesInvInputDframe(pop_size = 57800000, Dframe = GPS_Bayes_selec_SouthAfrica)

## Run Bayesian Invasiveness Model
SA_sero = Bim_Poisson(Bframe = s_pneumoniae_GPS_SA_df)
SA_null = Bim_Poisson(SeroAdjusted=FALSE, Bframe = s_pneumoniae_GPS_SA_df)


# -------------------------------------------------------
### Model Comparison
## Model comparison to avoid overfitting: LOO-CV & Bayes Factor

## LOO-CV uses likelihoods. It is not directly affected by prior distributions. (loo package in R)
## Bayes Factors are calculated using marginal likelihoods, which involves directly sampling from the prior. (bridgesampling package in Stan)

## LOO-CV  
progressionEstimation::compare_model_fits_with_loo(list(SA_sero$Model_Fit,
                                                        SA_null$Model_Fit)) %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling(latex_options = "scale_down")
## The output is a table. Each row is a different mode. The rows are ordered with the best-fitting model on the top row. 
## ELPD: expected log pointwise predictive density, negative value indicates worse fits
## One model is regarded as outperforming another when "elpd_diff" is greater than 4, and larger than the corresponding "se_diff"

## Bayes Factor
progressionEstimation::compare_model_fits_with_bf(list(SA_sero$Model_Fit,
                                                       SA_null$Model_Fit)) %>%
  dplyr::rename("log(Bayes Factor)" = log_Bayes_factor) %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling(latex_options = "scale_down")







## extract high freq serotype
{
  '''

high_freq_serotype = as.character(serotype_freq[which(serotype_freq$Freq > 10),][,1])
high_freq_serotype_freq = serotype_freq[which(serotype_freq$Var1 %in% high_freq_serotype), ]
GPS_Bayes_selec_SouthAfrica_high_freq_serotype = GPS_Bayes_selec_SouthAfrica[which(GPS_Bayes_selec_SouthAfrica$In_Silico_serotype %in% high_freq_serotype), ]

'''
  }  

## Generate Bayesian input table
{
  min_year = min(GPS_Bayes_selec_SouthAfrica$Year_collection)
  max_year = max(GPS_Bayes_selec_SouthAfrica$Year_collection)
  m = matrix(0, nrow=length(unique(GPS_Bayes_selec_SouthAfrica$In_Silico_serotype)), ncol=7)
  sero_count_df = as.data.frame(m)
  colnames(sero_count_df) = c("study", "type", "carriage", "disease","carriage_samples", "surveillance_population", "time_interval")
  sero_count_df$study = "SouthAfrica"
  sero_count_df$carriage_samples = nrow(GPS_Bayes_selec_SouthAfrica)
  sero_count_df$surveillance_population = 57800000
  sero_count_df$time_interval = max_year - min_year
  
  num_sero = 0
  for (sero in unique(GPS_Bayes_selec_SouthAfrica$In_Silico_serotype)){
    num_sero = num_sero + 1
    specific_serotype_df = GPS_Bayes_selec_SouthAfrica[which(GPS_Bayes_selec_SouthAfrica$In_Silico_serotype == sero),]
    print(specific_serotype_df)
    
    sero_count = nrow(specific_serotype_df)
    disease_count = nrow(specific_serotype_df[which(specific_serotype_df$Manifest_type == "IPD"),])
    carriage_count = sero_count - disease_count
    
    sero_count_df[num_sero, "type"] = sero
    sero_count_df[num_sero, "carriage"] = carriage_count
    sero_count_df[num_sero, "disease"] = disease_count
  }
  
  sero_count_df$study = as.factor(sero_count_df$study)
  sero_count_df$type = as.factor(sero_count_df$type)
}

## Bayesian Model
{
  ## 22 samples, 7 features
  s_pneumoniae_GPS_data <- progressionEstimation::process_input_data(s_pneumoniae_GPS_df) ## convert into a list
  
  
  ## Fit the serotype-specific model
  
  ## We can fit a model that assumes disease occurs as a Poisson process with a fixed progression rate (nu) per unit time, for each serotype j.
  ## The unit time will be number of year
  ## Hence, the expected number of disease isolates will be dependent only on 
  ## 1. the carriage prevalence of serotype j (rho)
  ## 2. the size of the population under surveillance for disease in the study (N)
  ## 3. the time interval over which disease surveillance was conducted
  ## 4. scaling factor accounts for variation in different populations (gamma). No scaling factor here since there is only one studied population. location_adjustment=FALSE
  
  s_pneumoniae_GPS_serotype_poisson_fit <- 
    progressionEstimation::fit_progression_rate_model(s_pneumoniae_GPS_data,
                                                      type_specific = TRUE, ## nu varies among serotypes
                                                      location_adjustment = FALSE,
                                                      num_chains = 2,
                                                      num_iter = 1e4)
  
  ## We can check whether MCMC converged in the specified number of iterations.
  ## We should specify the parameter of interest
  print(s_pneumoniae_GPS_serotype_poisson_fit@model_pars)
  
  ## The parameter nu_j refers to the progression rate estimates for each serotype j.
  ## We can assess the convergence between the two MCMCs run above for these progression rates
  rstan::traceplot(s_pneumoniae_GPS_serotype_poisson_fit,
                   pars= "nu_j")
  ## Both chains of MCMC have converged on similar values, suggesting we have identified genuine variation between serotypes.
  ## This can be formally tested by estimating rhat values across MCMCs
  plot(s_pneumoniae_sweden_poisson_fit, plotfun= "rhat", binwidth= 0.00005)
  ## All the parameters are estimated with rhat below the generally accepted threshold of 1.05.
  ## Most of the rhat values are close to 1.0
  ## It confirms convergence of the parameter estimates
  
  
  ## Predict outputs
  {
    ## We can now combine the results of the model fit with the original data to enable the MCMCs outputs to be interpreted.
    s_pneumoniae_GPS_poisson_model_output_df <-
      progressionEstimation::process_progression_rate_model_output(s_pneumoniae_GPS_poisson_fit,
                                                                   sero_count_df)
    progressionEstimation::plot_case_carrier_predictions(s_pneumoniae_GPS_poisson_model_output_df,
                                                         n_label = 3)
    ## Predicted and Observed case counts match well
    progressionEstimation::plot_progression_rates(s_pneumoniae_GPS_poisson_model_output_df,
                                                  unit_time= "year",
                                                  type_name= "Serotype")
    
  }
  
  ## Fit the Null (Alternative) model (no variation on progression rate, nu)
  {
    s_pneumoniae_GPS_null_fit <-
      progressionEstimation::fit_progression_rate_model(s_pneumoniae_GPS_data,
                                                        type_specific = FALSE,
                                                        location_adjustment = FALSE,
                                                        num_chains = 2,
                                                        num_iter = 1e4)
    
    s_pneumoniae_GPS_null_model_output_df <-
      progressionEstimation::process_progression_rate_model_output(s_pneumoniae_GPS_null_fit,
                                                                   sero_count_df)
    
    progressionEstimation::plot_case_carrier_predictions(s_pneumoniae_GPS_null_model_output_df,
                                                         n_label= 3)
  }  
  
  
}
