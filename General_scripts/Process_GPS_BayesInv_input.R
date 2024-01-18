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
library(ggpubr) # for multiple ggplot
library(gridExtra) # for multiple ggplot
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
Bim_Poisson = function(SeroAdjusted=TRUE, PopAdjusted=FALSE, 
                       NumChains=2, NumIter=1e4, 
                       Model="poisson",
                       GPSCMajor=FALSE,
                       GPSCMinor=FALSE,
                       Bframe){
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

## Select populations
{
## select South Africa
#GPS_Bayes_selec_SouthAfrica = GPS_Bayes_selec[which(GPS_Bayes_selec$"Country.x" == "SOUTH AFRICA"), ]
#country_freq = as.data.frame(table(GPS_Bayes_selec$Country.x))
#serotype_freq = as.data.frame(table(GPS_Bayes_selec_SouthAfrica$In_Silico_serotype))
}

### Generating input data
### Generating the dataframe with each row as a Population-Serotype pair
#---------------------------------------------------------------
{
### Population-serotype pair
for(country_sero in 1:nrow(GPS_Bayes_selec)){
  GPS_Bayes_selec[country_sero, "Pair"] = paste(GPS_Bayes_selec[country_sero, "Country.x"], GPS_Bayes_selec[country_sero, "In_Silico_serotype"], sep="_")
}
  
pair_freq = as.data.frame(table(GPS_Bayes_selec$Pair))
colnames(pair_freq) = c("name", "freq")
for(i in 1:nrow(pair_freq)){
  pop = strsplit(as.character(pair_freq[i,1]),"_")[[1]][1]
  pair_freq$pop[i] = pop
}
pair_freq$pop = as.factor(pair_freq$pop)

### Population-GPSC pair
for(country_gpsc in 1:nrow(GPS_Bayes_selec)){
  GPS_Bayes_selec[country_gpsc, "Pair"] = paste(GPS_Bayes_selec[country_gpsc, "Country.x"], GPS_Bayes_selec[country_gpsc, "GPSC_PoPUNK2"], sep="_")
}

pair_freq = as.data.frame(table(GPS_Bayes_selec$Pair))
colnames(pair_freq) = c("name", "freq")
for(i in 1:nrow(pair_freq)){
  pop = strsplit(as.character(pair_freq[i,1]),"_")[[1]][1]
  pair_freq$pop[i] = pop
}
pair_freq$pop = as.factor(pair_freq$pop)

### plot population_serotype frequency for all populations in grid
p <- list()
for(tpopi in 1:length(unique(pair_freq$pop))){
  tpop = unique(pair_freq$pop)[tpopi]
  tpop_df = pair_freq[which(pair_freq$pop == tpop),]
  p[[tpopi]] <- ggplot(data=tpop_df, aes(x=name, y=freq))+
    geom_bar(stat = "identity")+
    theme(axis.text.x = element_blank(), axis.ticks = element_blank(), 
          axis.title = element_blank(), plot.title = element_text(size=6))+
    ggtitle(tpop)
    #labs(x = "Population-Serotype Pair", y="Frequency", fill="Population")
}
population_serotype_frequency_plot = do.call(grid.arrange, p) ##par(mfrow) is not applicable to ggplot

freq_barplot <- ggplot(data=pair_freq, aes(x=name, y=freq, fill=pop))+
  geom_bar(stat = "identity")+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())+
  labs(x = "Population-Serotype Pair", y="Frequency", fill="Population")

### plot population_serotype frequency for all populations in a histogram
freq_barplot
ggsave("/Users/hc14/Documents/PhD_project/Stan_Bayesian/results/Popultion-Serotype_Freq.pdf", width=20, height=10)
}
#---------------------------------------------------------------

### generate input data for Bayesian Model
### study, type, carriage, disease, time_interval can be directly extracted from GPSC table
### carriage and disease information are extracted from manifest in GPSC table
### carriage samples and surveillance population information are from literature reviews
### pair_name is used to create population-serotype pair, and is not designed in the input format so it should be deleted before feeding into the model
dsize = matrix(0, nrow=nrow(pair_freq), ncol=8)
pop_sero_df = as.data.frame(dsize)
colnames(pop_sero_df) = c("study", "type", "carriage", "disease","carriage_samples", "surveillance_population", "time_interval", "pair_name")

for(i in 1:nrow(pair_freq)){
  pop = strsplit(as.character(pair_freq[i,1]),"_")[[1]][1]
  sero = strsplit(as.character(pair_freq[i,1]),"_")[[1]][2]
  pop_sero_df[i, "study"] = pop
  pop_sero_df[i, "type"] = sero
  pop_sero_df[i, "pair_name"] = as.character(pair_freq[i, "name"])
}

uniq_pair = unique(GPS_Bayes_selec$Pair)
## each row is a population-serotype pair
for(p in uniq_pair){
  q_pop = GPS_Bayes_selec[which(GPS_Bayes_selec$Pair == p),]
  t_count = nrow(q_pop)
  q_disease = q_pop[which(q_pop$Manifest_type == "IPD"),]
  disease_count = nrow(q_disease)
  q_carriage = q_pop[which(q_pop$Manifest_type != "IPD"),]
  #carriage_count = nrow(q_carriage)
  carriage_count = t_count - disease_count
  pop_sero_df[which(pop_sero_df$pair_name == p), "carriage"] = carriage_count
  pop_sero_df[which(pop_sero_df$pair_name == p), "disease"] = disease_count
}



### we do not know the number of samples we collected in the GPSC table
### carriage rate is collected from literature reviews
### surveillance population is retrieved from UN stats
carriage_rate_path = "/Users/hc14/Documents/PhD_project/PhD_datasets/GPS/GPS_carriage_rate.xlsx"
carriage_rate_file = read_excel(carriage_rate_path)
carriage_rate_df = as.data.frame(carriage_rate_file)

pop_select_sero_df = pop_sero_df[which(pop_sero_df$study %in% carriage_rate_df$Country), ]

for(study in unique(pop_select_sero_df$study)){
  GPS_subset = GPS_Bayes_selec[which(GPS_Bayes_selec$Country.x == study), ]
  GPS_time = max(GPS_subset$Year_collection, na.rm = TRUE) - min(GPS_subset$Year_collection, na.rm = TRUE) + 1
  study_subset = pop_select_sero_df[which(pop_select_sero_df$study == study), ]
  study_sample = sum(study_subset$carriage) + sum(study_subset$disease)
  study_ref = carriage_rate_df[which(carriage_rate_df$Country == study),]
  study_carriage_rate = study_ref$Carriage_rate
  study_pop_size = study_ref$Pop_size
  pop_select_sero_df[which(pop_select_sero_df$study == study), "carriage_samples"] = round(study_sample/study_carriage_rate)
  pop_select_sero_df[which(pop_select_sero_df$study == study), "surveillance_population"] = round(study_pop_size * 1000000)
  pop_select_sero_df[which(pop_select_sero_df$study == study), "time_interval"] = GPS_time
}

### Indonesia has no time data
{
'''
pop_select_sero_rmIndo_df = pop_select_sero_df[which(pop_select_sero_df$study != "INDONESIA"),]
pop_select_sero_rmIndo_df$study = as.factor(pop_select_sero_rmIndo_df$study)
pop_select_sero_rmIndo_df$type = as.factor(pop_select_sero_rmIndo_df$type)
'''
}

###try to run on a population (USA)
{
USA = pop_select_sero_rmIndo_df[which(pop_select_sero_rmIndo_df$study == "USA"),]
USA = USA[,1:7]
USA_pop_sero = Bim_Poisson(Bframe = USA)
USA$study = as.factor(USA$study)
USA$type = as.factor(USA$type)
}

## removing weird serotypes & small populaiton
weird_serotype = c("15BC", "6E(6A)", "6B/6D", "15A/15B/15C",
                   "UNTYPABLE", "POSSIBLE 6C", "POSSIBLE 6E",
                   "POSSIBLE 6D", "SEROGROUP 24", "COVERAGE TOO LOW",
                   "ALTERNATIVE_ALIB_NT", "33A/33F", "serogroup 24",
                   "ALTERNATIVE", "coverage too low", 
                   "POSSIBLE 06C", "POSSIBLE 06E")

small_pop_list = c("CAMBODIA", "ETHIOPIA", "IRELAND", "KENYA", 
                   "NETHERLANDS", "NEW ZEALAND", "SPAIN",
                   "THAILAND", "TURKEY", "INDONESIA")

bigpop_select_sero_df = pop_select_sero_df[which(!pop_select_sero_df$study %in% small_pop_list),]
bigpop_select_normal_sero_df = bigpop_select_sero_df[which(!bigpop_select_sero_df$type %in% weird_serotype),]
### make South Africa as reference population (make it as the first level)
bigpop_select_normal_sero_df$study = as.factor(bigpop_select_normal_sero_df$study)
bigpop_select_normal_sero_df_study_level = levels(bigpop_select_normal_sero_df$study)
bigpop_select_normal_sero_df_study_level = append(bigpop_select_normal_sero_df_study_level, "SOUTH AFRICA", after = 0)
bigpop_select_normal_sero_df_study_level = bigpop_select_normal_sero_df_study_level[-26]
bigpop_select_normal_sero_df$study = factor(bigpop_select_normal_sero_df$study, levels = bigpop_select_normal_sero_df_study_level)
### factorise character type column
bigpop_select_normal_sero_df$type = as.factor(bigpop_select_normal_sero_df$type)

### collecting lane id for GWAS analysis from GPS table
GPSC_country_select = GPS_merge[which(GPS_merge$Country.x %in% bigpop_select_normal_sero_df$study),]
## SWISS serotype 
GPSC_country_select[GPSC_country_select == "SWISS_NT"] = "SWISS"
GPSC_country_sero_select = GPSC_country_select[which(GPSC_country_select$In_Silico_serotype %in% bigpop_select_normal_sero_df$type),]
GPSC_country_sero_select_laneid = GPSC_country_sero_select$Lane_id
write.table(GPSC_country_sero_select_laneid, "/Users/hc14/Documents/PhD_project/invasiveness/Stan_Bayesian/results/GPSC_country_select_bigpop_normalsero_laneid.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE)


### run Poisson Bayesian model 
bigpop_normalsero_seroadjusted = Bim_Poisson(PopAdjusted = TRUE, Bframe = bigpop_select_normal_sero_df)
### checking convergence of each parameters
rstan::traceplot(bigpop_normalsero_seroadjusted$Model_Fit,
                 pars= "rho_ij")
rstan::traceplot(bigpop_normalsero_seroadjusted$Model_Fit,
                 pars= "log_nu_j")
rstan::traceplot(bigpop_normalsero_seroadjusted$Model_Fit,
                 pars= "gamma_i")
### generate output dataframe
output_path = "/Users/hc14/Documents/PhD_project/invasiveness/Stan_Bayesian/results/bigpop_normalsero_poisson_seroadjusted_output_df.txt"
bigpop_normalsero_poisson_seroadjusted_output_df = progressionEstimation::process_progression_rate_model_output(bigpop_normalsero_seroadjusted$Model_Fit,
                                                                                                        bigpop_select_normal_sero_df)
write.table(bigpop_normalsero_poisson_seroadjusted_output_df, output_path, sep="\t", quote=FALSE, row.names = FALSE)
### plot observed-predicted diagnosis plot
bigpop_normalsero_poisson_seroadjusted_prediction_plot = "/Users/hc14/Documents/PhD_project/invasiveness/Stan_Bayesian/results/bigpop_normalsero_poisson_seroadjusted.pdf"
pdf(file = bigpop_normalsero_poisson_seroadjusted_prediction_plot, width=16, height = 8)
progressionEstimation::plot_case_carrier_predictions(bigpop_normalsero_poisson_seroadjusted_output_df,
                                                     n_label = 5)
dev.off()
### plot progression rate
bigpop_normalsero_poisson_seroadjusted_progression_rate_plot = "/Users/hc14/Documents/PhD_project/invasiveness/Stan_Bayesian/results/bigpop_normalsero_poisson_seroadjusted_progression_rate.pdf"
pdf(file=bigpop_normalsero_poisson_seroadjusted_progression_rate_plot, width = 12, height = 8)
plot_ordered_progression_rates(bigpop_normalsero_poisson_seroadjusted_output_df, 
                               unit_time = "year",
                               type_name = "Serotype")
dev.off()

## select a result of null progression rate to include in the plot
bigpop_normalsero_poisson_null_output_df = progressionEstimation::process_progression_rate_model_output(bigpop_normalsero_poisson_null$Model_Fit,
                                                                                                                bigpop_select_normal_sero_df)
bigpop_normalsero_poisson_null_progression = bigpop_normalsero_poisson_null_output_df[1,]
bigpop_normalsero_poisson_null_progression$type = "Null"
bigpop_normalsero_poisson_null_progression_withNull = rbind(bigpop_normalsero_poisson_seroadjusted_output_df,bigpop_normalsero_poisson_null_progression)
plot_ordered_progression_rates(bigpop_normalsero_poisson_null_progression_withNull, 
                               unit_time = "year",
                               type_name = "Serotype")

### plot population scaling factor
bigpop_normalsero_poisson_seroadjusted_pop_factor_plot = "/Users/hc14/Documents/PhD_project/invasiveness/Stan_Bayesian/results/bigpop_normalsero_poisson_seroadjusted_pop_factor.pdf"
pdf(file=bigpop_normalsero_poisson_seroadjusted_pop_factor_plot, width = 12, height = 8)
plot_ordered_population_scale_factors(bigpop_normalsero_poisson_seroadjusted_output_df)
dev.off()

### generate sample-progression rate list for GWAS analysis
GPSC_country_sero_select_progression_rate = GPSC_country_sero_select
GPSC_country_sero_select_progression_rate$progression_rate = 0
for(seros in unique(bigpop_normalsero_poisson_seroadjusted_output_df$type)){
  progression_rate = bigpop_normalsero_poisson_seroadjusted_output_df[bigpop_normalsero_poisson_seroadjusted_output_df$type == seros, "nu"][1]
  GPSC_country_sero_select_progression_rate[GPSC_country_select_progression_rate$In_Silico_serotype == seros, "progression_rate"] = progression_rate 
}
GPSC_country_sero_select_pheno = GPSC_country_sero_select_progression_rate[,c("Lane_id", "progression_rate")]
write.table(GPSC_country_sero_select_pheno, "/Users/hc14/Documents/PhD_project/invasiveness/Stan_Bayesian/results/GPSC_country_sero_select_pheno.pheno", sep="\t", quote=FALSE, row.names = FALSE)


bigpop_normalsero_poisson_null = Bim_Poisson(PopAdjusted = TRUE, SeroAdjusted=FALSE, Bframe = bigpop_select_normal_sero_df)
bigpop_normalsero_negbio_null = Bim_Poisson(PopAdjusted=TRUE,SeroAdjusted=FALSE, Model = "negbin",Bframe = bigpop_select_normal_sero_df)
bigpop_normalsero_negbio_sero = Bim_Poisson(PopAdjusted=TRUE, Model = "negbin",Bframe = bigpop_select_normal_sero_df)




all_pop_sero = Bim_Poisson(PopAdjusted=TRUE,Bframe = pop_select_sero_rmIndo_df)
print(all_pop_sero$Model_Fit@model_pars)
rstan::traceplot(all_pop_sero$Model_Fit,
                 pars= "nu_j")
all_pop_sero_output_df <- progressionEstimation::process_progression_rate_model_output(all_pop_sero$Model_Fit,
                                                                                       pop_select_sero_rmIndo_df)
write.table(all_pop_sero_output_df, "/Users/hc14/Documents/PhD_project/invasiveness/Stan_Bayesian/results/allpop_binomial_sero_popadjusted_output.txt", sep="\t", quote=FALSE, row.names = FALSE)
all_pop_normal_sero_output_df = all_pop_sero_output_df[which(!all_pop_sero_output_df$type %in% weird_serotype), ]

pdf(file = "all_pop_normal_sero_prediction.pdf", width=16, height = 8)
progressionEstimation::plot_case_carrier_predictions(all_pop_normal_sero_output_df,
                                                     n_label = 5)
dev.off()

plot_ordered_progression_rates(all_pop_normal_sero_output_df, 
                               unit_time = "year",
                               type_name = "Serotype")
## select a result of null progression rate to include in the plot
allpop_null_poisson_progression = all_pop_null_output_df[1,]
allpop_null_poisson_progression$type = "Null"
all_pop_normal_sero_output_df_withNull = rbind(all_pop_normal_sero_output_df,allpop_null_poisson_progression)
plot_ordered_progression_rates(all_pop_normal_sero_output_df_withNull, 
                               unit_time = "year",
                               type_name = "Serotype")

progressionEstimation::plot_progression_rates(all_pop_sero_output_df, 
                                              unit_time = "year",
                                              type_name = "Serotype")
progressionEstimation::plot_study_scale_factors(all_pop_sero_output_df)


big_pop_sero_output_df = all_pop_sero_output_df[which(!all_pop_sero_output_df$study %in% small_pop_list),]
progressionEstimation::plot_study_scale_factors(big_pop_sero_output_df)
plot_population_scale_factors(big_pop_sero_output_df)

all_pop_null = Bim_Poisson(PopAdjusted=TRUE,SeroAdjusted=FALSE, Bframe = pop_select_sero_rmIndo_df)
rstan::traceplot(all_pop_sero$Model_Fit,
                 pars= "nu_j")
all_pop_null_output_df = progressionEstimation::process_progression_rate_model_output(all_pop_null$Model_Fit,
                                                                                      pop_select_sero_rmIndo_df)
write.table(all_pop_null_output_df, "/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/results/allpop_binomial_null_popadjusted_output.txt", sep="\t", quote=FALSE, row.names = FALSE)
progressionEstimation::plot_progression_rates(all_pop_null_output_df, 
                                              unit_time = "year",
                                              type_name = "Serotype")
progressionEstimation::plot_case_carrier_predictions(all_pop_null_output_df, 
                                                     n_label = 3)
save(all_pop_sero, all_pop_null, file = "all_pop_binomial_poisson.RData")


all_pop_sero_negbio = Bim_Poisson(PopAdjusted=TRUE, Model = "negbin",Bframe = pop_select_sero_rmIndo_df)
rstan::traceplot(all_pop_sero_negbio$Model_Fit,
                 pars= "nu_j")
all_pop_sero_negbio_output_df = progressionEstimation::process_progression_rate_model_output(all_pop_sero_negbio$Model_Fit,
                                                                                      pop_select_sero_rmIndo_df)
write.table(all_pop_sero_output_df, "/Users/hc14/Documents/PhD_project/Stan_Bayesian/results/allpop_negbin_sero_popadjusted_output.txt", sep="\t", quote=FALSE, row.names = FALSE)
progressionEstimation::plot_progression_rates(all_pop_sero_negbio_output_df, 
                                              unit_time = "year",
                                              type_name = "Serotype")
progressionEstimation::plot_case_carrier_predictions(all_pop_sero_negbio_output_df, 
                                                     n_label = 3)

all_pop_null_negbio = Bim_Poisson(PopAdjusted=TRUE,SeroAdjusted=FALSE, Model = "negbin",Bframe = pop_select_sero_rmIndo_df)
rstan::traceplot(all_pop_null_negbio$Model_Fit,
                 pars= "nu")
all_pop_null_negbio_output_df = progressionEstimation::process_progression_rate_model_output(all_pop_null_negbio$Model_Fit,
                                                                                      pop_select_sero_rmIndo_df)
write.table(all_pop_null_negbio_output_df, "/Users/hc14/Documents/PhD_project/Stan_Bayesian/results/allpop_null_negbio_popadjusted_output.txt", sep="\t", quote=FALSE, row.names = FALSE)
progressionEstimation::plot_progression_rates(all_pop_null_negbio_output_df, 
                                              unit_time = "year",
                                              type_name = "Serotype")
progressionEstimation::plot_case_carrier_predictions(all_pop_null_negbio_output_df, 
                                                     n_label = 3)




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

#------------------------------------
## generate binary outcome file
progression_rate_results = read.csv(file="/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/results/GPSC_country_sero_select_pheno.pheno", header=TRUE, sep="\t", quote="", row.names = NULL)
merged_rawdata = read.csv("/Users/hc14/Documents/PhD_project/PhD_datasets/GPS/GPS1_merge.csv", sep=",")
selected_laneid = progression_rate_results$Lane_id
merged_selected_rawdata = merged_rawdata[which(merged_rawdata$Lane_id %in% selected_laneid), c("Lane_id", "Manifest_type","Source")]
## check the coherence of "Manifest_type" and "Source"
merged_selected_rawdata[which(merged_selected_rawdata$Source %in% c("PERITONEAL FLUID","ASCETIC FLUID","PLEURAL ASPIRATE","LUNG ASPIRATE","CSF","BLOOD","PLEURAL FLUID","CEREBROSPINAL FLUID","OTHER")),"Manifest_type"] = "IPD"
merged_selected_rawdata[which(merged_selected_rawdata$Source %in% c("_","NP SWAB", "MIDDLE EAR FLUID", "NASOPHARYNGEAL SWAB", "EAR SWAB", "SPUTUM","NP SWAP","PUS")),"Manifest_type"] = "Carriage"
merged_selected_rawdata[which(merged_selected_rawdata$Manifest_type == "IPD"),"Disease_status"] = 1
merged_selected_rawdata[which(merged_selected_rawdata$Manifest_type == "Carriage"),"Disease_status"] = 0
write.table(merged_selected_rawdata, "GPSC_country_sero_select_pheno_diseasetable.txt", sep="\t", quote=FALSE, row.names = FALSE)
GPS1_binary_pheno = merged_selected_rawdata[,c("Lane_id", "Disease_status")]
write.table(GPS1_binary_pheno, "GPS1_binary_pheno.txt", sep="\t", quote=FALSE, row.names = FALSE)

merged_selected_covariate = merged_rawdata[which(merged_rawdata$Lane_id %in% selected_laneid), c("Lane_id","Country.x", "GPSC_PoPUNK2","Age_years", "HIV_status", "PCVType")]
merged_selected_covariate[which(merged_selected_covariate$HIV_status == "_"), "HIV_status"] = NA
merged_selected_covariate[which(merged_selected_covariate$PCVType == "_"), "PCVType"] = NA
write.table(merged_selected_covariate, "GPS1_covariate.txt", sep="\t", quote=FALSE, row.names = FALSE)

#------------------------------------------------
mainDir = "/Users/hc14/Documents/PhD_project/Invasiveness/GWAS/input/GPS1_phylo"
uniq_gpsc = unique(merged_selected_covariate$GPSC_PoPUNK2)

for(gpsc in uniq_gpsc){
  subDir = paste0("GPSC", gpsc)
  dir.create(file.path(mainDir, subDir))
  pathname = paste(mainDir,subDir,sep = "/")
  selected_gpsc_cov = merged_selected_covariate[which(merged_selected_covariate$GPSC_PoPUNK2 == gpsc),]
  filename_cov = paste0(pathname,"/GPSC",gpsc,"_covariates.txt")
  write.table(selected_gpsc_cov, filename_cov, sep="\t", quote=FALSE, row.names = FALSE)
  filename_lane = paste0(pathname,"/GPSC",gpsc,"_laneid.txt")
  selected_gpsc_lane = selected_gpsc_cov[,"Lane_id"]
  write.table(selected_gpsc_lane, filename_lane, sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
}


GPSC_freq_df = as.data.frame(table(merged_selected_covariate$GPSC_PoPUNK2), stringsAsFactors=FALSE)
colnames(GPSC_freq_df) = c("GPSC", "count")
GPSC_freq_phy_df = GPSC_freq_df[which(GPSC_freq_df$count >1),]
for(i in 1:nrow(GPSC_freq_phy_df)){
  GPSC_freq_phy_df[i,"GPSC"] = paste0("GPSC",i)
}
GPSC_freq_phy_df = GPSC_freq_phy_df$GPSC
write.table(GPSC_freq_phy_df, "/Users/hc14/Documents/PhD_project/Invasiveness/GWAS/input/GPSC_name_phy.txt", sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
