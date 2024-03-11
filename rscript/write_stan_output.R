### ---------------
BIM_sero_iga_input <- read.table("/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_input/BIM_sero_iga_input.txt", sep = "\t", header = TRUE, stringsAsFactors = TRUE)
BIM_iga_sero_input <- BIM_sero_iga_input
colnames(BIM_iga_sero_input) <- c("study", "feature2", "carriage", "disease", "carriage_samples", "surveillance_population", "time_interval", "feature1")

iga_based_serotype_adjusted_model_fit <- readRDS("/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_output_rdata/small_pop/iga_based_serotype_adjusted_model_fit.rds")

input_df$feature1 <- as.factor(input_df$feature1)
input_df <- BIM_iga_sero_input
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
