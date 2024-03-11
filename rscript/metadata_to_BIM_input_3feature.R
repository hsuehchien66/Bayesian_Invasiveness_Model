metadata_df <- GPS_curated_table
pop_name <- "Country.x"
pair_colname1 <- "In_Silico_serotype"
pair_colname2 <- "GATTATAATGTTACACCGAATTTTGTAGACC"
pair_colname3 <- "GATTTTCATTGCCGTTATGCCAAGCATAGCA"
BIM_colnames = c("study", "feature1", "carriage", "disease","carriage_samples", "surveillance_population", "time_interval", "feature2", "feature3", "pair_name")

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




