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
  

  
  
