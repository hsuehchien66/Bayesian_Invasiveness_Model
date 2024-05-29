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


process_progression_rate_model_output<-function(model_output,
                                                input_df,
                                                type = "type",
                                                strain_as_primary_type = FALSE,
                                                strain_as_secondary_type = FALSE,
                                                combined_strain_type = FALSE,
                                                condense = FALSE) {
  # Add model name
  input_df %<>%
    dplyr::mutate(model_name = model_output@model_name)
  # Process input data
  if (strain_as_primary_type | strain_as_secondary_type | combined_strain_type) {
    input_df %<>%
      tidyr::unite("combined",!!type,strain, sep='_', remove = FALSE) %>%
      dplyr::mutate(combined = factor(combined))
    if (combined_strain_type) {
      type = "combined"
    }
    input_df %<>%
      dplyr::select(model_name,
                    study,
                    !!type,
                    carriage,
                    disease,
                    carriage_samples,
                    surveillance_population,
                    time_interval,
                    strain)
  }
  # Remove unused rows
  if (condense) {
    if (strain_as_primary_type | strain_as_secondary_type) {
      input_df <- combine_rows(input_df %>%
                                 dplyr::select(model_name,
                                               study,
                                               !!type,
                                               carriage,
                                               disease,
                                               carriage_samples,
                                               surveillance_population,
                                               time_interval,
                                               strain)) %>%
        dplyr::distinct()
    } else {
      input_df <- combine_rows(input_df %>%
                                 dplyr::select(model_name,
                                               study,
                                               !!type,
                                               carriage,
                                               disease,
                                               carriage_samples,
                                               surveillance_population,
                                               time_interval),
                               col_name = type)
    }
  }
  # Extract factor levels
  i_levels = levels(input_df %>% dplyr::pull(study))
  j_levels = levels(input_df %>% dplyr::pull(!!type))
  if (strain_as_primary_type | strain_as_secondary_type) {
    k_levels = levels(input_df %>% dplyr::pull(strain))
  }
  # Carriage prevalence estimates
  carriage_df <- data.frame(
    "rho" = get_median("rho_ij",model_output),
    "rho_lower" = get_lower("rho_ij",model_output),
    "rho_upper" = get_upper("rho_ij",model_output)
  )
  input_df %<>% dplyr::bind_cols(carriage_df)
  # Variation by location
  scale_parameter <- 1
  if ("gamma_i" %in% model_output@model_pars) {
    location_parameters <- data.frame(
      "study" = i_levels,
      "gamma" = get_median("gamma_i",model_output),
      "gamma_lower" = get_lower("gamma_i",model_output),
      "gamma_upper" = get_upper("gamma_i",model_output)
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
  # Calculate invasiveness values
  nu_name = "nu"
  if ("nu_j" %in% model_output@model_pars) {
    nu_name = "nu_j"
  }
  progression_rates_df <- data.frame(
    "type" = j_levels,
    "nu" = as.numeric(get_median(nu_name,model_output)),
    "nu_lower" = as.numeric(get_lower(nu_name,model_output)),
    "nu_upper" = as.numeric(get_upper(nu_name,model_output))
  )
  input_df %<>% dplyr::left_join(progression_rates_df, by = setNames("type",type))
  
  if ("nu_k" %in% model_output@model_pars) {
    secondary_progression_rates_df <- data.frame(
      "type" = k_levels,
      "secondary_nu" = get_median("nu_k",model_output),
      "secondary_nu_lower" = get_lower("nu_k",model_output),
      "secondary_nu_upper" = get_upper("nu_k",model_output)
    )
    input_df %<>% dplyr::left_join(secondary_progression_rates_df, by = setNames("type","strain"))
  }
  
  if ("phi_nb" %in% model_output@model_pars) {
    precision_parameters_df <- data.frame(
      "phi" = get_median("phi_nb",model_output),
      "phi_lower" = get_lower("phi_nb",model_output),
      "phi_upper" = get_upper("phi_nb",model_output)
    )
    input_df %<>% dplyr::bind_cols(precision_parameters_df)
  }
  
  # Extract predictions and intervals
  input_df %<>%
    dplyr::mutate(carriage_prediction = get_median("c_ij_pred",model_output)) %>%
    dplyr::mutate(carriage_prediction_lower = get_lower("c_ij_pred",model_output)) %>%
    dplyr::mutate(carriage_prediction_upper =  get_upper("c_ij_pred",model_output)) %>%
    dplyr::mutate(disease_prediction = get_median("d_ij_pred",model_output)) %>%
    dplyr::mutate(disease_prediction_lower = get_lower("d_ij_pred",model_output)) %>%
    dplyr::mutate(disease_prediction_upper =  get_upper("d_ij_pred",model_output))
  
  # Add in absolute deviation
  input_df %<>%
    dplyr::mutate(carriage_abs_dev = abs(carriage - carriage_prediction)) %>%
    dplyr::mutate(disease_abs_dev = abs(disease - disease_prediction))
  
  return(input_df)
}