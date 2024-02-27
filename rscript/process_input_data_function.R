require(tidyverse)
require(magrittr)
require(rstan)
require(bridgesampling)
require(loo)
require(kableExtra)
require(xlsx)
require(cowplot)
require(ggrepel)

combine_rows <- function(df, col_name = "type") {
  df %<>%
    dplyr::group_by(study,!!! dplyr::syms(col_name)) %>%
    dplyr::mutate(carriage = sum(carriage)) %>%
    dplyr::mutate(disease = sum(disease)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct()
  return(df)
}


#' Generate model input from input data frame
#'
#' @param input_df Data frame containing details of case-carrier studies
#' @param type Name of column in input data frame used for typing information
#' @param use_strain Boolean specifying whether strain information should be used in addition to type information
#' @param combine_strain Boolean specifying whether strain information should be combined with type information
#' @param condense Boolean specifying whether repeated entries of the same type should be condensed into a single entry
#'
#' @return A list of lists used as an input to stan models
#' @export
#'
#' @importFrom rlang :=
#' @importFrom rlang !!
#'
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


