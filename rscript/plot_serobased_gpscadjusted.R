### need to reinstall packages after upgrading R
install.packages(c("readxl", "imputeTS", "tidyverse", "magrittr", "rstan", 
                   "bridgesampling", "loo", "cowplot", "ggrepel", "xlsx", 
                   "roxygen2", "ggpubr", "gridExtra", "dplyr"))
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

### function
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
### read BIM input
BIM_sero_gpsc_input <- read.table("/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_input/BIM_sero_gpsc_input.txt", 
                                  header = TRUE, sep = "\t", stringsAsFactors = TRUE)

BIM_sero_gpsc_input$strain <- as.factor(BIM_sero_gpsc_input$strain)

BIM_sero_variant_input <- read.table("/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_input/BIM_sero_variant_input.txt",
                                     header = TRUE, sep = "\t", stringsAsFactors = TRUE)
BIM_sero_variant_input$type <- as.factor(BIM_sero_variant_input$type)

### load model fit
setwd("/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_output_rdata")
load("serotype_gpsc_only_BIM.RData")
load("serotype_gpsc_BIM.RData")
load("gpsc_BIM.RData")
load("variantbased_nopopadjusted.RData")
load("variant_only_BIM.RData")
load("s_pneumoniae_serobased_v.RData")
load("serotype_variant_BIM.RData")
load("s_pneumoniae_variantbased_g.RData")
load("s_pneumoniae_gpscbased_v.RData")
load("s_pneumoniae_gpscbased_variantadjusted.RData")

### variant based output
s_pneumoniae_poisson_variantbased_nopopadjusted_output_df <- progressionEstimation::process_progression_rate_model_output(s_pneumoniae_poisson_variantbased_nopopadjusted_fit, 
                                                                                                         BIM_sero_variant_input,
                                                                                                         type = "type",
                                                                                                         strain_as_primary_type = FALSE,
                                                                                                         strain_as_secondary_type = FALSE, 
                                                                                                         condense = TRUE)
write.table(s_pneumoniae_poisson_variantbased_nopopadjusted_output_df, 
            file = "/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_output_results/05_variant_based/s_pneumoniae_poisson_variantbased_nopopadjusted_output_df.txt", 
            sep="\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

case_carrier_pred_variantbased_nopopadjusted = progressionEstimation::plot_case_carrier_predictions(s_pneumoniae_poisson_variantbased_nopopadjusted_output_df , n_label = 3)

variantbased_nopopadjusted_prate = plot_ordered_progression_rates(s_pneumoniae_poisson_variantbased_nopopadjusted_output_df,
                                                                           type="type",
                                                                           unit_time= "year",
                                                                           type_name= "variant")

pdf(file = "/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_output_results/05_variant_based/case_carrier_pred_variantbased_nopopadjusted.pdf", width=16, height = 8)
case_carrier_pred_variantbased_nopopadjusted
dev.off()

pdf(file = "/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_output_results/05_variant_based/variantbased_nopopadjusted_prate.pdf", width=16, height = 8)
variantbased_nopopadjusted_prate
dev.off()

### GPSC based output
s_pneumoniae_poisson_gpscbased_output_df <- progressionEstimation::process_progression_rate_model_output(s_pneumoniae_poisson_gpsc_fit, 
                                                                                                         BIM_sero_gpsc_input,
                                                                                                         type = "strain",
                                                                                                         strain_as_primary_type = FALSE,
                                                                                                         strain_as_secondary_type = FALSE, 
                                                                                                         condense = TRUE)
                                                                                                         
write.table(s_pneumoniae_poisson_gpscbased_output_df, 
            file = "/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_output_results/02_GPSC_based/s_pneumoniae_poisson_gpscbased_output_df.txt", 
            sep="\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

case_carrier_pred_gpscbased = progressionEstimation::plot_case_carrier_predictions(s_pneumoniae_poisson_gpscbased_output_df , n_label = 3, label_col = "strain")

gpscbased_prate = plot_ordered_progression_rates(s_pneumoniae_poisson_gpscbased_output_df,
                                                                  type="strain",
                                                                  unit_time= "year",
                                                                  type_name= "GPSC")

pdf(file = "/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_output_results/02_GPSC_based/case_carrier_pred_gpscbased.pdf", width=16, height = 8)
case_carrier_pred_gpscbased
dev.off()

pdf(file = "/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_output_results/02_GPSC_based/gpscbased_prate.pdf", width=16, height = 8)
gpscbased_prate
dev.off()



### serotypebased gpscadjusted

s_pneumoniae_poisson_serotypebased_gpscadjusted_output_df <- progressionEstimation::process_progression_rate_model_output(s_pneumoniae_poisson_serobased_gpsc_adjust_fit, 
                                                                                                                             BIM_sero_gpsc_input,
                                                                                                                             strain_as_secondary_type = TRUE)

write.table(s_pneumoniae_poisson_serotypebased_gpscadjusted_output_df, 
            file = "/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_output_results/03_serotype_based_GPSC_adjusted/s_pneumoniae_poisson_serotypebased_gpscadjusted_output_df.txt", 
            sep="\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

case_carrier_pred_serotypebased_gpscadjusted = progressionEstimation::plot_case_carrier_predictions(s_pneumoniae_poisson_serotypebased_gpscadjusted_output_df , n_label = 3)

serotype_prate_serotypebased_gpscadjusted = plot_ordered_progression_rates(s_pneumoniae_poisson_serotypebased_gpscadjusted_output_df,
                                                                                            type="type",
                                                                                            unit_time= "year",
                                                                                            type_name= "Serotype")

popfactor_serotypebased_gpscadjusted = plot_ordered_population_scale_factors(s_pneumoniae_poisson_serotypebased_gpscadjusted_output_df)

pdf(file = "/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_output_results/03_serotype_based_GPSC_adjusted/case_carrier_pred_serotypebased_gpscadjusted.pdf", width=16, height = 8)
case_carrier_pred_serotypebased_gpscadjusted
dev.off()

pdf(file = "/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_output_results/03_serotype_based_GPSC_adjusted/serotype_prate_serotypebased_gpscadjusted.pdf", width=16, height = 8)
serotype_prate_serotypebased_gpscadjusted
dev.off()

pdf(file = "/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_output_results/03_serotype_based_GPSC_adjusted/popfactor_serotypebased_gpscadjusted .pdf", width=16, height = 8)
popfactor_serotypebased_gpscadjusted 
dev.off()

### model selection

### LOO-CV can't be done using different input data
### Bayes factor can't be done using loading Rdata
# In this table, each row is a different model, ordered from the best-fitting at the top, to the worst-fitting at the bottom. This is based on the expected log pointwise predictive density (ELPD) for each model. The difference between  models is calculated in the *elpd_diff* column; negative values indicate worse fits. The standard error of the ELPD values across a model is given by the *se_diff* column; all of these values are relative to the best-fitting model. One model can be regarded as outperforming another when *elpd_diff* is greater than four, and larger than the corresponding *se_diff*. Here we can conclude that best-fitting model is that with type-specific progression rates; that is, there are significant differences in progression rates between serotypes.
# These processing commands also specify `condense = FALSE` (the default for the function). This means that the division of the population will replicate that in the rows of the input data. This is necessary for comparing models where progression rates are determined by different characteristics, because cross-validation requires the same number of data points in each compared model. However, `condense = TRUE` will collapse all entries for the same categorisation; e.g. if `type = 'serotype'`, then all duplicated entries of serotype belonging to the same study would be combined into a single entry. This would be necessary if appending this dataset to a wider dataset in which only serotype, and not strain, were known

s_pneumoniae_poisson_variantbased_fit@model_name <- "Variant-based"
s_pneumoniae_poisson_variantbased_nopopadjusted_fit@model_name <- "s_pneumoniae_poisson_variantbased_nopopadjusted_fit"
s_pneumoniae_poisson_gpsc_fit@model_name <- "s_pneumoniae_poisson_gpscbased_fit"
s_pneumoniae_poisson_serobased_v_fit@model_name <- "Serotype-based"
s_pneumoniae_poisson_serobased_variantadjusted_fit@model_name <- "Serotype-based Variant-adjusted"
s_pneumoniae_poisson_gpscbased_v_fit@model_name <- "GPSC-based"
s_pneumoniae_poisson_variantbased_g_fit@model_name <- "Variant_based"
s_pneumoniae_poisson_gpscbased_variantadjusted_fit@model_name <- "GPSC-based Variant-adjusted"
s_pneumoniae_poisson_gpsc_fit@model_name <- "GPSC-based"
s_pneumoniae_poisson_sero_fit@model_name <- "Serotype-based"
s_pneumoniae_poisson_gpscbased_seroadjust_fit@model_name <- "GPSC-based Serotype-adjusted"
s_pneumoniae_poisson_serobased_gpsc_adjust_fit@model_name <- "Serotype-based GPSC-adjusted"

# Serotype v.s. GPSC
## Cross-Validation
Sero_GPSC_loo_loglik <- progressionEstimation::compare_model_fits_with_loo(list(s_pneumoniae_poisson_gpsc_fit, s_pneumoniae_poisson_sero_fit
                                                                                ), 
                                                                          log_lik_param = "log_lik") 

# Serotype v.s. Variant
## Cross-Validation
Sero_Var_loo_loglik <- progressionEstimation::compare_model_fits_with_loo(list(s_pneumoniae_poisson_variantbased_fit, s_pneumoniae_poisson_serobased_v_fit, s_pneumoniae_poisson_serobased_variantadjusted_fit), 
                                                                       log_lik_param = "log_lik") 
Sero_Var_loo_loglik %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling(latex_options = "scale_down")

write.table(Sero_Var_loo_loglik, file = "/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_output_results/CV_Serotype_vs_Variant.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

Sero_Var_loo_disloglik <- progressionEstimation::compare_model_fits_with_loo(list(s_pneumoniae_poisson_variantbased_fit, s_pneumoniae_poisson_serobased_v_fit, s_pneumoniae_poisson_serobased_variantadjusted_fit), 
                                                                          log_lik_param = "disease_log_lik") 
Sero_Var_loo_disloglik %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling(latex_options = "scale_down")

write.table(Sero_Var_loo_loglik, file = "/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_output_results/Sero_Var_loo_disloglik.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


# GPSC v.s. Variant
## Cross-Validation

GPSC_Var_loo_loglik <- progressionEstimation::compare_model_fits_with_loo(list(s_pneumoniae_poisson_gpscbased_v_fit, s_pneumoniae_poisson_gpscbased_variantadjusted_fit, s_pneumoniae_poisson_variantbased_g_fit), 
                                                                          log_lik_param = "log_lik") 
GPSC_Var_loo_loglik %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling(latex_options = "scale_down")

write.table(GPSC_Var_loo_loglik, file = "/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_output_results/GPSC_Var_loo_loglik.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

GPSC_Var_loo_disloglik <- progressionEstimation::compare_model_fits_with_loo(list(s_pneumoniae_poisson_gpscbased_v_fit, s_pneumoniae_poisson_gpscbased_variantadjusted_fit, s_pneumoniae_poisson_variantbased_g_fit), 
                                                                          log_lik_param = "disease_log_lik") 

GPSC_Var_loo_disloglik %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling(latex_options = "scale_down")

write.table(GPSC_Var_loo_disloglik, file = "/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/BIM_output_results/GPSC_Var_loo_disloglik.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Serotype vs. GPSC
Sero_GPSC_loo_loglik <- progressionEstimation::compare_model_fits_with_loo(list(s_pneumoniae, s_pneumoniae_poisson_gpscbased_variantadjusted_fit, s_pneumoniae_poisson_variantbased_g_fit), 
                                                                          log_lik_param = "log_lik") 

### using bayes factor
# can't be done when loading
# check https://github.com/quentingronau/bridgesampling/issues/7

variant_vs_gpsc_comparison <- progressionEstimation::compare_model_fits_with_bf(list(s_pneumoniae_poisson_variantbased_fit, s_pneumoniae_poisson_variantbased_nopopadjusted_fit))
  dplyr::rename("log(Bayes Factor)" = log_Bayes_factor) %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling(latex_options = "scale_down")
