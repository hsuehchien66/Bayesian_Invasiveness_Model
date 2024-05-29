// The input data
data {

  int<lower=0> i_max; // study
  int<lower=0> j_max; // serotype
  int<lower=0> g1_max; // gene
  int<lower=0> n_obs; // number of combinations
  array[n_obs] int i_values;
  array[n_obs] int j_values;
  array[n_obs] int g1_values;
  array[n_obs] int c_ijg1;
  array[n_obs] int d_ijg1;
  array[n_obs] int n_i;
  array[n_obs] int N_i;

  // int<lower=0> i_values[n_obs]; // an integer array named "i_values" with a length of n_obs. The values in i_values must be greater than or equal to 0
  //int<lower=0> j_values[n_obs]; // an integer array named "j_values" with a length of n_obs. The values in j_values must be greater than or equal to 0
  //int<lower=0> g1_values[n_obs]; 
  //int<lower=0> c_ijg1[n_obs]; //carriage in certain serotype certain study
  //int<lower=0> d_ijg1[n_obs]; //disease in certain serotype certain study
  //int<lower=0> n_i[n_obs]; // number of all samples (including negative samples) in certain study
  //int<lower=0> N_i[n_obs]; // number of individuals in the certain study under surveillance
  vector<lower=0>[n_obs] t_i; // a vector named t_i with length n_obs
}

// The parameters accepted by the model
parameters {

  // carriage prevalence ~ U(0,1)
  vector<lower=0.0,upper=1.0>[n_obs] rho_ijg1; 

  // log serotype invasiveness ~ U(-6,1)
  vector<lower=-5.0,upper=-1.0>[j_max] log_nu_j;

  // log gene invasiveness ~ Cauchy
  vector<lower=-1.25, upper=1.25>[g1_max-1] log_nu_g1; // [g1_max-1] because relative to one genetype

  // dataset adjustment
  vector<lower=-pi()/2, upper=pi()/2>[i_max-1] gamma_varying; // [i_max-1] because relative to one population

}

// Transformed parameters
transformed parameters {

  // declare transformed parameters
  vector[i_max] gamma_i;
  real mu_pop = 0; // position parameter of Cauchy for gamma (inital mu)
  real tau_pop = 2; // scale parameter of Cauchy for gamma (inital tau)
  real mu_g1 = 0; // position parameter of Cauchy for gene invasiveness
  real tau_g1 = 0.5; // scale parameter of Cauchy for gene invasiveness
  vector<lower=0,upper=10.0>[j_max] nu_j;
  vector[g1_max] nu_g1;

  // calculate serotype invasiveness on a real scale (main characteristic, uniform distribution)
  for (j in 1:j_max) {
    nu_j[j] = pow(10, log_nu_j[j]);
  }

  // calculate gene invasiveness on a real scale (cauchy distribution, concentrate at 0)
  nu_g1[1] = 1;
  for (g1 in 2:g1_max) {
    nu_g1[g1] = pow(10, mu_g1 + tau_g1 * tan(log_nu_g1[g1-1]));
  }

  // add constant to gamma vector
  gamma_i[1] = 1;
  for (i in 2:i_max) {
    gamma_i[i] = pow(10, mu_pop + tau_pop * tan(gamma_varying[i-1]));
  }


}

// The model to be estimated
model {

  // Calculate prior probability for serotypes
  for (j in 1:j_max) {
    target += uniform_lpdf(log_nu_j[j] | -5, -1);
    // log_nu_j[j] ~ uniform(-6, 1)
  }


  // Calculate prior probability for genes
  for (g1 in 2:g1_max) {
    target += uniform_lpdf(log_nu_g1[g1-1] | -1.25, 1.25);
  }

  // Calculate prior probability for study adjustment
  for (i in 2:i_max) {
    target += uniform_lpdf(gamma_varying[i-1] | -pi()/2, pi()/2);
  }

  // iterate over datasets
  for (index in 1:n_obs) {

    // Get serotype
    int j = j_values[index];

    // Get gene
    int g1 = g1_values[index];

    // Get location adjustment
    int i = i_values[index];

    // Calculate prior probability for carriage frequency
    target += beta_lpdf(rho_ijg1[index] | 1, 1);

    // calculate likelihood given data
    target += binomial_lpmf(c_ijg1[index] | n_i[index], rho_ijg1[index]);
    target += poisson_lpmf(d_ijg1[index] | gamma_i[i]*nu_j[j]*nu_g1[g1]*rho_ijg1[index]*N_i[index]*t_i[index]);

  }
}

generated quantities {

  // Calculate and store log likelihood for loo
  vector[n_obs] carriage_log_lik;
  vector[n_obs] disease_log_lik;
  vector[n_obs] log_lik;

  // Calculate and store predictions for carriage
  vector[n_obs] c_ijg1_pred;

  // Calculate and store predictions for disease
  vector[n_obs] d_ijg1_pred;

  // iterate over datasets
  for (index in 1:n_obs) {

    // Get serotype
    int j = j_values[index];


    // Get gene
    int g1 = g1_values[index];

    // Get location adjustment
    int i = i_values[index];

    // Store predictions
    c_ijg1_pred[index] = n_i[index]*rho_ijg1[index];
    d_ijg1_pred[index] = gamma_i[i]*nu_j[j]*nu_g1[g1]*rho_ijg1[index]*N_i[index]*t_i[index];

    // Calculate likelihood given data
    carriage_log_lik[index] = binomial_lpmf( c_ijg1[index] | n_i[index], rho_ijg1[index] );
    disease_log_lik[index] = poisson_lpmf(  d_ijg1[index] | d_ijg1_pred[index] );
    log_lik[index] = carriage_log_lik[index] + disease_log_lik[index];

  }

}
