#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR(site_obs_ia);
  DATA_IVECTOR(site_index_i);
  DATA_SCALAR(log_sd_site);
  
  
  // Parameters
  PARAMETER_VECTOR(site_mean_a);
  PARAMETER(global_mean);
  PARAMETER(log_sd_global);
  
  
  // Objective function
  Type jnll = 0;
  
  Type sd_global = exp(log_sd_global);
  
  // Probability of site means
  int n_site = site_mean_a.size();
  for( int s = 0; s < n_site; s++) {
    jnll -= dnorm(site_mean_a(s),
                  global_mean, 
                  sd_global, true);
  }
  
  // Probability of observations
  Type sd_site = exp(log_sd_site);
  int n_site_obs = site_obs_ia.size();
  for(int i = 0; i < n_site_obs; i++) {
    jnll -= dnorm(site_obs_ia(i), 
                  site_mean_a(site_index_i(i)),
                  sd_site, true);
  }
  
  // Reporting
  ADREPORT(global_mean);
  ADREPORT(sd_global);
  ADREPORT(site_mean_a);
  
  
  return jnll;
}
