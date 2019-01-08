#include <TMB.hpp>
using namespace density;

template <class Type> // define function to bound phi between -1 and 1
Type bound(Type x){return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR(x);
  DATA_VECTOR(y);
  
  // Parameters
  PARAMETER(b0);
  PARAMETER(b1);
  PARAMETER(unbounded_phi);
  PARAMETER(log_ar_sd);
  
  // Transform variables
  Type ar_sd = exp(log_ar_sd);
  Type phi = bound(unbounded_phi);
  
 
  // Generate model fit
  int n_y = y.size(); // number of observations
  vector<Type> fit_y(n_y); // vector for fit
  fit_y = b0 + b1 * x; // make fit
  
  
  // Likelihood of observations
  Type jnll = 0;
  jnll += SCALE(AR1(phi), ar_sd)(y - fit_y); //already negative and logged 
  
  // Reporting
  ADREPORT(fit_y);
  ADREPORT(b0);
  ADREPORT(b1);
  ADREPORT(unbounded_phi);
  ADREPORT(log_ar_sd);
  
  return jnll;
  
}
