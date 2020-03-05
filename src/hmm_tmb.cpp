
#include "dist.hpp"

// Compute Negative log-likelihood for HMM
// DATA:
//   data: vector of observed counts
//   distname: name of observation distribution
//   n_states: number of states
// PARAMETERS:
//   wpar: vector of working parameters

template<class Type>
Type objective_function<Type>::operator() ()
{
  // DATA
  DATA_MATRIX(data); // data stream
  DATA_INTEGER(n_states); // number of states
  DATA_IVECTOR(distcode); // codes of observation distributions
  
  // PARAMETERS
  PARAMETER_VECTOR(ltpm);
  PARAMETER_VECTOR(wpar);
  
  int n_var = distcode.size();
  
  //======================//
  // Transform parameters //
  //======================//
  // Transition probability matrix
  matrix<Type> tpm(n_states, n_states);
  int cur = 0;
  for (int i = 0; i < n_states; ++i) {
    tpm(i, i) = 1;
    for (int j = 0; j < n_states; ++j) {
      if (i != j) {
        tpm(i, j) = exp(ltpm(cur));
        ++cur;
      }
    }
    tpm.row(i) /= tpm.row(i).sum();
  }
  
  //=================================//
  // Compute stationary distribution //
  //=================================//
  matrix<Type> delta(1, n_states);
  matrix<Type> I = matrix<Type>::Identity(n_states, n_states);
  matrix<Type> tpminv = I;
  tpminv -= tpm;
  tpminv = (tpminv.array() + 1).matrix();
  matrix<Type> ivec(1, n_states); 
  for (int i = 0; i < n_states; ++i) 
    ivec(0, i) = 1;
  
  // if tpm is ill-conditioned then just use uniform initial distribution
  try {
    tpminv = tpminv.inverse();
    delta = ivec * tpminv;
  } catch(...) {
    for (int i = 0; i < n_states; ++i) 
      delta(0, i) = 1.0 / n_states;
  }
  
  //===================================//  
  // Compute observation probabilities //
  //===================================//
  // Initialise matrix of probabilities to 1
  int n = data.rows();
  matrix<Type> prob(n, n_states);
  for(int i = 0; i < n; i++) {
    for(int s = 0; s < n_states; s++) {
      prob(i, s) = 1;
    }
  }
  
  // Counter to subset parameter vector
  int par_count = 0;
  
  // Loop over observed variables
  for(int var = 0; var < n_var; var++) {
    // Define observation distribution
    Dist <Type> obsdist(distcode(var));
    
    // Subset and transform observation parameters
    vector<Type> sub_wpar = wpar.segment(par_count, obsdist.npar() * n_states);
    par_count = par_count + obsdist.npar() * n_states;
    matrix<Type> par = obsdist.invlink(sub_wpar, n_states);
    
    // Loop over states (columns)
    for (int s = 0; s < n_states; ++s) {
      // Vector of parameters for state s
      vector<Type> subpar = par.row(s);
      
      // Loop over observations (rows)
      for (int i = 0; i < n; ++i) {
        prob(i, s) = prob(i, s) * obsdist.pdf(data(i, var), subpar, false);
      }
    }
  }
  
  //========================//
  // Compute log-likelihood //
  //========================//
  Type llk = 0;
  matrix<Type> phi(delta);
  Type sumphi = 0;
  for (int i = 0; i < n; ++i) {
    phi = (phi.array() * prob.row(i).array()).matrix();
    phi = phi * tpm;
    sumphi = phi.sum();
    llk += log(sumphi);
    phi /= sumphi;
  }
  Type nll = -llk;
  
  return nll;
}
