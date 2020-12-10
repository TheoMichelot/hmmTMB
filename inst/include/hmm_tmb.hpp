#include<memory>
#include "dist.hpp"

#ifndef _HMMTMB_
#define _HMMTMB_
#include <TMB.hpp>
#endif

//' Compute Negative log-likelihood for HMM
//' DATA:
//'   ID: vector of time series IDs
//'   data: matrix of response variables
//'   X_fe: block-diagonal design matrix for fixed effects
//'   X_re: block-diagonal design matrix for random effects
//'   S: block-diagonal penalty matrix
//'   ncol_re: number of columns of S and X_re for each random effect
//'   n_states: number of states
//'   distcode: vector of codes of observation distributions
//' PARAMETERS:
//'   ltpm: parameters for transition probabilities
//'   coeff_fe: observation parameters (fixed effects)
//'   coeff_re: observation parameters (random effects)
//'   log_lambda: vector of smoothness parameters

template<class Type>
Type objective_function<Type>::operator() ()
{
  // DATA
  DATA_VECTOR(ID); // vector of time series IDs
  DATA_MATRIX(data); // data stream
  DATA_INTEGER(n_states); // number of states
  DATA_IVECTOR(distcode); // codes of observation distributions
  // model matrices for observation process
  DATA_SPARSE_MATRIX(X_fe_obs); // design matrix for fixed effects
  DATA_SPARSE_MATRIX(X_re_obs); // design matrix for random effects
  DATA_SPARSE_MATRIX(S_obs); // penalty matrix
  DATA_IVECTOR(ncol_re_obs); // number of columns of S and X_re for each random effect
  // model matrices for hidden state process
  DATA_SPARSE_MATRIX(X_fe_hid); // design matrix for fixed effects
  DATA_SPARSE_MATRIX(X_re_hid); // design matrix for random effects
  DATA_SPARSE_MATRIX(S_hid); // penalty matrix
  DATA_IVECTOR(ncol_re_hid); // number of columns of S and X_re for each random effect
  // prior information 
  DATA_MATRIX(coeff_fe_obs_prior); // means, sds for prior on fixed effects for obs 
  DATA_MATRIX(coeff_fe_hid_prior); // means, sds for prior on fixed effects for hidden 
  DATA_MATRIX(log_lambda_obs_prior); // means, sds for prior on smoothing parameters for obs 
  DATA_MATRIX(log_lambda_hid_prior); // means, sds for prior on smoothing parameters for hidden 
  
  // PARAMETERS (fixed effects first, then random effects)
  PARAMETER_VECTOR(coeff_fe_obs); // observation parameters (fixed effects)
  PARAMETER_VECTOR(log_lambda_obs); // smoothness parameters
  PARAMETER_VECTOR(coeff_fe_hid); // state process parameters (fixed effects)
  PARAMETER_VECTOR(log_lambda_hid); // smoothness parameters
  PARAMETER_VECTOR(log_delta); // initial distribution
  PARAMETER_VECTOR(coeff_re_obs); // observation parameters (random effects)
  PARAMETER_VECTOR(coeff_re_hid); // state process parameters (random effects)

  // Number of observed variables
  int n_var = distcode.size();
  // Number of data rows
  int n = data.rows();
  
  //======================//
  // Transform parameters //
  //======================//
  // Observation parameters
  vector<Type> par_vec = X_fe_obs * coeff_fe_obs + X_re_obs * coeff_re_obs;
  matrix<Type> par_mat(n, par_vec.size()/n);
  for(int i = 0; i < par_mat.cols(); i++) {
    // Matrix with one row for each time step and
    // one column for each parameter
    par_mat.col(i) = par_vec.segment(i*n, n);
  }
  
  // Transition probabilities
  vector<Type> ltpm_vec = X_fe_hid * coeff_fe_hid + X_re_hid * coeff_re_hid;
  matrix<Type> ltpm_mat(n, ltpm_vec.size()/n);
  for(int i = 0; i < ltpm_mat.cols(); i++) {
    // Matrix with one row for each time step and
    // one column for each transition probability
    ltpm_mat.col(i) = ltpm_vec.segment(i*n, n);
  }

  // Create vector of transition probability matrices
  vector<matrix<Type> > tpm_array(n);
  for(int i = 0; i < n; i++) {
    matrix<Type> tpm(n_states, n_states);
    int cur = 0;
    for (int j = 0; j < n_states; j++) {
      tpm(j, j) = 1;
      for (int k = 0; k < n_states; k++) {
        if (j != k) {
          tpm(j, k) = exp(ltpm_mat(i, cur));
          cur++;
        }
      }
      tpm.row(j) = tpm.row(j)/tpm.row(j).sum();
    }
    tpm_array(i) = tpm;
  }
  
  // Initial distribution
  matrix<Type> delta(1, n_states);
  delta.setOnes();
  for(int i = 0; i < n_states - 1; i++)
    delta(0, i) = exp(log_delta(i));
  delta = delta/delta.sum();
  
  //===================================//  
  // Compute observation probabilities //
  //===================================//
  // Initialise matrix of probabilities to 1
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
    std::unique_ptr<Dist<Type>> obsdist = dist_generator<Type>(distcode(var));

    // Loop over observations (rows)
    for (int i = 0; i < n; ++i) {
      // Don't update likelihood if the observation is missing
      if(!R_IsNA(asDouble(data(i, var)))) {
        // Subset and transform observation parameters
        vector<Type> sub_wpar = 
          par_mat.row(i).segment(par_count, obsdist->npar() * n_states);
        matrix<Type> par = obsdist->invlink(sub_wpar, n_states);
        
        // Loop over states (columns)
        for (int s = 0; s < n_states; ++s) {
          // Vector of parameters for state s
          vector<Type> subpar = par.row(s);
          
          prob(i, s) = prob(i, s) * obsdist->pdf(data(i, var), subpar, false);
        }        
      }
    }

    par_count = par_count + obsdist->npar() * n_states;
  }
  
  //======================//
  // Priors               //
  //======================//
  Type llk = 0; 
  // fixed effects for observation 
  for (int i = 0; i < coeff_fe_obs.size(); ++i) {
    if (!R_IsNA(asDouble(coeff_fe_obs_prior(i, 0)))) {
      llk += dnorm(coeff_fe_obs(i), coeff_fe_obs_prior(i, 0), coeff_fe_obs_prior(i, 1), 1.0); 
    }
  }
  // fixed effects for hidden  
  for (int i = 0; i < coeff_fe_hid.size(); ++i) {
    if (!R_IsNA(asDouble(coeff_fe_hid_prior(i, 0)))) {
      llk += dnorm(coeff_fe_hid(i), coeff_fe_hid_prior(i, 0), coeff_fe_hid_prior(i, 1), 1.0); 
    }
  }
  // smoothing parameters for observation
  if (ncol_re_obs(0) > 0) {
    for (int i = 0; i < log_lambda_obs.size(); ++i) {
      if (!R_IsNA(asDouble(log_lambda_obs_prior(i, 0)))) {
        llk += dnorm(log_lambda_obs(i), log_lambda_obs_prior(i, 0), log_lambda_obs_prior(i, 1), 1.0); 
      }
    }
  }
  // smoothing parameters for hidden 
  if (ncol_re_hid(0) > 0) {
    for (int i = 0; i < log_lambda_hid.size(); ++i) {
      if (!R_IsNA(asDouble(log_lambda_hid_prior(i, 0)))) {
        llk += dnorm(log_lambda_hid(i), log_lambda_hid_prior(i, 0), log_lambda_hid_prior(i, 1), 1.0); 
      }
    }
  }
  
  //========================//
  // Compute log-likelihood //
  //========================//
  // Initialise log-likelihood
  matrix<Type> phi(delta);
  Type sumphi = 0;
  
  // Forward algorithm
  for (int i = 0; i < n; ++i) {
    // Re-initialise phi at first observation of each time series
    if(i == 0 || ID(i-1) != ID(i)) {
      phi = delta;
    }
    phi = (phi.array() * prob.row(i).array()).matrix();
    phi = phi * tpm_array(i);
    sumphi = phi.sum();
    llk = llk + log(sumphi);
    phi = phi / sumphi;
  }
  
  // Negative log-likelihood
  Type nllk = -llk;
  
  //===================//
  // Smoothing penalty //
  //===================//
  // Are there smooths in the observation model?
  if(ncol_re_obs(0) > 0) {
    // Index in matrix S
    int S_start = 0;

    // Loop over smooths
    for(int i = 0; i < ncol_re_obs.size(); i++) {
      // Size of penalty matrix for this smooth
      int Sn = ncol_re_obs(i);

      // Penalty matrix for this smooth
      Eigen::SparseMatrix<Type> this_S = S_obs.block(S_start, S_start, Sn, Sn);

      // Coefficients for this smooth
      vector<Type> this_coeff_re = coeff_re_obs.segment(S_start, Sn);

      // Add penalty
      nllk = nllk -
        Type(0.5) * Sn * log_lambda_obs(i) +
        Type(0.5) * exp(log_lambda_obs(i)) * density::GMRF(this_S).Quadform(this_coeff_re);

      // Increase index
      S_start = S_start + Sn;
    }
  }

  // Are there smooths in the hidden state model?
  if(ncol_re_hid(0) > 0) {
    // Index in matrix S
    int S_start = 0;

    // Loop over smooths
    for(int i = 0; i < ncol_re_hid.size(); i++) {
      // Size of penalty matrix for this smooth
      int Sn = ncol_re_hid(i);

      // Penalty matrix for this smooth
      Eigen::SparseMatrix<Type> this_S = S_hid.block(S_start, S_start, Sn, Sn);

      // Coefficients for this smooth
      vector<Type> this_coeff_re = coeff_re_hid.segment(S_start, Sn);

      // Add penalty
      nllk = nllk -
        Type(0.5) * Sn * log_lambda_hid(i) +
        Type(0.5) * exp(log_lambda_hid(i)) * density::GMRF(this_S).Quadform(this_coeff_re);

      // Increase index
      S_start = S_start + Sn;
    }
  }
  
  return nllk;
}
