
#include "dist.hpp"

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
//'   wpar_fe: observation parameters (fixed effects)
//'   wpar_re: observation parameters (random effects)
//'   log_lambda: vector of smoothness parameters

template<class Type>
Type objective_function<Type>::operator() ()
{
  // DATA
  DATA_VECTOR(ID); // vector of time series IDs
  DATA_MATRIX(data); // data stream
  DATA_SPARSE_MATRIX(X_fe); // design matrix for fixed effects
  DATA_SPARSE_MATRIX(X_re); // design matrix for random effects
  DATA_SPARSE_MATRIX(S); // Penalty matrix
  DATA_IVECTOR(ncol_re); // number of columns of S and X_re for each random effect
  DATA_INTEGER(n_states); // number of states
  DATA_IVECTOR(distcode); // codes of observation distributions
  
  // PARAMETERS
  PARAMETER_VECTOR(ltpm); // transition probabilities
  PARAMETER_VECTOR(wpar_fe); // observation parameters (fixed effects)
  PARAMETER_VECTOR(wpar_re); // observation parameters (random effects)
  PARAMETER_VECTOR(log_lambda); // smoothness parameters
  
  // Number of observed variables
  int n_var = distcode.size();
  // Number of data rows
  int n = data.rows();
  
  //======================//
  // Transform parameters //
  //======================//
  // Observation parameters
  vector<Type> par_vec = X_fe * wpar_fe + X_re * wpar_re;
  matrix<Type> par_mat(n, par_vec.size()/n);
  for(int i = 0; i < par_mat.cols(); i++) {
    // Matrix with one row for each time step and
    // one column for each parameter
    par_mat.col(i) = par_vec.segment(i*n, n);
  }
  
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
    tpm.row(i) = tpm.row(i)/tpm.row(i).sum();
  }
  
  //=================================//
  // Compute stationary distribution //
  //=================================//
  matrix<Type> delta(1, n_states);
  matrix<Type> I = matrix<Type>::Identity(n_states, n_states);
  matrix<Type> tpminv = I;
  tpminv = tpminv - tpm;
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
    
    // Loop over observations (rows)
    for (int i = 0; i < n; ++i) {
      // Subset and transform observation parameters
      vector<Type> sub_wpar = par_mat.row(i).segment(par_count, obsdist.npar() * n_states);
      matrix<Type> par = obsdist.invlink(sub_wpar, n_states);
      
      // Loop over states (columns)
      for (int s = 0; s < n_states; ++s) {
        // Vector of parameters for state s
        vector<Type> subpar = par.row(s);
        
        prob(i, s) = prob(i, s) * obsdist.pdf(data(i, var), subpar, false);
      }
    }
    
    par_count = par_count + obsdist.npar() * n_states;
  }
  
  //========================//
  // Compute log-likelihood //
  //========================//
  // Initialise log-likelihood
  Type llk = 0;
  matrix<Type> phi(delta);
  Type sumphi = 0;
  
  // Forward algorithm
  for (int i = 0; i < n; ++i) {
    // Re-initialise phi at first observation of each time series
    if(i == 0 || ID(i-1) != ID(i)) {
      phi = delta;
    }
    phi = (phi.array() * prob.row(i).array()).matrix();
    phi = phi * tpm;
    sumphi = phi.sum();
    llk = llk + log(sumphi);
    phi = phi / sumphi;
  }
  
  // Negative log-likelihood
  Type nllk = -llk;
  
  //===================//
  // Smoothing penalty //
  //===================//
  // Are there smooths?
  if(ncol_re(0) > 0) {
    // Index in matrix S
    int S_start = 0;
    
    // Loop over smooths
    for(int i = 0; i < ncol_re.size(); i++) {
      // Size of penalty matrix for this smooth
      int Sn = ncol_re(i);
      
      // Penalty matrix for this smooth
      Eigen::SparseMatrix<Type> this_S = S.block(S_start, S_start, Sn, Sn);
      
      // Coefficients for this smooth
      vector<Type> this_wpar_re = wpar_re.segment(S_start, Sn);
      
      // Add penalty
      nllk = nllk -
        Type(0.5) * Sn * log_lambda(i) +
        Type(0.5) * exp(log_lambda(i)) * density::GMRF(this_S).Quadform(this_wpar_re);    
      
      // Increase index
      S_start = S_start + Sn;
    }
  }
  
  return nllk;
}
