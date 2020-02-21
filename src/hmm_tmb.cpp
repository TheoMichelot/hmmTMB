
#include <TMB.hpp>

template <class Type>
class Dist {
  std::string name;
  
public:
  // Constructor
  Dist(std::string distname) {
    name = distname;
  };
  // Link function
  vector<Type> link(vector<Type> par, int n_states);
  // Inverse link function
  matrix<Type> invlink(vector<Type> wpar, int n_states);
  // Probability density/mass function
  Type pdf(Type x, vector<Type> par, bool logpdf);
};

template <class Type>
vector<Type> Dist<Type>::link(vector<Type> par, int n_states) {
  vector<Type> wpar(par.size());
  
  if(name == "pois") { // Poisson distribution
    wpar = log(par);
  } else if(name == "norm") { // Normal distribution
    // Mean
    for(int i = 0; i < n_states; i++)
      wpar(i) = par(i);
    
    // Standard deviation
    for(int i = n_states; i < 2*n_states; i++)
      wpar(i) = log(par(i));
  }
  
  return wpar;
}

template <class Type>
matrix<Type> Dist<Type>::invlink(vector<Type> wpar, int n_states) {
  // Number of parameters
  int n_par = wpar.size()/n_states;
  
  // Matrix of parameters
  matrix<Type> par(n_states, n_par);
  
  if(name == "pois") { // Poisson distribution
    for(int i = 0; i < n_states; i++) 
      par(i, 0) = exp(wpar(i));
  } else if(name == "norm") { // Normal distribution
    // Mean
    for(int i = 0; i < n_states; i++)
      par(i, 0) = wpar(i);
    
    // Standard deviation
    for(int i = 0; i < n_states; i++)
      par(i, 1) = exp(wpar(i + n_states));
  }
  
  return par;
}

template <class Type>
Type Dist<Type>::pdf(Type x, vector<Type> par, bool logpdf) {
  Type val = 0;
  
  if(name == "pois") { // Poisson distribution
    val = dpois(x, par(0), logpdf);
  } else if(name == "norm") { // Normal distribution
    val = dnorm(x, par(0), par(1), logpdf);
  }
  
  return val;
}

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
  DATA_VECTOR(data); // data stream
  DATA_INTEGER(n_states); // number of states
  DATA_STRING(distname); // name of observation distribution
  
  // PARAMETERS
  PARAMETER_VECTOR(ltpm);
  PARAMETER_VECTOR(wpar);
  
  // Define observation distribution
  Dist <Type> obsdist(distname);
  
  //======================//
  // Transform parameters //
  //======================//
  // Observation parameters
  matrix<Type> par = obsdist.invlink(wpar, n_states);

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
  // Initialise matrix of probabilities
  int n = data.rows();
  matrix<Type> prob(n, n_states);
  
  // Loop over states (columns)
  for (int s = 0; s < n_states; ++s) {
    // Vector of parameters for state s
    vector<Type> subpar = par.row(s);
    
    // Loop over observations (rows)
    for (int i = 0; i < n; ++i) {
      prob(i, s) = obsdist.pdf(data(i), subpar, false);
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
