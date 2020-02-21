
#ifndef _DIST_
#define _DIST_

#include <TMB.hpp>

//' Observation distribution class
//' 
//' This class comprises the link and inverse link functions for
//' the parameters of the distribution, as well as its probability
//' density (or mass) function.
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

//' Link function for distribution parameters
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

// Inverse link function for distribution parameters
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

// Probability density/mass function of the distribution
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

#endif
