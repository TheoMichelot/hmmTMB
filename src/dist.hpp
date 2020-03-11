
#ifndef _DIST_
#define _DIST_

#include <TMB.hpp>

// Include file for custom distribution if existing
#if __has_include("custom.hpp")
# include "custom.hpp"
#endif

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
  Dist(int distcode) {
    vector<std::string> distname(6);
    distname(0) = "pois";
    distname(1) = "norm";
    distname(2) = "gamma";
    distname(3) = "beta";
    distname(5) = "custom";
    name = distname(distcode);
  };
  // Link function
  vector<Type> link(vector<Type> par, int n_states);
  // Inverse link function
  matrix<Type> invlink(vector<Type> wpar, int n_states);
  // Probability density/mass function
  Type pdf(Type x, vector<Type> par, bool logpdf);
  // Number of parameters
  int npar();
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
  } else if(name == "gamma") { // Gamma distribution
    wpar = log(par);
  } else if(name == "beta") { // Beta distribution
    wpar = log(par);
  } else if(name == "custom") { // User-defined distribution
    wpar = custom_link(par, n_states);
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
  } else if(name == "gamma") { // Gamma distribution
    // Shape
    for(int i = 0; i < n_states; i++)
      par(i, 0) = exp(wpar(i));
    
    // Scale
    for(int i = 0; i < n_states; i++)
      par(i, 1) = exp(wpar(i + n_states));
  } else if(name == "beta") { // Beta distribution
    // Shape 1
    for(int i = 0; i < n_states; i++)
      par(i, 0) = exp(wpar(i));
    
    // Shape 2
    for(int i = 0; i < n_states; i++)
      par(i, 1) = exp(wpar(i + n_states));
  } else if(name == "custom") { // User-defined distribution
    par = custom_invlink(wpar, n_states);
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
  } else if(name == "gamma") { // Gamma distribution
    val = dgamma(x, par(0), par(1), logpdf);
  } else if(name == "beta") { // Beta distribution
    val = dbeta(x, par(0), par(1), logpdf);
  } else if(name == "custom") { // User-defined distribution
    val = custom_pdf(x, par, logpdf);
  }
  
  return val;
}

// Number of distribution parameters
template <class Type>
int Dist<Type>::npar() {
  int n = 0;
  
  if(name == "pois") { // Poisson distribution
    n = 1;
  } else if(name == "norm") { // Normal distribution
    n = 2;
  } else if(name == "gamma") { // Gamma distribution
    n = 2;
  } else if(name == "beta") { // Beta distribution
    n = 2;
  } else if(name == "custom") { // User-defined distribution
    n = custom_npar();
  }
  
  return n;
}

#endif
