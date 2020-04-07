#ifndef _DIST_
#define _DIST_

// Defines the abstract basic distribution Type and then a number of derived
// distribution types, e.g., Poisson, Normal, and Gamma. 
// Each distribution has a link function, an inverse link function, a pdf, and a 
// specified number of parameters. 

#ifndef _HMMTMB_
#define _HMMTMB_
#include <TMB.hpp>
#endif 
 
template <class Type>
class Dist {
public:
  // Constructor
  Dist() {};
  // Link function
  virtual vector<Type> link(const vector<Type>& par, const int& n_states) = 0; 
  // Inverse link function
  virtual matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) = 0; 
  // Probability density/mass function
  virtual Type pdf(const Type& x, const vector<Type>& par, const bool& logpdf) = 0; 
  // Number of parameters
  virtual int npar() = 0; 
};

template<class Type> 
class Poisson : public Dist<Type> {
public:
  // Constructor
  Poisson() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size()); 
    wpar = log(par); 
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    for(int i = 0; i < n_states; ++i) par(i, 0) = exp(wpar(i));
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const bool& logpdf) {
    Type val = dpois(x, par(0), logpdf);
    return(val); 
  }
  // Number of parameters 
  int npar() { return 1; }
};

template<class Type> 
class Normal : public Dist<Type> {
public:
  // Constructor
  Normal() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size()); 
    for (int i = 0; i < n_states; ++i) wpar(i) = par(i);
    for (int i = n_states; i < 2 * n_states; ++i) wpar(i) = log(par(i)); 
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    for (int i = 0; i < n_states; ++i) par(i, 0) = wpar(i); // mean 
    for (int i = 0; i < n_states; ++i) par(i, 1) = exp(wpar(i + n_states)); // sd
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const bool& logpdf) {
    Type val = dnorm(x, par(0), par(1), logpdf);
    return(val); 
  }
  // Number of parameters 
  int npar() { return 2; }
};

template<class Type> 
class Gamma : public Dist<Type> {
public:
  // Constructor
  Gamma() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size()); 
    wpar = log(par); 
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    for (int i = 0; i < n_states; ++i) par(i, 0) = exp(wpar(i)); // shape 
    for (int i = 0; i < n_states; ++i) par(i, 1) = exp(wpar(i + n_states)); // scale 
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const bool& logpdf) {
    Type val = dgamma(x, par(0), par(1), logpdf);
    return(val); 
  }
  // Number of parameters 
  int npar() { return 2; }
};

template<class Type> 
class Beta : public Dist<Type> {
public:
  // Constructor
  Beta() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size()); 
    wpar = log(par); 
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    for (int i = 0; i < n_states; ++i) par(i, 0) = exp(wpar(i)); // shape 
    for (int i = 0; i < n_states; ++i) par(i, 1) = exp(wpar(i + n_states)); // scale 
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const bool& logpdf) {
    Type val = dbeta(x, par(0), par(1), logpdf);
    return(val); 
  }
  // Number of parameters 
  int npar() { return 2; }
};

#endif
