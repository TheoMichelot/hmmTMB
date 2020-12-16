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
};

template<class Type> 
class ZeroInflatedPoisson : public Dist<Type> {
public:
  // Constructor
  ZeroInflatedPoisson() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size()); 
    for (int i = 0; i < n_states; ++i) wpar(i) = log(par(i) / (1.0 - par(i)));
    for (int i = n_states; i < 2 * n_states; ++i) wpar(i) = log(wpar(i)); 
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    for(int i = 0; i < n_states; ++i) par(i, 0) = 1.0 / (1.0 + exp(-wpar(i)));
    for(int i = 0; i < n_states; ++i) par(i, 1) = exp(wpar(i + n_states)); 
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const bool& logpdf) {
    Type val = dzipois(x, par(1), par(0), logpdf);
    return(val); 
  }
};

template<class Type> 
class ZeroTruncatedPoisson : public Dist<Type> {
public:
  // Constructor
  ZeroTruncatedPoisson() {}; 
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
    Type val = dpois(x, par(0)) / (1 - dpois(Type(0), par(0)));
    if (logpdf) val = log(val); 
    return(val); 
  }
};

template<class Type> 
class Binomial : public Dist<Type> {
public:
  // Constructor
  Binomial() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size());
    for (int i = 0; i < n_states; ++i) wpar(i) = par(i); 
    for (int i = n_states; i < 2 * n_states; ++i) wpar(i) = log(par(i) / (1.0 - par(i)));
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    for(int i = 0; i < n_states; ++i) par(i, 0) = wpar(i);
    for(int i = 0; i < n_states; ++i) par(i, 1) = 1.0 / (1.0 + exp(-wpar(i + n_states)));
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const bool& logpdf) {
    Type val = dbinom(x, par(0), par(1), logpdf);
    return(val); 
  }
};


template<class Type> 
class ZeroInflatedBinomial : public Dist<Type> {
public:
  // Constructor
  ZeroInflatedBinomial() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size()); 
    for (int i = 0; i < n_states; ++i) wpar(i) = log(par(i) / (1.0 - par(i)));
    for (int i = n_states; i < 2 * n_states; ++i) wpar(i) = par(i); 
    for (int i = 2 * n_states; i < 3 * n_states; ++i) wpar(i) = log(par(i) / (1.0 - par(i)));
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    for(int i = 0; i < n_states; ++i) par(i, 0) = 1.0 / (1.0 + exp(-wpar(i)));
    for(int i = 0; i < n_states; ++i) par(i, 1) = wpar(i + n_states); 
    for(int i = 0; i < n_states; ++i) par(i, 2) = 1.0 / (1.0 + exp(-wpar(i + 2 * n_states))); 
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const bool& logpdf) {
    Type val;
    if (x == Type(0)) val = par(0) + (1 - par(0)) * dbinom(x, par(1), par(2)); 
    else val = (1 - par(0)) * dbinom(x, par(1), par(2)); 
    if (logpdf) val = log(val); 
    return(val); 
  }
};

template<class Type> 
class ZeroInflatedNegativeBinomial : public Dist<Type> {
public:
  // Constructor
  ZeroInflatedNegativeBinomial() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size()); 
    for (int i = 0; i < n_states; ++i) wpar(i) = log(par(i) / (1.0 - par(i)));
    for (int i = n_states; i < 2 * n_states; ++i) wpar(i) = log(par(i)); 
    for (int i = 2 * n_states; i < 3 * n_states; ++i) wpar(i) = log(par(i) / (1.0 - par(i)));
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    for(int i = 0; i < n_states; ++i) par(i, 0) = 1.0 / (1.0 + exp(-wpar(i)));
    for(int i = 0; i < n_states; ++i) par(i, 1) = exp(wpar(i + n_states)); 
    for(int i = 0; i < n_states; ++i) par(i, 2) = 1.0 / (1.0 + exp(-wpar(i + 2 * n_states))); 
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const bool& logpdf) {
    Type val;
    if (x == Type(0)) val = par(0) + (1 - par(0)) * dnbinom(x, par(1), par(2)); 
    else val = (1 - par(0)) * dnbinom(x, par(1), par(2)); 
    if (logpdf) val = log(val); 
    return(val); 
  }
};

template<class Type> 
class ZeroTruncatedNegativeBinomial : public Dist<Type> {
public:
  // Constructor
  ZeroTruncatedNegativeBinomial() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size());
    for (int i = 0; i < n_states; ++i) wpar(i) = log(par(i)); 
    for (int i = n_states; i < 2 * n_states; ++i) wpar(i) = log(par(i) / (1.0 - par(i)));
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    for(int i = 0; i < n_states; ++i) par(i, 0) = exp(wpar(i));
    for(int i = 0; i < n_states; ++i) par(i, 1) = 1.0 / (1.0 + exp(-wpar(i + n_states)));
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const bool& logpdf) {
    Type val = dnbinom(x, par(0), par(1)) / (1.0 - dnbinom(Type(0.0), par(0), par(1)));
    if (logpdf) val = log(val); 
    return(val); 
  }
};

template<class Type> 
class NegativeBinomial : public Dist<Type> {
public:
  // Constructor
  NegativeBinomial() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size());
    for (int i = 0; i < n_states; ++i) wpar(i) = log(par(i)); 
    for (int i = n_states; i < 2 * n_states; ++i) wpar(i) = log(par(i) / (1.0 - par(i)));
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    for(int i = 0; i < n_states; ++i) par(i, 0) = exp(wpar(i));
    for(int i = 0; i < n_states; ++i) par(i, 1) = 1.0 / (1.0 + exp(-wpar(i + n_states)));
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const bool& logpdf) {
    Type val = dnbinom(x, par(0), par(1), logpdf);
    return(val); 
  }
};

template<class Type> 
class Categorical : public Dist<Type> {
public:
  // Constructor
  Categorical() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    int n_par = par.size() / n_states; 
    matrix<Type> wparmat(n_states, n_par);
    for (int i = 0; i < n_par; ++i) {
      wparmat.col(i) = par.segment(n_states * i, n_states * i + n_par); 
    }
    vector<Type> rowsums = wparmat.rowwise().sum(); 
    vector<Type> wpar(n_states * n_par); 
    for (int j = 1;  j < n_par; ++j) {
      for (int i = 0; i < n_states; ++i) {
        wpar(i + j * n_states) = log(wparmat(i, j) / (1 - rowsums(i))); 
      }
    }
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size() / n_states;
    matrix<Type> par(n_states, n_par);
    vector<Type> ewpar = exp(wpar);
    matrix<Type> wparmat(n_states, n_par);
    for (int i = 0; i < n_par; ++i) {
      wparmat.col(i) = ewpar.segment(n_states * i, n_states * i + n_par); 
    }
    vector<Type> etmp = wparmat.rowwise().sum(); 
    for (int i = 0; i < n_states; ++i) {
      Type s = 1.0 / (1.0 + etmp(i)); 
      for (int j = 0; j < n_par; ++j) {
        par(i, j) = exp(wpar(i + n_states * j)) * s;  
      }
    }
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const bool& logpdf) {
    int obs = int(asDouble(x)); 
    Type val; 
    if (obs == 0) {
      val = 1.0 - par.sum(); 
    } else {
      val = par(obs - 1);
    }
    if (logpdf) val = log(val); 
    return(val); 
  }
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
};

template<class Type> 
class VonMises : public Dist<Type> {
public:
  // Constructor
  VonMises() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size()); 
    // Mean in (-pi, pi]
    for(int i = 0; i < n_states; i++)
      wpar(i) = logit((par(i) + M_PI) / (2 * M_PI));
    // Concentration > 0
    for(int i = n_states; i < 2*n_states; i++)
      wpar(i) = log(par(i));
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    // Mean
    for (int i = 0; i < n_states; i++) 
      par(i, 0) = 2 * M_PI * invlogit(wpar(i)) - M_PI; 
    // Concentration
    for (int i = 0; i < n_states; i++) 
      par(i, 1) = exp(wpar(i + n_states));
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const bool& logpdf) {
    Type b = besselI(Type(par(1)), Type(0));
    Type val = 0;
    if(!logpdf)
      val = 1/(2 * M_PI * b) * exp(par(1) * cos(x - par(0)));
    else
      val = - log(2 * M_PI * b) + par(1) * cos(x - par(0));
    return(val); 
  }
};

#endif
