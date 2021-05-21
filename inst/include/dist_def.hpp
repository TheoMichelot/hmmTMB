#ifndef _DIST_
#define _DIST_

// Defines the abstract basic distribution Type and then a number of derived
// distribution types, e.g., Poisson, Normal, and Gamma. 
// Each distribution has a link function, an inverse link function, and a pdf.  

#ifndef _HMMTMB_
#define _HMMTMB_
#include <TMB.hpp>
#endif 

// Abstract Distribution Class 
template <class Type>
class Dist {
public:
  // Constructor
  Dist() {};
  // Link function
  virtual vector<Type> link(const vector<Type>& par, const int& n_states) = 0; 
  // Inverse link function
  virtual matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) = 0;
  
  // Both a univariate and vector input can be given to a pdf() function. 
  // Derived distributions can either have both options or only one. 
  
  // Probability density/mass function
  virtual Type pdf(const Type& x, const vector<Type>& par, const bool& logpdf) {
    return(0.0); 
  };
  // Vector input Probability density/mass function
  virtual Type pdf(const vector<Type>& x, const vector<Type>& par, const bool& logpdf) {
    return(0.0); 
  }
};

// DISCRETE DISTRIBUTIONS ----------------------

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

// CONTINUOUS DISTRIBUTIONS --------------------

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
class TruncatedNormal : public Dist<Type> {
public:
  // Constructor
  TruncatedNormal() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size()); 
    for (int i = 0; i < n_states; ++i) wpar(i) = par(i);
    for (int i = n_states; i < 2 * n_states; ++i) wpar(i) = log(par(i)); 
    for (int i = n_states * 2; i < 3 * n_states; ++i) wpar(i) = par(i);
    for (int i = n_states * 3; i < 4 * n_states; ++i) wpar(i) = par(i);
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    for (int i = 0; i < n_states; ++i) par(i, 0) = wpar(i); // mean 
    for (int i = 0; i < n_states; ++i) par(i, 1) = exp(wpar(i + n_states)); // sd
    for (int i = 0; i < n_states; ++i) par(i, 2) = wpar(i + 2 * n_states); // min
    for (int i = 0; i < n_states; ++i) par(i, 3) = wpar(i + 3 * n_states); // max
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const bool& logpdf) {
    Type left = pnorm(par(2), par(0), par(1)); 
    Type right = pnorm(par(3), par(0), par(1)); 
    Type val = dnorm(x, par(0), par(1)) / (right - left);
    if (logpdf) val = log(val); 
    return(val); 
  }
};

template<class Type> 
class FoldedNormal : public Dist<Type> {
public:
  // Constructor
  FoldedNormal() {}; 
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
    Type val = dnorm(x, par(0), par(1)) + dnorm(-x, par(0), par(1));
    if (logpdf) val = log(val); 
    return(val); 
  }
};

template<class Type> 
class Studentst : public Dist<Type> {
public:
  // Constructor
  Studentst() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size()); 
    for (int i = 0; i < n_states; ++i) wpar(i) = par(i);
    for (int i = n_states; i < 2 * n_states; ++i) wpar(i) = log(par(i)); ; 
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    for (int i = 0; i < n_states; ++i) par(i, 0) = wpar(i); // mean 
    for (int i = 0; i < n_states; ++i) par(i, 1) = exp(wpar(i + n_states)); // scale
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const bool& logpdf) {
    Type df = 2 * par(1) * par(1) / (par(1) * par(1) - 1); 
    Type val = dt(x - par(0), df, 0); 
    if (logpdf) val = log(val); 
    return(val); 
  }
};

template<class Type> 
class LogNormal : public Dist<Type> {
public:
  // Constructor
  LogNormal() {}; 
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
    for (int i = 0; i < n_states; ++i) par(i, 0) = wpar(i); // lmean 
    for (int i = 0; i < n_states; ++i) par(i, 1) = exp(wpar(i + n_states)); // lsd
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const bool& logpdf) {
    Type val = dnorm(log(x), par(0), par(1), logpdf) / x;
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
class AltGamma : public Dist<Type> {
public:
  // Constructor
  AltGamma() {}; 
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
    for (int i = 0; i < n_states; ++i) par(i, 0) = exp(wpar(i)); // mean 
    for (int i = 0; i < n_states; ++i) par(i, 1) = exp(wpar(i + n_states)); // sd 
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const bool& logpdf) {
    Type scale = par(1) * par(1) / par(0);
    Type shape = par(0) / scale; 
    Type val = dgamma(x, shape, scale, logpdf);
    return(val); 
  }
};

template<class Type> 
class Weibull : public Dist<Type> {
public:
  // Constructor
  Weibull() {}; 
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
    Type val = dweibull(x, par(0), par(1), logpdf);
    return(val); 
  }
};

template<class Type> 
class Exponential : public Dist<Type> {
public:
  // Constructor
  Exponential() {}; 
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
    for (int i = 0; i < n_states; ++i) par(i, 0) = exp(wpar(i)); // rate 
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const bool& logpdf) {
    Type val = dexp(x, par(0), logpdf);
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

// MIXED DISTRIBUTIONS -------------------------

template<class Type> 
class Tweedie : public Dist<Type> {
public:
  // Constructor
  Tweedie() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size()); 
    for (int i = 0; i < n_states; ++i) wpar(i) = par(i);
    for (int i = n_states; i < 2 * n_states; ++i) wpar(i) = log(par(i) / (1 - par(i))); 
    for (int i = 2 * n_states; i < 3 * n_states; ++i) wpar(i) = log(par(i)); 
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    for (int i = 0; i < n_states; ++i) par(i, 0) = wpar(i); // mean 
    for (int i = 0; i < n_states; ++i) par(i, 1) = 1 / (1 + exp(-wpar(i + n_states))); // power
    for (int i = 0; i < n_states; ++i) par(i, 2) = exp(wpar(i + 2 * n_states)); // dispersion 
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const bool& logpdf) {
    Type val = dtweedie(x, par(0), par(2), par(1) + 1.0, logpdf);
    return(val); 
  }
};

// ANGULAR DISTRIBUTIONS -----------------------
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

template<class Type> 
class WrpCauchy : public Dist<Type> {
public:
  // Constructor
  WrpCauchy() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size()); 
    // Mean in (-pi, pi]
    for(int i = 0; i < n_states; i++)
      wpar(i) = logit((par(i) + M_PI) / (2 * M_PI));
    // 0 < Concentration < 1
    for(int i = n_states; i < 2*n_states; i++)
      wpar(i) = log(par(i) / (1.0 - par(i)));
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
      par(i, 1) = 1.0 / (1.0 + exp(-wpar(i + n_states)));
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const bool& logpdf) {
    Type val = (1 - par(1)*par(1)) / (2 * M_PI * (1 + par(1)*par(1) - 2 * par(1) * cos(x - par(0)))); 
    if(logpdf) val = log(val); 
    return(val); 
  }
};

// MULTIVARIATE DISTRIBUTIONS ------------------

template<class Type> 
class MultivariateNormal : public Dist<Type> {
public:
  // Constructor
  MultivariateNormal() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size());
    int n_par = wpar.size()/n_states;
    int dim = this->dim(n_par); 
    int k = 0; 
    // means
    // means
    for (int d = 0; d < dim; ++d) {
      for (int i = 0; i < n_states; ++i) {
        wpar(k) = par(k); 
        ++k; 
      }
    }
    // sds 
    for (int d = 0; d < dim; ++d) {
      for (int i = 0; i < n_states; ++i) {
        wpar(k) = log(par(k)); 
        ++k; 
      }
    }
    // corr 
    for (int d = 0; d < 0.5*(dim * dim - dim); ++d) {
      for (int i = 0; i < n_states; ++i) {
       Type tmp = 0.5 * (par(k) + 1); 
       wpar(k) = log(tmp / (1 - tmp)); 
        ++k; 
      }
    }
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    int dim = this->dim(n_par);
    int k = 0; 
    // means
    for (int d = 0; d < dim; ++d) {
      for (int i = 0; i < n_states; ++i) {
        par(i, d) = wpar(k); 
        ++k; 
      }
    }
    // sds 
    for (int d = 0; d < dim; ++d) {
      for (int i = 0; i < n_states; ++i) {
        par(i, d + dim) = exp(wpar(k)); 
        ++k; 
      }
    }
    // corr 
    for (int d = 0; d < 0.5*(dim * dim - dim); ++d) {
      for (int i = 0; i < n_states; ++i) {
        par(i, d + 2 * dim) = 1 / (1 + exp(-wpar(k))); 
        par(i, d + 2 * dim) = 2 * par(i, d + 2 * dim) - 1; 
        ++k; 
      }
    }
    return(par); 
  }
  // Univariate Probability density function
  Type pdf(const Type& x, const vector<Type>& par, const bool& logpdf) {
    Type val = dnorm(x, par(0), par(1), logpdf); 
    return(val); 
  }
  
  // Multivariate Probability density function 
  Type pdf(const vector<Type>& x, const vector<Type>& par, const bool& logpdf) {
    int dim = this->dim(par.size()); 
    vector<Type> y(dim);
    for (int i = 0; i < dim; ++i) y(i) = x(i) - par(i); 
    vector<Type> sds(dim); 
    for (int i = 0; i < dim; ++i) sds(i) = par(i + dim); 
    vector<Type> corr((dim * dim - dim) / 2); 
    for (int i = 0; i < (dim * dim - dim) / 2; ++i) corr(i) = par(i + 2 * dim);
    matrix<Type> Sigma(dim, dim); 
    int k = 0; 
    for (int i = 0; i < dim; ++i) {
      for (int j = i; j < dim; ++j) {
        Sigma(j, i) = sds(j) * sds(i); 
        if (i != j) {
          Sigma(j, i) *= corr(k); 
          ++k; 
        }
      }
    }
    for (int i = 0; i < dim; ++i) {
      for (int j = 0; j < i; ++j) {
        Sigma(j, i) = Sigma(i, j); 
      }
    }
    Type val = density::MVNORM(Sigma)(y); 
    val = -val; 
    if (!logpdf) val = exp(val); 
    return(val); 
  }
  
  // Solve for dimension 
  int dim(const double& l) {
    double a = 1; 
    double b = 3; 
    double c = - 2 * l; 
    double d = b * b - 4 * a * c; 
    double root = (-b + sqrt(d)) / (2 * a); 
    return(int(root)); 
  }
};

template<class Type> 
class Dirichlet : public Dist<Type> {
public:
  // Constructor
  Dirichlet() {}; 
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
    for (int i = 0; i < n_par; ++i) {
      for (int j = 0; j < n_states; ++j) {
        par(j, i) = exp(wpar(i * n_states + j)); 
      }
    }
    return(par); 
  }

    // Multivariate Probability density function 
  Type pdf(const vector<Type>& x, const vector<Type>& par, const bool& logpdf) {
    Type val = 0; 
    for (int i = 0; i < x.size(); ++i) {
      val += (par(i) - 1) * log(x(i)); 
      val -= lgamma(par(i)); 
    }
    val += lgamma(par.sum()); 
    if (!logpdf) val = exp(val); 
    return(val); 
  }

};
#endif
