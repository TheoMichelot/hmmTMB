#ifndef _DISTCODES_
#define _DISTCODES_

#include<vector>
#include<string>
#include<memory>
#include "dist_def.hpp"
#include "added_dists.hpp"

// Function that takes a distribution code and produces a pointer to an object
// for that distribution class 
template <class Type>
std::unique_ptr<Dist<Type>> dist_generator(const int& code) {
  switch (code) {
  case 0: 
    return(std::unique_ptr<Dist<Type>>(new Poisson<Type>)); 
  case 1: 
    return(std::unique_ptr<Dist<Type>>(new Normal<Type>)); 
  case 2: 
    return(std::unique_ptr<Dist<Type>>(new Gamma<Type>)); 
  case 3: 
    return(std::unique_ptr<Dist<Type>>(new Beta<Type>)); 
  case 4:
    return(std::unique_ptr<Dist<Type>>(new VonMises<Type>));
  case 5: 
    return(std::unique_ptr<Dist<Type>>(new Binomial<Type>));
  case 6: 
    return(std::unique_ptr<Dist<Type>>(new NegativeBinomial<Type>));
  case 7: 
    return(std::unique_ptr<Dist<Type>>(new ZeroInflatedPoisson<Type>));
  case 8: 
    return(std::unique_ptr<Dist<Type>>(new ZeroInflatedBinomial<Type>));
  case 9: 
    return(std::unique_ptr<Dist<Type>>(new ZeroInflatedNegativeBinomial<Type>));
  case 10: 
    return(std::unique_ptr<Dist<Type>>(new ZeroTruncatedPoisson<Type>));
  case 11: 
    return(std::unique_ptr<Dist<Type>>(new ZeroTruncatedNegativeBinomial<Type>));
  case 12: 
    return(std::unique_ptr<Dist<Type>>(new Categorical<Type>)); 
  case 13: 
    return(std::unique_ptr<Dist<Type>>(new Exponential<Type>)); 
  case 14: 
    return(std::unique_ptr<Dist<Type>>(new LogNormal<Type>)); 
  case 15: 
    return(std::unique_ptr<Dist<Type>>(new Weibull<Type>)); 
  case 16: 
    return(std::unique_ptr<Dist<Type>>(new TruncatedNormal<Type>)); 
  case 17: 
    return(std::unique_ptr<Dist<Type>>(new FoldedNormal<Type>)); 
  case 18: 
    return(std::unique_ptr<Dist<Type>>(new Studentst<Type>)); 
  case 19: 
    return(std::unique_ptr<Dist<Type>>(new Tweedie<Type>)); 
  case 20: 
    return(std::unique_ptr<Dist<Type>>(new MultivariateNormal<Type>)); 
  case 21: 
    return(std::unique_ptr<Dist<Type>>(new Dirichlet<Type>));
  case 22: 
    return(std::unique_ptr<Dist<Type>>(new AltGamma<Type>));
  default: 
    return(std::unique_ptr<Dist<Type>>(new Poisson<Type>)); 
  }
}

#endif
