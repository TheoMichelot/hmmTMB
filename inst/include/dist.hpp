#ifndef _DISTCODES_
#define _DISTCODES_

#include<vector>
#include<string>
#include<memory>
#include "dist_def.hpp"
#include "added_dists.hpp"

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
  default: 
    return(std::unique_ptr<Dist<Type>>(new Poisson<Type>)); 
  }
}

#endif
