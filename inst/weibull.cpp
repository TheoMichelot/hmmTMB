
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
    for(int i = 0; i < n_states; ++i) par(i, 0) = exp(wpar(i));
    for(int i = 0; i < n_states; ++i) par(i, 1) = exp(wpar(i + n_states));
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const bool& logpdf) {
    Type val = dweibull(x, par(0), par(1), logpdf);
    return(val); 
  }
  // Number of parameters 
  int npar() { return 2; }
};


