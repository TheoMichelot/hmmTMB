
//====================//
// gamma distribution //
//====================//

//' Link function for distribution parameters
template <class Type>
vector<Type> custom_link(vector<Type> par, int n_states) {
  vector<Type> wpar(par.size());
  wpar = log(par);
  return wpar;
}

// Inverse link function for distribution parameters
template <class Type>
matrix<Type> custom_invlink(vector<Type> wpar, int n_states) {
  // Number of parameters
  int n_par = wpar.size()/n_states;
  
  // Matrix of parameters
  matrix<Type> par(n_states, n_par);
  
  // Shape
  for(int i = 0; i < n_states; i++)
    par(i, 0) = exp(wpar(i));
  
  // Scale
  for(int i = 0; i < n_states; i++)
    par(i, 1) = exp(wpar(i + n_states));
  
  return par;
}

// Probability density/mass function of the distribution
template <class Type>
Type custom_pdf(Type x, vector<Type> par, bool logpdf) {
  Type val = 0;
  val = dgamma(x, par(0), par(1), logpdf);
  return val;
}
