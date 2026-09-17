#ifndef _HMMTMB_LIKELIHOOD_
#define _HMMTMB_LIKELIHOOD_

#include <algorithm>

//' Banded forward algorithm
//'
//' The exact forward algorithm accumulates log-likelihood contributions that
//' are each conditional on the entire process history, through the scaled
//' forward variable. The Hessian with respect to anything entering those
//' contributions time point by time point -- the weights of a latent field,
//' say -- is therefore dense, which the automatic Laplace approximation cannot
//' afford.
//'
//' The banded version (Fischer, 2026) splits the series into blocks of length
//' bw and starts each block's forward variable from a fixed vector rho at the
//' beginning of the *previous* block. Each block then sees only the 2 * bw
//' observations nearest to it, cross-derivatives beyond that lag vanish, and
//' the Hessian is banded. Every block is traversed twice, so this costs about
//' twice the exact algorithm, and the error decays geometrically in bw because
//' an ergodic Markov chain forgets its initial condition exponentially fast.
//' This mirrors LaMa::forward_g(), so that the two can be checked against
//' each other.
//'
//' The exact algorithm stays where it was, in hmmTMB.hpp, and is what runs
//' when bw < 2. This function assumes bw >= 2.
//'
//' @param prob n x n_states matrix of state-dependent probabilities
//' @param tpm_array Transition probability matrices; entry i governs the
//' transition out of time step i
//' @param delta0 n_ID x n_states matrix of initial distributions
//' @param ID Vector of time series IDs, in blocks of equal value
//' @param bw Bandwidth, at least 2
//' @return Log-likelihood
template<class Type>
Type forward_alg_banded(const matrix<Type>& prob,
                        const vector<matrix<Type> >& tpm_array,
                        const matrix<Type>& delta0,
                        const vector<Type>& ID,
                        int bw) {
  int n = prob.rows();
  matrix<Type> phi(1, prob.cols()), rho(1, prob.cols());
  rho.setConstant(Type(1) / Type(prob.cols()));
  Type llk = 0;

  // One time series at a time
  int id = 0;
  for(int start = 0; start < n; id++) {
    int end = start;
    while(end + 1 < n && ID(end + 1) == ID(end)) end++;

    // First block, from the model's own initial distribution
    phi = (delta0.row(id).array() * prob.row(start).array()).matrix();
    llk += log(phi.sum());
    phi /= phi.sum();
    for(int i = start + 1; i <= std::min(start + bw - 1, end); i++) {
      phi = ((phi * tpm_array(i - 1)).array() * prob.row(i).array()).matrix();
      llk += log(phi.sum());
      phi /= phi.sum();
    }

    // Later blocks, each warmed up from rho over the preceding block
    for(int b = start + bw; b <= end; b += bw) {
      phi = (rho.array() * prob.row(b - bw).array()).matrix();
      phi /= phi.sum();
      for(int i = b - bw + 1; i < b; i++) {
        phi = ((phi * tpm_array(i - 1)).array() * prob.row(i).array()).matrix();
        phi /= phi.sum();
      }
      for(int i = b; i <= std::min(b + bw - 1, end); i++) {
        phi = ((phi * tpm_array(i - 1)).array() * prob.row(i).array()).matrix();
        llk += log(phi.sum());
        phi /= phi.sum();
      }
    }

    start = end + 1;
  }

  return llk;
}

#endif
