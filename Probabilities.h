//
// Created by Stephanos Tsoucas on 6/14/17.
//

#ifndef BOOSTPYTHONHELLOWORLD_PROBABILITIES_H
#define BOOSTPYTHONHELLOWORLD_PROBABILITIES_H

#include "Model.h"

using namespace arma;

// TODO change type of bin_data
struct Probabilities {

  const mat &bin_data;

  std::vector<double> bin_dec_data;

  // forward_probabilities(bin_state, domain_state, bin_position)
  cube &forward_probabilities;
  // backward_probabilities(bin_state, domain_state, bin_position)
  cube &backward_probabilities;
  std::vector<double> scaling_factors;
  const uword n_bin_positions;
  const uword n_histone_marks;

  double posterior_probability(unsigned long bin_position,
                               unsigned int bin_state,
                               unsigned int domain_state) const {
    return backward_probabilities(domain_state, bin_state, bin_position) *
           forward_probabilities(domain_state, bin_state, bin_position);
  }

  double log_likelihood() const {
    double log_lik = 0.0;
    for (double s : scaling_factors) {
      log_lik += log(s);
    }

    return log_lik;
  }

  bool operator<(const Probabilities &other) {
    return n_bin_positions < other.n_bin_positions;
  }

  Probabilities(const mat &bin_data, cube &forward_probabilities,
                cube &backward_probabilities,
                std::vector<double> scaling_factors)
      : bin_data(bin_data), forward_probabilities(forward_probabilities),
        backward_probabilities(backward_probabilities),
        scaling_factors(scaling_factors), n_bin_positions(bin_data.n_cols),
        n_histone_marks(bin_data.n_rows) {

    for (uword i = 0; i < bin_data.n_cols; i++) {
      // vec to int
      bin_dec_data.push_back(
          (double)Emissions_Dec::vec_to_int(bin_data.col(i)));
    }
  }
};

#endif // BOOSTPYTHONHELLOWORLD_PROBABILITIES_H
