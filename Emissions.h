//
// Created by Stephanos Tsoucas on 6/14/17.
//

#ifndef BOOSTPYTHONHELLOWORLD_EMISSIONS_H
#define BOOSTPYTHONHELLOWORLD_EMISSIONS_H

#include <armadillo>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdio.h>
//#include <unordered_map>

using namespace arma;

class Emissions {

  const mat emissions_probabilities;
  const mat emissions_probabilities_complement;
  const uword n_bin_states;

  // Memoization table for emissions probabilities
  mutable std::vector<std::vector<double>> emissions_memo;

  uword vec_to_int(const vec &v) const {

    uword ret = 0;
    uword s = v.n_elem;
#pragma omp simd
    for (int i = 0; i < s; i++) {
      ret += (((int)v(i)) << i);
    }

    return ret;
  }

  double emission_probability_from_index(const unsigned int bin_state,
                                         const uword idx_of_observation,
                                         const vec &observation) const;

  double emission_probability_direct(const unsigned int bin_state,
                                     const vec &observation) const {

    const vec &log_obs_prob =
        log((abs(1 - observation) %
             emissions_probabilities_complement.col(bin_state)) +
            (observation % emissions_probabilities.col(bin_state)));

    const double prob = exp(accu(log_obs_prob));
    return prob == NAN ? 0.0 : prob;
  }

public:
  Emissions(mat emissions_probabilities)
      : emissions_probabilities(emissions_probabilities),
        emissions_probabilities_complement(1 - emissions_probabilities),
        n_bin_states(emissions_probabilities.n_cols) {

    const uword n_histone_marks = emissions_probabilities.n_rows;
    const uword total_possibilities = (int)pow(2, n_histone_marks);

    const uword n_bin_states = emissions_probabilities.n_cols;

    for (int i = 0; i <= total_possibilities; i++) {
      emissions_memo.emplace_back(n_bin_states, -1);
    }
  }

  double emission_probability(const unsigned int bin_state,
                              const vec &observation) const;

  colvec emission_probabilities_col(const vec &observation) const;

  rowvec emission_probabilities_row(const vec &observation) const;

  const mat &emissions_probabilities_m() const {
    return emissions_probabilities;
  }
};

Emissions emissions_from_states(
    const unsigned int n_bin_states, const int n_histone_marks,
    const std::vector<unsigned int> &bins_to_bin_states, const mat &bin_data);

#endif // BOOSTPYTHONHELLOWORLD_EMISSIONS_H
