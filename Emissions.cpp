//
// Created by Stephanos Tsoucas on 6/14/17.
//

#include "Emissions.h"

double
Emissions::emission_probability_from_index(const unsigned int bin_state,
                                           const uword idx_of_observation,
                                           const vec &observation) const {

  double val = emissions_memo[idx_of_observation][bin_state];

  if (val != -1) {
    return val;
  } else {
    double ret = emission_probability_direct(bin_state, observation);
    emissions_memo[idx_of_observation][bin_state] = ret;
    return ret;
  }
}

double Emissions::emission_probability(const unsigned int bin_state,
                                       const vec &observation) const {

  uword idx = vec_to_int(observation);

  return emission_probability_from_index(bin_state, idx, observation);
}

colvec Emissions::emission_probabilities_col(const vec &observation) const {

  uword idx = vec_to_int(observation);

  colvec emissions(n_bin_states);

  for (int i = 0; i < n_bin_states; i++) {
    emissions(i) = emission_probability_from_index(i, idx, observation);
  }

  return emissions;
}

rowvec Emissions::emission_probabilities_row(const vec &observation) const {

  uword idx = vec_to_int(observation);

  rowvec emissions(n_bin_states);

  for (int i = 0; i < n_bin_states; i++) {
    emissions(i) = emission_probability_from_index(i, idx, observation);
  }

  return emissions;
}

// TODO update this separately for multiple chromosomes (?) may not need to,
// since initial model  is gotten using concatenation of all bin data.  Move this
// into emissions_dec
Emissions emissions_from_states(
    const unsigned int n_bin_states, const int n_histone_marks,
    const std::vector<unsigned int> &bins_to_bin_states, const mat &bin_data) {

  std::set<std::vector<unsigned int>> observations_seen;

  mat em = zeros(n_histone_marks, n_bin_states);
  std::vector<double> state_counts(n_bin_states, 0);

  for (int bin_position = 0; bin_position < bins_to_bin_states.size();
       bin_position++) {
    const int bin_state = bins_to_bin_states.at(bin_position);

    state_counts.at(bin_state) += 1;

    const vec &obs_vec = bin_data.col(bin_position);

    std::vector<unsigned int> observation =
        conv_to<std::vector<unsigned int>>::from(obs_vec);
    observations_seen.insert(observation);

    em.col(bin_state) += obs_vec;
  }

  for (int bin_state = 0; bin_state < n_bin_states; bin_state++) {
    em.col(bin_state) /= state_counts.at(bin_state);
  }

  em.replace(NAN, 0.0);

  const double emissions_kick = 0.1;

  em *= (1 - emissions_kick);
  em += emissions_kick;

  return Emissions(em);
}
