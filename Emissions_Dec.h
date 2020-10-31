#ifndef BOOSTPYTHONHELLOWORLD_EMISSIONS_DEC_H
#define BOOSTPYTHONHELLOWORLD_EMISSIONS_DEC_H

#include <armadillo>
#include <set>

using namespace arma;

class Emissions_Dec {

  const mat emissions_probabilities;
  const mat emissions_probabilities_complement;
  const uword n_bin_states;
  mutable std::vector<std::vector<double>> emissions_memo;

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
  Emissions_Dec(mat emissions_probabilities)
      : emissions_probabilities(emissions_probabilities),
        emissions_probabilities_complement(1 - emissions_probabilities),
        n_bin_states(emissions_probabilities.n_cols) {

    const uword n_histone_marks = emissions_probabilities.n_rows;
    const uword total_possibilities = (int)pow(2, n_histone_marks);

    const uword n_bin_states = emissions_probabilities.n_cols;

    for (int i = 0; i <= total_possibilities; i++) {
      emissions_memo.emplace_back(n_bin_states, 0);
    }

    for (int i = 0; i <= total_possibilities; i++) {

      vec observation = int_to_vec(i, n_histone_marks);

      for (int j = 0; j <= n_bin_states; j++) {
        emissions_memo[i][j] = emission_probability_direct(j, observation);
      }
    }
  }

  static vec int_to_vec(const unsigned int x, const unsigned int num_places);

  static uword vec_to_int(const vec &v);

  double emission_probability(const unsigned int bin_state,
                              unsigned int observation) const {
    return emissions_memo[observation][bin_state];
  }

  colvec emission_probabilities_col(const unsigned int observation) const {
    return colvec(emissions_memo[observation]);
  }

  rowvec emission_probabilities_row(const unsigned int observation) const {
    return rowvec(emissions_memo[observation]);
  }

  const mat &emissions_probabilities_m() const {
    return emissions_probabilities;
  }

  static Emissions_Dec
  emissions_from_states(const unsigned int n_bin_states,
                        const int n_histone_marks,
                        const std::vector<unsigned int> &bins_to_bin_states,
                        const mat &bin_data) {

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

    return Emissions_Dec(em);
  }
};
#endif
