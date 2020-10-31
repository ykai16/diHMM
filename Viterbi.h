//
// Created by Stephanos Tsoucas on 6/26/17.
//

#ifndef BOOSTPYTHONHELLOWORLD_VITERBI_H
#define BOOSTPYTHONHELLOWORLD_VITERBI_H

#include "Annotation.h"
#include "Model.h"
#include <future>
#include <utility>

void print_trellis(
    const std::vector<std::vector<std::vector<std::pair<uword, uword>>>>
        &sub_solution_states) {

  std::cout << "Sub solution states size: " << sub_solution_states.size()
            << std::endl;

  for (auto &v : sub_solution_states) {
    for (auto &v1 : v) {
      for (auto &p : v1) {
        std::cout << "first: " << p.first << " second: " << p.second
                  << std::endl;
      }
      std::cout << "--" << std::endl;
    }
    std::cout << "---------" << std::endl;
  }
}

void initialize_bin_states(
    const uword n_bin_positions, const uword n_bin_states,
    const uword n_domain_states,
    std::vector<std::vector<std::vector<std::pair<uword, uword>>>>
        &sub_solution_states) {

  for (uword i = 0; i < n_bin_positions; i++) {

    sub_solution_states.push_back(
        std::vector<std::vector<std::pair<uword, uword>>>());

    for (uword bin_state = 0; bin_state < n_bin_states; bin_state++) {

      sub_solution_states[i].push_back(std::vector<std::pair<uword, uword>>());

      for (uword domain_state = 0; domain_state < n_domain_states;
           domain_state++) {
        sub_solution_states[i][bin_state].push_back(
            std::pair<uword, uword>(0, 0));
      }
    }
  }
}

// TODO: reduce code duplication
// Given annotations, find the total numbers
vec bin_state_distributions(
    const Model &m,
    const std::vector<std::pair<uword, uword>> &state_annotations) {

  vec bin_state_counts = zeros(m.n_bin_states);

  for (const auto &p : state_annotations) {
    bin_state_counts(p.first) += 1.0;
  }

  double sum = accu(bin_state_counts);

  bin_state_counts /= sum;

  return bin_state_counts;
}

vec domain_state_distributions(
    const Model &m,
    const std::vector<std::pair<uword, uword>> &state_annotations) {

  vec domain_state_counts = zeros(m.n_domain_states);

  for (const auto &p : state_annotations) {
    domain_state_counts.at(p.second) += 1.0;
  }

  double sum = accu(domain_state_counts);

  domain_state_counts /= sum;

  return domain_state_counts;
}

// Returns a list of pairs- first is bin state, second domain state.
Annotation annotate_data(const Model &m, const mat &bin_data) {

  const uword bin_positions = bin_data.n_cols;

  cube sub_solutions(m.n_bin_states, m.n_domain_states, bin_positions);

  std::vector<std::vector<std::vector<std::pair<uword, uword>>>>
      sub_solutions_states;
  initialize_bin_states(bin_positions, m.n_bin_states, m.n_domain_states,
                        sub_solutions_states);

  std::vector<std::pair<uword, uword>> annotations(
      bin_positions, std::pair<uword, uword>(0, 0));

  // Base case
  vec emissions_0 = log(
      m.emission_probabilities_col(Emissions_Dec::vec_to_int(bin_data.col(0))));

  mat initial_probabilities = log(m.initial_probabilities);

  cube bin_transitions = log(m.bin_transition_probabilities_t());

  mat domain_transitions = log(m.domain_transition_probabilities.t());

  sub_solutions.slice(0) = initial_probabilities.each_col() + emissions_0;

  for (uword bin_state = 0; bin_state < m.n_bin_states; bin_state++) {
    for (uword domain_state = 0; domain_state < m.n_domain_states;
         domain_state++) {
      sub_solutions_states.at(0).at(bin_state).at(domain_state) =
          std::pair<uword, uword>(bin_state, domain_state);
    }
  }

  for (uword i = 1; i < bin_positions; i++) {

    const vec &observation = bin_data.col(i);
    const mat &prev_solution = sub_solutions.slice(i - 1);
    mat solution_i(m.n_bin_states, m.n_domain_states);

    // TODO: emission probabilities
    if (i % m.domain_size == 0) {

      const mat &solution_prev = sub_solutions.slice(i - 1);

      for (uword dest_domain_state = 0; dest_domain_state < m.n_domain_states;
           dest_domain_state++) {

        const mat &b = bin_transitions.slice(dest_domain_state);

        vec emissions = log(m.emission_probabilities_col(
            Emissions_Dec::vec_to_int(observation)));

        const mat &prev_solutions_domains =
            solution_prev.each_row() +
            domain_transitions.row(dest_domain_state);

        mat interm(m.n_bin_states, m.n_domain_states);

        // interm_states(dest_bin_state, src_domain_state) gives the argmax bin
        // state originating in  src_domain_state  interm tells me bin
        umat interm_states(m.n_bin_states, m.n_domain_states);

        for (uword src_domain_state = 0; src_domain_state < m.n_domain_states;
             src_domain_state++) {

          for (uword dest_bin_state = 0; dest_bin_state < m.n_bin_states;
               dest_bin_state++) {

            const rowvec &b_p_s =
                b.row(dest_bin_state) +
                prev_solutions_domains.col(src_domain_state).t();

            // TODO wrong indexing?
            interm_states(dest_bin_state, src_domain_state) = index_max(b_p_s);
            interm(dest_bin_state, src_domain_state) = max(b_p_s);
          }
        }

        // Take max over rows of interm;
        solution_i.col(dest_domain_state) = max(interm, 1);
        // this tells me domain -> put together bin, domain pair
        const uvec &max_domain_indices = index_max(interm, 1);
        // Get the domain index. Look up the bin state
        for (uword dest_bin_state = 0; dest_bin_state < m.n_bin_states;
             dest_bin_state++) {

          // TODO: verify
          const uword max_src_domain_state = max_domain_indices(dest_bin_state);
          const uword max_src_bin_state =
              interm_states(dest_bin_state, max_src_domain_state);

          sub_solutions_states.at(i).at(dest_bin_state).at(dest_domain_state) =
              std::pair<uword, uword>(max_src_bin_state, max_src_domain_state);
        }
      }

      const vec &emissions = log(
          m.emission_probabilities_col(Emissions_Dec::vec_to_int(observation)));
      solution_i.each_col() += emissions;

    } else {

      // Inter-domain
      for (uword domain_state = 0; domain_state < m.n_domain_states;
           domain_state++) {

        const mat &b = bin_transitions.slice(domain_state);

        const rowvec &prev = prev_solution.col(domain_state).t();

        const mat &b_prev = b.each_row() + prev;

        // Take max over bin states
        solution_i.col(domain_state) = max(b_prev, 1);
        const uvec &max_bin_states = index_max(b_prev, 1);

        for (uword dest_bin_state = 0; dest_bin_state < m.n_bin_states;
             dest_bin_state++) {
          // assign to solution_states
          // make pair of bin-, domain-state
          sub_solutions_states.at(i).at(dest_bin_state).at(domain_state) =
              std::pair<uword, uword>(max_bin_states(dest_bin_state),
                                      domain_state);
        }
      }

      const vec &emissions = log(
          m.emission_probabilities_col(Emissions_Dec::vec_to_int(observation)));
      solution_i.each_col() += emissions;
    }

    // TODO populate the sub_solutions.slice(i) directly?
    sub_solutions.slice(i) = solution_i;
  }

  // Find the (bin, domain) state pair that produces the highest likelihood
  const uvec max_likelihood =
      ind2sub(size(sub_solutions.slice(bin_positions - 1)),
              sub_solutions.slice(bin_positions - 1).index_max());

  uword max_bin_state = max_likelihood(0);
  uword max_domain_state = max_likelihood(1);

  for (sword i = bin_positions - 1; i >= 0; i--) {

    const std::pair<uword, uword> &prev_state =
        sub_solutions_states.at(i).at(max_bin_state).at(max_domain_state);

    annotations.at(i) = prev_state;
    max_bin_state = prev_state.first;
    max_domain_state = prev_state.second;
  }

  Annotation a(bin_state_distributions(m, annotations),
               domain_state_distributions(m, annotations), bin_data,
               annotations);
  return a;
}

std::vector<Annotation> annotate(const Model &m,
                                 const std::vector<mat> &bin_data) {

  std::vector<std::future<Annotation>> a_f;
  std::vector<Annotation> as;

  for (const auto &b : bin_data) {
    a_f.push_back(std::async(std::launch::async, annotate_data, m, b));
  }

  for (auto &a : a_f) {
    as.push_back(a.get());
  }

  return as;
}

#endif // BOOSTPYTHONHELLOWORLD_VITERBI_H
