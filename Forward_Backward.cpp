//
// Created by Stephanos Tsoucas on 6/14/17.
//

#include "Forward_Backward.h"

std::pair<mat, vec> interm_emission_probabilities(const Model &m,
                                                  const Probabilities &p) {
  vec posteriors_bin = zeros(m.n_bin_states);

  const unsigned long n_bin_positions = p.forward_probabilities.n_slices;

  mat updated_emissions = zeros(m.n_bin_states, p.n_histone_marks);

  for (int i = 0; i < n_bin_positions; i++) {

    // Sum over the domain
    const vec &bin_state_probabilities = sum(
        p.forward_probabilities.slice(i) % p.backward_probabilities.slice(i),
        1);

    const rowvec &observation = p.bin_data.col(i).t();

    posteriors_bin += bin_state_probabilities;

    mat intermediate = ones(m.n_bin_states, p.n_histone_marks);
    intermediate.each_col() %= bin_state_probabilities;
    intermediate.each_row() %= observation;

    updated_emissions += intermediate;
  }

  return std::pair<mat, vec>(updated_emissions.t(), posteriors_bin);
}

mat Forward_Backward::forward_vars_position_0(const Model &m,
                                              const Probabilities &p) {

  mat forward_probabilities(m.n_bin_states, m.n_domain_states);

  const uword observation_0 = p.bin_dec_data[0];

  for (int bin_state = 0; bin_state < m.n_bin_states; bin_state++) {
    for (int domain_state = 0; domain_state < m.n_domain_states;
         domain_state++) {
      const double state_probability =
          m.initial_probabilities(bin_state, domain_state);
      const double emission_probability =
          m.emission_probability(bin_state, observation_0);

      double forward_var = state_probability * emission_probability;
      forward_probabilities(bin_state, domain_state) = forward_var;
    }
  }

  return forward_probabilities;
}

// Updates p.forward_probabilities in place
bool Forward_Backward::forward_probabilities(const Model &m, Probabilities &p) {
  const unsigned long n_bins = p.n_bin_positions;
  const mat &domain_transitions = m.domain_transition_probabilities.t();
  const cube &bin_transitions = m.bin_transition_probabilities_t();

  mat forward_vars_0 = forward_vars_position_0(m, p);
  double scaling_factor_0 = accu(forward_vars_0);

  p.forward_probabilities.slice(0) = forward_vars_0 / scaling_factor_0;

  p.scaling_factors.at(0) = scaling_factor_0;

  for (uword bin_position = 1; bin_position < n_bins; bin_position++) {
    const uword observation = p.bin_dec_data[bin_position];
    const uword prev_bin_position = bin_position - 1;
    p.forward_probabilities.slice(bin_position) =
        mat(m.n_bin_states, m.n_domain_states);

    if (bin_position % m.domain_size == 0) {
      // Domain boundary (bin_position is at the beginning of a domain)
      const mat &forward_var_prev_position =
          p.forward_probabilities.slice(prev_bin_position);

      mat &forward_vars = p.forward_probabilities.slice(bin_position);

      for (int dest_domain_state = 0; dest_domain_state < m.n_domain_states;
           dest_domain_state++) {

        const mat &b = bin_transitions.slice(dest_domain_state);

        const mat &forward_vars_times_domain_transitions =
            forward_var_prev_position.each_row() %
            domain_transitions.row(dest_domain_state);
        mat forward_vars_dest_domain =
            b * forward_vars_times_domain_transitions;
        forward_vars.col(dest_domain_state) = sum(forward_vars_dest_domain, 1);
      }

      vec emissions_col = m.emission_probabilities_col(observation);
      forward_vars.each_col() %= emissions_col;

      double scaling_factor = accu(forward_vars);
      p.scaling_factors[bin_position] = scaling_factor;

      forward_vars /= scaling_factor;

    } else {
      // Intra-domain position

      mat &forward_vars = p.forward_probabilities.slice(bin_position);
      const mat &forward_var_prev_position =
          p.forward_probabilities.slice(prev_bin_position);

      for (int domain_state = 0; domain_state < m.n_domain_states;
           domain_state++) {
        const mat &b = bin_transitions.slice(domain_state);
        forward_vars.col(domain_state) =
            b * forward_var_prev_position.col(domain_state);
      }

      const vec &emissions_col = m.emission_probabilities_col(observation);
      forward_vars.each_col() %= emissions_col;

      double scaling_factor = accu(forward_vars);
      p.scaling_factors[bin_position] = scaling_factor;

      forward_vars /= scaling_factor;
    }
  }
  std::cout << "Finished forward probabilities" << std::endl;
  return true;
}

// Updates p.backward_probabilities in place. Note that p.forward_probabilities
// must have already been computed with
// the given model to produce the correct results.
bool Forward_Backward::backward_probabilities(const Model &m,
                                              Probabilities &p) {
  const unsigned long n_bins = p.n_bin_positions;

  p.backward_probabilities.slice(n_bins - 1) =
      ones(m.n_bin_states, m.n_domain_states);

  const mat &domain_transitions = m.domain_transition_probabilities;
  const cube &bin_transition_probabilities = m.bin_transition_probabilities;

  for (int bin_position = n_bins - 2; bin_position >= 0; bin_position--) {

    const long bin_position_ahead = bin_position + 1;
    const uword observation = p.bin_dec_data[bin_position_ahead];
    p.backward_probabilities.slice(bin_position) =
        zeros(m.n_bin_states, m.n_domain_states);

    if (bin_position_ahead % m.domain_size == 0) {

      mat &backward_vars = p.backward_probabilities.slice(bin_position);

      for (int dest_domain_state = 0; dest_domain_state < m.n_domain_states;
           dest_domain_state++) {

        const vec &b_i_1 = p.backward_probabilities.slice(bin_position_ahead)
                               .col(dest_domain_state);

        mat intermediate = ones(m.n_bin_states, m.n_domain_states);

        intermediate.each_col() %= b_i_1;

        const vec &emissions = m.emission_probabilities_col(observation);
        intermediate.each_col() %= emissions;

        // *= ?
        intermediate = bin_transition_probabilities.slice(dest_domain_state) *
                       intermediate;

        const rowvec &domain_transition =
            domain_transitions.col(dest_domain_state).t();
        intermediate.each_row() %= domain_transition;

        //#pragma critical
        backward_vars += intermediate;
      }

      backward_vars /= p.scaling_factors.at(bin_position_ahead);

    } else {

      mat &backward_vars = p.backward_probabilities.slice(bin_position);

      const vec &emissions = m.emission_probabilities_col(observation);
      const mat &backward_var_next_position =
          p.backward_probabilities.slice(bin_position_ahead).each_col() %
          emissions;

      for (int src_domain_state = 0; src_domain_state < m.n_domain_states;
           src_domain_state++) {
        const mat &bin_transitions =
            bin_transition_probabilities.slice(src_domain_state);
        backward_vars.col(src_domain_state) =
            bin_transitions * backward_var_next_position.col(src_domain_state);
      }

      backward_vars /= p.scaling_factors.at(bin_position_ahead);
    }
  }
  std::cout << "Finished backwards probabilities" << std::endl;
  return true;
}

// Compute this asynchronously
Intermediate_Model intermediate_model_update(const Model &m,
                                             const Probabilities &p) {
  std::cout << "Starting intermediate model update" << std::endl;
  const long n_bins = p.n_bin_positions;
  mat updated_domain_transitions = zeros(m.n_domain_states, m.n_domain_states);
  cube updated_bin_transitions =
      zeros(m.n_bin_states, m.n_bin_states, m.n_domain_states);
  mat flattened_bin_transition_probabilities(m.n_bin_states * m.n_bin_states,
                                             m.n_domain_states);

  for (int dest_domain = 0; dest_domain < m.n_domain_states; dest_domain++) {
    flattened_bin_transition_probabilities.col(dest_domain) =
        vectorise(m.bin_transition_probabilities.slice(dest_domain));
  }

  mat domain_transitions = m.domain_transition_probabilities;
  cube bin_transitions = m.bin_transition_probabilities;

  const mat domain_transitions_t = m.domain_transition_probabilities.t();

  mat inter_domain_base = ones(m.n_bin_states * m.n_domain_states,
                               m.n_bin_states * m.n_domain_states);
  for (int dest_domain = 0; dest_domain < m.n_domain_states; dest_domain++) {
    // Bin transition probabilities
    const int submat_start = dest_domain * m.n_domain_states;
    inter_domain_base
        .submat(0, submat_start,
                size(inter_domain_base.n_rows, m.n_domain_states))
        .each_col() %= flattened_bin_transition_probabilities.col(dest_domain);

    // Domain transition probabilities
    inter_domain_base
        .submat(0, submat_start,
                size(inter_domain_base.n_rows, m.n_domain_states))
        .each_row() %= domain_transitions_t.row(dest_domain);
  }

  int num_threads_for_update =
      n_bins < 1000000 ? 1 : ((n_bins - 600000) / 100000) + 1;
  std::cout << "bin positions " << n_bins << std::endl;
  std::cout << "num threads: " << num_threads_for_update << std::endl;
  //	#pragma omp parallel for schedule(static, 1000)
  //num_threads(num_threads_for_update)
  for (int bin_position = 0; bin_position < n_bins - 1; bin_position++) {
    const int bin_position_ahead = bin_position + 1;
    const mat &f = p.forward_probabilities.slice(bin_position);
    const mat &b = p.backward_probabilities.slice(bin_position_ahead);
    const uword observation = p.bin_dec_data[bin_position_ahead];
    const double scaling_factor = p.scaling_factors[bin_position_ahead];

    if (bin_position_ahead % m.domain_size == 0) {
      // Inter-domain position

      const vec &emission = m.emission_probabilities_col(observation);
      mat b_e = b.each_col() % emission;
      mat b_f = kron(b_e, f);
      b_f %= inter_domain_base;

      // Add up results for both domains and bins
      mat d_vec = sum(b_f, 0);
      mat domain_transitions_i =
          reshape(d_vec, size(m.n_domain_states, m.n_domain_states));
      domain_transitions_i /= scaling_factor;

      updated_domain_transitions += domain_transitions_i;

      for (int dest_domain = 0; dest_domain < m.n_domain_states;
           dest_domain++) {
        const int submat_start = dest_domain * m.n_domain_states;

        mat b_vec = sum(b_f.submat(0, submat_start,
                                   size(m.n_bin_states * m.n_bin_states,
                                        m.n_domain_states)),
                        1);
        mat bin_transitions_i = reshape(b_vec, m.n_bin_states, m.n_bin_states);
        bin_transitions_i /= scaling_factor;
        updated_bin_transitions.slice(dest_domain) += bin_transitions_i;
      }
    } else {
      const rowvec emissions = m.emission_probabilities_row(observation);
      const mat b_e = (b.each_col() % emissions);

      for (int dest_domain = 0; dest_domain < m.n_domain_states;
           dest_domain++) {

        mat bin_transitions =
            kron(f.col(dest_domain), b_e.col(dest_domain).t());
        bin_transitions %= m.bin_transition_probabilities.slice(dest_domain);
        // bin_transitions.each_row() %= emissions;
        bin_transitions /= scaling_factor;
        //		#pragma omp critical
        { updated_bin_transitions.slice(dest_domain) += bin_transitions; }
      }
    }
  }

  std::pair<mat, vec> updated_emission_probs =
      interm_emission_probabilities(m, p);

  std::cout << "Finished intermediate model update" << std::endl;
  return Intermediate_Model(updated_domain_transitions, updated_bin_transitions,
                            updated_emission_probs.first,
                            updated_emission_probs.second);
}

Model aggregate_intermediate_models(
    const Model &prev_m, const std::vector<Intermediate_Model> &models,
    const std::vector<Probabilities> &probabilities) {

  mat updated_domain_transitions =
      zeros(prev_m.n_domain_states, prev_m.n_domain_states);

  cube updated_bin_transitions =
      zeros(prev_m.n_bin_states, prev_m.n_bin_states, prev_m.n_domain_states);

  mat updated_emission_probs =
      zeros(prev_m.n_histone_marks, prev_m.n_bin_states);

  vec expected_bin_states = zeros(prev_m.n_bin_states);

  // Compute initial probabilities
  mat updated_initial_probabilities =
      zeros(prev_m.n_bin_states, prev_m.n_domain_states);

  for (auto &p : probabilities) {
    updated_initial_probabilities +=
        p.forward_probabilities.slice(0) % p.backward_probabilities.slice(0);
  }

  // Average over all chromosomes
  updated_initial_probabilities /= probabilities.size();

  for (auto &i_m : models) {

    // Domain transitions
    updated_domain_transitions += i_m.interm_domain_transitions;

    // Bin transitions
    for (int i = 0; i < prev_m.n_domain_states; i++) {
      updated_bin_transitions.slice(i) += i_m.interm_bin_transitions.slice(i);
    }

    // Emission probabilities
    updated_emission_probs += i_m.interm_emissions;
    expected_bin_states += i_m.interm_expected_bin_states;
  }

  // Domain transitions
  const vec &domain_transition_totals = sum(updated_domain_transitions, 1);
  updated_domain_transitions.each_col() /= domain_transition_totals;
  updated_domain_transitions.replace(NAN, 0.0);

  // Bin transitions
  for (int i = 0; i < prev_m.n_domain_states; i++) {
    mat totals = sum(updated_bin_transitions.slice(i), 1);
    updated_bin_transitions.slice(i).each_col() /= totals;
  }

  updated_bin_transitions.replace(NAN, 0.0);

  // Emissions
  // Sum across histone marks, take average over the observation in the bin
  // state
  const rowvec &emission_totals = expected_bin_states.t();

  updated_emission_probs.each_row() /= emission_totals;
  updated_emission_probs.replace(NAN, 0.0);

  return Model(prev_m.n_bin_states, prev_m.n_domain_states, prev_m.domain_size,
               updated_domain_transitions, updated_bin_transitions,
               updated_emission_probs, updated_initial_probabilities);
}

Model Forward_Backward::update_model_parameters(
    const Model &m, std::vector<Probabilities> &ps) {

  std::vector<Intermediate_Model> interm_ms;
  interm_ms.reserve(ps.size());

#pragma omp parallel for schedule(dynamic, 1) num_threads(56)
  for (int i = 0; i < ps.size(); i++) {

    forward_probabilities(m, ps.at(i));

    backward_probabilities(m, ps.at(i));

    interm_ms.push_back(intermediate_model_update(m, ps.at(i)));
  }

  return aggregate_intermediate_models(m, interm_ms, ps);
}
