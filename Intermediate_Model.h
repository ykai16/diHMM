
#ifndef BOOSTPYTHONHELLOWORLD_INTERMEDIATE_MODEL_H
#define BOOSTPYTHONHELLOWORLD_INTERMEDIATE_MODEL_H

// TODO: make these matrices into Future<mat>
struct Intermediate_Model {

  const mat interm_domain_transitions;

  const cube interm_bin_transitions;

  const mat interm_emissions;

  // Used to divide interm_emissions by
  const vec interm_expected_bin_states;

  Intermediate_Model(const mat interm_domain_transitions,
                     const cube interm_bin_transitions,
                     const mat interm_emissions,
                     const vec interm_expected_bin_states)
      : interm_domain_transitions(interm_domain_transitions),
        interm_bin_transitions(interm_bin_transitions),
        interm_emissions(interm_emissions),
        interm_expected_bin_states(interm_expected_bin_states) {}
};

#endif // BOOSTPYTHONHELLOWORLD_INTERMEDIATE_MODEL_H