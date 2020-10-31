//
// Created by Stephanos Tsoucas on 7/5/17.
//

#ifndef BOOSTPYTHONHELLOWORLD_ANNOTATION_H
#define BOOSTPYTHONHELLOWORLD_ANNOTATION_H

#include <armadillo>

using namespace arma;

class Annotation {

  const vec bins;
  const vec domains;
  const mat &bin_data;
  const std::vector<std::pair<uword, uword>> anno;

public:
  Annotation(vec bin_state_distributions, vec domain_state_distributions,
             const mat &bin_data,
             const std::vector<std::pair<uword, uword>> &annotation_data)
      : bins(bin_state_distributions), domains(domain_state_distributions),
        bin_data(bin_data), anno(annotation_data) {}

  const vec bin_state_coverage() const { return bins; }

  const vec domain_state_coverage() const { return domains; }

  const std::vector<std::pair<uword, uword>> &annotations() const {
    return anno;
  }
};

#endif // BOOSTPYTHONHELLOWORLD_ANNOTATION_H
