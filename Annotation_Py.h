//
// Created by Stephanos Tsoucas on 7/5/17.
//

#ifndef BOOSTPYTHONHELLOWORLD_ANNOTATION_PY_H
#define BOOSTPYTHONHELLOWORLD_ANNOTATION_PY_H

#include "Annotation.h"
#include "Arma_Numpy_Conversions.h"

class Annotation_Py {

  const Annotation a;

public:
  Annotation_Py(const Annotation annotation) : a(annotation) {}

  np::ndarray bin_state_distributions() {
    return from_vec(a.bin_state_coverage());
  }

  np::ndarray domain_state_distributions() {
    return from_vec(a.domain_state_coverage());
  }

  np::ndarray annotations() const { return from_vec_of_pairs(a.annotations()); }
};

#endif // BOOSTPYTHONHELLOWORLD_ANNOTATION_PY_H
