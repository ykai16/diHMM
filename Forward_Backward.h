//
// Created by Stephanos Tsoucas on 6/14/17.
//

#ifndef BOOSTPYTHONHELLOWORLD_FORWARD_BACKWARD_H
#define BOOSTPYTHONHELLOWORLD_FORWARD_BACKWARD_H

#include "DiHMM_Trainer.h"
#include "Intermediate_Model.h"
#include <cblas.h>
#include <chrono>
#include <ctime>
using namespace arma;

class Forward_Backward {

  static mat forward_vars_position_0(const Model &m, const Probabilities &p);

public:
  static bool forward_probabilities(const Model &m, Probabilities &p);
  static bool backward_probabilities(const Model &m, Probabilities &p);
  static Model update_model_parameters(const Model &model,
                                       std::vector<Probabilities> &p);
};

#endif // BOOSTPYTHONHELLOWORLD_FORWARD_BACKWARD_H
