//
// Created by Stephanos Tsoucas on 6/14/17.
//
#ifndef BOOSTPYTHONHELLOWORLD_DIHMM_TRAINER_H
#define BOOSTPYTHONHELLOWORLD_DIHMM_TRAINER_H

#include "Probabilities.h"
#include <boost/filesystem.hpp>
#include <chrono>
#include <ctime>
#include <future>
#include <memory>
#include <omp.h>

extern "C" void openblas_set_num_threads(int num_threads);

/*
    Computes the forward and backward probabilities based on Model m. p stores
   the results.
*/
template <typename T>
void forward_backward_probabilities(const Model &m,
                                    std::vector<Probabilities> &ps) {

  //    #pragma omp parallel for num_threads(56)
  for (int i = 0; i < ps.size(); i++) {
    T::forward_probabilities(m, ps.at(i));
    T::backward_probabilities(m, ps.at(i));
  }
}

template <typename T>
const Model run_dihmm(Model model, const std::vector<mat> &all_bin_data,
                      const int max_iterations, const double tolerance) {

  std::shared_ptr<Model> m = std::make_shared<Model>(model);

  std::vector<cube> backward_ps_cubes;
  backward_ps_cubes.reserve(all_bin_data.size());
  std::vector<cube> forward_ps_cubes;
  forward_ps_cubes.reserve(all_bin_data.size());
  std::vector<std::vector<double>> scaling_factors;
  std::vector<Probabilities> ps;
  ps.reserve(all_bin_data.size());

  bool converged = false;
  for (uword i = 0; i < all_bin_data.size(); i++) {

    const mat &bin_data_i = all_bin_data.at(i);
    const uword data_size = bin_data_i.n_cols;
    const uword n_bin_positions = data_size - (data_size % m->domain_size);
    std::cout << "BIN POSITIONS " << n_bin_positions << std::endl;

    backward_ps_cubes.emplace_back(m->n_bin_states, m->n_domain_states,
                                   n_bin_positions);
    forward_ps_cubes.emplace_back(m->n_bin_states, m->n_domain_states,
                                  n_bin_positions);
    scaling_factors.emplace_back(n_bin_positions, 0);
  }

  for (int i = 0; i < all_bin_data.size(); i++) {
    ps.emplace_back(all_bin_data.at(i), forward_ps_cubes.at(i),
                    backward_ps_cubes.at(i), scaling_factors.at(i));
  }

  double initial_log_lik = 0;
  for (auto &p : ps) {
    initial_log_lik += p.log_likelihood();
  }

  const int interval = 20;

  std::cout << "Initial log likelihood: " << initial_log_lik << std::endl;

  for (int current_iter = 0; !converged && current_iter < max_iterations;
       current_iter++) {

    double prev_loglikelihood = 0;

    // TODO make a helper for this: it's all over the place
    for (auto &p : ps) {
      prev_loglikelihood += p.log_likelihood();
    }

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    m = std::make_shared<Model>(T::update_model_parameters(*m, ps));
    end = std::chrono::system_clock::now();
    //  forward_backward_probabilities<T>(*m, ps);
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    std::cout << "finished computation at " << std::ctime(&end_time)
              << "elapsed time: " << elapsed_seconds.count() << "s\n";
    std::cout << "Current iteration: " << current_iter << std::endl;
    // Add up sum of all log_likelihoods across all Probabilities
    double current_log_lik = 0;

    for (const auto &p : ps) {
      current_log_lik += p.log_likelihood();
    }

    double delta = std::abs(current_log_lik - prev_loglikelihood) /
                   (1.0 + std::abs(prev_loglikelihood));
    std::cout << "Delta in likelihood: "
              << (current_log_lik - prev_loglikelihood) << std::endl;
    converged = delta < tolerance;
    std::cout << "Converged: " << converged << std::endl;
    std::cout << "Old ll: " << prev_loglikelihood << std::endl;
    std::cout << "New ll: " << current_log_lik << std::endl;

    if (current_iter % interval == 0) {
      std::string path =
          "./partially_trained_" + std::to_string(current_iter) + "/";
      boost::filesystem::create_directory(path);
      save_model(*m, path);
    }
  }

  if (!converged) {
    std::cout << "Model did not converge " << std::endl;
  }

  return *m;
}

template <typename T>
const Model run_dihmm(const int n_bin_states, const int n_domain_states,
                      const int domain_size, const int max_iterations,
                      const double tolerance,
                      const std::vector<mat> &all_bin_data) {
  openblas_set_num_threads(1);
  const uword num_histone_marks = all_bin_data.at(0).n_rows;
  bool converged = false;
  // TODO: concatenate all matrices into one big matrix to pass into
  // initial_model
  // TODO: do this in a different scope so it gets collected
  mat all_bin_data_concat = all_bin_data.at(0);

  for (int i = 1; i < all_bin_data.size(); i++) {
    all_bin_data_concat = join_rows(all_bin_data_concat, all_bin_data.at(i));
  }

  std::shared_ptr<Model> m = std::make_shared<Model>(
      initial_model(all_bin_data_concat, n_bin_states, n_domain_states,
                    num_histone_marks, domain_size));

  std::vector<cube> backward_ps_cubes;
  backward_ps_cubes.reserve(all_bin_data.size());
  std::vector<cube> forward_ps_cubes;
  forward_ps_cubes.reserve(all_bin_data.size());
  std::vector<std::vector<double>> scaling_factors;
  std::vector<Probabilities> ps;
  ps.reserve(all_bin_data.size());

  for (uword i = 0; i < all_bin_data.size(); i++) {

    const mat &bin_data_i = all_bin_data.at(i);
    const uword data_size = bin_data_i.n_cols;
    const uword n_bin_positions = data_size - (data_size % domain_size);
    std::cout << "BIN POSITIONS " << n_bin_positions << std::endl;

    backward_ps_cubes.emplace_back(m->n_bin_states, m->n_domain_states,
                                   n_bin_positions);
    forward_ps_cubes.emplace_back(m->n_bin_states, m->n_domain_states,
                                  n_bin_positions);
    scaling_factors.emplace_back(n_bin_positions, 0);
  }

  for (int i = 0; i < all_bin_data.size(); i++) {
    ps.emplace_back(all_bin_data.at(i), forward_ps_cubes.at(i),
                    backward_ps_cubes.at(i), scaling_factors.at(i));
  }

  // forward_backward_probabilities<T>(*m, ps);

  // Add up sum of all log_likelihoods across all Probabilities
  double initial_log_lik = 0;

  for (auto &p : ps) {
    initial_log_lik += p.log_likelihood();
  }

  const int interval = 40;

  std::cout << "Initial log likelihood: " << initial_log_lik << std::endl;

  for (int current_iter = 0; !converged && current_iter < max_iterations;
       current_iter++) {

    double prev_loglikelihood = 0;

    // TODO make a helper for this: it's all over the place
    for (auto &p : ps) {
      prev_loglikelihood += p.log_likelihood();
    }

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    m = std::make_shared<Model>(T::update_model_parameters(*m, ps));
    end = std::chrono::system_clock::now();
    //  forward_backward_probabilities<T>(*m, ps);
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    std::cout << "finished computation at " << std::ctime(&end_time)
              << "elapsed time: " << elapsed_seconds.count() << "s\n";
    std::cout << "Current iteration: " << current_iter << std::endl;
    // Add up sum of all log_likelihoods across all Probabilities
    double current_log_lik = 0;

    for (const auto &p : ps) {
      current_log_lik += p.log_likelihood();
    }

    double delta = std::abs(current_log_lik - prev_loglikelihood) /
                   (1.0 + std::abs(prev_loglikelihood));
    std::cout << "Delta in likelihood: "
              << (current_log_lik - prev_loglikelihood) << std::endl;
    converged = delta < tolerance;
    std::cout << "Converged: " << converged << std::endl;
    std::cout << "Old ll: " << prev_loglikelihood << std::endl;
    std::cout << "New ll: " << current_log_lik << std::endl;

    if (current_iter % interval == 0) {
      std::string path =
          "./partially_trained" + std::to_string(current_iter) + "/";
      boost::filesystem::create_directory(path);
      save_model(*m, path);
    }
  }

  if (!converged) {
    std::cout << "Model did not converge " << std::endl;
  }

  return *m;
}

#endif // BOOSTPYTHONHELLOWORLD_DIHMM_TRAINER_H
