//
// Created by Stephanos Tsoucas on 6/14/17.
//

#ifndef BOOSTPYTHONHELLOWORLD_ARMA_NUMPY_CONVERSIONS_H
#define BOOSTPYTHONHELLOWORLD_ARMA_NUMPY_CONVERSIONS_H

#include <armadillo>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace p = boost::python;
namespace np = boost::python::numpy;
namespace a = arma;

np::ndarray from_vec(const a::vec &v) {
  p::tuple shape = p::make_tuple(v.size());
  np::dtype dtype = np::dtype::get_builtin<double>();
  np::ndarray n = np::zeros(shape, dtype);

  for (int i = 0; i < v.size(); i++) {
    n[i] = v(i);
  }

  return n;
}

np::ndarray from_mat(const a::mat &m) {
  p::tuple shape = p::make_tuple(m.n_rows, m.n_cols);
  np::dtype dtype = np::dtype::get_builtin<double>();
  np::ndarray n = np::zeros(shape, dtype);

  for (int i = 0; i < m.n_rows; i++) {
    for (int j = 0; j < m.n_cols; j++) {
      n[i][j] = m(i, j);
    }
  }
  return n;
}

np::ndarray from_cube(const a::cube &c) {
  p::tuple shape = p::make_tuple(c.n_slices, c.n_rows, c.n_cols);
  np::dtype dtype = np::dtype::get_builtin<double>();
  np::ndarray n = np::zeros(shape, dtype);

  for (int slice = 0; slice < c.n_slices; slice++) {
    for (int i = 0; i < c.n_rows; i++) {
      for (int j = 0; j < c.n_cols; j++) {
        n[slice][i][j] = c(i, j, slice);
      }
    }
  }

  return n;
}

np::ndarray
from_vec_of_pairs(const std::vector<std::pair<a::uword, a::uword>> &v) {
  p::tuple shape = p::make_tuple(v.size(), 2);
  np::dtype dtype = np::dtype::get_builtin<int>();
  np::ndarray n = np::zeros(shape, dtype);

  for (int i = 0; i < v.size(); i++) {
    n[i][0] = v.at(i).first;
    n[i][1] = v.at(i).second;
  }

  return n;
}

#endif // BOOSTPYTHONHELLOWORLD_ARMA_NUMPY_CONVERSIONS_H
