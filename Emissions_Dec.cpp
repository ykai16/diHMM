#include "Emissions_Dec.h"

vec Emissions_Dec::int_to_vec(const unsigned int x,
                              const unsigned int num_places) {

  vec ret(num_places);

  for (int i = 0; i < num_places; i++) {
    ret(i) = (x >> i) & 1;
  }

  return ret;
}

uword Emissions_Dec::vec_to_int(const vec &v) {

  uword ret = 0;
  uword s = v.n_elem;
#pragma omp simd
  for (int i = 0; i < s; i++) {
    ret += (((int)v(i)) << i);
  }

  return ret;
}
