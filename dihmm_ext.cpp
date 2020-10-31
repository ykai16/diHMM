#include "Model_Py.h"

namespace p = boost::python;
namespace np = boost::python::numpy;

BOOST_PYTHON_MODULE(dihmm_ext) {
  using namespace boost::python;

  Py_Initialize();
  np::initialize();
  def("test_gpu", test_gpu);
  def("run_dihmm_existing_model", run_dihmm_existing_py);
  def("run_dihmm", run_dihmm_py);
  def("annotate", annotate_py);
  def("save_model", save_model_py);
  def("load_model", load_model_py);

  class_<Model_Py>("Model", no_init)
      .def_readonly("domain_transition_probabilities",
                    &Model_Py::domain_transition_probabilities)
      .def_readonly("bin_transition_probabilities",
                    &Model_Py::bin_transition_probabilities)
      .def_readonly("emission_probabilities", &Model_Py::emission_probabilities)
      .def_readonly("n_bin_states", &Model_Py::n_bin_states)
      .def_readonly("n_domain_states", &Model_Py::n_domain_states);

  class_<Annotation_Py>("Annotation", no_init)
      .def_readonly("bin_state_distributions",
                    &Annotation_Py::bin_state_distributions)
      .def_readonly("domain_state_distributions",
                    &Annotation_Py::domain_state_distributions)
      .def_readonly("annotations", &Annotation_Py::annotations);
}
