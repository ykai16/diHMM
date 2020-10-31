cmake_minimum_required( VERSION 3.3 )

project( BoostPythonHelloWorld )

#set(Boost_USE_STATIC_LIBS   ON)

if(APPLE)
set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")
endif(APPLE)

#set(Armadillo_DIR /homes10/ext/tsoucas/workspace/usr/share/Armadillo/CMake/Armadillo )

#include_directories( BEFORE /homes10/ext/tsoucas/workspace/usr/include /homes10/ext/tsoucas/workspace/usr/share /homes10/ext/tsoucas/workspace/usr/lib )
#set(Boost_NO_BOOST_CMAKE ON )
#set(Boost_NO_SYSTEM_PATHS ON)
#set(BOOST_ROOT /homes10/ext/tsoucas/workspace/boost_1_64_0 )
#set(BOOST_LIBRARYDIR /homes10/ext/tsoucas/workspace/usr/lib)

set( ENV{NVBLAS_CONFIG_FILE} "~/diHMM/dihmm_py/nvblas.conf" )

set( ENV{ARMA_DONT_USE_WRAPPER} )
find_package( Armadillo REQUIRED )
#find_package( Armadillo REQUIRED PATHS /homes10/ext/tsoucas/workspace/usr/include /homes10/ext/tsoucas/workspace/usr/share NO_DEFAULT_PATHS ) 
include_directories(${ARMADILLO_INCLUDE_DIRS})

# Find necessary packages
find_package( PythonLibs 2.7 REQUIRED )
include_directories( ${PYTHON_INCLUDE_DIRS} )
 
find_package( Boost COMPONENTS python numpy filesystem REQUIRED  )
include_directories( ${Boost_INCLUDE_DIR} )

set( ENV{BLA_VENDOR} "OpenBLAS" )
find_package( BLAS REQUIRED ) 
 
# Build our library
add_library( dihmm SHARED Model.cpp Model.h Emissions.cpp Emissions.h Emissions_Dec.cpp Emissions_Dec.h Viterbi.h Probabilities.h DiHMM_Trainer.h Forward_Backward.cpp Forward_Backward.h Annotation.h Annotation_Py.h Intermediate_Model.h)

#set (CMAKE_CXX_FLAGS "-std=c++0x")
set (CMAKE_CXX_FLAGS "-std=c++0x -fopenmp -DARMA_DONT_USE_WRAPPER -DOPENBLAS_NUM_THREADS=4 -O3")
set (CMAKE_SHARED_LINKER_FLAGS "-L/usr/lib -lopenblas -lboost_python -lboost_system -lgomp -lpython2.7 ")


target_link_libraries( dihmm ${Boost_LIBRARIES} ${BOOST_NUMPY_LIBRARY} ${PYTHON_LIBRARIES} ${ARMADILLO_LIBRARIES} openblas lapack) 
# Define the wrapper library that wraps our library
add_library( dihmm_ext SHARED dihmm_ext.cpp Arma_Numpy_Conversions.h Model_Py.h)
target_link_libraries( dihmm_ext dihmm ${Boost_LIBRARIES} ${BOOST_NUMPY_LIBRARY} ${PYTHON_LIBRARIES} ${ARMADILLO_LIBRARIES} openblas lapack)
# don't prepend wrapper library name with lib
set_target_properties( dihmm_ext PROPERTIES PREFIX "" )
