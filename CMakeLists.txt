cmake_minimum_required(VERSION 3.6...3.29)

project (ACTIONet)
list(APPEND CMAKE_MODULE_PATH "${ACTIONet_SOURCE_DIR}/cmake")

############ Compiler options
# C++11 is not supported, as the code uses C++14 features
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_compile_options(-w)
add_compile_options(-O3)

find_package(Python REQUIRED COMPONENTS Interpreter Development NumPy)
include_directories(${Python_NumPy_INCLUDE_DIRS})

# Check for pybind11
find_package(pybind11 CONFIG REQUIRED)

# Check for CHOLMOD
find_package(SuiteSparse REQUIRED)
include_directories(${SUITESPARSE_INCLUDE_DIRS})



# Print out the discovered paths
include(CMakePrintHelpers)
cmake_print_variables(Python_INCLUDE_DIRS)
cmake_print_variables(Python_NumPy_INCLUDE_DIRS)

############ Headers
include_directories(${ACTIONet_SOURCE_DIR}/_ACTIONet/include)
include_directories(${ACTIONet_SOURCE_DIR}/lib/ACTIONet/include)
include_directories(${ACTIONet_SOURCE_DIR}/lib/ACTIONet/include/math/arma)
include_directories(${ACTIONet_SOURCE_DIR}/lib/ACTIONet/include/math/arma/armadillo_bits)
include_directories(${ACTIONet_SOURCE_DIR}/lib/ACTIONet/include/clustering/HDBSCAN)
include_directories(${ACTIONet_SOURCE_DIR}/lib/ACTIONet/include/clustering/leiden)
include_directories(${ACTIONet_SOURCE_DIR}/lib/ACTIONet/include/clustering/leiden/igraph)
include_directories(${ACTIONet_SOURCE_DIR}/lib/ACTIONet/include/math)
include_directories(${ACTIONet_SOURCE_DIR}/lib/ACTIONet/include/math/StatsLib)
include_directories(${ACTIONet_SOURCE_DIR}/lib/ACTIONet/include/math/StatsLib/gcem)
include_directories(${ACTIONet_SOURCE_DIR}/lib/ACTIONet/include/math/pcg)
include_directories(${ACTIONet_SOURCE_DIR}/lib/ACTIONet/include/layout/uwot)
include_directories(${ACTIONet_SOURCE_DIR}/lib/ACTIONet/include/layout/s_gd2)
include_directories(${ACTIONet_SOURCE_DIR}/lib/ACTIONet/include/nn_search/hnsw)
include_directories(${ACTIONet_SOURCE_DIR}/lib/ACTIONet/include/nn_search/nmslib)
include_directories(${ACTIONet_SOURCE_DIR}/lib/ACTIONet/include/misc)
include_directories(${ACTIONet_SOURCE_DIR}/lib/pybind11/include)
#include_directories(${Python3_INCLUDE_DIRS})

############ Compile
add_subdirectory(lib/ACTIONet)

pybind11_add_module(_ACTIONet MODULE _ACTIONet/wrapper.cc)
target_link_libraries(_ACTIONet PUBLIC ACTIONet_core)
