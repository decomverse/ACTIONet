cmake_minimum_required (VERSION 3.1)
project (ACTIONet)
list(APPEND CMAKE_MODULE_PATH "${ACTIONet_SOURCE_DIR}/cmake")

############ Compiler options
# C++11 is not supported, as the code uses C++14 features
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_compile_options(-w)

############ Third-party dependencies
find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})

if (DEFINED NUMPY_INCLUDE_DIRS)
    include_directories(${NUMPY_INCLUDE_DIRS})
    find_package(Python REQUIRED)
else()
    find_package(Python REQUIRED COMPONENTS NumPy)
    include_directories(${Python_NumPy_INCLUDE_DIRS})
endif()
include_directories(${Python_INCLUDE_DIRS})

# Check for MKL and depending on whether or not MKL is present
# check for lapack and blas with suitesparse
find_package(MKL)
include_directories(${MKL_INCLUDE_DIRS})
if(NOT ${MKL_FOUND})
    set(SUITESPARSE_USE_LAPACK_BLAS ON)
endif()
find_package(SuiteSparse REQUIRED)
include_directories(${SUITESPARSE_INCLUDE_DIRS})

############ Headers
include_directories(${ACTIONet_SOURCE_DIR}/include)
include_directories(${ACTIONet_SOURCE_DIR}/include/arma)
include_directories(${ACTIONet_SOURCE_DIR}/include/arma/armadillo_bits)
include_directories(${ACTIONet_SOURCE_DIR}/include/carma)
include_directories(${ACTIONet_SOURCE_DIR}/include/ACTIONet)
include_directories(${ACTIONet_SOURCE_DIR}/include/ACTIONet/boost)
include_directories(${ACTIONet_SOURCE_DIR}/include/ACTIONet/clixo)
include_directories(${ACTIONet_SOURCE_DIR}/include/ACTIONet/HDBSCAN)
include_directories(${ACTIONet_SOURCE_DIR}/include/ACTIONet/hnsw)
include_directories(${ACTIONet_SOURCE_DIR}/include/ACTIONet/leiden)
include_directories(${ACTIONet_SOURCE_DIR}/include/ACTIONet/leiden/igraph)
include_directories(${ACTIONet_SOURCE_DIR}/include/ACTIONet/s_gd2)
include_directories(${ACTIONet_SOURCE_DIR}/include/ACTIONet/uwot)

############ Compile
add_subdirectory(src)
add_subdirectory(pybind11)

pybind11_add_module(_ACTIONet python_interface/wrapper.cc)
target_link_libraries(_ACTIONet PUBLIC ACTIONet_src)
