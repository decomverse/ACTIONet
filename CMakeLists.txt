cmake_minimum_required (VERSION 3.1)
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

# find python
execute_process(COMMAND python-config --prefix OUTPUT_VARIABLE PYTHON_SEARCH_PATH)
string(REGEX REPLACE "\n$" "" PYTHON_SEARCH_PATH "${PYTHON_SEARCH_PATH}")
file(GLOB_RECURSE PYTHON_DY_LIBS ${PYTHON_SEARCH_PATH}/lib/libpython*.dylib ${PYTHON_SEARCH_PATH}/lib/libpython*.so)
if (PYTHON_DY_LIBS)
    list(GET PYTHON_DY_LIBS 0 PYTHON_LIBRARY)
    message("-- Find shared libpython: ${PYTHON_LIBRARY}")
else()
    message(WARNING "Cannot find shared libpython, try find_package")
endif()

find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(Python REQUIRED Development)
find_package(Python REQUIRED Interpreter)

find_package(PythonInterp)
find_package(PythonLibs ${PYTHON_VERSION_STRING} EXACT)

############ Third-party dependencies
find_package(PythonLibs REQUIRED)
find_package(PythonPyEnv REQUIRED)
include_directories(${Python_INCLUDE_DIRS})



# Check for CHOLMOD
find_package(SuiteSparse REQUIRED)
include_directories(${SUITESPARSE_INCLUDE_DIRS})

if (DEFINED NUMPY_INCLUDE_DIRS)
    include_directories(${NUMPY_INCLUDE_DIRS})
    find_package(Python REQUIRED)
else()
    find_package(Python REQUIRED COMPONENTS NumPy)
    include_directories(${Python_NumPy_INCLUDE_DIRS})
endif()


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
#include_directories(${PYTHON_INCLUDE_DIRS})

############ Compile
add_subdirectory(lib/ACTIONet)
add_subdirectory(lib/pybind11)

pybind11_add_module(_ACTIONet MODULE _ACTIONet/wrapper.cc)
target_link_libraries(_ACTIONet PUBLIC ACTIONet_core)
