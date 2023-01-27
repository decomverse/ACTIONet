cmake_minimum_required(VERSION 3.1)
project(ACTIONet)
list(APPEND CMAKE_MODULE_PATH "${ACTIONet_SOURCE_DIR}/cmake")

# ########### Compiler options
# C++11 is not supported, as the code uses C++17 features
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_compile_options(-w)
add_compile_options(-O3)

# ########### Third-party dependencies
find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})

# Check for BLAS/LAPACK
find_package(BLAS REQUIRED)
message(STATUS "BLAS_COMPILER_FLAGS: ${BLAS_COMPILER_FLAGS}")
add_definitions("-DEIGEN_USE_BLAS") # is adding  ${BLAS_COMPILER_FLAGS} necessary?
list(APPEND EXTERNAL_LIBS "${BLAS_LIBRARIES}")

# Check for CHOLMOD
find_package(SuiteSparse REQUIRED)
include_directories(${SUITESPARSE_INCLUDE_DIRS})

# ########### Headers
include_directories(${ACTIONET_HEADER_DIRS} ${ACTIONet_SOURCE_DIR}/include ${ACTIONet_SOURCE_DIR}/include/math/aarand )

# ########### Compile
file(GLOB_RECURSE SOURCES *.cpp *.cc *.c)
add_library(ACTIONet_core STATIC ${SOURCES})
target_link_libraries(
    ACTIONet_core PUBLIC
    ${BLAS_LIBRARIES}
    ${SUITESPARSE_LIBRARIES}
    ${HDF5_LIBRARIES}
)
