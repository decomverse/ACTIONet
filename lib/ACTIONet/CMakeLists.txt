file(GLOB_RECURSE SOURCES *.cpp *.cc *.c)

add_library(ACTIONet_src STATIC ${SOURCES})
target_link_libraries(
    ACTIONet_src PUBLIC
    ${BLAS_LIBRARIES}
    ${SUITESPARSE_LIBRARIES}
    ${HDF5_LIBRARIES}
)
