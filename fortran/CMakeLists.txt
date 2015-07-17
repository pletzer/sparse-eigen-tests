cmake_minimum_required (VERSION 2.8)
project(sparse_tests)

enable_language(Fortran)

find_package(MPI)
set(CMAKE_FORTRAN_COMPILER ${MPI_FORTRAN_COMPILER})

add_executable(test_read test_read.f90 sparse_matrix.f90)

add_executable(computeEigen computeEigen.f90 sparse_matrix.f90)
target_link_libraries(computeEigen ${ARPACK_LIBS} ${LAPACK_LIB} ${BLAS_LIB})


