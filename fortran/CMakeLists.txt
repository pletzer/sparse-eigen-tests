cmake_minimum_required (VERSION 2.8)
project(sparse_tests)

enable_language(Fortran)

find_package(MPI)
set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})

add_executable(test_read test_read.f90 sparse_matrix.f90)

set(ARPACK_LIBRARY_DIR
    /opt/aci/sw/arpack-ng/3.1.5_gcc-4.9.2-1/lib)
set(LAPACK_LIB 
   /opt/aci/sw/lapack/3.5.0_gcc-4.9.2-1/lib/liblapack.so)
set(BLAS_LIB
   /opt/aci/sw/blas/3.5.0_gcc-4.9.2-1/libblas.so)
set(ARPACK_LIBS 
    ${ARPACK_LIBRARY_DIR}/libparpack.so 
    ${ARPACK_LIBRARY_DIR}/libarpack.so)
add_executable(computeEigen computeEigen.f90 sparse_matrix.f90)
target_link_libraries(computeEigen ${ARPACK_LIBS} ${LAPACK_LIB} ${BLAS_LIB})


