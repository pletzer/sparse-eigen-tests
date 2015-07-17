cmake_minimum_required (VERSION 2.8)
project(sparse_tests)

enable_language(Fortran)
add_executable(test_read test_read.f90 sparse_matrix.f90)


