program test 
  use sparse_matrix_mod
  implicit none
  type(sparse_matrix_type) :: mat
  
  call sparse_new_from_file(mat, '../data/10x10.dat')
  call sparse_print(mat)
  call sparse_del(mat)

end program test
