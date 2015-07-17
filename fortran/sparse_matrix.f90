module sparse_matrix_mod
  implicit none

  type sparse_matrix_type
    integer, allocatable :: irows(:)
    integer, allocatable :: jcols(:)
    real*8, allocatable  :: values(:)
    integer              :: nrows, ncols, nnz
  end type sparse_matrix_type

  contains

  subroutine sparse_new_from_file(self, filename)
    type(sparse_matrix_type), intent(inout) :: self
    character*(*), intent(in)               :: filename
    integer :: i, j, k, iu
    real*8 :: val
    integer, parameter :: buff_size = 1000000
    integer, dimension(buff_size) :: ibuff, jbuff
    real*8, dimension(buff_size) :: buff

    iu = 10
    open(iu, file=filename, status='old', action='read') 
    k = 1
    do 
      read(iu, *, end=1) i, j, val
      jbuff(k) = j
      ibuff(k) = i
      buff(k) = val
      k = k + 1
      if (k > buff_size) then 
        exit
      endif 
    end do
    1 close(iu)

    allocate(self%irows(k))
    allocate(self%jcols(k))
    allocate(self%values(k))
    self%irows(:) = ibuff(1:k)
    self%jcols(:) = jbuff(1:k)
    self%values(:) = buff(1:k)
    self%nrows = maxval(self%irows)
    self%ncols = maxval(self%jcols)
    self%nnz = k
    ! should add some error code when buff_size is too small
    ! just write an error message for the time being
    if (k == buff_size) then
       print *, 'ERROR: need to increase buff_size!!!'
    endif

  end subroutine sparse_new_from_file

  subroutine sparse_del(self)
    type(sparse_matrix_type), intent(inout) :: self

    deallocate(self%irows)
    deallocate(self%jcols)
    deallocate(self%values)
  end subroutine sparse_del

  subroutine sparse_matmult(self, vec, res)
    type(sparse_matrix_type), intent(inout) :: self
    real*8, intent(in)                      :: vec(*)
    real*8, intent(out)                     :: res(*)

    integer :: i, j, k

    do j = 1, self%ncols
      res(j) = 0
    end do

    do k = 1, self%nnz
      i = self%irows(k)
      j = self%jcols(k)
      res(i) = res(i) + self%values(k) * vec(j)
    enddo
    
  end subroutine sparse_matmult

end module sparse_matrix_mod
