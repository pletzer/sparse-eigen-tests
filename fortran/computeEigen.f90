program computeEigen
  use sparse_matrix_mod
  implicit none

  ! blas/lapack
  real*8 :: dlapy2, dnrm2
  external dlapy2, dnrm2, daxpy

  character*32 :: filename = '../data/10x10.dat'
  integer, parameter :: num_eig = 5
  character*2, parameter :: which = 'LA' ! largest eigenvalues
  character*1, parameter :: bmat = 'I'
  type(sparse_matrix_type) :: mat
  integer :: info, ido, n
  real*8, parameter :: tol = 0 ! use machine precision
  real*8, allocatable :: resid(:)
  integer, parameter :: ncv = 30 ! max number of basis vectors used by the Implicitely Restarted Arnoldi Process
  real*8, allocatable :: v(:,:), workd(:), workev(:), workl(:)
  integer :: iparam(11)
  integer :: ipntr(14)
  integer :: lworkl
 
  ! read the data
  call sparse_new_from_file(mat, filename)

  ! initialize
  ido = 0
  n = mat%nrows * mat%ncols
  allocate(resid(mat%ncols))
  allocate(v(mat%nrows, ncv))
  allocate(workd(3*mat%nrows))
  allocate(workev(3*ncv))
  lworkl = 3*ncv**2 + 6*ncv
  allocate(workl(lworkl))
  
  ! compute the eigenvalues/vectors
  10 continue

!        %---------------------------------------------%
!        | Repeatedly call the routine DNAUPD and take |
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%

         call dnaupd(ido, bmat, n, which, num_eig, tol, resid, ncv, & 
     &                 v, mat%nrows, iparam, ipntr, workd, workl, lworkl, &
     &                 info )

         if (ido /= -1 .or. ido == 1) then

!           %-------------------------------------------%
!           | Perform matrix vector multiplication      |
!           |                y <--- Op*x                |
!           | The user should supply his/her own        |
!           | matrix vector multiplication routine here |
!           | that takes workd(ipntr(1)) as the input   |
!           | vector, and return the matrix vector      |
!           | product to workd(ipntr(2)).               | 
!           %-------------------------------------------%

            call sparse_matmult(mat, workd(ipntr(1)), workd(ipntr(2)))

!           %-----------------------------------------%
!           | L O O P   B ACK call DNAUPD again. |
!           %-----------------------------------------%

            go to 10

         endif

!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%

      if (info < 0) then
      endif
  

  ! clean up
  call sparse_del(mat)

end program computeEigen
