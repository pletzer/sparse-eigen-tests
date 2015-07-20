      program dssimp

      implicit none


!     %------------------------------------------------------%
!     | Storage Declarations:                                |
!     |                                                      |
!     | The maximum dimensions for all arrays are            |
!     | set here to accommodate a problem size of            |
!     | N .le. MAXN                                          |
!     |                                                      |
!     | NEV is the number of eigenvalues requested.          |
!     |     See specifications for ARPACK usage below.       |
!     |                                                      |
!     | NCV is the largest number of basis vectors that will |
!     |     be used in the Implicitly Restarted Arnoldi      |
!     |     Process.  Work per major iteration is            |
!     |     proportional to N*NCV*NCV.                       |
!     |                                                      |
!     | You must set:                                        |
!     |                                                      |
!     | MAXN:   Maximum dimension of the A allowed.          |
!     | MAXNEV: Maximum NEV allowed.                         |
!     | MAXNCV: Maximum NCV allowed.                         |
!     %------------------------------------------------------%
!
      integer :: maxn, maxnev, maxncv, ldv
!
!     %--------------%
!     | Local Arrays |
!     %--------------%
!
      real*8, allocatable :: &
     &                 v(:, :), workl(:), &
     &                 workd(:), d(:, :), resid(:), &
     &                 ax(:)
      logical, allocatable :: slct(:)
      integer          iparam(11), ipntr(11)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      character        bmat*1, which*2
      integer          ido, n, nev, ncv, lworkl, info, ierr, &
     &                 j, nx, ishfts, maxitr, mode1, nconv
      logical          rvec
      real*8 tol, sigma
!
!     %------------%
!     | Parameters |
!     %------------%
!
      real*8, parameter :: zero = 0
!  
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%
!
      real*8 dnrm2
      external         dnrm2, daxpy
!
!     %--------------------%
!     | Intrinsic function |
!     %--------------------%
!
      intrinsic        abs
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!     %-------------------------------------------------%
!     | The following include statement and assignments |
!     | initiate trace output from the internal         |
!     | actions of ARPACK.  See debug.doc in the        |
!     | DOCUMENTS directory for usage.  Initially, the  |
!     | most useful information will be a breakdown of  |
!     | time spent in the various stages of computation |
!     | given by setting msaupd = 1.                    |
!     %-------------------------------------------------%
!
      integer  logfil, ndigit, mgetv0, &
     &         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,&
     &         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,&
     &         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
      common /debug/ &
     &         logfil, ndigit, mgetv0, &
     &         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,&
     &         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,&
     &         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd

      ndigit = -3
      logfil = 6
      msgets = 0
      msaitr = 0 
      msapps = 0
      msaupd = 1
      msaup2 = 0
      mseigt = 0
      mseupd = 0

!
!     %-----------------------------------------------%
!     |                                               | 
!     | Specifications for ARPACK usage are set       | 
!     | below:                                        |
!     |                                               |
!     |    1) NEV = 4  asks for 4 eigenvalues to be   |  
!     |       computed.                               | 
!     |                                               |
!     |    2) NCV = 20 sets the length of the Arnoldi |
!     |       factorization                           |
!     |                                               |
!     |    3) This is a standard problem              |
!     |         (indicated by bmat  = 'I')            |
!     |                                               |
!     |    4) Ask for the NEV eigenvalues of          |
!     |       largest magnitude                       |
!     |         (indicated by which = 'LM')           |
!     |       See documentation in DSAUPD for the     |
!     |       other options SM, LA, SA, LI, SI.       | 
!     |                                               |
!     | Note: NEV and NCV must satisfy the following  |
!     | conditions:                                   |
!     |              NEV <= MAXNEV                    |
!     |          NEV + 1 <= NCV <= MAXNCV             |
!     %-----------------------------------------------%
!
      nev   = 4
      ncv   = 20 
!     
!     %-------------------------------------------------%
!     | The following sets dimensions for this problem. |
!     %-------------------------------------------------%
!
      nx = 10
      n = nx*nx
      maxnev = nev
      maxncv = ncv
      maxn = n
      ldv=maxn
      
!     %-------------------------------------------------%
!     | Allocate data                                   |
!     %-------------------------------------------------%

      allocate(v(ldv,maxncv))
      allocate(workl(maxncv*(maxncv+8)))
      allocate(workd(3*maxn))
      allocate(d(maxncv,2))
      allocate(resid(maxn))
      allocate(ax(maxn))
      allocate(slct(maxncv))

      bmat  = 'I'
      which = 'LM'
!
      if ( n .gt. maxn ) then
         print *, ' ERROR with _SSIMP: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _SSIMP: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _SSIMP: NCV is greater than MAXNCV '
         go to 9000
      end if
!
!     %-----------------------------------------------------%
!     |                                                     |
!     | Specification of stopping rules and initial         |
!     | conditions before calling DSAUPD                    |
!     |                                                     |
!     | TOL  determines the stopping criterion.             |
!     |                                                     |
!     |      Expect                                         |
!     |           abs(lambdaC - lambdaT) < TOL*abs(lambdaC) |
!     |               computed   true                       |
!     |                                                     |
!     |      If TOL .le. 0,  then TOL <- macheps            |
!     |           (machine precision) is used.              |
!     |                                                     |
!     | IDO  is the REVERSE COMMUNICATION parameter         |
!     |      used to specify actions to be taken on return  |
!     |      from DSAUPD. (See usage below.)                |
!     |                                                     |
!     |      It MUST initially be set to 0 before the first |
!     |      call to DSAUPD.                                | 
!     |                                                     |
!     | INFO on entry specifies starting vector information |
!     |      and on return indicates error codes            |
!     |                                                     |
!     |      Initially, setting INFO=0 indicates that a     | 
!     |      random starting vector is requested to         |
!     |      start the ARNOLDI iteration.  Setting INFO to  |
!     |      a nonzero value on the initial call is used    |
!     |      if you want to specify your own starting       |
!     |      vector (This vector must be placed in RESID.)  | 
!     |                                                     |
!     | The work array WORKL is used in DSAUPD as           | 
!     | workspace.  Its dimension LWORKL is set as          |
!     | illustrated below.                                  |
!     |                                                     |
!     %-----------------------------------------------------%
!
      lworkl = ncv*(ncv+8)
      tol = zero 
      info = 0
      ido = 0
!
!     %---------------------------------------------------%
!     | Specification of Algorithm Mode:                  |
!     |                                                   |
!     | This program uses the exact shift strategy        |
!     | (indicated by setting PARAM(1) = 1).              |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 1 of DSAUPD is used     |
!     | (IPARAM(7) = 1). All these options can be changed |
!     | by the user. For details see the documentation in |
!     | DSAUPD.                                           |
!     %---------------------------------------------------%
!
      ishfts = 1
      maxitr = 300 
      mode1 = 1
!
      iparam(1) = ishfts
!                
      iparam(3) = maxitr
!                  
      iparam(7) = mode1
!
!     %------------------------------------------------%
!     | M A I N   L O O P (Reverse communication loop) |
!     %------------------------------------------------%
!
 10   continue
!
!        %---------------------------------------------%
!        | Repeatedly call the routine DSAUPD and take | 
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
!
         call dsaupd(ido, bmat, n, which, nev, tol, resid, &
     &               ncv, v, ldv, iparam, ipntr, workd, workl, &
     &               lworkl, info)
!
         if (ido .eq. -1 .or. ido .eq. 1) then
!
!           %--------------------------------------%
!           | Perform matrix vector multiplication |
!           |              y <--- OP*x             |
!           | The user should supply his/her own   |
!           | matrix vector multiplication routine |
!           | here that takes workd(ipntr(1)) as   |
!           | the input, and return the result to  |
!           | workd(ipntr(2)).                     |
!           %--------------------------------------%
!
            call av(nx, workd(ipntr(1)), workd(ipntr(2)))
!
!           %-----------------------------------------%
!           | L O O P   B A C K to call DSAUPD again. |
!           %-----------------------------------------%
!
            go to 10
!
         end if 
!
!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%
!
      if ( info .lt. 0 ) then
!
!        %--------------------------%
!        | Error message. Check the |
!        | documentation in DSAUPD. |
!        %--------------------------%
!
         print *, ' '
         print *, ' Error with _saupd, info = ', info
         print *, ' Check documentation in _saupd '
         print *, ' '
!
      else 
!
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using DSEUPD.                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |  
!        |                                           |
!        | Eigenvectors may be also computed now if  |
!        | desired.  (indicated by rvec = .true.)    | 
!        |                                           |
!        | The routine DSEUPD now called to do this  |
!        | post processing (Other modes may require  |
!        | more complicated post processing than     |
!        | mode1.)                                   |
!        |                                           |
!        %-------------------------------------------%
!           
          rvec = .true.
!
          call dseupd (rvec, 'All', slct, d, v, ldv, sigma, & 
     &         bmat, n, which, nev, tol, resid, ncv, v, ldv, &
     &         iparam, ipntr, workd, workl, lworkl, ierr)
!
!         %----------------------------------------------%
!         | Eigenvalues are returned in the first column |
!         | of the two dimensional array D and the       |
!         | corresponding eigenvectors are returned in   |
!         | the first NCONV (=IPARAM(5)) columns of the  |
!         | two dimensional array V if requested.        |
!         | Otherwise, an orthogonal basis for the       |
!         | invariant subspace corresponding to the      |
!         | eigenvalues in D is returned in V.           |
!         %----------------------------------------------%
!
          if (ierr .ne. 0) then
!
!            %------------------------------------%
!            | Error condition:                   |
!            | Check the documentation of DSEUPD. |
!            %------------------------------------%
!
             print *, ' '
             print *, ' Error with _seupd, info = ', ierr
             print *, ' Check the documentation of _seupd. '
             print *, ' '
!
          else
!
             nconv =  iparam(5)
             do 20 j=1, nconv
!
!               %---------------------------%
!               | Compute the residual norm |
!               |                           |
!               |   ||  A*x - lambda*x ||   |
!               |                           |
!               | for the NCONV accurately  |
!               | computed eigenvalues and  |
!               | eigenvectors.  (iparam(5) |
!               | indicates how many are    |
!               | accurate to the requested |
!               | tolerance)                |
!               %---------------------------%
!
                call av(nx, v(1,j), ax)
                call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                d(j,2) = dnrm2(n, ax, 1)
                d(j,2) = d(j,2) / abs(d(j,1))
!
 20          continue
!
!            %-----------------------------%
!            | Display computed residuals. |
!            %-----------------------------%
!
             call dmout(6, nconv, 2, d, maxncv, -6, &
     &            'Ritz values and relative residuals')
          end if
!
!         %-------------------------------------------%
!         | Print additional convergence information. |
!         %-------------------------------------------%
!
          if (info .eq. 1) then
             print *, ' '
             print *, ' Maximum number of iterations reached.'
             print *, ' '
          else if (info .eq. 3) then
             print *, ' ' 
             print *, ' No shifts could be applied during implicit', &
     &                ' Arnoldi update, try increasing NCV.'
             print *, ' '
          end if      
!
          print *, ' '
          print *, ' _SSIMP '
          print *, ' ====== '
          print *, ' '
          print *, ' Size of the matrix is ', n
          print *, ' The number of Ritz values requested is ', nev
          print *, ' The number of Arnoldi vectors generated', &
     &             ' (NCV) is ', ncv
          print *, ' What portion of the spectrum: ', which
          print *, ' The number of converged Ritz values is ', &
     &               nconv 
          print *, ' The number of Implicit Arnoldi update', &
     &             ' iterations taken is ', iparam(3)
          print *, ' The number of OP*x is ', iparam(9)
          print *, ' The convergence criterion is ', tol
          print *, ' '
!
      end if
!
!     %---------------------------%
!     | Done with program dssimp. |
!     %---------------------------%
!
 9000 continue
!
      end
! 
! ------------------------------------------------------------------
!     matrix vector subroutine
!
!     The matrix used is the 2 dimensional discrete Laplacian on unit
!     square with zero Dirichlet boundary condition.
!
!     Computes w <--- OP*v, where OP is the nx*nx by nx*nx block 
!     tridiagonal matrix
!
!                  | T -I          | 
!                  |-I  T -I       |
!             OP = |   -I  T       |
!                  |        ...  -I|
!                  |           -I T|
!
!     The subroutine TV is called to computed y<---T*x.
!
      subroutine av(nx, v, w)
      integer           nx, j, lo, n2
      real*8 v(nx*nx), w(nx*nx), h2
      integer, parameter :: one = 1
!
      call tv(nx,v(1),w(1))
      call daxpy(nx, -one, v(nx+1), 1, w(1), 1)
!
      do 10 j = 2, nx-1
         lo = (j-1)*nx
         call tv(nx, v(lo+1), w(lo+1))
         call daxpy(nx, -one, v(lo-nx+1), 1, w(lo+1), 1)
         call daxpy(nx, -one, v(lo+nx+1), 1, w(lo+1), 1)
  10  continue 
!
      lo = (nx-1)*nx
      call tv(nx, v(lo+1), w(lo+1))
      call daxpy(nx, -one, v(lo-nx+1), 1, w(lo+1), 1)
!
!     Scale the vector w by (1/h^2), where h is the mesh size
!
      n2 = nx*nx
      h2 = one / dble((nx+1)*(nx+1))
      call dscal(n2, one/h2, w, 1)
      return
      end
!
!-------------------------------------------------------------------
      subroutine tv(nx, x, y)
!
      integer           nx, j 
      real*8 &
     &                  x(nx), y(nx), dd, dl, du
!
      real*8, parameter :: one = 1, four = 4
!
!     Compute the matrix vector multiplication y<---T*x
!     where T is a nx by nx tridiagonal matrix with DD on the 
!     diagonal, DL on the subdiagonal, and DU on the superdiagonal.
!     
!
      dd  = four 
      dl  = -one 
      du  = -one
! 
      y(1) =  dd*x(1) + du*x(2)
      do 10 j = 2,nx-1
         y(j) = dl*x(j-1) + dd*x(j) + du*x(j+1) 
 10   continue 
      y(nx) =  dl*x(nx-1) + dd*x(nx) 
      return
      end

