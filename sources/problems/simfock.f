C     =================================================================
C     File: simfock.f
C     =================================================================

C     =================================================================
C     Module: Subroutines that define the problem
C     =================================================================

C     Last update of any of the component of this module:

C     March 4, 2008.

C     Users are encouraged to download periodically updated versions of
C     this code at the COLLECTION home page:
C
C     www.ime.usp.br/~egbirgin/collection/
C
C     and periodically updated versions of the TANGO Project solvers at
C     the TANGO home page:
C
C     www.ime.usp.br/~egbirgin/tango/

C     =================================================================

C     Simulation of Hartree-Fock problem
C     ----------------------------------

C     Given kdim, ndim in 1, 2, 3, ..., ndim <= kdim, n = kdim * kdim,
C     x in R^n, we define mat(x) in R^{kdim x kdim} as the matrix whose
C     columnwise representation is x. The problem is:
C
C                  Minimize 0.5 x^T A x + b^T x
C
C                  subject to
C
C                          mat(x) mat(x) = mat(x)
C
C                          trace[mat(x)] = ndim
C
C                          mat(x) = mat(x)^T
C
C     Therefore, the problem has kdim^2 variables, kdim^2 nonlinear
C     equality constraints and kdim * ( kdim - 1 ) / 2 + 1 linear
C     equality constraints.

C     ******************************************************************
C     ******************************************************************

      subroutine inip(n,x,l,u,m,lambda,equatn,linear,coded,checkder)

      implicit none

C     SCALAR ARGUMENTS
      integer m,n
      logical checkder

C     ARRAY ARGUMENTS
      logical coded(10),equatn(*),linear(*)
      double precision l(*),lambda(*),u(*),x(*)

C     COMMON SCALARS
      integer kdim,ndim

C     LOCAL SCALARS
      integer i,j

C     COMMON BLOCKS
      common /probdata/ kdim,ndim

C     Set problem data

      write(*,*) 'Dimension of the space (kdim) maximum 100: '
      read(*,*) kdim

      write(*,*) 'Rank (ndim) of the desired idempotent solution: ',
     +           '(ndim must be less than or equal to kdim)'
      read(*,*) ndim

C     Number of variables

      n = kdim ** 2

C     Initial point
      do i = 1,kdim
          do j = 1,kdim
              x((j-1)*kdim+i) = 1.0d0
          end do
      end do

C     Lower and upper bounds
      do i = 1,n
          l(i) = - 1.0d+05
          u(i) =   1.0d+05
      end do

C     Number of constraints (equalities plus inequalities)

      m = 2 * kdim ** 2 + 1

C     Lagrange multipliers approximation.

      do i = 1,m
          lambda(i) =  0.0d0
      end do

C     For each constraint i, set equatn(i) = .true. if it is an equality
C     constraint of the form c_i(x) = 0, and set equatn(i) = .false. if
C     it is an inequality constraint of the form c_i(x) <= 0.

      do i = 1,m
          equatn(i) = .true.
      end do

C     For each constraint i, set linear(i) = .true. if it is a linear
C     constraint, otherwise set linear(i) = .false.

      do i = 1,m
          linear(i) = .false.
      end do

C     Indicate which subroutines did you code.

      coded( 1) = .true.  ! evalf
      coded( 2) = .true.  ! evalg
      coded( 3) = .false. ! evalh
      coded( 4) = .true.  ! evalc
      coded( 5) = .true.  ! evaljac
      coded( 6) = .false. ! evalhc
      coded( 7) = .false. ! evalfc
      coded( 8) = .false. ! evalgjac
      coded( 9) = .false. ! evalhl
      coded(10) = .false. ! evalhlp

C     Set checkder = .true. if you code some derivatives and you would
C     like them to be tested by finite differences. It is highly
C     recommended.

      checkder = .false.

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalf(n,x,f,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,n
      double precision f

C     ARRAY ARGUMENTS
      double precision x(n)

C     COMMON SCALARS
      integer kdim,ndim

C     LOCAL SCALARS
      integer i,j

C     LOCAL ARRAYS
      double precision g(10000000)

C     COMMON BLOCKS
      common /probdata/ kdim,ndim

      flag = 0

      f = 0.0d0

      do i = 2,kdim - 1
          f = f - x((i-2)*kdim+i) + 2.0d0 * x((i-1)*kdim+i) -
     +        x(i*kdim+i)
      end do

      f = f + 2.0d0 * x(1) - x(kdim+1) - x((kdim-2)*kdim+kdim) +
     +    2.0d0 * x((kdim-1)*kdim+kdim)

      f = 2.0d0 * f

      call gedede(n,x,g)

      do i = 1,kdim
          do j = 1,kdim
              f = f + g((i-1)*kdim+j) * x((i-1)*kdim+j)
          end do
      end do

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalg(n,x,g,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,n

C     ARRAY ARGUMENTS
      double precision g(n),x(n)

C     COMMON SCALARS
      integer kdim,ndim

C     LOCAL SCALARS
      integer i

C     COMMON BLOCKS
      common /probdata/ kdim,ndim

      flag = 0

      call gedede (n,x,g)

      do i = 2,kdim - 1
          g((i-2)*kdim+i) = g((i-2)*kdim+i) - 1.0d0
          g((i-1)*kdim+i) = g((i-1)*kdim+i) + 2.0d0
          g(i*kdim+i)     = g(i*kdim+i)     - 1.0d0
      end do

      g(1) = g(1) + 2.0d0

      g((kdim-2)*kdim+kdim) = g((kdim-2)*kdim+kdim) - 1.0d0
      g((kdim-1)*kdim+kdim) = g((kdim-1)*kdim+kdim) + 2.0d0
      g(kdim+1)             = g(kdim+1)             - 1.0d0

      do i = 1,n
          g(i) = 2.0d0 * g(i)
      end do

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalh(n,x,hlin,hcol,hval,hnnz,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,hnnz,n

C     ARRAY ARGUMENTS
      integer hcol(*),hlin(*)
      double precision hval(*),x(n)

      flag = - 1

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalc(n,x,ind,c,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer ind,flag,n
      double precision c

C     ARRAY ARGUMENTS
      double precision x(n)

C     COMMON SCALARS
      integer kdim,ndim

C     LOCAL SCALARS
      integer i,j,k,iminus

C     COMMON BLOCKS
      common /probdata/ kdim,ndim

      flag = 0

      if ( ind .le. kdim ** 2 ) then

          j = mod(ind,kdim)
          if ( j .eq. 0 ) then
              j = kdim
          end if

          i = ( ind - j ) / kdim + 1

          c = - x((j-1)*kdim+i)
          do k = 1,kdim
              c = c + x((k-1)*kdim+i) * x((j-1)*kdim+k)
          end do

          return

      end if

      if ( ind .le. 2 * kdim ** 2 ) then

          iminus = ind - kdim ** 2

          j = mod(iminus,kdim)
          if( j .eq. 0 ) then
              j = kdim
          end if

          i = ( iminus - j ) / kdim + 1

          c =  x((j-1)*kdim+i) - x((i-1)*kdim+j)

          return

      end if

      c = - dfloat( ndim )

      do i = 1,kdim
          c = c + x((i-1)*kdim+i)
      end do

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evaljac(n,x,ind,jcvar,jcval,jcnnz,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,ind,jcnnz,n

C     ARRAY ARGUMENTS
      integer jcvar(n)
      double precision x(n),jcval(n)

C     COMMON SCALARS
      integer kdim,ndim

C     LOCAL SCALARS
      integer i,j,k,iminus

C     COMMON BLOCKS
      common /probdata/ kdim,ndim

      flag = 0

      if ( ind .le. kdim ** 2 ) then

          jcnnz = 2 * kdim + 1

          j = mod(ind,kdim)
          if ( j .eq. 0 ) then
              j = kdim
          end if

          i = ( ind - j ) / kdim + 1

          do k = 1,kdim
              jcvar(k) = ( j - 1 ) * kdim + k
              jcval(k) = x((k-1)*kdim+i)
              jcvar(kdim+k) = ( k - 1 ) * kdim + i
              jcval(kdim+k) = x((j-1)*kdim+k)
          end do

          jcvar(2*kdim+1) = ( j - 1 ) * kdim + i
          jcval(2*kdim+1) = - 1.0d0

          return

      end if

      if ( ind .le. 2 * kdim ** 2 ) then

          jcnnz = 2

          iminus = ind - kdim ** 2

          j = mod(iminus,kdim)
          if ( j .eq. 0 ) then
              j = kdim
          end if

          i = ( iminus - j ) / kdim + 1

          jcvar(1) = ( j - 1 ) * kdim + i
          jcval(1) =   1.0d0
          jcvar(2) = ( i - 1 ) * kdim + j
          jcval(2) = - 1.0d0

          return

      end  if

      jcnnz = kdim

      do i = 1,kdim
          jcvar(i) = ( i - 1 ) * kdim + i
          jcval(i) = 1.0d0
      end do

      end

C     ******************************************************************
C     ******************************************************************

      subroutine gedede(n,x,g)

      implicit none

C     SCALAR ARGUMENTS
      integer n

C     ARRAY ARGUMENTS
      double precision x(n),g(n)

C     LOCAL SCALARS
      integer i,j,jj

      do i = 1,n

          g(i) = x(i)

          do jj = 1,10

              j = i + jj

              if ( j .le. n ) then
                  g(i) = g(i) + 1.0d0 / dfloat(i+j) * x(j)
              end if

              j = i - jj

              if ( j .ge. 1 ) then
                  g(i) = g(i) + 1.0d0 / dfloat(i+j) * x(j)
              end if

          end do

      end do

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalhc(n,x,ind,hclin,hccol,hcval,hcnnz,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,hcnnz,ind,n

C     ARRAY ARGUMENTS
      integer hccol(*),hclin(*)
      double precision hcval(*),x(n)

      flag = - 1

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalfc(n,x,f,m,c,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,m,n
      double precision f

C     ARRAY ARGUMENTS
      double precision c(m),x(n)

      flag = - 1

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,jcnnz,m,n

C     ARRAY ARGUMENTS
      integer jcfun(*),jcvar(*)
      double precision g(n),jcval(*),x(n)

      flag = - 1

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalhl(n,x,m,lambda,scalef,scalec,hllin,hlcol,hlval,
     +hlnnz,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,hlnnz,m,n
      double precision scalef

C     ARRAY ARGUMENTS
      integer hlcol(*),hllin(*)
      double precision hlval(*),lambda(m),scalec(m),x(n)

      flag = - 1

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalhlp(n,x,m,lambda,scalef,scalec,p,hp,goth,flag)

      implicit none

C     SCALAR ARGUMENTS
      logical goth
      integer flag,m,n
      double precision scalef

C     ARRAY ARGUMENTS
      double precision hp(n),lambda(m),p(n),scalec(m),x(n)

      flag = - 1

      end

C     ******************************************************************
C     ******************************************************************

      subroutine endp(n,x,l,u,m,lambda,equatn,linear)

      implicit none

C     SCALAR ARGUMENTS
      integer m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision l(n),lambda(m),u(n),x(n)

      end
