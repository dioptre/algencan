C     =================================================================
C     File: simfock2.f
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

C     Simulation of Hartree-Fock problem (version 2)
C     ----------------------------------------------

C     Description:

C     Given kdim and ndim satisfying ndim <= kdim, the problem is:
C
C                  Minimize 0.5 x^T A x + b^T x
C
C                  subject to
C
C                          M orthonornal,
C
C     where x is the columnwise representation of the matrix X = M M^T,
C     and M in R^{kdim \times ndim}. Therefore, the problem has
C     ndim * kdim variables (the elements of matrix M) and
C     ndim + ndim * ( ndim - 1 ) / 2 nonlinear equality constraints.

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

C     PARAMETERS
      integer ndimax
      parameter ( ndimax = 1000 )

C     COMMON SCALARS
      integer kdim,ndim

C     COMMON ARRAYS
      integer pair(2,ndimax * ( ndimax - 1 ) / 2)

C     LOCAL SCALARS
      integer i,j,k
      double precision seed

C     COMMON BLOCKS
      common /probdata/ pair,kdim,ndim

C     EXTERNAL FUNCTIONS
      double precision drand

C     Set problem data

      write(*,*) 'Dimension of the space (kdim) maximum 100: '
      read(*,*) kdim

      write(*,*) 'Rank (ndim) of the desired idempotent solution: ',
     +           '(ndim must be less than or equal to kdim)'
      read(*,*) ndim

C     Set pairs
      k = 0
      do i = 1,ndim
          do j = i + 1,ndim
              k = k + 1
              pair(1,k) = i
              pair(2,k) = j
          end do
      end do

C     Number of variables

      n = ndim * kdim

C     Initial point
      seed = 1324.0d0

      do j = 1,ndim
          do i = 1,kdim
              x((j-1)*kdim+i) = drand(seed)
          end do
      end do

C     Lower and upper bounds
      do i = 1,n
          l(i) = - 1.0d+05
          u(i) =   1.0d+05
      end do

C     Number of constraints (equalities plus inequalities)

      m = ndim + ndim * ( ndim - 1 ) / 2

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

C     PARAMETERS
      integer ndimax
      parameter ( ndimax = 1000 )

C     COMMON SCALARS
      integer kdim,ndim

C     COMMON ARRAYS
      integer pair(2,ndimax * ( ndimax - 1 ) / 2)

C     LOCAL SCALARS
      integer i,j,k,ny
      double precision tmp

C     LOCAL ARRAYS
      double precision y(10000000),gy(1000000)

C     COMMON BLOCKS
      common /probdata/ pair,kdim,ndim

      flag = 0

C     Y = M M^t

C     y is vector which saves Y (kdim \times kdim) in the columnwise
C     representation and x is the same thing for M (kdim \times ndim)

      do j = 1,kdim
          do i = 1,kdim
C             Y(i,j) = row i by row j of matrix M
              tmp = 0.0d0
              do k = 1,ndim
                  tmp = tmp + x((k-1)*kdim+i) * x((k-1)*kdim+j)
              end do
              y((j-1)*kdim+i) = tmp
          end do
      end do

C     Now the quadratic objective function in y

      f = 0.0d0

      do i = 2,kdim - 1
          f = f - y((i-2)*kdim+i) + 2.0d0 * y((i-1)*kdim+i) -
     +        y(i*kdim+i)
      end do

      f = f + 2.0d0 * y(1) - y(kdim+1) - y((kdim-2)*kdim+kdim) +
     +    2.0d0 * y((kdim-1)*kdim+kdim)

      f = 2.0d0 * f

      ny = kdim ** 2
      call gedede(ny,y,gy)

      do i = 1,kdim
          do j = 1,kdim
              f = f + gy((i-1)*kdim+j) * y((i-1)*kdim+j)
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

C     PARAMETERS
      integer ndimax
      parameter ( ndimax = 1000 )

C     COMMON SCALARS
      integer kdim,ndim

C     COMMON ARRAYS
      integer pair(2,ndimax * ( ndimax - 1 ) / 2)

C     LOCAL SCALARS
      integer i,j,k,ny
      double precision tmp

C     LOCAL ARRAYS
      double precision y(1000000),gy(1000000)

C     COMMON BLOCKS
      common /probdata/ pair,kdim,ndim

      flag = 0

C     Y = M M^t

C     y is vector which saves Y (kdim \times kdim) in the columnwise
C     representation and x is the same thing for M (kdim \times ndim)

      do j = 1,kdim
          do i = 1,kdim
C             Y(i,j) = row i by row j of matrix M
              tmp = 0.0d0
              do k = 1,ndim
                  tmp = tmp + x((k-1)*kdim+i) * x((k-1)*kdim+j)
              end do
              y((j-1)*kdim+i) = tmp
          end do
      end do

C     Now the derivatives of the quadratic objective function
C     with respect to y

      ny = kdim ** 2
      call gedede (ny,y,gy)

      do i = 2,kdim - 1
          gy((i-2)*kdim+i) = gy((i-2)*kdim+i) - 1.0d0
          gy((i-1)*kdim+i) = gy((i-1)*kdim+i) + 2.0d0
          gy(i*kdim+i)     = gy(i*kdim+i)     - 1.0d0
      end do

      gy(1) = gy(1) + 2.0d0

      gy((kdim-2)*kdim+kdim) = gy((kdim-2)*kdim+kdim) - 1.0d0
      gy((kdim-1)*kdim+kdim) = gy((kdim-1)*kdim+kdim) + 2.0d0
      gy(kdim+1)             = gy(kdim+1)             - 1.0d0

      do i = 1,ny
          gy(i) = 2.0d0 * gy(i)
      end do

C     Finally, the chain rule to compute the derivatives with respect
C     to the elements of matrix M

      do j = 1,ndim
          do i = 1,kdim
              g((j-1)*kdim+i) = 0.0d0
          end do
      end do

      do j = 1,kdim
          do i = 1,kdim
              tmp = gy((j-1)*kdim+i)
              do k = 1,ndim
                  g((k-1)*kdim+i) =
     +            g((k-1)*kdim+i) + tmp * x((k-1)*kdim+j)
                  g((k-1)*kdim+j) =
     +            g((k-1)*kdim+j) + tmp * x((k-1)*kdim+i)
              end do
          end do
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

C     PARAMETERS
      integer ndimax
      parameter ( ndimax = 1000 )

C     COMMON SCALARS
      integer kdim,ndim

C     COMMON ARRAYS
      integer pair(2,ndimax * ( ndimax - 1 ) / 2)

C     LOCAL SCALARS
      integer i,i1,i2

C     COMMON BLOCKS
      common /probdata/ pair,kdim,ndim

      flag = 0

      if ( ind .le. ndim ) then

C         Euclidian norm of column ind of matrix M equal to 1

          c = - 1.0d0
          do i = 1,kdim
              c = c + x((ind-1)*kdim+i) ** 2
          end do

          return

      end if

      if ( ind .le. ndim + ndim * ( ndim - 1 ) / 2 ) then

C         Given the number of the constraint, identify the indices
C         i1 and i2 of both columns of matrix M

          i  = ind - ndim
          i1 = pair(1,i)
          i2 = pair(2,i)

C         Orthogonality between columns i1 and i2 of matrix M

          c = 0.0d0
          do i = 1,kdim
              c = c + x((i1-1)*kdim+i) * x((i2-1)*kdim+i)
          end do

          return

      end if

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

C     PARAMETERS
      integer ndimax
      parameter ( ndimax = 1000 )

C     COMMON SCALARS
      integer kdim,ndim

C     COMMON ARRAYS
      integer pair(2,ndimax * ( ndimax - 1 ) / 2)

C     LOCAL SCALARS
      integer i,i1,i2

C     COMMON BLOCKS
      common /probdata/ pair,kdim,ndim

      flag = 0

      if ( ind .le. ndim ) then

C         Euclidian norm of column ind of matrix M equal to 1

          jcnnz = kdim

          do i = 1,kdim
              jcvar(i) = ( ind - 1 ) * kdim + i
              jcval(i) = 2.0d0 * x((ind-1)*kdim+i)
          end do

          return

      end if

      if ( ind .le. ndim + ndim * ( ndim - 1 ) / 2 ) then

C         Given the number of the constraint, identify the indices
C         i1 and i2 of both columns of matrix M

          i  = ind - ndim
          i1 = pair(1,i)
          i2 = pair(2,i)

C         Orthogonality between columns i1 and i2 of matrix M

          jcnnz = 2 * kdim

          do i = 1,kdim
              jcvar(i) = ( i1 - 1 ) * kdim + i
              jcval(i) = x((i2-1)*kdim+i)
              jcvar(kdim+i) = ( i2 - 1 ) * kdim + i
              jcval(kdim+i) = x((i1-1)*kdim+i)
          end do

          return

      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine gedede(ny,y,gy)

      implicit none

C     SCALAR ARGUMENTS
      integer ny

C     ARRAY ARGUMENTS
      double precision y(ny),gy(ny)

C     LOCAL SCALARS
      integer i,j,jj

      do i = 1,ny

          gy(i) = y(i)

          do jj = 1,10

              j = i + jj

              if ( j .le. ny ) then
                  gy(i) = gy(i) + 1.0d0 / dfloat(i+j) * y(j)
              end if

              j = i - jj

              if ( j .ge. 1 ) then
                  gy(i) = gy(i) + 1.0d0 / dfloat(i+j) * y(j)
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
