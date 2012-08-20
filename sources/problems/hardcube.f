C     =================================================================
C     File: hardcube.f
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

C     Hardcube problem
C     ----------------

C     We consider the cube [-1,1]^{ndim} intersected with the
C     complement of the unitary sphere. We want to find npun points in
C     this set such that the minimum distance between them is maximal.

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
      integer npumax
      parameter ( npumax = 1000 )

C     COMMON SCALARS
      integer ndim,npun

C     COMMON ARRAYS
      integer pair(2,npumax * ( npumax - 1 ) / 2)

C     LOCAL SCALARS
      integer i,j,k
      double precision seed

C     COMMON BLOCKS
      common /probdata/ pair,ndim,npun

C     EXTERNAL FUNCTIONS
      double precision drand

C     Set problem data

      write(*,*) 'Dimension of the space (ndim) and number of ',
     +           'points (npun): '
      read(*,*) ndim,npun

      write(*,*) 'Seed for the initial point random generation: '
      read(*,*) seed

C     Set pairs

      k = 0
      do i = 1,npun
          do j = i + 1,npun
              k = k + 1
              pair(1,k) = i
              pair(2,k) = j
          end do
      end do

C     Number of variables

      n = ndim * npun + 1

C     Initial point

C     seed = 120927.0d0

      do i = 1,n
          x(i) = 2.0d0 * drand(seed) - 1.0d0
      end do

C     Lower and upper bounds
      do i = 1,n - 1
          l(i) = - 1.0d0
          u(i) =   1.0d0
      end do

      l(n) = - 1.0d+20
      u(n) =   1.0d+20

C     Number of constraints (equalities plus inequalities)

      m = npun + npun * ( npun - 1 ) / 2

C     Lagrange multipliers approximation.

      do i = 1,m
          lambda(i) =  0.0d0
      end do

C     For each constraint i, set equatn(i) = .true. if it is an equality
C     constraint of the form c_i(x) = 0, and set equatn(i) = .false. if
C     it is an inequality constraint of the form c_i(x) <= 0.

      do i = 1,m
          equatn(i) = .false.
      end do

C     For each constraint i, set linear(i) = .true. if it is a linear
C     constraint, otherwise set linear(i) = .false.

      do i = 1,m
          linear(i) = .false.
      end do

C     Indicate which subroutines did you code.

      coded( 1) = .true.  ! evalf
      coded( 2) = .true.  ! evalg
      coded( 3) = .true.  ! evalh
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

      flag = 0

      f = x(n)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalg(n,x,g,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,n

C     ARRAY ARGUMENTS
      double precision g(n),x(n)

C     LOCAL SCALARS
      integer i

      flag = 0

      do i = 1,n - 1
          g(i) = 0.0d0
      end do

      g(n) = 1.0d0

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

      flag = 0

      hnnz = 0

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
      integer npumax
      parameter ( npumax = 1000 )

C     COMMON SCALARS
      integer ndim,npun

C     COMMON ARRAYS
      integer pair(2,npumax * ( npumax - 1 ) / 2)

C     LOCAL SCALARS
      integer i,i1,i2

C     COMMON BLOCKS
      common /probdata/ pair,ndim,npun

      flag = 0

      if ( 0 .lt. ind .and. ind .le. npun ) then

          c = 1.0d0

          do i = 1,ndim
              c = c - x(ndim * ( ind - 1 ) + i) ** 2
          end do

          return

      end if

      if ( npun .lt. ind .and.
     +     ind .le. npun + npun * (npun - 1) / 2 ) then

          i1 = pair(1,ind - npun)
          i2 = pair(2,ind - npun)

          c = - x(n)
          do i = 1,ndim
              c = c - ( x(ndim*(i1-1)+i) - x(ndim*(i2-1)+i) ) ** 2
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
      integer npumax
      parameter ( npumax = 1000 )

C     COMMON SCALARS
      integer ndim,npun

C     COMMON ARRAYS
      integer pair(2,npumax * ( npumax - 1 ) / 2)

C     LOCAL SCALARS
      integer i,i1,i2
      double precision tmp

C     COMMON BLOCKS
      common /probdata/ pair,ndim,npun

      flag = 0

      if ( 0 .lt. ind .and. ind .le. npun ) then

          jcnnz = ndim

          do i = 1,ndim
              jcvar(i) = ndim * ( ind - 1 ) + i
              jcval(i) = - 2.0d0 * x(ndim * ( ind - 1 ) + i)
          end do

          return

      end if

      if ( npun .lt. ind .and.
     +     ind .le. npun + npun * (npun - 1) / 2 ) then

          i1 = pair(1,ind - npun)
          i2 = pair(2,ind - npun)

          jcnnz = 2 * ndim + 1

          do i = 1,ndim
              tmp = 2.0d0 * ( x(ndim*(i1-1)+i) - x(ndim*(i2-1)+i) )
              jcvar(i)      = ndim * ( i1 - 1 ) + i
              jcval(i)      = - tmp
              jcvar(ndim+i) = ndim * ( i2 - 1 ) + i
              jcval(ndim+i) =   tmp
          end do

          jcvar(jcnnz) = n
          jcval(jcnnz) = - 1.0d0

          return

      end if

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

C     PARAMETERS
      integer npumax
      parameter ( npumax = 1000 )

C     COMMON SCALARS
      integer ndim,npun

C     COMMON ARRAYS
      integer pair(2,npumax * ( npumax - 1 ) / 2)

C     LOCAL SCALARS
      integer i,j,k
      double precision z,zmin

C     COMMON BLOCKS
      common /probdata/ pair,ndim,npun

      write(*,*) 'Solution: '
      do i = 1,npun
          write(*,*) i,(x((i-1)*ndim+j),j=1,ndim)
      end do

      zmin = 1.0d+99
      do i = 1,npun - 1
          do j = i + 1,npun
              z = 0.0d0
              do k = 1,ndim
                  z = z + ( x((i-1)*ndim+k) - x((j-1)*ndim+k) ) ** 2
              end do
              z = sqrt(z)
              zmin = min( zmin, z )
          end do
      end do

      write(*,*) 'Minimum distance between point = ',zmin

      end
