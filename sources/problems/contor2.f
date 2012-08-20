C     =================================================================
C     File: contor2.f
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

C     Contor problem 2
C     ----------------
C
C     Find the solution of the differential equation
C     y'' = a y^2 + b y + c which is closest to a set of given data
C     points, considering y, a, b, and c as unknowns. The data (ti, yi)
C     are such that yi is a random [-1, 1] perturbation of a true
C     discrete solution. The true discrete solution is such that
C     y(0) = y'(0) = 1.

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
      integer npunmax
      parameter ( npunmax = 10000 )

C     COMMON SCALARS
      integer npun
      double precision h

C     COMMON ARRAYS
      double precision y(npunmax)

C     LOCAL SCALARS
      integer i
      double precision deri,seed,tmp

C     COMMON BLOCKS
      common /probdata/ y,h,npun

C     EXTERNAL FUNCTIONS
      double precision drand

C     Set problem data

      write(*,*) 'Number of points, including extremes: '
      read(*,*) npun

      h = 1.0d0 / dfloat(npun - 1)

C     Generation of the solution
      x(1) = 1.0d0
      deri = 1.0d0
      x(2) = x(1) + deri * h
      do i = 3,npun
          deri = x(i-1) ** 2
          x(i) = h ** 2 * deri + 2.0d0 * x(i-1) - x(i-2)
      end do

C     Invention of data
      seed = 173437.0d0
      do i = 1,npun
          tmp = 2.0d0 * drand(seed) - 1.0d0
          y(i) = x(i) + tmp
      end do

C     Number of variables
      n = npun + 3

C     Initial point
      do i = 1,npun
          x(i) = y(i)
      end do
      x(npun+1) = drand(seed)
      x(npun+2) = drand(seed)
      x(npun+3) = drand(seed)

C     Lower and upper bounds
      do i = 1,n
          l(i) = - 1.0d+05
          u(i) =   1.0d+05
      end do

C     Number of constraints (equalities plus inequalities)

      m = npun - 2

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
      integer npunmax
      parameter ( npunmax = 10000 )

C     COMMON SCALARS
      integer npun
      double precision h

C     COMMON ARRAYS
      double precision y(npunmax)

C     LOCAL SCALARS
      integer i

C     COMMON BLOCKS
      common /probdata/ y,h,npun

      flag = 0

      f = 0.0d0
      do i = 1,npun
          f = f + ( y(i) - x(i) ) ** 2
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
      integer npunmax
      parameter ( npunmax = 10000 )

C     COMMON SCALARS
      integer npun
      double precision h

C     COMMON ARRAYS
      double precision y(npunmax)

C     LOCAL SCALARS
      integer i

C     COMMON BLOCKS
      common /probdata/ y,h,npun

      flag = 0

      do i = 1,npun
          g(i) = 2.0d0 * ( x(i) - y(i) )
      end do

      do i = npun + 1,n
          g(i) = 0.0d0
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
      integer npunmax
      parameter ( npunmax = 10000 )

C     COMMON SCALARS
      integer npun
      double precision h

C     COMMON ARRAYS
      double precision y(npunmax)

C     COMMON BLOCKS
      common /probdata/ y,h,npun

      flag = 0

      c = ( x(ind+2) - 2.0d0 * x(ind+1) + x(ind) ) / h ** 2 -
     +    ( x(npun+1) * x(ind+1) ** 2 + x(npun+2) * x(ind+1) +
     +      x(npun+3) )

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
      integer npunmax
      parameter ( npunmax = 10000 )

C     COMMON SCALARS
      integer npun
      double precision h

C     COMMON ARRAYS
      double precision y(npunmax)

C     COMMON BLOCKS
      common /probdata/ y,h,npun

      flag = 0

      jcnnz = 6

      jcvar(1) = ind + 2
      jcval(1) = 1.0d0 / h ** 2

      jcvar(2) = ind + 1
      jcval(2) = - 2.0d0 / h ** 2 -
     +            ( 2.0d0 * x(ind+1) * x(npun+1) + x(npun+2) )

      jcvar(3) = ind
      jcval(3) = 1.0d0 / h ** 2

      jcvar(4) = npun + 1
      jcval(4) = - x(ind+1) ** 2

      jcvar(5) = npun + 2
      jcval(5) = - x(ind+1)

      jcvar(6) = npun + 3
      jcval(6) = - 1.0d0

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
