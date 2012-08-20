C     =================================================================
C     File: mountain1.f
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

C     Mountain Pass problem
C     ---------------------

C     Given a surface S(x,y) and initial and final points pi, pf in R^2
C     the problem consists on finding a path pi, p1, p2, ..., pN, pf
C     from pi to pf such that max_{1 <= k <= N} S(pk) is as small as
C     possible. Moreover, the distance between consecutive points in the
C     path must be less than or equal to a prescribed tolerance.
C
C     The nonlinear programming formulation of the problem follows:
C
C     min z
C
C     subject to
C
C     d(pi     ,p1)^2 <= dmax^2
C     d(p_{k-1},pk)^2 <= dmax^2, k = 2, ..., N
C     d(pN     ,pf)^2 <= dmax^2
C
C     S(pk) <= z, k = 1, 2, ..., N
C
C     where dmax is the prescribed maximum distance and d(.,.) is the
C     Euclidian distance.
C
C     The problem has n = 2 N + 1 variables and m = 2 N + 1 inequality
C     constraints, where N is the number of intermediate points in the
C     path.
C
C     In the implementation, pi = (-10,-10), pf = (10,10). N is a user
C     defined parameter. dmax is also defined by the user and must
C     satisfy dmax >= d(pi,pf) / (N + 1). The initial values for p1, p2,
C     ..., pN correspond to a perturbation of the equally spaced points
C     in the segment [pi,pf]. The intial value of z is z = max { S(pi),
C     S(p1), S(p2), ..., S(pN), S(pf) }. The considered surfaces are
C
C     Mountain 1: S(x,y) = sin( x * y ) + sin( x + y )
C
C     Mountain 2: S(x,y) = sin(x) + cos(y)

C     References:
C
C     [1] J. Horák, Constrained mountain pass algorithm for the
C     numerical solution of semilinear elliptic problems, Numerische
C     Mathematik 98, pp. 251-276, 2004.
C
C     [2] J. J. Moré and T. S. Munson, Computing mountain passes and
C     transition states, Mathematical Programming 100, pp. 151-182, 2004.

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
      integer npmax
      parameter ( npmax = 100000 )

C     COMMON SCALARS
      integer np
      double precision dmax2

C     COMON ARRAYS
      double precision pi(npmax),pf(npmax)

C     LOCAL SCALARS
      integer b,i,j
      double precision dist,drand,fmoun,seed,tmp

C     COMMON BLOCKS
      common /probdata/ pi,pf,dmax2,np

C     Set problem data

      pi(1) = - 10.0d0
      pi(2) = - 10.0d0

      pf(1) =   10.0d0
      pf(2) =   10.0d0

      dist = sqrt( ( pi(1) - pf(1) ) ** 2 + ( pi(2) - pf(2) ) ** 2 )

      write(*,*) 'The distance between the extremes of the path is ',
     +            dist,'.'

      write(*,*) 'Number of points in the path excluding extremes ',
     +           '(maximum 100,000): '
      read(*,*) np

      dist = dist / ( np + 1 )

      write(*,*) 'The minimum distance between consecutive points in ',
     +           'the path is ',dist,'.'

      write(*,100)

      read(*,*) dmax2

      if ( dmax2 .lt. 0.0d0 ) then
          dmax2 = 2.0d0 * dist
      end if

      dmax2 = dmax2 ** 2

      write(*,*) 'Seed for the initial point random generation: '
      read(*,*) seed

C     Number of variables

      n = 2 * np + 1

C     Initial point
      x(n) = max( fmoun(pi(1),pi(2)), fmoun(pf(1),pf(2)) )
      do i = 1,np
          b = 2 * ( i - 1 )
          do j = 1,2
              tmp    = pi(j) + i * ( pf(j) - pi(j) ) / ( np + 1 )
              x(b+j) = tmp * ( 1.0d0 + 0.1d0 * drand(seed) )
          end do
          x(n) = max( x(n), fmoun(x(b+1),x(b+2)) )
      end do

C     Lower and upper bounds
      do i = 1,n
          l(i) = - 1.0d+20
          u(i) =   1.0d+20
      end do

C     Number of constraints (equalities plus inequalities)

      m = np + 1 + np

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
      coded( 6) = .true.  ! evalhc
      coded( 7) = .false. ! evalfc
      coded( 8) = .false. ! evalgjac
      coded( 9) = .false. ! evalhl
      coded(10) = .false. ! evalhlp

C     Set checkder = .true. if you code some derivatives and you would
C     like them to be tested by finite differences. It is highly
C     recommended.

      checkder = .false.

 100  format(/,1X,'Enter the maximum allowed distance between ',
     +         'consecutive points in the path',/,1X,'(It must be ',
     +         'greater than the minimum distance. Enter a negative ',
     +         'number to',/,1X,'use twice the minimum distance): ')

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
      integer npmax
      parameter ( npmax = 100000 )

C     COMMON SCALARS
      integer np
      double precision dmax2

C     COMON ARRAYS
      double precision pi(npmax),pf(npmax)

C     COMMON BLOCKS
      common /probdata/ pi,pf,dmax2,np

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

C     PARAMETERS
      integer npmax
      parameter ( npmax = 100000 )

C     COMMON SCALARS
      integer np
      double precision dmax2

C     COMON ARRAYS
      double precision pi(npmax),pf(npmax)

C     COMMON BLOCKS
      common /probdata/ pi,pf,dmax2,np

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
      integer npmax
      parameter ( npmax = 100000 )

C     COMMON SCALARS
      integer np
      double precision dmax2

C     COMON ARRAYS
      double precision pi(npmax),pf(npmax)

C     LOCAL SCALARS
      integer b,b1,b2
      double precision fmoun

C     COMMON BLOCKS
      common /probdata/ pi,pf,dmax2,np

      flag = 0

      if ( ind .eq. 1 ) then

          b1 = 2 * ( ind - 1 )
          c  = (x(b1+1) - pi(1))**2   + (x(b1+2) - pi(2))**2   - dmax2

      else if ( ind .le. np ) then

          b1 = 2 * ( ind - 1 )
          b2 = 2 * ( ind - 2 )
          c  = (x(b1+1) - x(b2+1))**2 + (x(b1+2) - x(b2+2))**2 - dmax2

      else if ( ind .eq. np + 1 ) then

          b2 =  2 * ( ind - 2 )
          c  = (pf(1) - x(b2+1))**2   + (pf(2) - x(b2+2))**2   - dmax2

      else if ( ind .le. np + 1 + np ) then

          b = 2 * ( ind - np - 2 )
          c = fmoun(x(b+1),x(b+2)) - x(n)

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
      integer npmax
      parameter ( npmax = 100000 )

C     COMMON SCALARS
      integer np
      double precision dmax2

C     COMON ARRAYS
      double precision pi(npmax),pf(npmax)

C     LOCAL SCALARS
      integer b,b1,b2,j

C     COMMON BLOCKS
      common /probdata/ pi,pf,dmax2,np

      flag = 0

      if ( ind .eq. 1 ) then

          jcnnz = 2

          b1 = 2 * ( ind - 1 )

          do j = 1,2
              jcvar(j) = b1 + j
              jcval(j) = 2.0d0 * ( x(b1 + j) - pi(j) )
          end do

      else if ( ind .le. np ) then

          jcnnz = 4

          b1 = 2 * ( ind - 1 )
          b2 = 2 * ( ind - 2 )

          do j = 1,2
              jcvar(j  ) = b1 + j
              jcval(j  ) =   2.0d0 * ( x(b1 + j) - x(b2 + j) )
              jcvar(2+j) = b2 + j
              jcval(2+j) = - 2.0d0 * ( x(b1 + j) - x(b2 + j) )
          end do

      else if ( ind .eq. np + 1 ) then

          jcnnz = 2

          b2 = 2 * ( ind - 2 )

          do j = 1,2
              jcvar(j) = b2 + j
              jcval(j) = - 2.0d0 * ( pf(j) - x(b2 + j) )
          end do

      else if ( ind .le. np + 1 + np ) then

          jcnnz = 3

          b = 2 * ( ind - np - 2 )

          call gmoun(x(b+1),x(b+2),jcval(1),jcval(2))

          jcvar(1) = b + 1
          jcvar(2) = b + 2

          jcvar(3) = n
          jcval(3) = - 1.0d0

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

C     PARAMETERS
      integer npmax
      parameter ( npmax = 100000 )

C     COMMON SCALARS
      integer np
      double precision dmax2

C     COMON ARRAYS
      double precision pi(npmax),pf(npmax)

C     LOCAL SCALARS
      integer b,b1,b2,j

C     COMMON BLOCKS
      common /probdata/ pi,pf,dmax2,np

      flag = 0

      if ( ind .eq. 1 ) then

          hcnnz = 2

          b1 = 2 * ( ind - 1 )

          do j = 1,2
              hclin(j) = b1 + j
              hccol(j) = b1 + j
              hcval(j) = 2.0d0
          end do

      else if ( ind .le. np ) then

          hcnnz = 6

          b1 = 2 * ( ind - 1 )
          b2 = 2 * ( ind - 2 )

          do j = 1,2
              hclin(j)   = b1 + j
              hccol(j)   = b1 + j
              hcval(j)   = 2.0d0

              hclin(2+j) = b2 + j
              hccol(2+j) = b2 + j
              hcval(2+j) = 2.0d0

              hclin(4+j) = b1 + j
              hccol(4+j) = b2 + j
              hcval(4+j) = - 2.0d0
          end do

      else if ( ind .eq. np + 1 ) then

          hcnnz = 2

          b2 = 2 * ( ind - 2 )

          do j = 1,2
              hclin(j) = b2 + j
              hccol(j) = b2 + j
              hcval(j) = 2.0d0
          end do

      else if ( ind .le. np + 1 + np ) then

          hcnnz = 3

          b = 2 * ( ind - np - 2 )

          call hmoun(x(b+1),x(b+2),hcval(1),hcval(2),hcval(3))

          hclin(1) = b + 1
          hccol(1) = b + 1

          hclin(2) = b + 2
          hccol(2) = b + 2

          hclin(3) = b + 2
          hccol(3) = b + 1

      end if

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
      integer npmax
      parameter ( npmax = 100000 )

C     COMMON SCALARS
      integer np
      double precision dmax2

C     COMON ARRAYS
      double precision pi(npmax),pf(npmax)

C     LOCAL SCALARS
      integer b,i,imax
      double precision f,fmax,fmoun

C     LOCAL ARRAYS
      double precision g(2),h(2,2)

C     COMMON BLOCKS
      common /probdata/ pi,pf,dmax2,np

      write(*,*)

      fmax = - 1.0d+99
      do i = 1,np
          b = 2 * ( i - 1 )
          f = fmoun(x(b+1),x(b+2))
          if ( f .gt. fmax ) then
              fmax = f
              imax = i
          end if
      end do

      b = 2 * ( imax - 1 )

      write(*,*) 'Maximizer on optimal path: ',x(b+1),x(b+2)
      write(*,*) 'Maximum: ',fmax

      call gmoun(x(b+1),x(b+2),g(1),g(2))
      write(*,*) 'Gradient of the energy: ',g(1),g(2)

      call hmoun(x(b+1),x(b+2),h(1,1),h(2,2),h(1,2))
      write(*,*) 'Hessian of the energy: '

      write(*,*) 'Determinant: ',h(1,1) * h(2,2) - h(1,2) ** 2

      if ( imax .ne. 1 ) then
          b = 2 * ( imax - 2 )
          write(*,*) 'Previous point in the path: ',x(b+1),x(b+2)
          write(*,*) 'Energy at previous point in the path: ',
     +                fmoun(x(b+1),x(b+2))
      end if

      if ( imax .ne. np ) then
          b = 2 * imax
          write(*,*) 'Following point in the path: ',x(b+1),x(b+2)
          write(*,*) 'Energy at following point in the path: ',
     +                fmoun(x(b+1),x(b+2))
      end if

      end

C     ******************************************************************
C     ******************************************************************

      double precision function fmoun(x,y)

      implicit none

C     SCALAR ARGUMENTS
      double precision x,y

      fmoun = sin( x * y ) + sin( x + y )

      return

      end

C     ******************************************************************
C     ******************************************************************

      subroutine gmoun(x,y,dfdx,dfdy)

      implicit none

C     SCALAR ARGUMENTS
      double precision x,y,dfdx,dfdy

      dfdx = y * cos( x * y ) + cos( x + y )
      dfdy = x * cos( x * y ) + cos( x + y )

      end

C     ******************************************************************
C     ******************************************************************

      subroutine hmoun(x,y,dfdxx,dfdyy,dfdxy)

      implicit none

C     SCALAR ARGUMENTS
      double precision x,y,dfdxx,dfdyy,dfdxy

      dfdxx = - y ** 2 * sin( x * y ) - sin( x + y )
      dfdyy = - x ** 2 * sin( x * y ) - sin( x + y )
      dfdxy = cos( x * y ) - y * x * sin( x * y ) - sin( x + y )

      end
