C     =================================================================
C     File: pedazos4.f
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

C     Pedazos4 problem: Fit a continuous piecewise polynomial
C     -------------------------------------------------------

C     We generate a table (t_i, y_i) with 30 equally spaced points 
C     between 1 and 30. The first 10 points are on a line, the second 
C     10 on a parabola and the third 10 on a cubic. The overall 
C     function is continuous at 10 and 20. We fit a line to the first 
C     10, a quadratic to the second 10 and a cubic to the third ten, 
C     with the conditions that we preserve continuity at 10 and 20. So,
C     we have 9 unknowns and 2 equality constraints. The objective 
C     function is quadratic and the constraints are linear. 

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

C     COMMON ARRAYS
      double precision t(30),y(30)

C     LOCAL SCALARS
      integer i
      double precision drand,seed

C     COMMON BLOCKS
      common /probdata/ t,y
      
C     Number of variables

      n = 9

C     Adjust a line, a parabola and a cubic

      do i = 1,10
          t(i) = dfloat(i)
          y(i) = t(i)
      end do

      do i = 11,20
          t(i) = dfloat(i)
          y(i) = 0.4d0 * ( t(i) - 15.0d0 ) ** 2
      end do

      do i = 21,30
          t(i) = dfloat(i)
          y(i) = - 10.0d0 * ( t(i) - 25.0d0 ) ** 3 / 125.0d0
      end do

C     Perturbation

      seed = 17172937.0d0

      do i = 1,30
          y(i) = y(i) +
     +           y(i) * 0.20d0 * ( drand(seed) * 2.0d0 - 1.0d0 )
      end do


C     Initial point

      do i = 1,n
          x(i) = 0.0d0
      end do
      
C     Lower and upper bounds

      do i = 1,n
          l(i) = - 1.0d+20
          u(i) =   1.0d+20
      end do

C     Number of constraints (equalities plus inequalities)

      m = 2

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
          linear(i) = .true.
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

C     COMMON ARRAYS
      double precision t(30),y(30)

C     LOCAL SCALARS
      integer i
      double precision r

C     COMMON BLOCKS
      common /probdata/ t,y
      
      flag = 0

      f = 0.0d0

      do i = 1,30

          if ( t(i) .le. 10.0d0 ) then
              r = ( x(1) + x(2) * t(i) - y(i) ) ** 2

          else if ( t(i) .le. 20.0d0 ) then
              r = ( x(3) + x(4) * t(i) + x(5) * t(i) ** 2 - y(i) ) ** 2

          else
              r = ( x(6) + x(7) * t(i) + x(8) * t(i) ** 2 + 
     +              x(9) * t(i) ** 3 - y(i) ) ** 2
          end if

          f = f + r

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

C     COMMON ARRAYS
      double precision t(30),y(30)

C     LOCAL SCALARS
      integer i
      double precision r

C     COMMON BLOCKS
      common /probdata/ t,y
      
      flag = 0

      do i = 1,n
          g(i) = 0.0d0
      end do

      do i = 1,30

          if ( t(i) .le. 10.0d0 ) then
              r    = x(1) + x(2) * t(i) - y(i)
              g(1) = g(1) + 2.0d0 * r
              g(2) = g(2) + 2.0d0 * r * t(i)

          else if ( t(i) .le. 20.0d0 ) then
              r    = x(3) + x(4) * t(i) + x(5) * t(i) ** 2 - y(i)
              g(3) = g(3) + 2.0d0 * r
              g(4) = g(4) + 2.0d0 * r * t(i)
              g(5) = g(5) + 2.0d0 * r * t(i) ** 2

          else
              r    = x(6) + x(7) * t(i) + x(8) * t(i) ** 2 + 
     +               x(9) * t(i) ** 3 - y(i)
              g(6) = g(6) + 2.0d0 * r
              g(7) = g(7) + 2.0d0 * r * t(i)
              g(8) = g(8) + 2.0d0 * r * t(i) ** 2
              g(9) = g(9) + 2.0d0 * r * t(i) ** 3
          end if

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

      flag = 0

      hnnz = 19 

      hlin(1)  = 1
      hcol(1)  = 1
      hval(1)  = 20.0d0  

      hlin(2)  = 2
      hcol(2)  = 1
      hval(2)  = 110.0d0 

      hlin(3)  = 2
      hcol(3)  = 2
      hval(3)  = 770.0d0  

      hlin(4)  = 3
      hcol(4)  = 3
      hval(4)  = 20.0d0  

      hlin(5)  = 4
      hcol(5)  = 3
      hval(5)  = 310.0d0

      hlin(6)  = 4    
      hcol(6)  = 4
      hval(6)  = 4970.0d0  

      hlin(7)  = 5
      hcol(7)  = 3
      hval(7)  = 4970.0d0

      hlin(8)  = 5     
      hcol(8)  = 4
      hval(8)  = 82150.0d0  

      hlin(9)  = 5
      hcol(9)  = 5
      hval(9)  = 1394666.0d0 

      hlin(10) = 6
      hcol(10) = 6 
      hval(10) = 20.0d0  

      hlin(11) = 7    
      hcol(11) = 6
      hval(11) = 510.0d0  

      hlin(12) = 7 
      hcol(12) = 7 
      hval(12) = 13170.0d0

      hlin(13) = 8
      hcol(13) = 6
      hval(13) = 13170.0d0  

      hlin(14) = 8
      hcol(14) = 7
      hval(14) = 344250.0d0  

      hlin(15) = 8
      hcol(15) = 8
      hval(15) = 9102666.0d0 

      hlin(16) = 9
      hcol(16) = 6
      hval(16) = 344250.0d0  

      hlin(17) = 9
      hcol(17) = 7
      hval(17) = 9102666.0d0  

      hlin(18) = 9
      hcol(18) = 8
      hval(18) = 243308250.0d0  

      hlin(19) = 9
      hcol(19) = 9
      hval(19) = 6.56895081d+09

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

      flag = 0

      if ( ind .eq. 1 ) then
          c = ( x(1) + x(2) * 10.d0 ) - 
     +        ( x(3) + x(4) * 10.d0 + x(5) * 100.d0 )

      else if( ind .eq. 2 ) then
          c = ( x(3) + x(4) * 20.d0 + x(5) * 400.d0 ) - 
     +        ( x(6) + x(7) * 20.d0 + x(8) * 400.d0 + x(9) * 8000.d0 )
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

      flag = 0

      if ( ind .eq. 1 ) then

          jcnnz = 5

          jcvar(1) =        1
          jcval(1) =     1.d0

          jcvar(2) =        2
          jcval(2) =    10.d0

          jcvar(3) =        3
          jcval(3) = -   1.d0

          jcvar(4) =        4
          jcval(4) = -  10.d0

          jcvar(5) =        5
          jcval(5) = - 100.d0

      else if( ind .eq. 2 ) then

          jcnnz = 7

          jcvar(1) =         3
          jcval(1) =      1.d0

          jcvar(2) =         4
          jcval(2) =     20.d0

          jcvar(3) =         5
          jcval(3) =    400.d0

          jcvar(4) =         6
          jcval(4) = -    1.d0

          jcvar(5) =         7
          jcval(5) = -   20.d0

          jcvar(6) =         8
          jcval(6) = -  400.d0

          jcvar(7) =         9
          jcval(7) = - 8000.d0

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

      flag = 0

      hcnnz = 0

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

C     COMMON ARRAYS
      double precision t(30),y(30)

C     LOCAL SCALARS
      integer i
      double precision z,tt

C     COMMON BLOCKS
      common /probdata/ t,y
      
      write(*,*) 't, observed y, computed y'

      do i = 1,30

          if ( t(i) .le. 10.0d0 ) then
              z = x(1) + x(2) * t(i)

          else if ( t(i) .le. 20.0d0 ) then
              z = x(3) + x(4) * t(i) + x(5) * t(i) ** 2

          else
              z = x(6) + x(7) * t(i) + x(8) * t(i) ** 2 + 
     +            x(9) * t(i) ** 3
          end if

          write(*,*) t(i),y(i),z

      end do

      tt = 1.0d0
 10   if ( tt .le. 30.01d0 ) then
          if ( tt .le. 10.0d0 ) then
              z = x(1) + x(2) * tt

          else if ( tt .le. 20.0d0 ) then
              z = x(3) + x(4) * tt + x(5) * tt ** 2

          else
              z = x(6) + x(7) * tt + x(8) * tt ** 2 + 
     +            x(9) * tt ** 3
          end if

          write(*,*) tt,z

          tt = tt + 0.01d0
          go to 10

      end if

      end
