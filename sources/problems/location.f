
C     =================================================================
C     File: location.f
C     =================================================================

C     =================================================================
C     Module: Subroutines that define the problem
C     =================================================================

C     Last update of any of the component of this module: 
 
C     September 4, 2007.

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

C     Location problem
C     ----------------

C     This is a variant of the family of location problems introduced 
C     in [1]. In the original problem, given a set of np disjoint 
C     polygons P1, P2, ..., Pnp in R^2 one wishes to find the point z1 
C     in P1 that minimizes the sum of the distances to the other 
C     polygons. In this variant [2], we have, in addition to the np 
C     polygons, nc circles. Moreover, there is an ellipse which has a 
C     non empty intersection with P1 and such that z1 must be inside 
C     the ellipse and zi, i = 2, ..., np+nc must be outside. The 
C     detailed formulation can be found in [2].
C
C     [1] E. G. Birgin, J. M. Martínez and M. Raydan, "Algorithm 813: 
C     SPG - software for convex-constrained optimization", ACM 
C     Transactions on Mathematical Software 27, pp. 340-349, 2001.
C
C     [2] R. Andreani, E. G. Birgin, J. M. Martínez and M. L. Schuverdt, 
C     "On Augmented Lagrangian methods with general lower-level 
C     constraints", submitted, 2005.

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

C     This subroutine must set some problem data. For achieving this 
C     objective YOU MUST MODIFY it according to your problem. See below 
C     where your modifications must be inserted.
C     
C     Parameters of the subroutine:
C
C     On Entry:
C
C     This subroutine has no input parameters.
C
C     On Return
C
C     n        integer,
C              number of variables,
C
C     x        double precision x(n),
C              initial point,
C
C     l        double precision l(n),
C              lower bounds on x,
C
C     u        double precision u(n),
C              upper bounds on x,
C
C     m        integer,
C              number of constraints (excluding the bounds),
C
C     lambda   double precision lambda(m),
C              initial estimation of the Lagrange multipliers,
C
C     equatn   logical equatn(m)
C              for each constraint j, set equatn(j) = .true. if it is an 
C              equality constraint of the form c_j(x) = 0, and set 
C              equatn(j) = .false. if it is an inequality constraint of 
C              the form c_j(x) <= 0,
C
C     linear   logical linear(m)
C              for each constraint j, set linear(j) = .true. if it is a 
C              linear constraint, and set linear(j) = .false. if it is a
C              nonlinear constraint.

C     PARAMETERS
      integer mmax,ncmax,nmax,npmax,nvsmax
      parameter ( ncmax    =                  1000 )
      parameter ( npmax    =                  1000 )
      parameter ( nvsmax   =            13 * npmax )
      parameter ( nmax     = 2 * ( npmax + ncmax ) )
      parameter ( mmax     =         npmax + ncmax )

C     COMMON SCALARS
      integer nc,np,nx,ny,totnvs
      double precision rmax,xdisp,xscal,xstep,ydisp,yscal,ystep

C     COMMON ARRAYS
      integer nvs(npmax),pol(nvsmax)
      double precision ccent(2*ncmax),edges(3*nvsmax),pcent(2*npmax),
     +        radii(ncmax),vert(2*nvsmax)

C     LOCAL SCALARS
      integer i,inform,nvmapp,nvmipp,pnum
      double precision procit,propol,rmin,secmar,step

C     EXTERNAL SUBROUTINES
      external genpro

C     COMMON BLOCKS
      common /const/ nx,ny,xstep,ystep,rmax
      common /circl/ ccent,radii,nc
      common /poly1/ nvs,np,totnvs
      common /poly2/ vert
      common /poly3/ edges
      common /poly4/ pcent
      common /poly5/ pol
      common /ellip/ xdisp,xscal,ydisp,yscal

C     DEFINE SOME PROBLEM PARAMETERS

C     GRID STEP (THIS IS JUST FOR SCALING THE PROBLEM)

      step   = 1.0d0

C     SECURITY MARGIN (BETWEEN 0 AND 1) TO AVOID CITIES OVERLAPPING

      secmar = 0.1d0

C     DIFFERENT STEPS MAY BE USED IN THE GRID

      xstep  =  step
      ystep  =  step

C     MINIMUM AND MAXIMUM RADII FOR THE CIRCLE- AND POLYGON-CITIES 
C     GENERATION

      rmin   =   0.5d0            * ( step / 2.0d0 )
      rmax   = ( 1.0d0 - secmar ) * ( step / 2.0d0 )

C     READ PROBLEM VARIABLE DATA

C     THIS PROBLEM DATA DEFINES DIFFERENT INSTANCIES OF THE PROBLEM AND 
C     IS RELATED TO THE NUMBER OF VARIABLES AND CONSTRAINTS WHICH 
C     APPROXIMATELY ARE:
C
C     - THE NUMBER OF POLYGON-CITIES IS APPROXIMATELY:
C       NPC = NX * NY * PROCIT * PROPOL
C
C     - THE NUMBER OF CIRCLE-CITIES IS APPROXIMATELY:
C       NCC = MX 8 NY * PROCIT * ( 1 - PROPOL )
C
C     - THE TOTAL NUMBER IF CITIES IS APPROXIMATELY:
C       NC = NPC + NCC
C
C     - NUMBER OF VARIABLES: 
C       2 * ( NC + 1 )
C
C     - NUMBER OF LINEAR CONSTRAINTS:
C       NPC * (NVMAPP - NVMIPP ) / 2 + 4
C
C     - NUMBER OF NONLINEAR CONSTRAINTS:
C       NCC + NC + 1 
C
C     WHERE:
C
C     - PNUM IS A NUMBER THAT IDENTIFIES THE PROBLEM INSTANCE,
C     - NX IS THE NUMBER OF POINTS IN THE GRID ABSCISSA,
C     - NY IS THE NUMBER OF POINTS IN THE GRID ORDINATE,
C     - PROCIT IS THE PROBABILITY OF HAVING A CITY AT A GRID POINT,
C     - PROPOL IS THE PROBABILITY OF A CITY TO BE A POLYGON (THE CITIES
C       WHICH ARE NOT POLYGONS ARE CIRCLES),
C     - NVMIPP IS THE MINIMUM NUNBER OF VERTICES OF A POLYGON,
C     - NVMAPP IS THE MAXIMUM NUMBER OF VERTICES OF A POLYGON.

      read (*,fmt=*) pnum,nx,ny,procit,propol,nvmipp,nvmapp

      if ( procit .lt. 0.0d0 .or. procit .gt. 1.0d0 ) then
          write (*,fmt=*) 'PROCIT MUST BE BETWEEN ZERO AND ONE'
          stop
      end if

      if (propol .lt. 0.0d0 .or. propol .gt. 1.0d0 ) then
          write (*,fmt=*) 'PROPOL MUST BE BETWEEN ZERO AND ONE'
          stop
      end if

      if ( nvmipp .lt. 3 ) then
          write (*,fmt=*) 'NVMIPP MUST BE GREATER THAN OR EQUAL TO 3'
          stop
      end if

      if ( nvmapp .lt. nvmipp ) then
          write (*,fmt=*) 'NVMAPP MUST BE GREATER THAN OR EQUAL TO'
          write (*,fmt=*) 'NVMIPP. INCREASE NVMAPP OR REDUCE NVMIPP.'
          stop
      end if

C     GENERATE THE PROBLEM

      call genpro(nx,ny,xstep,ystep,procit,propol,nvmipp,nvmapp,rmin,
     +rmax,inform)
      if ( inform .ne. 0 ) then
          write (*,fmt=*) 'THERE WAS AN ERROR IN THE PROBLEM GENERATION'
          write (*,fmt=*) '(PROBABLY DUE TO MEMORY SPACE AVAILABILITY).'
          stop
      end if

C     Number of variables

      n = 2 * ( np + nc )

C     Initial point

      do i = 1,np
          x(2*i-1) = pcent(2*i-1)
          x(2*i  ) = pcent(2*i  )
      end do

      do i = 1,nc
          x(2*(np+i)-1) = ccent(2*i-1)
          x(2*(np+i)  ) = ccent(2*i  )
      end do
      
C     Lower and upper bounds

      do i = 1,n
          l(i) = - 1.0d+20
          u(i) =   1.0d+20
      end do

C     Number of constraints (equalities plus inequalities)

      m = np + nc + totnvs + nc

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

C     Ellipse in and out
      do i = 1,np + nc
          linear(i) = .false.
      end do

C     Polygons
      do i = np + nc + 1,np + nc + totnvs
          linear(i) = .true.
      end do

C     Circles
      do i = np + nc + totnvs + 1,np + nc + totnvs + nc
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

C     This subroutine must compute the objective function. For achieving 
C     this objective YOU MUST MODIFY it according to your problem. See 
C     below where your modifications must be inserted.
C     
C     Parameters of the subroutine:
C
C     On Entry:
C
C     n        integer,
C              number of variables,
C
C     x        double precision x(n),
C              current point,
C
C     On Return
C
C     f        double precision,
C              objective function value at x,
C
C     flag     integer,
C              You must set it to any number different of 0 (zero) if 
C              some error ocurred during the evaluation of the objective 
C              function. (For example, trying to compute the square root 
C              of a negative number, dividing by zero or a very small 
C              number, etc.) If everything was o.k. you must set it 
C              equal to zero.

C     LOCAL SCALARS
      integer i,ndist
      double precision diff1,diff2,dist

C     INTRINSIC FUNCTIONS
      intrinsic sqrt

      flag = 0

      f = 0.0d0

      ndist = n / 2 - 1

      do i = 1,ndist
          diff1 = x(1) - x(2 * i + 1)
          diff2 = x(2) - x(2 * i + 2)
          dist = sqrt( diff1 ** 2 + diff2 ** 2)
          if ( dist .le. 1.0d-4 ) then
              write (*,fmt=*)
     +          'ERROR IN PROBLEM DEFINITION (DIST TOO SMALL)'
              flag = - 1
          end if
          f = f + dist
      end do

      f = f / ndist

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine evalg(n,x,g,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,n

C     ARRAY ARGUMENTS
      double precision g(n),x(n)

C     This subroutine must compute the gradient vector of the objective 
C     function. For achieving these objective YOU MUST MODIFY it in the 
C     way specified below. However, if you decide to use numerical 
C     derivatives (we dont encourage this option at all!) you dont need
C     to modify evalg.
C
C     Parameters of the subroutine:
C
C     On Entry:
C
C     n        integer,
C              number of variables,
C
C     x        double precision x(n),
C              current point,
C
C     On Return
C
C     g        double precision g(n),
C              gradient vector of the objective function evaluated at x,
C
C     flag     integer,
C              You must set it to any number different of 0 (zero) if 
C              some error ocurred during the evaluation of any component 
C              of the gradient vector. (For example, trying to compute 
C              the square root of a negative number, dividing by zero or 
C              a very small number, etc.) If everything was o.k. you 
C              must set it equal to zero.

C     LOCAL SCALARS
      integer i,ndist
      double precision diff1,diff2,dist

C     INTRINSIC FUNCTIONS
      intrinsic sqrt

      flag = 0

      g(1) = 0.0d0
      g(2) = 0.0d0

      ndist = n / 2 - 1

      do i = 1,ndist
          diff1 = x(1) - x(2 * i + 1)
          diff2 = x(2) - x(2 * i + 2)
          dist = sqrt( diff1 ** 2 + diff2 ** 2)
          g(2*i+1) = - ( diff1 / dist ) / ndist
          g(2*i+2) = - ( diff2 / dist ) / ndist
          g(1) = g(1) - g(2 * i + 1)
          g(2) = g(2) - g(2 * i + 2)
      end do

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine evalh(n,x,hlin,hcol,hval,hnnz,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,n,hnnz

C     ARRAY ARGUMENTS
      integer hcol(*),hlin(*)
      double precision hval(*),x(n)

C     This subroutine might compute the Hessian matrix of the objective 
C     function. For achieving this objective YOU MAY MODIFY it according 
C     to your problem. To modify this subroutine IS NOT MANDATORY. See 
C     below where your modifications must be inserted.
C     
C     Parameters of the subroutine:
C
C     On Entry:
C
C     n        integer,
C              number of variables,
C
C     x        double precision x(n),
C              current point,
C
C     On Return
C
C     hnnz     integer,
C              number of perhaps-non-null elements of the computed 
C              Hessian,
C
C     hlin     integer hlin(hnnz),
C              see below,
C
C     hcol     integer hcol(hnnz),
C              see below,
C
C     hval     double precision hval(hnnz),
C              the non-null value of the (hlin(k),hcol(k)) position 
C              of the Hessian matrix of the objective function must 
C              be saved at hval(k). Just the lower triangular part of
C              Hessian matrix must be computed,
C
C     flag     integer,
C              You must set it to any number different of 0 (zero) if 
C              some error ocurred during the evaluation of the Hessian
C              matrix of the objective funtion. (For example, trying 
C              to compute the square root of a negative number, 
C              dividing by zero or a very small number, etc.) If 
C              everything was o.k. you must set it equal to zero.

C     LOCAL SCALARS
      integer i,ndist
      double precision diff1,diff2,dist,g11,g22,g21

C     INTRINSIC FUNCTIONS
      intrinsic sqrt

      flag = 0

      ndist = n / 2 - 1

      hnnz = 0

      do i = 1,ndist
          diff1 = x(1) - x(2 * i + 1)
          diff2 = x(2) - x(2 * i + 2)
          dist  = sqrt( diff1 ** 2 + diff2 ** 2 )

          g11 = - ( - 1.0d0 + diff1 * diff1 / dist ** 2 ) / dist / ndist
          g22 = - ( - 1.0d0 + diff2 * diff2 / dist ** 2 ) / dist / ndist
          g21 = - (         + diff1 * diff2 / dist ** 2 ) / dist / ndist

          hlin(hnnz+1)  = 1
          hcol(hnnz+1)  = 1
          hval(hnnz+1)  =   g11 

          hlin(hnnz+2)  = 2 * i + 1
          hcol(hnnz+2)  = 2 * i + 1
          hval(hnnz+2)  =   g11

          hlin(hnnz+3)  = 2 * i + 1
          hcol(hnnz+3)  = 1
          hval(hnnz+3)  = - g11

          hlin(hnnz+4)  = 2
          hcol(hnnz+4)  = 2
          hval(hnnz+4)  =   g22

          hlin(hnnz+5)  = 2 * i + 2
          hcol(hnnz+5)  = 2 * i + 2
          hval(hnnz+5)  =   g22

          hlin(hnnz+6)  = 2 * i + 2
          hcol(hnnz+6)  = 2
          hval(hnnz+6)  = - g22

          hlin(hnnz+7)  = 2
          hcol(hnnz+7)  = 1
          hval(hnnz+7)  =   g21

          hlin(hnnz+8)  = 2 * i + 2
          hcol(hnnz+8)  = 2 * i + 1
          hval(hnnz+8)  =   g21

          hlin(hnnz+9)  = 2 * i + 1
          hcol(hnnz+9)  = 2
          hval(hnnz+9)  = - g21

          hlin(hnnz+10) = 2 * i + 2
          hcol(hnnz+10) = 1
          hval(hnnz+10) = - g21

          hnnz = hnnz + 10
      end do

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

C     This subroutine computes the ind-th constraint.
C
C     Parameters of the subroutine:
C
C     On Entry:
C
C     n        integer,
C              number of variables,
C
C     x        double precision x(n),
C              current point,
C
C     ind      integer,
C              index of the constraint to be computed,
C
C     On Return
C
C     c        double precision,
C              ind-th constraint evaluated at x,
C
C     flag     integer
C              You must set it to any number different of 0 (zero) if 
C              some error ocurred during the evaluation of the 
C              constraint. (For example, trying to compute the square 
C              root of a negative number, dividing by zero or a very 
C              small number, etc.) If everything was o.k. you must set 
C              it equal to zero.
 
C     PARAMETERS
      integer mmax,ncmax,nmax,npmax,nvsmax
      parameter ( ncmax    =                  1000 )
      parameter ( npmax    =                  1000 )
      parameter ( nvsmax   =            13 * npmax )
      parameter ( nmax     = 2 * ( npmax + ncmax ) )
      parameter ( mmax     =         npmax + ncmax )

C     COMMON SCALARS
      integer nc,np,totnvs
      double precision xdisp,xscal,ydisp,yscal

C     COMMON ARRAYS
      integer nvs(npmax),pol(nvsmax)
      double precision ccent(2*ncmax),edges(3*nvsmax),pcent(2*npmax),
     +        radii(ncmax),vert(2*nvsmax)

C     LOCAL SCALARS
      integer i,j

C     LOCAL ARRAYS
      double precision p(2)

C     COMMON BLOCKS
      common /circl/ ccent,radii,nc
      common /poly1/ nvs,np,totnvs
      common /poly2/ vert
      common /poly3/ edges
      common /poly4/ pcent
      common /poly5/ pol
      common /ellip/ xdisp,xscal,ydisp,yscal

      flag = 0

C     Ellipse_in
      if ( ind .eq. 1 ) then

          p(1) = ( x(2 * ind - 1) - xdisp ) / xscal
          p(2) = ( x(2 * ind)     - ydisp ) / yscal
          c = p(1) ** 2 + p(2) ** 2 - 1.0d0

C     Ellipse_out
      else if ( ind .le. np + nc ) then

          p(1) = ( x(2 * ind - 1) - xdisp ) / xscal
          p(2) = ( x(2 * ind)     - ydisp ) / yscal
          c = 1.0d0 - p(1) ** 2 - p(2) ** 2

C     Polygons
      else if ( ind .le. np + nc + totnvs ) then

          j = ind - ( np + nc )
          i = pol(j)
          c = edges(3*j-2) * x(2*i-1) + edges(3*j-1) * x(2*i) 
     +      + edges(3*j)

C     Circles
      else ! if ( ind .le. no + nc + totnvs + nc ) then

          j = ind - ( np + nc + totnvs )
          i = j + np
          c = ( x(2*i-1) - ccent(2*j-1) ) ** 2 +
     +        ( x(2*i  ) - ccent(2*j  ) ) ** 2 - radii(j) ** 2

      end if

      end
 
C     ******************************************************************
C     ******************************************************************
 
      subroutine evaljac(n,x,ind,jcvar,jcval,jcnnz,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,ind,n,jcnnz

C     ARRAY ARGUMENTS
      integer jcvar(n)
      double precision x(n),jcval(n)

C     This subroutine computes the gradient of the ind-th constraint. 
C
C     Parameters of the subroutine:
C
C     On Entry:
C
C     n        integer,
C              number of variables,
C
C     x        double precision x(n),
C              current point,
C
C     ind      integer,
C              index of the constraint whose gradient will be computed,
C
C     On Return
C
C     jcnnz   integer,
C              number of perhaps-non-null elements of the computed 
C              gradient,
C
C     jcvar   integer jcvar(jcnnz),
C              see below,
C
C     jcval   double precision jcval(jcnnz),
C              the non-null value of the partial derivative of the 
C              ind-th constraint with respect to the jcvar(k)-th 
C              variable must be saved at jcval(k).
C
C     flag     integer
C              You must set it to any number different of 0 (zero) if 
C              some error ocurred during the evaluation of the 
C              constraint. (For example, trying to compute the square 
C              root of a negative number, dividing by zero or a very 
C              small number, etc.) If everything was o.k. you must set 
C              it equal to zero.

C     PARAMETERS
      integer mmax,ncmax,nmax,npmax,nvsmax
      parameter ( ncmax    =                  1000 )
      parameter ( npmax    =                  1000 )
      parameter ( nvsmax   =            13 * npmax )
      parameter ( nmax     = 2 * ( npmax + ncmax ) )
      parameter ( mmax     =         npmax + ncmax )

C     COMMON SCALARS
      integer nc,np,totnvs
      double precision xdisp,xscal,ydisp,yscal

C     COMMON ARRAYS
      integer nvs(npmax),pol(nvsmax)
      double precision ccent(2*ncmax),edges(3*nvsmax),pcent(2*npmax),
     +        radii(ncmax),vert(2*nvsmax)

C     LOCAL SCALARS
      integer i,j

C     LOCAL ARRAYS
      double precision p(2)

C     COMMON BLOCKS
      common /circl/ ccent,radii,nc
      common /poly1/ nvs,np,totnvs
      common /poly2/ vert
      common /poly3/ edges
      common /poly4/ pcent
      common /poly5/ pol
      common /ellip/ xdisp,xscal,ydisp,yscal

      flag = 0

C     Ellipse_in
      if ( ind .eq. 1 ) then

          p(1) = ( x(2 * ind - 1) - xdisp ) / xscal
          p(2) = ( x(2 * ind)     - ydisp ) / yscal

          jcnnz = 2

          jcvar(1) = 2 * ind - 1
          jcval(1) = 2.0d0 * p(1) / xscal

          jcvar(2) = 2 * ind
          jcval(2) = 2.0d0 * p(2) / yscal

C     Ellipse_out
      else if ( ind .le. np + nc ) then

          p(1) = ( x(2 * ind - 1) - xdisp ) / xscal
          p(2) = ( x(2 * ind)     - ydisp ) / yscal

          jcnnz = 2

          jcvar(1) = 2 * ind - 1
          jcval(1) = - 2.0d0 * p(1) / xscal

          jcvar(2) = 2 * ind
          jcval(2) = - 2.0d0 * p(2) / yscal

C     Polygons
      else if ( ind .le. np + nc + totnvs ) then

          j = ind - ( np + nc )
          i = pol(j)

          jcnnz = 2

          jcvar(1) = 2 * i - 1
          jcval(1) = edges(3*j-2)

          jcvar(2) = 2 * i
          jcval(2) = edges(3*j-1)

C     Circles
      else ! if ( ind .le. no + nc + totnvs + nc ) then

          j = ind - ( np + nc + totnvs )
          i = j + np

          jcnnz = 2

          jcvar(1) = 2 * i - 1
          jcval(1) = 2.0d0 * ( x(2*i-1) - ccent(2*j-1) )

          jcvar(2) = 2 * i
          jcval(2) = 2.0d0 * ( x(2*i  ) - ccent(2*j  ) )

      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine evalhc(n,x,ind,hclin,hccol,hcval,hcnnz,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,ind,n,hcnnz

C     ARRAY ARGUMENTS
      integer hccol(*),hclin(*)
      double precision hcval(*),x(n)

C     This subroutine might compute the Hessian matrix of the ind-th
C     constraint. For achieving this objective YOU MAY MODIFY it 
C     according to your problem. To modify this subroutine IS NOT 
C     MANDATORY. See below where your modifications must be inserted.
C     
C     Parameters of the subroutine:
C
C     On Entry:
C
C     n        integer,
C              number of variables,
C
C     x        double precision x(n),
C              current point,
C
C     ind      integer,
C              index of the constraint whose Hessian will be computed,
C
C     On Return
C
C     hcnnz    integer,
C              number of perhaps-non-null elements of the computed 
C              Hessian,
C
C     hclin    integer hclin(hcnnz),
C              see below,
C
C     hccol    integer hccol(hcnnz),
C              see below,
C
C     hcval    double precision hcval(hcnnz),
C              the non-null value of the (hclin(k),hccol(k)) position 
C              of the Hessian matrix of the ind-th constraint must 
C              be saved at hcval(k). Just the lower triangular part of
C              Hessian matrix must be computed,
C
C     flag     integer,
C              You must set it to any number different of 0 (zero) if 
C              some error ocurred during the evaluation of the Hessian
C              matrix of the ind-th constraint. (For example, trying 
C              to compute the square root of a negative number, 
C              dividing by zero or a very small number, etc.) If 
C              everything was o.k. you must set it equal to zero.

C     PARAMETERS
      integer mmax,ncmax,nmax,npmax,nvsmax
      parameter ( ncmax    =                  1000 )
      parameter ( npmax    =                  1000 )
      parameter ( nvsmax   =            13 * npmax )
      parameter ( nmax     = 2 * ( npmax + ncmax ) )
      parameter ( mmax     =         npmax + ncmax )

C     COMMON SCALARS
      integer nc,np,totnvs
      double precision xdisp,xscal,ydisp,yscal

C     COMMON ARRAYS
      integer nvs(npmax),pol(nvsmax)
      double precision ccent(2*ncmax),edges(3*nvsmax),pcent(2*npmax),
     +        radii(ncmax),vert(2*nvsmax)

C     LOCAL SCALARS
      integer i,j

C     COMMON BLOCKS
      common /circl/ ccent,radii,nc
      common /poly1/ nvs,np,totnvs
      common /poly2/ vert
      common /poly3/ edges
      common /poly4/ pcent
      common /poly5/ pol
      common /ellip/ xdisp,xscal,ydisp,yscal

      flag = 0

C     Ellipse_in
      if ( ind .eq. 1 ) then

          hcnnz = 2

          hclin(1) = 2 * ind - 1
          hccol(1) = 2 * ind - 1
          hcval(1) = 2.0d0 / xscal ** 2

          hclin(2) = 2 * ind
          hccol(2) = 2 * ind
          hcval(2) = 2.0d0 / yscal ** 2

C     Ellipse_out
      else if ( ind .le. np + nc ) then

          hcnnz = 2

          hclin(1) = 2 * ind - 1
          hccol(1) = 2 * ind - 1
          hcval(1) = - 2.0d0 / xscal ** 2

          hclin(2) = 2 * ind
          hccol(2) = 2 * ind
          hcval(2) = - 2.0d0 / yscal ** 2

C     Polygons
      else if ( ind .le. np + nc + totnvs ) then

          hcnnz = 0

C     Circles
      else ! if ( ind .le. no + nc + totnvs + nc ) then

          j = ind - ( np + nc + totnvs )
          i = j + np

          hcnnz = 2

          hclin(1) = 2 * i - 1
          hccol(1) = 2 * i - 1
          hcval(1) = 2.0d0

          hclin(2) = 2 * i
          hccol(2) = 2 * i
          hcval(2) = 2.0d0

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

C     This subroutine might compute the product of the Hessian of the
C     Lagrangian times vector p (just the Hessian of the objective
C     function in the unconstrained or bound-constrained case).
C
C     Parameters of the subroutine:
C
C     On Entry:
C
C     n        integer,
C              number of variables,
C
C     x        double precision x(n),
C              current point,
C
C     m        integer,
C              number of constraints,
C
C     lambda   double precision lambda(m),
C              vector of Lagrange multipliers,
C
C     p        double precision p(n),
C              vector of the matrix-vector product,
C
C     goth     logical,
C              can be used to indicate if the Hessian matrices were
C              computed at the current point. It is set to .false.
C              by the optimization method every time the current
C              point is modified. Sugestion: if its value is .false.
C              then compute the Hessians, save them in a common
C              structure and set goth to .true.. Otherwise, just use
C              the Hessians saved in the common block structure,
C
C     On Return
C
C     hp       double precision hp(n),
C              Hessian-vector product,
C
C     goth     logical,
C              see above,
C
C     flag     integer,
C              You must set it to any number different of 0 (zero) if
C              some error ocurred during the evaluation of the
C              Hessian-vector product. (For example, trying to compute
C              the square root of a negative number, dividing by zero
C              or a very small number, etc.) If everything was o.k. you
C              must set it equal to zero.

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

C     This subroutine can be used to do some extra job after the solver
C     has found the solution,like some extra statistics, or to save the
C     solution in some special format or to draw some graphical
C     representation of the solution. If the information given by the
C     solver is enough for you then leave the body of this subroutine
C     empty.
C     
C     Parameters of the subroutine:
C
C     The paraemters of this subroutine are the same parameters of
C     subroutine inip. But in this subroutine there are not output
C     parameter. All the parameters are input parameters.

      call drawsol(n,x) 

      end

C     *****************************************************************
C     *****************************************************************

      subroutine genpro(nx,ny,xstep,ystep,procit,propol,nvmipp,nvmapp,
     +rmin,rmax,inform)

C     This subroutine generates a location problem (see the Figure).
C     First, a regular grid with nx horizontal points and ny vertical 
C     points in the positive orthant is considered. The points of the 
C     grid start at the origin with an horizontal distance of xstep and 
C     a vertical distance of ystept. This grid will be the space at 
C     which the cities (represented by polygons) will be distributed. 
C     Before building the cities, an area (rectangle) of preservation 
C     where almost nothing can be done, is defined. This area of 
C     preservation will receive, after the construction of the cities, 
C     an hydraulic plant of energy generation (to supply the energy to 
C     the cities). Then, in the rest of the space, the cities are built. 
C     At each point of the grid (out of the central region) a city 
C     (represented by a polygon) will be built with probability procit. 
C     The definition of the polygon uses variables nvmipp, nvmapp, rmin 
C     and rmax in a way described in the genpol (generate polygon) 
C     subroutine. To transmit the energy from the plant to the cities, a 
C     tower inside each city and a tower inside the central region must 
C     be built. The objective of the problem is to determine the 
C     location of this towers in order to minimize the sum of the 
C     distances from each city tower to the central one. 
C
C     On Entry:
C
C     nx    integer,
C           number of horizontal points in the grid,
C
C     ny    integer,
C           number of vertical points in the grid,
C
C     xstep double precision,
C           horizontal distance between points of the grid,
C
C     ystep double precision,
C           vertical distance between points of the grid,
C
C     procit  double precision, 
C           probability of defining a city at point of the grid
C           (0 <= procit <= 1),
C
C     propol double precision, 
C           probability of the city been defined by a polygon (0 <= 
C           propol <= 1), if the city is not defined by a polygon then 
C           it is defined by a circle,
C
C     nvmipp integer,
C           parameter for the polygon generation (described in genpol 
C           subroutine),
C
C     nvmapp integer,
C           parameter for the polygon generation (described in genpol 
C           subroutine),
C
C     rmin  double precision,
C           parameter for the polygon generation (described in genpol 
C           subroutine),
C
C     rmax  double precision,
C           parameter for the polygon generation (described in genpol 
C           subroutine).
C
C     On output:
C
C     As described in the genpol subroutine, the output is saved in the 
C     polyg common block.


C     PARAMETERS
      integer ncmax
      parameter ( ncmax  =                  1000 )
      integer npmax
      parameter ( npmax  =                  1000 )
      integer nvsmax
      parameter ( nvsmax =            13 * npmax )

C     SCALAR ARGUMENTS
      integer inform,nvmapp,nvmipp,nx,ny
      double precision procit,propol,rmax,rmin,xstep,ystep

C     COMMON SCALARS
      integer nc,np,totnvs
      double precision xdisp,xscal,ydisp,yscal

C     COMMON ARRAYS
      integer nvs(npmax),pol(nvsmax)
      double precision ccent(2*ncmax),edges(3*nvsmax),pcent(2*npmax),
     +        radii(ncmax),vert(2*nvsmax)

C     LOCAL SCALARS
      integer i,j,k
      logical inscre,outcre,insrec
      double precision c,cx,cy,seed,vx,vy,lx,ly,ux,uy

C     EXTERNAL FUNCTIONS
      double precision drand
      external drand

C     EXTERNAL SUBROUTINES
      external genpol

C     COMMON BLOCKS
      common /circl/ ccent,radii,nc
      common /poly1/ nvs,np,totnvs
      common /poly2/ vert
      common /poly3/ edges
      common /poly4/ pcent
      common /poly5/ pol
      common /ellip/ xdisp,xscal,ydisp,yscal

      inform = 0

C     SEED FOR THE RANDOM GENERATION

      seed = 760013.0d0

C     DEFINE CENTRAL REGION

      xdisp = 0.20d0 * ( nx - 1 ) * xstep
      ydisp = 0.50d0 * ( ny - 1 ) * ystep

      xscal = 0.15d0 * ( nx - 1 ) * xstep
      yscal = 0.50d0 * ( ny - 1 ) * ystep

C     DEFINE CENTRAL POINT POLYGON

      nvs(1)    =       4

      lx        =   xdisp + 0.50d0 * xscal 
      ux        =   xdisp + 1.50d0 * xscal
      ly        =   ydisp - 0.75d0 * yscal
      uy        =   ydisp + 0.75d0 * yscal

      vert(1)   =      lx
      vert(2)   =      ly
      vert(3)   =      lx
      vert(4)   =      uy
      vert(5)   =      ux
      vert(6)   =      uy
      vert(7)   =      ux
      vert(8)   =      ly

      edges(1)  = - 1.0d0
      edges(2)  =   0.0d0
      edges(3)  =      lx
      edges(4)  =   0.0d0
      edges(5)  =   1.0d0
      edges(6)  = -    uy
      edges(7)  =   1.0d0
      edges(8)  =   0.0d0
      edges(9)  = -    ux
      edges(10) =   0.0d0
      edges(11) = - 1.0d0
      edges(12) =      ly

      np        =       1
      totnvs    =       4

      pcent(1)  = ( lx + ux ) / 2.0d0
      pcent(2)  = ( ly + uy ) / 2.0d0

      do i = 1,4
          pol(i) = 1
      end do

C     THESE THREE LINES BELOW SHOULD BE USED INSTEAD OF THE 27 LINES
C     ABOVE IF THE CENTRAL RECTANGLE WOULD NOT TO BE CONSIDERED.

C     nvs(1)    =       0
C     np        =       1
C     totnvs    =       0

      nc        =       0

C     DEFINE CITY-POLYGONS AND CITY-CIRCLES CENTERED AT THE GRID 
C     POINTS. THE CITIES MUST BE OUTSIDE THE CENTRAL REGION AND NOT
C     COMPLETELY INSIDE THE CASSINIAN CURVE

      do i = 0,nx - 1
          do j = 0,ny - 1

              cx = i * xstep
              cy = j * ystep

              if ( drand(seed) .le. propol ) then

                ! GENERATE A CITY-POLYGON

                  call genpol(cx,cy,nvmipp,nvmapp,rmin,rmax,seed,inform)

                  if ( inform .ne. 0 ) then
                      return
                  end if

                ! TEST WHETHER THE POLYGON IS IN THE INTERIOR, IT IS 
                ! OUTSIDE OR IT CUTS CENTRAL REGION

                ! 1) CITIES IN THE INTERIOR ARE FORBIDDEN TO AVOID 
                ! INFEASIBLE PROBLEMS

                ! 2) CITIES THAT TOUCH THE CENTRAL REGION ARE DESIRED
                ! WITH THE AIM OF MAKING THE CENTRAL-REGION CONSTRAINT
                ! ACTIVE AT THE SOLUTION AS MANY TIMES AS POSSIBLE

                ! 3) CITIES THAT DO NOT TOUCH THE CENTRAL REGION ARE
                ! ARE INTRODUCED INTO THE PROBLEM WITH PROBABILITY PROCIT

                ! IT IS NOT EASY TO DETERMINE THE PREVIOUS RELATIONS (INSIDE,
                ! OUTSIDE, HALF AND A HALF) BETWEEN THE POLYGON AND THE 
                ! NON-C0NVEX CENTRAL REGION. SO WE WILL CONSIDER THAT:

                ! A POLYGON IS INSIDE THE CENTRAL REGION IF ALL ITS VERTICES
                ! ARE INSIDE OR IN THE BORDER OF THE CENTRAL REGION.

                ! A POLYGON IS OUTSIDE THE CENTRAL REGION IF ALL ITS VERTICES
                ! ARE OUTSIDE OR IN THE BORDER OF THE CENTRAL REGION.

                ! A POLYGON CUTS THE CENTRAL REGION IF IT HAS AT LEAST A 
                ! VERTEX INSIDE OR IN THE BORDER OF THE CENTRAL REGION AND
                ! AT LEAST A VERTEX OUTSIDE OR IN THE BORDER OF THE CENTRAL
                ! REGION. SO, TO BE HALF AND HALF IS EQUIVALENT TO NOT TO
                ! INSIDE, NOR OUTISDE.
                 
                  inscre = .true.
                  outcre = .true.

                  do k = 1,nvs(np+1)

                      vx = ( vert(2*(totnvs+k)-1) - xdisp ) / xscal
                      vy = ( vert(2*(totnvs+k))   - ydisp ) / yscal 

                      c = vx ** 2 + vy ** 2 - 1.0d0

                      if ( c .gt. 0.0d0 ) then
                          inscre = .false.
                      end if

                      if ( c .lt. 0.0d0 ) then
                          outcre = .false.
                      end if

                  end do

                  if ( cx .ge. lx - rmax .and. cx .le. ux + rmax .and.
     +                 cy .ge. ly - rmax .and. cy .le. uy + rmax ) then
                      insrec = .true.
                  else
                      insrec = .false.
                  end if

                ! ADD THE CITY-POLYGON TO THE PROBLEM
 
                  if ( .not. insrec .and.
     +               ( ( .not. inscre .and. .not. outcre ) .or.
     +                 ( outcre .and. drand(seed) .le. procit ) ) ) then

                      np = np + 1
                      totnvs = totnvs + nvs(np)

                  end if

              else

                ! GENERATE A CITY-CIRCLE

                  call gencir(cx,cy,rmin,rmax,seed,inform)

                  if ( inform .ne. 0 ) then
                      return
                  end if

                ! VERIFY WHERE IT IS

                  inscre = .true.
                  outcre = .true.

                  do k = 1,4

                      if ( k .eq. 1 ) then
                          vx = ( cx - radii(nc + 1) - xdisp ) / xscal
                          vy = ( cy                 - ydisp ) / yscal 
                      else if ( k .eq. 2 ) then
                          vx = ( cx + radii(nc + 1) - xdisp ) / xscal
                          vy = ( cy                 - ydisp ) / yscal 
                      else if ( k .eq. 3 ) then
                          vx = ( cx                 - xdisp ) / xscal
                          vy = ( cy - radii(nc + 1) - ydisp ) / yscal 
                      else if ( k .eq. 4 ) then
                          vx = ( cx                 - xdisp ) / xscal
                          vy = ( cy + radii(nc + 1) - ydisp ) / yscal 
                      end if

                      c = vx ** 2 + vy ** 2 - 1.0d0

                      if ( c .gt. 0.0d0 ) then
                          inscre = .false.
                      end if

                      if ( c .lt. 0.0d0 ) then
                          outcre = .false.
                      end if

                 end do

                 if ( cx .ge. lx - rmax .and. cx .le. ux + rmax .and.
     +                cy .ge. ly - rmax .and. cy .le. uy + rmax ) then
                     insrec = .true.
                 else
                     insrec = .false.
                 end if

                ! ADD THE CITY-CIRCLE TO THE PROBLEM

                  if ( .not. insrec .and.
     +               ( ( .not. inscre .and. .not. outcre ) .or.
     +                 ( outcre .and. drand(seed) .le. procit ) ) ) then

                      nc = nc + 1

                  end if

              end if

          end do
      end do

      end

C     *****************************************************************
C     *****************************************************************

      subroutine gencir(cx,cy,rmin,rmax,seed,inform)

C     PARAMETERS
      integer ncmax
      parameter ( ncmax  = 1000 )

C     SCALAR ARGUMENTS
      integer inform
      double precision cx,cy,rmax,rmin,seed

C     COMMON SCALARS
      integer nc

C     COMMON ARRAYS
      double precision ccent(2*ncmax),radii(ncmax)

C     EXTERNAL FUNCTIONS
      double precision drand
      external drand

C     COMMON BLOCKS
      common /circl/ ccent,radii,nc

      inform = 0

C     VERIFY SPACE AVAILABILITY FOR CIRCLES

      if ( nc .eq. ncmax ) then
          write (*,fmt=*) 'THE MAXIMUM NUMBER OF CIRCLES WAS ACHIEVED.'
          write (*,fmt=*) 'INCREASE NCMAX OR REDUCE THE NUMBER OF GRID '
          write (*,fmt=*) 'POINTS (NX*NY) OR THE PROBABILITY OF HAVING '
          write (*,fmt=*) 'A CITY-CIRCLE AT A GRID POINT'
          write (*,fmt=*) '(PROCIT*(1-PROPOL)).'
          inform = - 1
          return
      end if

C     SAVE THE CENTER 

      ccent(2 * nc + 1) = cx
      ccent(2 * nc + 2) = cy

C     GENERATE THE RADIUS

      radii(nc + 1) = rmin + ( rmax - rmin ) * drand(seed)

      end

C     *****************************************************************
C     *****************************************************************

      subroutine genpol(cx,cy,nvmipp,nvmapp,rmin,rmax,seed,inform)

C     This subroutine generates a polygon in R^2 with its 
C     vertices in a sphere centered at point (cx,cy). The
C     number of vertices is randomly generated 
C     satisfying nvmipp <= number of vertices  <= 
C     nvmapp. The ratio of the sphere is also randomly 
C     generated satisfying rmin <= ratio < rmax. The 
C     generated polygon is stored in the common block "polyg".
C 
C     On Entry:
C
C     cx    double precision,
C           first coordinate of the center of the sphere,
C
C     cy    double precision,
C           second coordinate of the center of the sphere,
C
C     nvmipp integer,
C           minimum number of vertices,
C
C     nvmapp integer,
C           maximum number of vertices,
C
C     rmin double precision,
C           minimum ratio of the sphere,
C
C     rmax double precision,
C           maximum ratio of the sphere,
C
C     seed double precision,
C           seed for the random generation.
C
C     On Output:
C
C     inform integer
C            termination parameter
C
C            0 = everything was ok
C
C          < 0 = there was not enough memory space.
C
C     The generated polygon is stored in the polyg common block
C     described below.
C
C     Common block polyg:
C
C     common /polyg/nvs,vert,edges,np,totnvs
C 
C     This structure represents, at any time, np polygons.
C     Position i of array nvs indicates the number of vertices
C     of polygon i. Arrays vert and edges store the 
C     vertices and edges of the polygons.
C
C     For example, if nvs(1) = 3 it indicates that the first
C     polygon has 3 vertices (edges). Then, if the vertices
C     are (x1,y1), (x2,y2) and (x3,y3), we have that vert(1) 
C     = x1, vert(2) = y1, vert(3) = x2, vert(4) = y2, vert(5) 
C     = x3, and vert(6) = y3. And, if the edges (written as
C     ax + by + c = 0) are a1 x + b1 y + c1 = 0, a2 x + b2 y 
C     + c2 = 0, and a3 x + b3 y + c3 = 0 then edges(1) = a1,
C     edges(2) = b1, edges(3) = c1, edges(4) = a2, edges(5) = 
C     b2, edges(6) = c2, edges(7) = a3, edges(8) = b3 and
C     edges(9) = c3.
C
C     totnvs indicates the total number of vertices  
C     of the set of polygons. This information is used when
C     a new polygon is created to know the first free position
C     of arrays vert and edges at which the vertices and edges 
C     of the new polygon will be saved.
C
C     Two additional details: 
C
C     1) For each polygon, the vertices are ordered clockwise 
C     and edge i corresponds to the edge between vertices 
C     i and i+1 (0 if i=n).
C
C     2) For each edge of the form ax + bx + c = 0, constants
C     a, b and c are chosen in such a way that 
C     (|a| = 1 or |b| = 1) and (a cx + b cy + c <= 0).

C     PARAMETERS
      integer npmax
      parameter ( npmax  =                  1000 )
      integer nvsmax
      parameter ( nvsmax =            13 * npmax )
      double precision pi
      parameter ( pi     =     3.1415926535898d0 )

C     SCALAR ARGUMENTS
      integer inform,nvmapp,nvmipp
      double precision cx,cy,rmax,rmin,seed

C     COMMON SCALARS
      integer np,totnvs

C     COMMON ARRAYS
      integer nvs(npmax),pol(nvsmax)
      double precision edges(3*nvsmax),pcent(2*npmax),vert(2*nvsmax)

C     LOCAL SCALARS
      integer i
      double precision dist,lseed,mindist,r

C     LOCAL ARRAYS
      double precision angl(nvsmax)

C     EXTERNAL FUNCTIONS
      double precision drand
      external drand

C     EXTERNAL SUBROUTINES
      external class,constr

C     INTRINSIC FUNCTIONS
      intrinsic cos,int,sin

C     COMMON BLOCKS
      common /poly1/ nvs,np,totnvs
      common /poly2/ vert
      common /poly3/ edges
      common /poly4/ pcent
      common /poly5/ pol

      inform = 0

C     VERIFY SPACE AVAILABILITY FOR POLYGONS

      if ( np .eq. npmax ) then
          write (*,fmt=*) 'THE MAXIMUM NUMBER OF POLYGONS WAS ACHIEVED.'
          write (*,fmt=*) 'INCREASE NPMAX OR REDUCE THE NUMBER OF GRID '
          write (*,fmt=*) 'POINTS (NX*NY) OR THE PROBABILITY OF HAVING '
          write (*,fmt=*) 'A CITY-POLYGON AT A GRID POINT'
          write (*,fmt=*) '(PROCIT*PROPOL).'
          inform = - 1
          return
      end if

C     GENERATE THE NUMBER OF VERTICES 

      lseed = 157318.0d0 + seed

      nvs(np+1) = nvmipp + int((nvmapp-nvmipp+1)*drand(lseed))

C     VERIFY SPACE AVAILABILITY FOR VERTICES AND SIDES

      if ( totnvs + nvs(np+1) .gt. nvsmax ) then
          write (*,fmt=*) 'THE MAXIMUM NUMBER OF POLYGONS VERTICES WAS'
          write (*,fmt=*) 'ACHIEVED. INCREASE NVSMAX OR REDUCE THE'
          write (*,fmt=*) 'AVERAGE NUMBER OF VERTICES (NVMAPP-NVMIPP)/2'
          write (*,fmt=*) 'OR THE NUMBER OF POLYGONS. TO REDUCE THE'
          write (*,fmt=*) 'NUMBER OF POLYGONS, REDUCE THE NUMBER OF'
          write (*,fmt=*) 'GRID POINTS (NX*NY) OR THE PROBABILITY OF'
          write (*,fmt=*) 'HAVING A CITY-POLYGONS AT A GRID POINT'
          write (*,fmt=*) '(PROCIT*PROPOL).'
          inform = - 1
          return
      end if

C     SAVE THE "CENTER" TO BE USED AS INITIAL POINT

      pcent(2*np+1)  = cx
      pcent(2*np+2)  = cy

C     IDENTIFY THE CONSTRAINTS WITH THE POLYGON

      do i = totnvs + 1,totnvs + nvs(np+1)
          pol(i) = np + 1
      end do

C     GENERATE THE RADIUS OF THE SPHERE

      r = rmin + ( rmax - rmin ) * drand( seed )

C     GENERATE ALL ANGLES SATISFYING 0 <= ANGLE_I < 2*PI

 10   continue

      do i = 1,nvs(np + 1)
          angl(i) = 2 * pi * drand( lseed )
      end do

C     CLASSIFY THE ANGLES IN DECREASING ORDER

      call class(nvs(np + 1),angl)

C     CONSTRUCT THE VERTICES

      do i = 1,nvs(np + 1)
          vert(2 * ( totnvs + i ) - 1) = cx + r * cos( angl(i) )
          vert(2 * ( totnvs + i )    ) = cy + r * sin( angl(i) )
      end do

C     FOR NUMERICAL STABILITY IN THE CONSTRAINT CALCULATION,
C     AVOID TOO SIMILAR ANGLES (EQUIVALENT TO TOO NEAR POINTS)

      mindist = 1.0D+99

      do i = totnvs + 1,totnvs + nvs(np+1) - 1

          dist = ( vert(2*i-1) - vert(2*i+1) ) ** 2 +
     +           ( vert(2*i)   - vert(2*i+2) ) ** 2

          if ( dist .lt. mindist ) then
              mindist = dist
          end if

      end do

      i = totnvs + nvs(np+1)

      dist = ( vert(2*i-1) - vert(2*totnvs+1) ) ** 2 +
     +       ( vert(2*i)   - vert(2*totnvs+2) ) ** 2

      if ( dist .lt. mindist ) then
          mindist = dist
      end if

      if ( mindist .lt. 1.0D-8 ) then
C         write (*,fmt=*) 'DISCARDING ANGLES.'
          go to 10
      end if

C     CONSTRUCT THE EDGES

      do i = totnvs + 1,totnvs + nvs(np+1) - 2
          call constr(vert(2*i-1),vert(2*i),vert(2*i+1),vert(2*i+2),
     +                vert(2*i+3),vert(2*i+4),edges(3*i-2),edges(3*i-1),
     +                edges(3*i),inform)
          if ( inform .ne. 0 ) then
              return
          end if
      end do

      i = totnvs + nvs(np+1) - 1

      call constr(vert(2*i-1),vert(2*i),vert(2*i+1),vert(2*i+2),
     +            vert(2*totnvs+1),vert(2*totnvs+2),edges(3*i-2),
     +            edges(3*i-1),edges(3*i),inform)
      if ( inform .ne. 0 ) then
          return
      end if

      i = totnvs + nvs(np+1)

      call constr(vert(2*i-1),vert(2*i),vert(2*totnvs+1),
     +            vert(2*totnvs+2),vert(2*totnvs+3),vert(2*totnvs+4),
     +            edges(3*i-2),edges(3*i-1),edges(3*i),inform)
      if ( inform .ne. 0 ) then
          return
      end if

      end

C     *****************************************************************
C     *****************************************************************

      subroutine constr(x1,y1,x2,y2,x3,y3,a,b,c,inform)

C     This subroutine computes the real constants a, b and c of 
C     the straight line ax + by + c = 0 in R^2 defined by the 
C     points (x1,y1) and (x2,y2); such that the point (x3,y3) 
C     satisfies the constraint ax + by + c <= 0.
C
C     On Entry:
C
C     x1    double precision,
C           first coordinate of point (x1,y1),
C
C     y1    double precision,
C           second coordinate of point (x1,y1),
C
C     x2    double precision,
C           first coordinate of point (x2,y2),
C
C     y2    double precision,
C           second coordinate of point (x2,y2),
C
C     x3    double precision,
C           first coordinate of point (x3,y3),
C
C     y3    double precision,
C           second coordinate of point (x3,y3).
C
C     On Return
C
C     a,b,c double precision
C           the desired constants.

C     SCALAR ARGUMENTS
      integer inform
      double precision a,b,c,x1,x2,x3,y1,y2,y3

      if (x1.eq.x2 .and. y1.eq.y2) then
          write (*,fmt=*)
     +      'ERROR IN FUNCTION CONSTRAINT: X1=X2 AND Y1=Y2'
          inform = - 1
          return
      end if

      if (y1.ne.y2) then
          a = 1.0d0
          b = - (x2-x1)/ (y2-y1)
          c = - (x1+b*y1)

      else
          a = 0.0d0
          b = 1.0d0
          c = -y1
      end if

      if (a*x3+b*y3+c.gt.0.0d0) then
          a = -a
          b = -b
          c = -c
      end if

      end

C     *****************************************************************
C     *****************************************************************

      subroutine class(n,x)

C     This subroutine classifies the elements of a vector in
C     decreasing order, i.e., on output: x(1) >= x(2) >= 
C     ... >= x(n).
C
C     On Entry:
C
C     n     integer,
C           number of elements of the vector to be classified,
C
C     x     double precision x(n),
C           vector to be classified.
C
C     On Return
C
C     x     double precision x(n),
C           classified vector.

C     SCALAR ARGUMENTS
      integer n

C     ARRAY ARGUMENTS
      double precision x(n)

C     LOCAL SCALARS
      double precision aux,xmax
      integer i,j,pos

      do i = 1,n

          xmax = x(i)
          pos = i

          do j = i + 1,n
              if (x(j).gt.xmax) then
                  xmax = x(j)
                  pos = j
              end if
          end do

          if (pos.ne.i) then
              aux = x(i)
              x(i) = x(pos)
              x(pos) = aux
          end if

      end do

      end

C     *****************************************************************
C     *****************************************************************

      subroutine drawsol(n,x)

C     SCALAR ARGUMENTS
      integer n

C     ARRAY ARGUMENTS
      double precision x(n)

C     This subroutine generate a metapost file with the
C     graphical representation of the problem.

C     On Entry:
C
C     n     integer,
C           number of elements of array x,
C
C     x     double precision x(n)
C           current approximation to the solution.

C     PARAMETERS
      integer mmax,ncmax,nmax,npmax,nvsmax
      parameter ( ncmax    =                  1000 )
      parameter ( npmax    =                  1000 )
      parameter ( nvsmax   =            13 * npmax )
      parameter ( nmax     = 2 * ( npmax + ncmax ) )
      parameter ( mmax     =         npmax + ncmax )

C     COMMON SCALARS
      integer nc,np,nx,ny,totnvs
      double precision rmax,xdisp,xscal,xstep,ydisp,yscal,ystep

C     COMMON ARRAYS
      integer nvs(npmax),pol(nvsmax)
      double precision ccent(2*ncmax),edges(3*nvsmax),pcent(2*npmax),
     +        radii(ncmax),vert(2*nvsmax)

C     LOCAL SCALARS
      integer base,i,j
      double precision scale,vx,vy

C     COMMON BLOCKS
      common /const/ nx,ny,xstep,ystep,rmax
      common /circl/ ccent,radii,nc
      common /poly1/ nvs,np,totnvs
      common /poly2/ vert
      common /poly3/ edges
      common /poly4/ pcent
      common /poly5/ pol
      common /ellip/ xdisp,xscal,ydisp,yscal

      scale = min( 21.6d0 / ( ( nx - 1 ) * xstep + 2.0d0 * rmax ), 
     +             27.9d0 / ( ( ny - 1 ) * ystep + 2.0d0 * rmax ) )

C     DRAW THE PROBLEM

      open(unit=10,file='solution.mp')

      write(10,10) scale

C     ELLIPSE

      write(10,50) 2.0d0 * xscal, 2.0d0 * yscal,xdisp,ydisp

C     CIRCLES

      do i = 1,nc
          write(10,50) 2.0d0 * radii(i),2.0d0 * radii(i),
     +                 ccent(2 * i - 1),ccent(2 * i)
      end do

C     POLYGONS

      base = 0

      do i = 1,np

          if ( nvs(i) .ne. 0 ) then

              write(10,20)

              do j = base + 1,base + nvs(i)
                  vx = vert(2 * j - 1)
                  vy = vert(2 * j)
                  write(10,60) vx,vy
              end do

              write(10,40)

              base = base + nvs(i)

          end if

      end do

C     LINES

      do i = 1,np + nc
          write(10,70) x(1),x(2),x(2 * i - 1),x(2 * i)
      end do

C     DOTS WITH LABELS

C     write(10,80)

C     do i = 1,np + nc
C         write(10,90) i,x(2 * i - 1),x(2 * i)
C     end do

C     DOTS WITHOUT LABELS

      write(10,80)

      do i = 1,np + nc
          write(10,95) x(2 * i - 1),x(2 * i)
      end do

C     END

      write(10,100)

      close(10)

C     NON-EXECUTABLE STATEMENTS

  10  format('beginfig(1);'/,'u = ',f20.10,' cm;') 
  20  format('draw')
  30  format(5X,'(',f20.10,'u,',f20.10,'u)..')
  40  format(5X,'cycle;')
  50  format('draw fullcircle',
     +       /,5X,'xscaled   ',f20.10,'u',
     +       /,5X,'yscaled   ',f20.10,'u',
     +       /,5X,'shifted (',f20.10,'u,',f20.10,'u);')
  60  format(5X,'(',f20.10,'u,',f20.10,'u)--')
  70  format('draw (',f20.10,'u,',f20.10,'u)--',
     +       /,5X,'(',f20.10,'u,',f20.10,'u) dashed evenly;')
  80  format('pickup pencircle scaled 4pt;')
  90  format('dotlabel.top (btex $',i4,'$ etex, (',f20.10,'u,',
     +        f20.10,'u) );')
  95  format('draw (',f20.10,'u,',f20.10,'u);')
 100  format('endfig;'/,'end;') 

      end

C     ******************************************************************
C     ******************************************************************

c     double precision function drand(ix)
c     double precision ix

c     double precision a,p,b15,b16,xhi,xalo,leftlo,fhi,k
c     data a/16807.d0/,b15/32768.d0/,b16/65536.d0/,p/2147483647.d0/

c     xhi= ix/b16
c     xhi= xhi - dmod(xhi,1.d0)
c     xalo= (ix-xhi*b16)*a
c     leftlo= xalo/b16
c     leftlo= leftlo - dmod(leftlo,1.d0)
c     fhi= xhi*a + leftlo
c     k= fhi/b15
c     k= k - dmod(k,1.d0)
c     ix= (((xalo-leftlo*b16)-p)+(fhi-k*b15)*b16)+k
c     if (ix.lt.0) ix= ix + p
c     drand= ix*4.656612875d-10

c     return
c     end
 
