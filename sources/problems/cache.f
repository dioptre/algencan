C     =================================================================
C     File: cache.f
C     =================================================================

C     =================================================================
C     Module: Subroutines that define the problem
C     =================================================================

C     Last update of any of the component of this module:

C     January 30, 2007.

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

C     Cache problem
C     -------------

C     Mathematical programming model for fixing content location cache
C     expiration dates at aggregation nodes in a content network, in
C     order to maximize the information discovered by the nodes of the
C     networks, while taking into account bandwidth constraints at the
C     aggregation nodes.

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
      integer nmax
      parameter ( nmax      =  500000 )

C     COMMON SCALARS
      double precision alphaS,alphaB,betaS,betaB,BWin,BWout,t

C     COMMON ARRAYS
      double precision ll(nmax),ff(nmax),mu(nmax),llambda(nmax),w(nmax)

C     COMMON BLOCKS
      common /probdata/ ll,ff,mu,llambda,w,alphaS,alphaB,betaS,betaB,
     +                  BWin,BWout,t

C     LOCAL SCALARS
      integer i,k
      double precision tmp

C     Set problem data

      open(10,file='cache.dat')

      read(10,*) k

      read(10,*) alphaS,alphaB,betaS,betaB,BWin,BWout

      do i = 1,k
          read(10,*) ff(i)
      end do

      do i = 1,k
          read(10,*) llambda(i)
      end do

      do i = 1,k
          read(10,*) mu(i)
      end do

      do i = 1,k
          read(10,*) ll(i)
      end do

      close(10)

      t = 0.0d0
      do i = 1,k
          tmp  = ll(i) * llambda(i) / mu(i)
          t    = t + ff(i) * tmp
          w(i) = tmp / mu(i)
      end do

C     Number of variables

      n = 3 * k

C     Initial point

      do i = 1,k
          x(2*k+i) = 1.0d-01
          x(    i) = ( 1.0d0 - exp( - ff(i) * x(2*k+i) ) ) / x(2*k+i)
          x(  k+i) = ( 1.0d0 - exp( - mu(i) * x(2*k+i) ) ) / x(2*k+i)
      end do

C     Lower and upper bounds

      do i = 1,k
          l(    i) =   0.0d0
          u(    i) =   ff(i)
          l(  k+i) =   0.0d0
          u(  k+i) =   mu(i)
          l(2*k+i) =   0.0d0
          u(2*k+i) = 1.0d+20
      end do

C     Number of constraints (equalities plus inequalities)

      m = 2 * k + 2

C     Lagrange multipliers approximation.

      do i = 1,m
          lambda(i) =  0.0d0
      end do

C     For each constraint i, set equatn(i) = .true. if it is an equality
C     constraint of the form c_i(x) = 0, and set equatn(i) = .false. if
C     it is an inequality constraint of the form c_i(x) <= 0.

      do i = 1,2 * k
          equatn(i) = .true.
      end do
      equatn(2 * k + 1) = .false.
      equatn(2 * k + 2) = .false.

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

C     PARAMETERS
      integer nmax
      parameter ( nmax      =  500000 )

C     COMMON SCALARS
      double precision alphaS,alphaB,betaS,betaB,BWin,BWout,t

C     COMMON ARRAYS
      double precision ll(nmax),ff(nmax),mu(nmax),llambda(nmax),w(nmax)

C     LOCAL SCALARS
      integer i,k

C     COMMON BLOCKS
      common /probdata/ ll,ff,mu,llambda,w,alphaS,alphaB,betaS,betaB,
     +                  BWin,BWout,t

      flag = 0

      k = n / 3

      f = 0.0d0
      do i = 1,k
          f = f - w(i) * (mu(i) * x(i) + ff(i) * x(k+i) - x(i) * x(k+i))
      end do

      f = f / t

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

C     PARAMETERS
      integer nmax
      parameter ( nmax      =  500000 )

C     COMMON SCALARS
      double precision alphaS,alphaB,betaS,betaB,BWin,BWout,t

C     COMMON ARRAYS
      double precision ll(nmax),ff(nmax),mu(nmax),llambda(nmax),w(nmax)

C     LOCAL SCALARS
      integer i,k

C     COMMON BLOCKS
      common /probdata/ ll,ff,mu,llambda,w,alphaS,alphaB,betaS,betaB,
     +                  BWin,BWout,t

      flag = 0

      k = n / 3

      do i = 1,k
          g(  i)   = - w(i) * ( mu(i) - x(k+i) ) / t
          g(k+i)   = - w(i) * ( ff(i) -   x(i) ) / t
          g(2*k+i) = 0.0d0
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

C     This subroutine must compute the ind-th constraint of your
C     problem. For achieving this objective YOU MUST MOFIFY it
C     according to your problem. See below the places where your
C     modifications must be inserted.
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
      integer nmax
      parameter ( nmax      =  500000 )

C     COMMON SCALARS
      double precision alphaS,alphaB,betaS,betaB,BWin,BWout,t

C     COMMON ARRAYS
      double precision ll(nmax),ff(nmax),mu(nmax),llambda(nmax),w(nmax)

C     LOCAL SCALARS
      integer i,k

C     COMMON BLOCKS
      common /probdata/ ll,ff,mu,llambda,w,alphaS,alphaB,betaS,betaB,
     +                  BWin,BWout,t

      flag = 0

      k = n / 3

      if ( ind .le. k ) then

          i = ind
          c = ( 1.0d0 - exp( - ff(i) * x(2*k+i) ) ) - x(i) * x(2*k+i)

      else if ( ind .le. 2 * k ) then

          i = ind - k
          c = ( 1.0d0 - exp( - mu(i) * x(2*k+i) ) ) - x(k+i) * x(2*k+i)

      else if ( ind .eq. 2 * k + 1 ) then

          c = 0.0d0
          do i = 1,k
              c = c + betaS * ll(i) * ff(i) + alphaB * ll(i) *
     +        llambda(i) / mu(i) * ff(i) / ( 1.0d0 + x(2*k+i) * ff(i) )
          end do
          c = c - BWin

      else if ( ind .eq. 2 * k + 2 ) then

    	  c = 0.0d0
	  do i = 1,k
	     c = c + alphaS * ll(i) * llambda(i) / mu(i) * ff(i) +
     +       betaB * ll(i) * ff(i) / ( 1.0d0 + x(2*k+i) * ff(i) )
          end do
    	  c = c - BWout

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

C     This subroutine must compute the gradient of the ind-th constraint.
C     For achieving these objective YOU MUST MODIFY it in the way
C     specified below.
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
C     jcnnz    integer,
C              number of perhaps-non-null elements of the computed
C              gradient,
C
C     jcvar    integer jcvar(jcnnz),
C              see below,
C
C     jcval    double precision jcval(jcnnz),
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
      integer nmax
      parameter ( nmax      =  500000 )

C     COMMON SCALARS
      double precision alphaS,alphaB,betaS,betaB,BWin,BWout,t

C     COMMON ARRAYS
      double precision ll(nmax),ff(nmax),mu(nmax),llambda(nmax),w(nmax)

C     LOCAL SCALARS
      integer i,k

C     COMMON BLOCKS
      common /probdata/ ll,ff,mu,llambda,w,alphaS,alphaB,betaS,betaB,
     +                  BWin,BWout,t

      flag = 0

      k = n / 3

      if ( ind .le. k ) then

          i = ind

          jcnnz = 2

          jcvar(1) = 2 * k + i
          jcval(1) = ff(i) *  exp( - ff(i) * x(2*k+i) ) - x(i)

          jcvar(2) = i
          jcval(2) = - x(2*k+i)

      else if ( ind .le. 2 * k ) then

          i = ind - k

          jcnnz = 2

          jcvar(1) = 2 * k + i
          jcval(1) = mu(i) * exp( - mu(i) * x(2*k+i) ) - x(k+i)

          jcvar(2) = k + i
          jcval(2) = - x(2*k+i)

      else if ( ind .eq. 2 * k + 1 ) then

         jcnnz = k
         do i = 1,k
             jcvar(i) = 2 * k + i
             jcval(i) = - ff(i) ** 2 * alphaB * ll(i) * llambda(i) /
     +                     mu(i) / ( 1.0d0 + x(2*k+i) * ff(i) ) ** 2
         end do

      else if ( ind .eq. 2 * k + 2 ) then

         jcnnz = k
         do i = 1,k
             jcvar(i) = 2 * k + i
             jcval(i) = - betaB * ll(i) * ff(i) ** 2 /
     +                     ( 1.0d0 + x(2*k+i) * ff(i) ) ** 2
         end do

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

      end
