C     =================================================================
C     File: toyprob2.f
C     =================================================================

C     =================================================================
C     Module: Subroutines that define the problem
C     =================================================================

C     Last update of any of the component of this module:

C     January 20, 2010.

C     Users are encouraged to download periodically updated versions of
C     this code at the TANGO home page:

C     www.ime.usp.br/~egbirgin/tango/

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

C     LOCAL SCALARS
      integer i

C     ******************************************************************
C     FROM HERE ON YOU MUST MODIFY THE SUBROUTINE TO SET YOUR PROBLEM
C     DATA:
C     ******************************************************************

C     Number of variables

      n = 2

C     Initial point

      do i = 1,n - 1
          x(i) = 0.0d0
      end do

      x(n) = 0.0d0

C     Lower and upper bounds

      do i = 1,n - 1
          l(i) = - 10.0d0
          u(i) =   10.0d0
      end do

      l(n) = - 1.0d+20
      u(n) =   1.0d+20

C     Number of constraints (equalities plus inequalities)

      m = 2

C     Lagrange multipliers approximation. Most users prefer to use the
C     null initial Lagrange multipliers estimates. However, if the
C     problem that you are solving is "slightly different" from a
C     previously solved problem of which you know the correct Lagrange
C     multipliers, we encourage you to set these multipliers as initial
C     estimates. Of course, in this case you are also encouraged to use
C     the solution of the previous problem as initial estimate of the
C     solution. Similarly, most users prefer to use rho = 10 as initial
C     penalty parameters. But in the case mentioned above (good
C     estimates of solution and Lagrange multipliers) larger values of
C     the penalty parameters (say, rho = 1000) may be more useful. More
C     warm-start procedures are being elaborated.

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

      linear(1) = .false.
      linear(2) =  .true.

C     Indicate which subroutines did you code.

      coded( 1) = .false. ! evalf
      coded( 2) = .false. ! evalg
      coded( 3) = .false. ! evalh
      coded( 4) = .false. ! evalc
      coded( 5) = .false. ! evaljac
      coded( 6) = .false. ! evalhc
      coded( 7) = .true.  ! evalfc
      coded( 8) = .true.  ! evalgjac
      coded( 9) = .true.  ! evalhl
      coded(10) = .false. ! evalhlp

C     Set checkder = TRUE if you code some derivatives and you would
C     like them to be tested by finite differences. It is highly
C     recommended.

      checkder = .true.

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

C     Objective function

C     ******************************************************************
C     FROM HERE ON YOU MUST MODIFY THE SUBROUTINE TO SET YOUR OBJECTIVE
C     FUNCTION:
C     ******************************************************************

      flag = - 1

C     ******************************************************************
C     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALF.
C     ******************************************************************

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

C     Gradient vector

C     ******************************************************************
C     FROM HERE ON YOU MUST MODIFY THE SUBROUTINE TO SET THE GRADIENT
C     VECTOR OF YOUR OBJECTIVE FUNCTION:
C     ******************************************************************

      flag = - 1

C     ******************************************************************
C     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALG.
C     ******************************************************************

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

C     ******************************************************************
C     FROM HERE ON YOU MAY (OPTIONALY) MODIFY THE SUBROUTINE TO SET THE
C     HESSIAN MATRIX OF YOUR OBJECTIVE FUNCTION:
C     ******************************************************************

      flag = - 1

C     ******************************************************************
C     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALH.
C     ******************************************************************

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

C     Constraints

C     ******************************************************************
C     FROM HERE ON YOU MUST MODIFY THE SUBROUTINE TO SET YOUR CONSTRAINTS
C     ******************************************************************

      flag = - 1

C     ******************************************************************
C     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALC.
C     ******************************************************************

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

C     Sparse gradient vector of the ind-th constraint

C     ******************************************************************
C     FROM HERE ON YOU MUST MODIFY THE SUBROUTINE TO SET THE GRADIENTS
C     OF YOUR CONSTRAINTS:
C     ******************************************************************

      flag = - 1

C     ******************************************************************
C     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALJAC.
C     ******************************************************************

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

C     ******************************************************************
C     FROM HERE ON YOU MAY (OPTIONALY) MODIFY THE SUBROUTINE TO SET THE
C     HESSIANS OF YOUR CONSTRAINTS:
C     ******************************************************************

      flag = - 1

C     ******************************************************************
C     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALHC.
C     ******************************************************************

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

      flag = 0

C     Objective function

      f = x(2)

C     Constraints

      c(1) = x(1) ** 2 + 1 - x(n)
      c(2) = 2.0d0 - x(1) - x(n)

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

      flag = 0

C     Gradient of the objective function

      g(1) = 0.0d0
      g(2) = 1.0d0

C     Gradient of constraints

      jcnnz = 0

C     Constraint 1

      jcfun(jcnnz + 1) = 1
      jcvar(jcnnz + 1) = 1
      jcval(jcnnz + 1) = 2.0d0 * x(1)

      jcfun(jcnnz + 2) = 1
      jcvar(jcnnz + 2) = 2
      jcval(jcnnz + 2) = - 1.0d0

      jcnnz = jcnnz + 2

C     Constraint 2

      jcfun(jcnnz + 1) = 2
      jcvar(jcnnz + 1) = 1
      jcval(jcnnz + 1) = - 1.0d0

      jcfun(jcnnz + 2) = 2
      jcvar(jcnnz + 2) = 2
      jcval(jcnnz + 2) = - 1.0d0

      jcnnz = jcnnz + 2

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

      flag = 0

      hlnnz = 1

      hllin(1) = 1
      hlcol(1) = 1
      hlval(1) = lambda(1) * scalec(1) * 2.0d0

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalhlp(n,x,m,lambda,sf,sc,p,hp,goth,flag)

      implicit none

C     SCALAR ARGUMENTS
      logical goth
      integer flag,m,n
      double precision sf

C     ARRAY ARGUMENTS
      double precision hp(n),lambda(m),p(n),sc(m),x(n)

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

