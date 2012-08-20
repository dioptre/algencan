C     =================================================================
C     File: hardspheres.f
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

C     Hard spheres problem
C     --------------------

C     This problem is associated to the family of Hard-Spheres
C     problem. It belongs to a family of sphere packing
C     problems, a class of challenging problems dating from the
C     beginning of the seventeenth century. In the tradition of famous
C     problems in mathematics, the statements of these problems are
C     elusively simple, and have withstood the attacks of many worthy
C     mathematicians (e.g. Newton, Hilbert, Gregory), while most of its
C     instances remain open problems. Furthermore, it is related to
C     practical problems in chemistry, biology and physics.
C     The Hard-Spheres Problem is to maximize the minimum pairwise
C     distance between np points on a sphere in R^ndim. This problem
C     may be reduced to a nonlinear optimization problem that turns
C     out, as might be expected from the mentioned history, to be a
C     particularly hard, nonconvex problem, with a potentially large
C     number of (nonoptimal) points satisfying KKT conditions. We have
C     thus a class of problems indexed by the parameters ndim and np,
C     that provides a suitable set of test problems for evaluating
C     Nonlinear Programming codes.

C     After some algebric manipulations, we can formulate this problem
C     as
C
C               Minimize z
C
C               subject to
C
C               z \geq <x_i,x_j> for all different pair of indices i,j
C
C               ||x_i||^2 = 1    for all i = 1,...,NP
C
C     The goal is to find an objective value less than 0.5 (This means
C     that the NP points stored belong to the sphere and every distance
C     between two of them is greater than 1.0).
C
C     Obs: the starting point is aleatorally chosen although each
C     variable belongs to [-1,1].
C
C     References:
C
C     [1] "Validation of an Augmented Lagrangian algorithm with a
C         Gauss-Newton Hessian approximation using a set of
C         Hard-Spheres problems", N. Krejic, J. M. Martinez, M. Mello
C         and E. A. Pilotta, Tech. Report RP 29/98, IMECC-UNICAMP,
C         Campinas, 1998.
C
C     [2] "Inexact-Restoration Algorithm for Constrained Optimization",
C         J. M. Martinez and E. A. Pilotta, Tech. Report,
C         IMECC-UNICAMP, Campinas, 1998.
C
C     [3] "Sphere Packings, Lattices and Groups", J. H. Conway and
C         N. J. C. Sloane, Springer-Verlag, NY, 1988.

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
      integer npmax
      parameter ( npmax = 1000 )

C     COMMON SCALARS
      integer ndim,np

C     COMMON ARRAYS
      integer pair(2,npmax * ( npmax - 1 ) / 2)

C     LOCAL SCALARS
      integer i,j,k
      double precision seed

C     COMMON BLOCKS
      common /probdata/ pair,ndim,np

C     EXTERNAL FUNCTIONS
      double precision drand

C     Set problem data

      write(*,*) 'Dimension of the space (ndim) and number of ',
     +           'points (np): '
      read(*,*) ndim,np

      write(*,*) 'Seed for the initial point random generation: '
      read(*,*) seed

C     Alternatively:
C     read(*,*) ndim,np,seed
C     seed = 123456.0d0 * seed

C     Set pairs

      k = 0
      do i = 1,np
          do j = i + 1,np
              k = k + 1
              pair(1,k) = i
              pair(2,k) = j
          end do
      end do

C     Number of variables

      n = ndim * np + 1

C     Initial point

      do i = 1,n - 1
          x(i) = 2.0d0 * drand(seed) - 1.0d0
      end do

      x(n) = drand(seed)

C     Lower and upper bounds
      do i = 1,n
          l(i) = - 1.0d+20
          u(i) =   1.0d+20
      end do

C     Number of constraints (equalities plus inequalities)

      m = np + np * ( np - 1 ) / 2

C     Lagrange multipliers approximation.

      do i = 1,m
          lambda(i) =  0.0d0
      end do

C     For each constraint i, set equatn(i) = .true. if it is an equality
C     constraint of the form c_i(x) = 0, and set equatn(i) = .false. if
C     it is an inequality constraint of the form c_i(x) <= 0.

      do i = 1,np
          equatn(i) = .true.
      end do

      do i = np + 1, np + np * ( np - 1 ) / 2
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
      integer npmax
      parameter ( npmax = 1000 )

C     COMMON SCALARS
      integer ndim,np

C     COMMON ARRAYS
      integer pair(2,npmax * ( npmax - 1 ) / 2)

C     LOCAL SCALARS
      integer i,i1,i2

C     COMMON BLOCKS
      common /probdata/ pair,ndim,np

      flag = 0

      if ( 0 .lt. ind .and. ind .le. np ) then

          c = - 1.0d0

          do i = 1,ndim
              c = c + x(ndim * ( ind - 1 ) + i) ** 2
          end do

          return

      end if

      if ( np .lt. ind .and.
     +     ind .le. np + np * (np - 1) / 2 ) then

          i1 = pair(1,ind - np)
          i2 = pair(2,ind - np)

          c = - x(n)
          do i = 1,ndim
              c = c + x(ndim * (i1 - 1) + i) * x(ndim * (i2 - 1) + i)
          end do

          return

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

C     This subroutine must compute the gradient of the constraint ind.
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
      integer npmax
      parameter ( npmax = 1000 )

C     COMMON SCALARS
      integer ndim,np

C     COMMON ARRAYS
      integer pair(2,npmax * ( npmax - 1 ) / 2)

C     LOCAL SCALARS
      integer i,i1,i2

C     COMMON BLOCKS
      common /probdata/ pair,ndim,np

      flag = 0

      if ( 0 .lt. ind .and. ind .le. np ) then

          jcnnz = ndim

          do i = 1,ndim
              jcvar(i) = ndim * ( ind - 1 ) + i
              jcval(i) = 2.0d0 * x(ndim * ( ind - 1 ) + i)
          end do

          return

      end if

      if ( np .lt. ind .and.
     +     ind .le. np + np * (np - 1) / 2 ) then

          i1 = pair(1,ind - np)
          i2 = pair(2,ind - np)

          jcnnz = 2 * ndim + 1

          do i = 1,ndim
              jcvar(i) = ndim * ( i1 - 1 ) + i
              jcval(i) = x(ndim * ( i2 - 1 ) + i)
              jcvar(ndim+i) = ndim * ( i2 - 1 ) + i
              jcval(ndim+i) = x(ndim * ( i1 - 1 ) + i)
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
      integer npmax
      parameter ( npmax = 1000 )

C     COMMON SCALARS
      integer ndim,np

C     COMMON ARRAYS
      integer pair(2,npmax * ( npmax - 1 ) / 2)

C     LOCAL SCALARS
      integer i,i1,i2,tmp

C     COMMON BLOCKS
      common /probdata/ pair,ndim,np

      flag = 0

      if ( 1 .le. ind .and. ind .le. np ) then

          hcnnz = ndim

          do i = 1,ndim
              hclin(i) = ndim * ( ind - 1 ) + i
              hccol(i) = ndim * ( ind - 1 ) + i
              hcval(i) = 2.0d0
          end do

          return

      end if

      if ( ind .ge. np + 1 .and.
     +     ind .le. np + np * (np - 1) / 2 ) then

          i1 = pair(1,ind - np)
          i2 = pair(2,ind - np)

          if ( i1 .lt. i2 ) then
              tmp = i1
              i1  = i2
              i2  = tmp
          end if

          hcnnz = ndim

          do i = 1,ndim
              hclin(i)      = ndim * ( i1 - 1 ) + i
              hccol(i)      = ndim * ( i2 - 1 ) + i
              hcval(i)      = 1.0d0
          end do

          return

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

      end
