C     =================================================================
C     File: balls.f
C     =================================================================

C     =================================================================
C     Module: Subroutines that define the problem
C     =================================================================

C     Last update of any of the component of this module:

C     March 13, 2007.

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

C     Balls problem
C     -------------
C
C     This problem consists of finding np points in R^ndim such that
C     the distance between any pair of them is not less than 1 and the
C     maximum distance is as small as possible:
C
C     Minimize z
C
C     subject to 1 <= ||x_i - x_j||_2^2 <= leq z,
C
C     for all i not equal to j, i,j =1, ..., np.

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
      integer i,ind1,ind2,j,k
      double precision dist,mdist,pnorm,seed

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

      do i = 1,np
          pnorm = 0.0d0
          do j = 1,ndim
              x(ndim*(i-1)+j) = 2.0d0 * drand(seed) - 1.0d0
              pnorm = pnorm + x(ndim*(i-1)+j) ** 2
          end do
          pnorm = sqrt( pnorm )
          do j = 1,ndim
              x(ndim*(i-1)+j) = x(ndim*(i-1)+j) / pnorm
          end do
      end do

      mdist = 0.0d0
      do i = 1,np
          do j = i + 1,np
              dist = 0.0d0
              do k = 1,ndim
                  ind1 = ndim * ( i - 1 ) + k
                  ind2 = ndim * ( j - 1 ) + k
                  dist = dist + ( x(ind1) - x(ind2) ) ** 2
              end do
              mdist = max( mdist, dist )
          end do
      end do

      x(n) = mdist

C     Lower and upper bounds
      do i = 1,n
          l(i) = - 1.0d+20
          u(i) =   1.0d+20
      end do

C     Number of constraints (equalities plus inequalities)

      m = np * ( np - 1 )

C     Lagrange multipliers approximation.

      do i = 1,m
          lambda(i) = 0.0d0
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
      integer i,i1,i2,ind1,ind2

C     COMMON BLOCKS
      common /probdata/ pair,ndim,np

      flag = 0

      if ( 1 .le. ind .and. ind .le. np * ( np - 1 ) / 2 ) then

          i1 = pair(1,ind)
          i2 = pair(2,ind)

          c = 1.0d0

          do i = 1,ndim
              ind1 = ndim * ( i1 - 1 ) + i
              ind2 = ndim * ( i2 - 1 ) + i
              c = c - ( x(ind1) - x(ind2) ) ** 2
          end do

      else if ( np * ( np - 1 ) / 2 + 1 .le. ind .and.
     +          ind .le. np * (np - 1) ) then

          i1 = pair(1,ind - np * ( np - 1 ) / 2)
          i2 = pair(2,ind - np * ( np - 1 ) / 2)

          c = - x(n)
          do i = 1,ndim
              ind1 = ndim * ( i1 - 1 ) + i
              ind2 = ndim * ( i2 - 1 ) + i
              c = c + ( x(ind1) - x(ind2) ) ** 2
          end do

      else
          flag = - 1
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
      integer i,i1,i2,ind1,ind2
      double precision diff

C     COMMON BLOCKS
      common /probdata/ pair,ndim,np

      flag = 0

      if ( 1 .le. ind .and. ind .le. np * ( np - 1 ) / 2 ) then

          i1 = pair(1,ind)
          i2 = pair(2,ind)

          jcnnz = 2 * ndim

          do i = 1,ndim
              ind1 = ndim * ( i1 - 1 ) + i
              ind2 = ndim * ( i2 - 1 ) + i
              diff = 2.0d0 * ( x(ind1) - x(ind2) )
              jcvar(i)      =   ind1
              jcval(i)      = - diff
              jcvar(ndim+i) =   ind2
              jcval(ndim+i) =   diff
          end do

      else if ( np * ( np - 1 ) / 2 + 1 .le. ind .and.
     +          ind .le. np * (np - 1) ) then

          i1 = pair(1,ind - np * ( np - 1 ) / 2)
          i2 = pair(2,ind - np * ( np - 1 ) / 2)

          jcnnz = 2 * ndim + 1

          do i = 1,ndim
              ind1 = ndim * ( i1 - 1 ) + i
              ind2 = ndim * ( i2 - 1 ) + i
              diff = 2.0d0 * ( x(ind1) - x(ind2) )
              jcvar(i)      =   ind1
              jcval(i)      =   diff
              jcvar(ndim+i) =   ind2
              jcval(ndim+i) = - diff
          end do

          jcvar(jcnnz) = n
          jcval(jcnnz) = - 1.0d0

      else
          flag = - 1
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
      integer i,i1,i2,ind1,ind2

C     COMMON BLOCKS
      common /probdata/ pair,ndim,np

      flag = 0

      if ( 1 .le. ind .and. ind .le. np * ( np - 1 ) / 2 ) then

          i1 = pair(1,ind)
          i2 = pair(2,ind)

          hcnnz = 3 * ndim

          do i = 1,ndim
              ind1 = ndim * ( i1 - 1 ) + i
              ind2 = ndim * ( i2 - 1 ) + i

              hclin(3*(i-1)+1) = ind1
              hccol(3*(i-1)+1) = ind1
              hcval(3*(i-1)+1) = - 2.0d0

              hclin(3*(i-1)+2) = ind2
              hccol(3*(i-1)+2) = ind1
              hcval(3*(i-1)+2) = 2.0d0

              hclin(3*(i-1)+3) = ind2
              hccol(3*(i-1)+3) = ind2
              hcval(3*(i-1)+3) = - 2.0d0
          end do

      else if ( np * ( np - 1 ) / 2 + 1 .le. ind .and.
     +          ind .le. np * (np - 1) ) then

          i1 = pair(1,ind - np * ( np - 1 ) / 2)
          i2 = pair(2,ind - np * ( np - 1 ) / 2)

          hcnnz = 3 * ndim

          do i = 1,ndim
              ind1 = ndim * ( i1 - 1 ) + i
              ind2 = ndim * ( i2 - 1 ) + i

              hclin(3*(i-1)+1) = ind1
              hccol(3*(i-1)+1) = ind1
              hcval(3*(i-1)+1) = 2.0d0

              hclin(3*(i-1)+2) = ind2
              hccol(3*(i-1)+2) = ind1
              hcval(3*(i-1)+2) = - 2.0d0

              hclin(3*(i-1)+3) = ind2
              hccol(3*(i-1)+3) = ind2
              hcval(3*(i-1)+3) = 2.0d0
          end do

      else
          flag = - 1
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

C     PARAMETERS
      integer npmax
      parameter ( npmax = 1000 )

C     COMMON SCALARS
      integer ndim,np

C     COMMON ARRAYS
      integer pair(2,npmax * ( npmax - 1 ) / 2)

C     COMMON BLOCKS
      common /probdata/ pair,ndim,np

      call drawsol(ndim,np,n,x,'solution.mp')

      end

C     ******************************************************
C     ******************************************************

      subroutine drawsol(ndim,np,n,x,solfile)

C     This subroutine generate a metapost file with the
C     graphical representation of the problem.

C     SCALAR ARGUMENTS
      integer n,ndim,np

C     ARRAY ARGUMENTS
      double precision x(n)
      character * 11 solfile

C     LOCAL SCALARS
      integer i,i1,i2,j,j1,j2
      double precision scale

C     JSUT 2D PICTURES
      if ( ndim .ne. 2 ) then
          return
      end if

C     OPENING OUTPUT FILE
      open(unit=10,file=solfile)

C     SCALING
      scale = 10.0d0 / x(n)
      write(10,10) scale

C     LINKS
      do i = 1,np
          i1 = ndim * ( i - 1 ) + 1
          i2 = ndim * ( i - 1 ) + 2
          do j = i + 1,np
              j1 = ndim * ( j - 1 ) + 1
              j2 = ndim * ( j - 1 ) + 2
              write(10,20) x(i1),x(i2),x(j1),x(j2)
          end do
      end do

C     POINTS
      write(10,30)
      do i = 1,np
          i1 = ndim * ( i - 1 ) + 1
          i2 = ndim * ( i - 1 ) + 2
          write(10,40) x(i1),x(i2)
      end do

      write(10,50)

      close(10)

C     NON-EXECUTABLE STATEMENTS

  10  format('beginfig(1);'/,'u = ',f20.10,' cm;')
  20  format('draw (',f20.10,'u,',f20.10,'u)--',
     +            '(',f20.10,'u,',f20.10,'u);')
  30  format('pickup pencircle scaled 4pt;')
  40  format('drawdot (',f20.10,'u,',f20.10,'u);')
  50  format('endfig;',/,'end;')

      end
