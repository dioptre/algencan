C     =================================================================
C     File: ellipsoid.f
C     =================================================================

C     =================================================================
C     Module: Subroutines that define the problem
C     =================================================================

C     Last update of any of the component of this module:

C     January 31, 2007.

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

C     Smallest enclosing ellipsoid problem
C     ------------------------------------

C     Description:

C     Given npun fixed points pi in R^ndim, minimize the volume of the
C     ellipsoid, centered at the origin, that fits them. The model can
C     be written as:
C
C                  Minimize - \sum_{i=1,ndim} log ( l_{ii} )
C
C                  subject to
C
C                          pi^T L L^T pi <= 1, i = 1, ..., npun,
C
C                          l_{ii} >= 0, i = 1, ..., ndim,
C
C     where L is ndim x ndim lower-triangular matrix.
C
C     In this version we aim to find the smallest encolsing ellipsoid
C     centered at the origin when considering the 4278 points in R^3
C     representing the atoms of a protein of the Isoform alpha of the
C     Thyroid Hormone Receptor.
C
C     References:
C
C     M.J. Todd and E.A. Yildirim, ``On Khachiyan's algorithm for the
C     computation of minimum volume enclosing ellipsoids'' (9/05).

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
      integer ndimmax,npunmax
      parameter ( ndimmax =    10 )
      parameter ( npunmax = 20000 )

C     COMMON SCALARS
      integer ndim,npun

C     COMMON ARRAYS
      integer E(ndimmax,ndimmax)
      double precision p(ndimmax,npunmax)

C     LOCAL SCALARS
      integer i,j,k

C     LOCAL ARRAYS
      double precision c(ndimmax)

C     COMMON BLOCKS
      common /probdata/ p,E,ndim,npun

C     Set problem data (using ellipsoid.dat)

      read(*,*) ndim,npun

      do j = 1,npun
          read(*,*) (p(i,j),i=1,ndim)
      end do

C     or asking to the user and generating the point by
C     Cauchy distribution

C     write(*,*) 'Dimension of the space (ndim) maximum 10: '
C     read(*,*) ndim

C     write(*,*) 'Number of points (npun) maximum 20,000: '
C     read(*,*) npun

C     write(*,*) 'Seed for the random generation of point: '
C     read(*,*) seed

C     Alternatively
C     read(*,*) ndim,npun,seed
C     seed = 123456.0d0 * seed

C     Generate the npun points in R^ndim using the Cauchy distribution

C     do j = 1,npun
C         do i = 1,ndim
C             call cauchyd(seed,p(i,j))
C         end do
C     end do

C     Compute the center of mass and substract from each point

      do i = 1,ndim
          c(i) = 0.0d0
      end do

      do j = 1,npun
          do i = 1,ndim
              c(i) = c(i) + p(i,j)
          end do
      end do

      do i = 1,ndim
          c(i) = c(i) / npun
      end do

      do j = 1,npun
          do i = 1,ndim
              p(i,j) = p(i,j) - c(i)
          end do
      end do

C     Set the mapping between the elements of the lower-triagular matrix
C     L and its column-wise representation at vector x.

      k = 0
      do j = 1,ndim
          do i = j,ndim
              k = k + 1
              E(i,j) = k
          end do
      end do

C     Number of variables (elements of a lower-triangular ndim x ndim matrix)
      n = ndim * ( ndim + 1 ) / 2

C     Lower and upper bounds
      do j= 1,ndim
          l(E(j,j)) = 1.0d-16
          u(E(j,j)) = 1.0d+20
          do i = j + 1,ndim
              l(E(i,j)) = - 1.0d+20
              u(E(i,j)) =   1.0d+20
          end do
      end do

C     Initial point

      do j = 1,ndim
          x(E(j,j)) = 1.0d0
          do i = j + 1,ndim
              x(E(i,j)) = 0.0d0
          end do
      end do

C     Number of constraints (equalities plus inequalities)

      m = npun

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
      integer ndimmax,npunmax
      parameter ( ndimmax =    10 )
      parameter ( npunmax = 20000 )

C     COMMON SCALARS
      integer ndim,npun

C     COMMON ARRAYS
      integer E(ndimmax,ndimmax)
      double precision p(ndimmax,npunmax)

C     LOCAL SCALARS
      integer i

C     COMMON BLOCKS
      common /probdata/ p,E,ndim,npun

      flag = 0

      f = 0.0d0
      do i = 1,ndim
          f = f - log( x(E(i,i)) )
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
      integer ndimmax,npunmax
      parameter ( ndimmax =    10 )
      parameter ( npunmax = 20000 )

C     COMMON SCALARS
      integer ndim,npun

C     COMMON ARRAYS
      integer E(ndimmax,ndimmax)
      double precision p(ndimmax,npunmax)

C     LOCAL SCALARS
      integer i,j

C     COMMON BLOCKS
      common /probdata/ p,E,ndim,npun

      flag = 0

      do j = 1,ndim
          g(E(j,j)) = - 1.0d0 / x(E(j,j))
          do i = j + 1,ndim
              g(E(i,j)) = 0.0d0
          end do
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

C     PARAMETERS
      integer ndimmax,npunmax
      parameter ( ndimmax =    10 )
      parameter ( npunmax = 20000 )

C     COMMON SCALARS
      integer ndim,npun

C     COMMON ARRAYS
      integer E(ndimmax,ndimmax)
      double precision p(ndimmax,npunmax)

C     LOCAL SCALARS
      integer i

C     COMMON BLOCKS
      common /probdata/ p,E,ndim,npun

      flag = 0

      hnnz = ndim

      do i = 1,ndim
          hlin(i) = E(i,i)
          hcol(i) = E(i,i)
          hval(i) = 1.0d0 / x(E(i,i)) ** 2
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
      integer ndimmax,npunmax
      parameter ( ndimmax =    10 )
      parameter ( npunmax = 20000 )

C     COMMON SCALARS
      integer ndim,npun

C     COMMON ARRAYS
      integer E(ndimmax,ndimmax)
      double precision p(ndimmax,npunmax)

C     LOCAL SCALARS
      integer i,j

C     LOCAL ARRAYS
      double precision y(ndimmax)

C     COMMON BLOCKS
      common /probdata/ p,E,ndim,npun

      flag = 0

C     Compute y = L^T p_{ind}

      do j = 1,ndim
          y(j) = 0.0d0
          do i = j,ndim
              y(j) = y(j) + x(E(i,j)) * p(i,ind)
          end do
      end do

C     Compute c = y^t y - 1.0d0

      c = - 1.0d0
      do j = 1,ndim
          c = c + y(j) ** 2
      end do

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
      integer ndimmax,npunmax
      parameter ( ndimmax =    10 )
      parameter ( npunmax = 20000 )

C     COMMON SCALARS
      integer ndim,npun

C     COMMON ARRAYS
      integer E(ndimmax,ndimmax)
      double precision p(ndimmax,npunmax)

C     LOCAL SCALARS
      integer i,j

C     LOCAL ARRAYS
      double precision y(ndimmax)

C     COMMON BLOCKS
      common /probdata/ p,E,ndim,npun

      flag = 0

      do j = 1,ndim
          y(j) = 0.0d0
          do i = j,ndim
              y(j) = y(j) + x(E(i,j)) * p(i,ind)
          end do
      end do

      jcnnz = ndim * ( ndim + 1 ) / 2

      do j = 1,ndim
          do i = j,ndim
              jcvar(E(i,j)) = E(i,j)
              jcval(E(i,j)) = 2.0d0 * y(j) * p(i,ind)
          end do
      end do

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
      integer ndimmax,npunmax
      parameter ( ndimmax =    10 )
      parameter ( npunmax = 20000 )

C     COMMON SCALARS
      integer ndim,npun

C     COMMON ARRAYS
      integer E(ndimmax,ndimmax)
      double precision p(ndimmax,npunmax)

C     LOCAL SCALARS
      integer i,j,k

C     COMMON BLOCKS
      common /probdata/ p,E,ndim,npun

      flag = 0

      hcnnz = 0
      do j = 1,ndim
          do i = j,ndim
              do k = i,ndim
                  hcnnz = hcnnz + 1
                  hclin(hcnnz) = E(k,j)
                  hccol(hcnnz) = E(i,j)
                  hcval(hcnnz) = 2.0d0 * p(k,ind) * p(i,ind)
              end do
          end do
      end do

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

      subroutine drawsol(n,x)

C     SCALAR ARGUMENTS
      integer n

C     ARRAY ARGUMENTS
      double precision x(n)

C     This subroutine generate a metapost file with the graphical
C     representation of the problem (just for ndim=2).

C     PARAMETERS
      integer ndimmax,npunmax
      double precision pi
      parameter ( ndimmax =                10 )
      parameter ( npunmax =             20000 )
      parameter ( pi      = 3.1415926535898d0 )

C     COMMON SCALARS
      integer ndim,npun

C     COMMON ARRAYS
      integer E(ndimmax,ndimmax)
      double precision p(ndimmax,npunmax)

C     LOCAL SCALARS
      integer i,j,nprob
      double precision xmin,xmax,ymin,ymax,scale,aa,bb,cc,dd,alpha

C     LOCAL ARRAYS
      double precision lambda(2),M(2,2),v(2)

C     COMMON BLOCKS
      common /probdata/ p,E,ndim,npun

C     JUST FOR NDIM=2!

      if ( ndim .ne. 2 ) then
          return
      end if

C     PROBLEM ID

      nprob = 1

C     SCALING

      xmin =   1.0d+99
      xmax = - 1.0d+99

      do j = 1,npun
          xmin = min( xmin, p(1,j) )
          xmax = max( xmax, p(1,j) )
      end do

      ymin =   1.0d+99
      ymax = - 1.0d+99

      do j = 1,npun
          ymin = min( ymin, p(2,j) )
          ymax = max( ymax, p(2,j) )
      end do

      scale = min( 10.0d0 / ( xmax - xmin ), 10.0d0 / ( ymax - ymin ) )

C     DRAW THE PROBLEM

      open(unit=10,file='sol.mp')

      write(10,10) nprob,1.0d0

C     ELLIPSE

C     COMPUTE M = L L^T

      M(1,1) = x(E(1,1)) ** 2
      M(2,1) = x(E(1,1)) * x(E(2,1))
      M(2,2) = x(E(2,1)) ** 2 + x(E(2,2)) ** 2

C     COMPUTE M EIGENVALUES

      aa = 1.0d0
      bb = - ( M(1,1) + M(2,2) )
      cc = M(1,1) * M(2,2) - M(2,1) ** 2

      dd = sqrt( bb ** 2 - 4.0d0 * aa * cc )

      lambda(1) = ( - bb + dd ) / ( 2.0d0 * aa )
      lambda(2) = ( - bb - dd ) / ( 2.0d0 * aa )

C     COMPUTE ANGLE (USING THE EIGENVECTOR RELATED TO LAMBDA1)

      if ( abs( M(2,1) ) .gt. 1.0d-08 ) then
          v(1) = 1.0d0
          v(2) = - ( M(1,1) - lambda(1) ) / M(2,1)
      else
          v(2) = 1.0d0
          v(1) = - M(2,1) / ( M(1,1) - lambda(1) )
      end if

      alpha = v(1) / sqrt( v(1) ** 2 + v(2) ** 2 )
      alpha = acos(alpha) / pi * 180.0d0
      alpha = sign(alpha,v(2))

C     COMPUTE ELLIPSE SEMI-AXIS

      lambda(1) = sqrt( 1.0d0 / lambda(1) )
      lambda(2) = sqrt( 1.0d0 / lambda(2) )

      write(10,20) 2.0d0*lambda(1)*scale,
     +             2.0d0*lambda(2)*scale,0.0d0,0.0d0,alpha

C     DOTS WITHOUT LABELS

      write(10,30)

      do i = 1,npun
          write(10,40) scale*p(1,i),scale*p(2,i)
      end do

C     END

      write(10,100)

      close(10)

C     NON-EXECUTABLE STATEMENTS

  10  format('beginfig(',I2,');'/,'u = ',f20.10,' cm;')
  20  format('draw fullcircle',
     +       /,5X,'xscaled  ',f20.10,'u',
     +       /,5X,'yscaled  ',f20.10,'u',
     +       /,5X,'shifted (',f20.10,'u,',f20.10,'u) ',
     +       /,5X,'rotated  ',f20.10,';')
  30  format('pickup pencircle scaled 4pt;')
  40  format('draw (',f20.10,'u,',f20.10,'u);')
 100  format('endfig;'/,'end;')

      end
