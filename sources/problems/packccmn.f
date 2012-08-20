C     =================================================================
C     File: packccmn.f
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

C     Packing fixed-dimension circular items within a fixed-dimension
C     ---------------------------------------------------------------
C     circular object maximizing the number of items
C     ----------------------------------------------
C
C     We wish to place k circles of radii ri (i=1, ..., k) into a
C     circular object with radius R in such a way that the circles are
C     not overlapped. Therefore, given k, the radii ri (i=1, ..., k)
C     and R, the goal is to solve the problem:
C
C     Minimize \sum_{i not equal j} max(0, (ri+rj)^2 - d(pi,pj)^2 )^2
C
C     subject to d(0,pi)^2 <= (R - ri)^2, i = 1, ..., k,
C
C     where d(.,.) is the Euclidian distance. If the objective function
C     value at the global minimizer of this problem is zero then the
C     answer of the decision problem: "Given k circular items of radii
C     ri (i=1, ..., k) and a circular object of radius R, whether is it
C     possible to locate all the items within the object or not." is YES,
C     otherwise, the answer is NO.

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

C     COMMON SCALARS
      integer nite,ndim
      double precision iterad,objrad

C     LOCAL SCALARS
      integer i
      double precision drand,seed,r,t

C     COMMON BLOCKS
      common /packdata/ nite,ndim,iterad,objrad

C     SET PROBLEM DATA

C     Dimension of the space
      ndim   =      2

C     Number of identical circular items to be packed
      nite   =     80

C     Radius of the circular items to be packed
      iterad =  1.0d0

C     Radius of the circular object within which the items will
C     be packed
      objrad = 10.0d0

C     Number of variables

      n = 2 * nite

C     Initial point

      seed = 12337.0d0

      do i = 1,nite
          r = ( objrad - iterad ) * drand(seed)
          t = 2.0d0 * 3.14159d0 * drand(seed)
          x(2*i-1) = r * cos( t )
          x(2*i  ) = r * sin( t )
      end do

C     Lower and upper bounds
      do i = 1,n
          l(i) = - 1.0d+20
          u(i) =   1.0d+20
      end do

C     Number of constraints (equalities plus inequalities)

      m = nite

C     Lagrange multipliers approximation

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

C     COMMON SCALARS
      integer ndim,nite
      double precision iterad,objrad

C     LOCAL SCALARS
      integer i,j,k,ind1,ind2
      double precision fparc,dist

C     COMMON BLOCKS
      common /packdata/ nite,ndim,iterad,objrad

      flag = 0

C     COMPUTE DENSE OVERLAPPING

      f = 0.0d0
      do i = 1,nite
          do j = i + 1,nite
              dist = 0.0d0
              do k = 1,ndim
                  ind1 = ( i - 1 ) * ndim + k
                  ind2 = ( j - 1 ) * ndim + k
                  dist = dist + ( x(ind1) - x(ind2) ) ** 2
              end do
              fparc = max( 0.0d0, ( 2.0d0 * iterad ) ** 2 - dist )
              f = f + fparc ** 2
          end do
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

C     COMMON SCALARS
      integer ndim,nite
      double precision iterad,objrad

C     LOCAL SCALARS
      integer i,ind1,ind2,j,k
      double precision dist,fparc,tmp

C     COMMON BLOCKS
      common /packdata/ nite,ndim,iterad,objrad

      flag = 0

C     COMPUTE DENSE OVERLAPPING

      do i = 1,n
          g(i) = 0.0d0
      end do

      do i = 1,nite
          do j = i + 1,nite

              dist = 0.0d0
              do k = 1,ndim
                  ind1 = ( i - 1 ) * ndim + k
                  ind2 = ( j - 1 ) * ndim + k
                  dist = dist + ( x(ind1) - x(ind2) ) ** 2
              end do
              fparc = max( 0.0d0, ( 2.0d0 * iterad ) ** 2 - dist )

              if ( fparc .ne. 0.d0 ) then
                  do k = 1,ndim
                      ind1 = ( i - 1 ) * ndim + k
                      ind2 = ( j - 1 ) * ndim + k
                      tmp = 4.0d0 * fparc * ( x(ind1) - x(ind2) )
                      g(ind1) = g(ind1) - tmp
                      g(ind2) = g(ind2) + tmp
                  end do
              end if

          end do
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

C     PARAMETERS
      integer nitemax,ndimmax
      parameter ( nitemax =        10000 )
      parameter ( ndimmax =            2 )

C     COMMON SCALARS
      integer ndim,nite
      double precision iterad,objrad

C     LOCAL SCALARS
      integer i,j,k,l,ind1,ind2,ind3,ind4
      double precision dist,tmp

C     LOCAL ARRAYS
      double precision diagb(ndimmax,ndimmax,nitemax)

C     COMMON BLOCKS
      common /packdata/ nite,ndim,iterad,objrad

      flag = 0

C     INITALIZE DIAGONAL BLOCKS

      do k = 1,nite
          do j = 1,ndim
              do i = j,ndim
                  diagb(i,j,k) = 0.0d0
              end do
          end do
      end do

C     COMPUTE DENSE OVERLAPPING SECOND DERIVATIVES

      hnnz = 0

      do i = 1,nite
          do j = i + 1,nite

              dist = 0.0d0
              do k = 1,ndim
                  ind1 = ( i - 1 ) * ndim + k
                  ind2 = ( j - 1 ) * ndim + k
                  dist = dist + ( x(ind1) - x(ind2) ) ** 2
              end do

              if ( dist .le. ( 2.0d0 * iterad ) ** 2 ) then
                  do k = 1,ndim
                      ind1 = ( i - 1 ) * ndim + k
                      ind2 = ( j - 1 ) * ndim + k
                      tmp = 8.0d0 * ( x(ind1) - x(ind2) ) ** 2
     +                    - 4.0d0 * ( ( 2.0d0 * iterad ) ** 2 - dist )
                      if ( tmp .ne. 0.0d0 ) then
C                         H(ind1,ind1) = H(ind1,ind1) + tmp
                          diagb(k,k,i) = diagb(k,k,i) + tmp
C                         H(ind2,ind2) = H(ind2,ind2) + tmp
                          diagb(k,k,j) = diagb(k,k,j) + tmp
C                         H(ind2,ind1) = H(ind2,ind1) - tmp
                          hnnz = hnnz + 1
                          hlin(hnnz) = ind2
                          hcol(hnnz) = ind1
                          hval(hnnz) = - tmp
                      end if
                      do l = 1,k - 1
                          ind3 = ( i - 1 ) * ndim + l
                          ind4 = ( j - 1 ) * ndim + l
                          tmp = 8.0d0 * ( x(ind3) - x(ind4) )
     +                                * ( x(ind1) - x(ind2) )
                          if ( tmp .ne. 0.0d0 ) then
C                             H(ind1,ind3) = H(ind1,ind3) + tmp
                              diagb(k,l,i) = diagb(k,l,i) + tmp
C                             H(ind2,ind4) = H(ind2,ind4) + tmp
                              diagb(k,l,j) = diagb(k,l,j) + tmp
C                             H(ind2,ind3) = H(ind2,ind3) - tmp
                              hnnz = hnnz + 1
                              hlin(hnnz) = ind2
                              hcol(hnnz) = ind3
                              hval(hnnz) = - tmp
                          end if
                      end do
                      do l = k + 1,ndim
                          ind3 = ( i - 1 ) * ndim + l
                          ind4 = ( j - 1 ) * ndim + l
                          tmp = 8.0d0 * ( x(ind3) - x(ind4) )
     +                                * ( x(ind1) - x(ind2) )
                          if ( tmp .ne. 0.0d0 ) then
C                             H(ind2,ind3) = H(ind2,ind3) - tmp
                              hnnz = hnnz + 1
                              hlin(hnnz) = ind2
                              hcol(hnnz) = ind3
                              hval(hnnz) = - tmp
                          end if
                      end do
                  end do
              end if

          end do
      end do

C     ADD DIAGONAL BLOCKS TO THE SPARSE STRUCTURE

      do k = 1,nite
          do j = 1,ndim
              do i = j,ndim
                  if ( diagb(i,j,k) .ne. 0.0d0 ) then
                      hnnz = hnnz + 1
                      hlin(hnnz) = ( k - 1 ) * ndim + i
                      hcol(hnnz) = ( k - 1 ) * ndim + j
                      hval(hnnz) = diagb(i,j,k)
                  end if
              end do
          end do
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

C     COMMON SCALARS
      integer nite,ndim
      double precision iterad,objrad

C     LOCAL SCALARS
      integer i

C     COMMON BLOCKS
      common /packdata/ nite,ndim,iterad,objrad

      flag = 0

      if ( 0 .lt. ind .and. ind .le. nite ) then

          c = - ( objrad -iterad ) ** 2

          do i = 1,ndim
              c = c + x(ndim * ( ind - 1 ) + i) ** 2
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

C     COMMON SCALARS
      integer nite,ndim
      double precision iterad,objrad

C     LOCAL SCALARS
      integer i

C     COMMON BLOCKS
      common /packdata/ nite,ndim,iterad,objrad

      flag = 0

      if ( 0 .lt. ind .and. ind .le. nite ) then

          jcnnz = ndim

          do i = 1,ndim
              jcvar(i) = ndim * ( ind - 1 ) + i
              jcval(i) = 2.0d0 * x(ndim * ( ind - 1 ) + i)
          end do

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

C     COMMON SCALARS
      integer nite,ndim
      double precision iterad,objrad

C     LOCAL SCALARS
      integer k

C     COMMON BLOCKS
      common /packdata/ nite,ndim,iterad,objrad

      flag = 0

      if ( 0 .lt. ind .and. ind .le. nite ) then

          hcnnz = ndim

          do k = 1,ndim
              hclin(k) = ndim * ( ind - 1 ) + k
              hccol(k) = ndim * ( ind - 1 ) + k
              hcval(k) = 2.0d0
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

