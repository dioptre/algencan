C     =================================================================
C     File: piecefit.f
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

C     Piecefit problem
C     ----------------

C     We defined usol, a piecewise linear function defined on [0,4] 
C     with kinks in 1 and 3 and discotinuity in 2. We discretized [0,4] 
C     with n points, including extremes. We defined uobs, a 10% 
C     perturbation of usol. The objective function is the sum of the p 
C     smaller squares of the discretized second derivatives of x. p, 
C     the OVO parameter, is given by the user. The restriction is that 
C     the root mean squared deviation of the solution found with 
C     respect to uobs is smaller than 0.2. In other words, using LOVO, 
C     this code aims to find a piecewise harmonic one-dimensional 
C     function. 

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
      integer p
      double precision h

C     COMMON ARRAYS
      double precision uobs(2000),usol(2000)

C     LOCAL SCALARS
      integer i
      double precision drand,seed,t

C     COMMON BLOCKS
      common /ovo/ uobs,usol,h,p
      
      write(*,*) 'Number of variables (points in [0,4]): '
      read(*,*) n

      write(*,*) 'p (integer close to n-2), perhaps n-3: '
      read(*,*) p 

C     Generation of uobs

      seed = 17172937.0d0

      h = 4.0d0 / dfloat(n - 1)

      do i = 1,n
          t = dfloat(i-1) * h

          if ( t .le. 1.0d0 ) then
              uobs(i) = 1.0d0 + t

          else if ( t .le. 2.0d0 ) then
              uobs(i) = 1.0d0 - ( t - 2.0d0 )

          else if ( t .le. 3.0d0 ) then
              uobs(i) = 2.0d0

          else 
              uobs(i) = 2.0d0 - ( t - 3.0d0 )
          end if

          usol(i) = uobs(i)
          uobs(i) = uobs(i) + 
     +              uobs(i) * 0.1d0 * ( 2.0d0 * drand(seed) - 1.0d0 )
      end do

C     write(*,*) 'uobs: ',(uobs(i),i=1,n)

C     Initial point

      do i = 1,n
          x(i) = uobs(i)
      end do
      
C     Lower and upper bounds

      do i = 1,n
          l(i) = - 1.0d+20
          u(i) =   1.0d+20
      end do

C     Number of constraints (equalities plus inequalities)

      m = 1

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

C     COMMON SCALARS
      integer p
      double precision h

C     COMMON ARRAYS
      double precision uobs(2000),usol(2000)

C     LOCAL SCALARS
      integer i,imax,j
      double precision zmax

C     LOCAL ARRAYS
      double precision resid(2000)

C     COMMON BLOCKS
      common /ovo/ uobs,usol,h,p
      
      flag = 0

      f = 0.0d0

      do i = 1,n - 2
          resid(i) = ( x(i+2) - 2.0d0 * x(i+1) + x(i) ) / h ** 2
      end do

      do j = 1,p
          zmax = 1.0d+99
          do i = 1,n - 2
              if ( abs( resid(i) ) .lt. zmax ) then
                  zmax = abs( resid(i) )
                  imax = i
              end if
          end do

          f = f + resid(imax) ** 2
          resid(imax) = 1.0d+99
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
      integer p
      double precision h

C     COMMON ARRAYS
      double precision uobs(2000),usol(2000)

C     LOCAL SCALARS
      integer i,imax,j
      double precision zmax

C     LOCAL ARRAYS
      double precision resid(2000)

C     COMMON BLOCKS
      common /ovo/ uobs,usol,h,p
      
      flag = 0

      do i = 1,n
          g(i) = 0.0d0
      end do

      do i = 1,n - 2
          resid(i) = ( x(i+2) - 2.0d0 * x(i+1) + x(i) ) / h ** 2
      end do

      do j = 1,p
          zmax = 1.0d+99
          do i = 1,n - 2
              if ( abs( resid(i) ) .lt. zmax ) then
                  zmax = abs( resid(i) )
                  imax = i
              end if
          end do

          g(imax+2) = g(imax+2) + 2.0d0 * resid(imax) / h ** 2
          g(imax+1) = g(imax+1) - 4.0d0 * resid(imax) / h ** 2
          g(imax)   = g(imax)   + 2.0d0 * resid(imax) / h ** 2   

          resid(imax) = 1.0d+99
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

C     COMMON SCALARS
      integer p
      double precision h

C     COMMON ARRAYS
      double precision uobs(2000),usol(2000)

C     LOCAL SCALARS
      integer i
      double precision z

C     COMMON BLOCKS
      common /ovo/ uobs,usol,h,p
      
      flag = 0

      z = 0.0d0
      do i = 1,n
          z = z + ( x(i) - uobs(i) ) ** 2
      end do

      z = z / dfloat(n)

      c = z - 0.01d0 

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
      integer p
      double precision h

C     COMMON ARRAYS
      double precision uobs(2000),usol(2000)

C     LOCAL SCALARS
      integer i

C     COMMON BLOCKS
      common /ovo/ uobs,usol,h,p
      
      flag = 0

      jcnnz = n

      do i = 1,n
          jcvar(i) = i
          jcval(i) = 2.0d0 * ( x(i) - uobs(i) ) / dfloat(n)
      end do

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

C     COMMON SCALARS
      integer p
      double precision h

C     COMMON ARRAYS
      double precision uobs(2000),usol(2000)

C     LOCAL SCALARS
      integer i
      double precision z,zmax

C     COMMON BLOCKS
      common /ovo/ uobs,usol,h,p
      
      zmax = 0.0d0

      write(*,*) 'i, t, sol-found, true-sol, second-derivative'

      do i = 1,n
          if ( i .eq. 1 .or. i .eq. n ) then
              z = 0.0d0
          else
              z = ( x(i+1) - 2.0d0 * x(i) + x(i-1) ) / h ** 2
          end if

          write(*,*) i,(i-1)*h,x(i),usol(i),z

          zmax = max( zmax, z ) 
      end do

      write(*, *)' Maximum second derivative: ',zmax

      end
