C     =================================================================
C     File: condor.f
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

C     Condor problem
C     --------------

C     We generate (using subroutine evalu) a function called usol as a 
C     Fourier polynomial with up to q = 20 terms. A grid in [0, 2 pi] 
C     with up to 1000 points is established. We generate function fi(t) 
C     in such a way that usol is a solution of the discretization of 
C     the differential equation u''(t) + [u'(t)]^2 + u(t) = fi(t). We
C     select up to 1000 observation points equally spaced in [-2, 18]
C     and observe the solution usol at these observation points, 
C     obtaining the observed solution uobs. Finally, uint, the 
C     interpolation of uobs in the grid points is computed. The 
C     optimization problem is to obtain the solution of the 
C     discretized differential equation which is closest, in the 
C     least-squares sense, to uint. Instances with q = 20, 500 
C     observation points and 1000 grid points are challenging.

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
      double precision h

C     COMMON ARRAYS
      double precision uint(1000),fon(1000),usol(1000)

C     LOCAL SCALARS
      integer i,j,p,q
      double precision a0,z,pen,aux,u1,u2

C     LOCAL ARRAYS
      double precision a(20),b(20),t(1000),uobs(1000)

C     COMMON BLOCKS
      common /objective/ uint
      common /constraints/ fon,h
      common /solution/ usol
      
      write(*, *) 'Number of points in the grid (nmax=1000): '
      read(*,*) n

      write(*, *) 'Number of observed point (pmax<=nmax): '
      read(*,*) p
      write(*,*) 'Number of terms of the Fourier polynomial ',
     +           '(qmax=20): '
      read(*,*) q

C     Fourier polynomial coefficients

      do i = 1,q
          a(i) = cos( dfloat(i) )      / dfloat(i)
          b(i) = sin( dfloat(i) ) ** 2 / dfloat(i)
      end do

      a0 = 2.0d0
      z  = 0.0d0

      do i = 1,p
          t(i) = -2.0d0 + ( 20.0d0 / dfloat(p) ) * i
      end do

      do i = 1,p
          call evalu(a0,a,b,q,t(i),uobs(i),u1,u2,z)
      end do

C     write(*,*) 't = ',(t(i),i=1,p)
C     write(*,*) 'uobs = ',(uobs(i),i=1,p)

C     Independent term of the differential equation (fon)

      h = atan(1.0d0) * 8.0d0 / dfloat(n - 1)

      do i = 1,n - 2
          z = i * h
          call evalu(a0,a,b,q,z,aux,u1,u2,fon(i))
      end do

C     Analytic solution at each grid point (usol)

      do i = 1,n
          z = ( i - 1 ) * h
          call evalu(a0,a,b,q,z,usol(i),u1,u2,aux)
      end do
 
C     Interpolation

      do i = 1,n
          z = dfloat(i - 1) * h

          if( z .le. t(2) ) then
              pen = ( uobs(2) - uobs(1) ) / ( t(2) - t(1) )
              uint(i) = uobs(1) + pen * ( z - t(1) )

          else if( z .ge. t(p-1) )then
              pen = ( uobs(p) - uobs(p-1) ) / ( t(p) - t(p-1) )
              uint(i) = uobs(p-1) + pen * ( z - t(p-1) )

          else
              do j = 2,p - 2
                  if ( z .ge. t(j) .and. z .le. t(j+1) ) then
                      pen = ( uobs(j+1) - uobs(j) ) / ( t(j+1) - t(j) )
                      uint(i) = uobs(j) + pen * ( z - t(j) )
                  end if
              end do
          end if
      end do      

C     write(*,*) 'Interpolated points: '
C     do i = 1,n
C         write(*,*) i,(i-1)*h,uint(i)
C     end do
      
C     Initial point

      do i = 1,n
          x(i) = uint(i)
      end do
      
C     Lower and upper bounds

      do i = 1,n
          l(i) = - 1.0d+20
          u(i) =   1.0d+20
      end do

C     Number of constraints (equalities plus inequalities)

      m = n - 2

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

      subroutine evalu(a0,a,b,m,x,u,u1,u2,fi)

      implicit none

C     SCALAR ARGUMENTS
      integer m
      double precision a0,fi,u,u1,u2,x

C     ARRAY ARGUMENTS
      double precision a(20),b(20)

C     LOCAL SCALARS
      integer k

      u  = a0
      u1 = 0.0d0
      u2 = 0.0d0

      do k = 1,m
          u  = u  + a(k) * cos(k*x)         + b(k) * sin(k*x)
          u1 = u1 - a(k) * k * sin(k*x)     + b(k) * k * cos(k*x)
          u2 = u2 - a(k) * k * k * cos(k*x) - b(k) * k * k * sin(k*x)
      end do

      fi = u2 + u1 ** 2 + u

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
      double precision uint(1000)

C     LOCAL SCALAR
      integer i

C     COMMON BLOCKS
      common /objective/ uint

      flag = 0

      f = 0.0d0
      do i = 1,n
          f = f + ( x(i) - uint(i) ) ** 2
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
      double precision uint(1000)

C     LOCAL SCALAR
      integer i

C     COMMON BLOCKS
      common /objective/ uint

      flag = 0

      do i = 1,n 
          g(i) = 2.0d0 * ( x(i) - uint(i) ) 
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
      double precision h

C     COMMON ARRAYS
      double precision fon(1000)

C     COMMON BLOCKS
      common /constraints/ fon,h
      
      flag = 0

      c = ( x(ind+2) - 2.0d0 * x(ind+1) + x(ind) ) / h ** 2 + 
     +    ( ( x(ind+2) - x(ind) ) / ( 2.0d0 * h ) ) ** 2 + 
     +    x(ind+1) - fon(ind)

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

C     COMMON SCALARS
      double precision h

C     COMMON ARRAYS
      double precision fon(1000)

C     COMMON BLOCKS
      common /constraints/ fon,h
      
C     Sparse gradient vector of the ind-th constraint

      flag = 0

      jcnnz = 3

      jcvar(1) = ind
      jcval(1) = 1.0d0 / h ** 2 - 
     +           2.0d0 * ( x(ind+2) - x(ind) ) / ( 2.0d0 * h ) ** 2 

      jcvar(2) = ind + 1
      jcval(2) = - 2.0d0 / h ** 2 + 1.0d0

      jcvar(3) = ind + 2

      jcval(3) = 1.0d0 / h ** 2 + 
     +           2.0d0 * ( x(ind+2) - x(ind) ) / ( 2.0d0 * h ) ** 2 

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

C     COMMON ARRAYS
      double precision usol(1000)

C     LOCAL SCALARS
      integer i
      double precision zmax

C     COMMON BLOCKS
      common /solution/ usol
      
      zmax = 0.0d0
      do i = 1,n
          zmax = max( zmax, abs( x(i) - usol(i) ) )
      end do

      write(*,*)
      write(*,*) 'Comparison between analytic and computed solution: '

      write(*,*) 'Analytic Computed Difference'
      do i = 1,n
          write(*,*) usol(i),x(i),usol(i) - x(i)
      end do

      write(*, *) 'Maximum error: ',zmax        

      end
