C     =================================================================
C     File: vi8EA.f
C     =================================================================

C     =================================================================
C     Module: Subroutines that define the problem
C     =================================================================

C     Last update of any of the component of this module:

C     March 3, 2008.

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

C     Elastic beam problem
C     --------------------

C     This problem consists in finding the thickness distribution
C     minimizing the compliance of an elastic beam [1]. The
C     optimization model is based on the formulation in [2,3]. Finite
C     elements method is used to obtain a suitable discretized problem.
C     When the beam with endpoints fixed-simply supported is discretized
C     by using 8 elements the optimization problem has 23 variables and
C     30 constraints.
C
C     [1] Haslinger, J. and Mäkinen, R.A.E. (2003) Introduction to Shape
C     Optimization. Theory, Approximation, and Computation. SIAM,
C     Philadelphia.
C
C     [2] Maciel, M.C., Pilotta, E.A. and Sottosanto, G.N. (2004)
C     Aplicación de un modelo de optimización al diseño de una viga
C     elástica, en G. Buscaglia, E. Dari, O. Zamonsky (Eds.), Mecánica
C     Computacional, Vol. XXIII, pp. 2831-2844.
C
C     [3] Maciel, M.C., Pilotta, E.A. and Sottosanto, G.N. (2007)
C     Thickness Optimization of an Elastic Beam, To appear in
C     MAT-Serie A.

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

C     LOCAL SCALARS
      integer i

C     Number of variables

      n = 23

C     Initial point

      do i = 1,8
          x(i) = 5.0d-02
      end do

      x(9)  = -0.06515503841146d+02
      x(10) = -0.92244479687500d+02
      x(11) = -0.20385745312500d+02
      x(12) = -1.20279968750000d+02
      x(13) = -0.34561163085938d+02
      x(14) = -0.99731467187500d+02
      x(15) = -0.43945320833333d+02
      x(16) = -0.46223975000000d+02
      x(17) = -0.45394907226563d+02
      x(18) =  0.24617507812500d+02
      x(19) = -0.37719735937500d+02
      x(20) =  0.97167981250000d+02
      x(21) = -0.21682745638021d+02
      x(22) =  1.55802445312500d+02
      x(23) =  1.84895900000001d+02

C     Lower and upper bounds

      do i = 1,8
          l(i) = 0.01d0
          u(i) = 0.10d0
      end do

      do i = 9,21,2
          l(i) = - 1.0d+09
          u(i) =   0.0d0
      end do

      do i=10,22,2
          l(i) = - 1.0d+09
          u(i) =   1.0d+09
      end do

      l(23) = - 1.0d+09
      u(23) =   1.0d+09

C     Number of constraints (equalities plus inequalities)

      m = 30

C     Lagrange multipliers approximation.

      do i = 1,m
          lambda(i) = 0.0d0
      end do

C     For each constraint i, set equatn(i) = .true. if it is an equality
C     constraint of the form c_i(x) = 0, and set equatn(i) = .false. if
C     it is an inequality constraint of the form c_i(x) <= 0.

      do i = 1,16
          equatn(i) = .true.
      end do

      do i = 17,30
          equatn(i) = .false.
      end do

C     For each constraint i, set linear(i) = .true. if it is a linear
C     constraint, otherwise set linear(i) = .false.

      do i = 1,15
          linear(i) = .false.
      end do

      do i = 16,30
          linear(i) = .true.
      end do

C     Indicate which subroutines did you code.

      coded( 1) = .true.  ! evalf
      coded( 2) = .true.  ! evalg
      coded( 3) = .true.  ! evalh
      coded( 4) = .true.  ! evalc
      coded( 5) = .true.  ! evaljac
      coded( 6) = .false. ! evalhc
      coded( 7) = .false. ! evalfc
      coded( 8) = .false. ! evalgjac
      coded( 9) = .false. ! evalhl
      coded(10) = .false. ! evalhlp

C     Set checkder = TRUE if you code some derivatives and you would
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

      flag = 0

      f = - 0.125d0 *
     +      ( x(9) + x(11) + x(13) + x(15) + x(17) + x(19) + x(21) )

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalg(n,x,g,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,n

C     ARRAY ARGUMENTS
      double precision g(n),x(n)

C     LOCAL SCALARS
      integer i

      flag = 0

      do i = 1,n
          g(i) = 0.0d0
      end do

      do i = 9,21,2
          g(i) = - 0.125d0
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

C     LOCAL SCALARS
      integer j

      flag = 0

      if ( ind .eq. 1 ) then
          c =  ( 6144.0d0 * x(1) ** 3 + 6144.0d0 * x(2) ** 3 ) * x(9)
     +      + ( -384.0d0 * x(1) ** 3 + 384.0d0 * x(2) ** 3 ) * x(10)
     +      - 6144.0d0 * x(2) ** 3 * x(11) + 384.0d0 * x(2) ** 3
     +      * x(12) + 0.125d0

      else if ( ind .eq. 2 ) then
          c = ( - 384.0d0 * x(1) ** 3 + 384.0d0 * x(2) ** 3 ) * x(9)
     +      + ( 32.0d0 * x(1) ** 3 + 32.0d0 * x(2) ** 3 ) * x(10)
     +      - 384.0d0 * x(2) ** 3 * x(11) + x(2) ** 3 * 16.0d0 * x(12)

      else if ( ind .eq. 3 ) then
          c = - 6144.0d0 * x(2) ** 3 * x(9) - 384.0d0 * x(2) ** 3
     +      * x(10) + ( 6144.0d0 * x(2) ** 3 + 6144.0d0 * x(3) ** 3 )
     +      * x(11) + ( - 384.0d0 * x(2) ** 3 + 384.0d0 * x(3) ** 3 )
     +      * x(12) - 6144.0d0 * x(3) ** 3 * x(13) + 384.0d0 * x(3) ** 3
     +      * x(14) + 0.125d0

      else if ( ind .eq. 4 ) then
           c = 384.0d0 * x(2) ** 3 * x(9) + 16.0d0 * x(2) ** 3 * x(10)
     +       + ( - 384.0d0 * x(2) ** 3 + 384.0d0 * x(3) ** 3 ) * x(11)
     +       + ( 32.0d0 * x(2) ** 3 + 32.0d0 * x(3) ** 3 ) * x(12)
     +       - 384.0d0 * x(3) ** 3 * x(13) + 16.0d0 * x(3) ** 3 * x(14)

      else if ( ind .eq. 5 ) then
          c = - 6144.0d0 * x(3) ** 3 * x(11) - 384.0d0 * x(3) ** 3
     +      * x(12) + ( 6144.0d0 * x(3) ** 3 + 6144.0d0 * x(4) ** 3 )
     +      * x(13) + ( -384.0d0 * x(3) ** 3 + 384.0d0 * x(4) ** 3 )
     +      * x(14) - 6144.0d0 * x(4) ** 3 * x(15) + 384.0d0 * x(4) ** 3
     +      * x(16) + 0.125d0

      else if ( ind .eq. 6 ) then
          c = 384.0d0 * x(3) ** 3 * x(11) + 16.0d0 * x(3) **3 * x(12)
     +      + ( - 384.0d0 * x(3) ** 3 + 384.0d0 * x(4) ** 3 ) * x(13)
     +      + ( 32.0d0 * x(3) ** 3 + 32.0d0 * x(4) ** 3 ) * x(14)
     +      - 384.0d0 * x(4) ** 3 * x(15) + 16.0d0 * x(4) ** 3 * x(16)

      else if ( ind .eq. 7 ) then
          c = - 6144.0d0 * x(4) ** 3 * x(13) - 384.0d0 * x(4) ** 3
     +      * x(14) + ( 6144.0d0 * x(4) ** 3 + 6144.0d0 * x(5) ** 3 )
     +      * x(15) + ( - 384.0d0 * x(4) ** 3 + 384.0d0 * x(5) ** 3 )
     +      * x(16) - 6144.0d0 * x(5) ** 3 * x(17) + 384.0d0 * x(5) ** 3
     +      * x(18) + 0.125d0

      else if ( ind .eq. 8 ) then
          c = 384.0d0 * x(4) ** 3 * x(13) + 16.0d0 * x(4) ** 3 * x(14)
     +      + ( - 384.0d0 * x(4) ** 3 + 384.0d0 * x(5) ** 3 ) * x(15)
     +      + ( 32.0d0 * x(4) ** 3 + 32.0d0 * x(5) ** 3 ) * x(16)
     +      - 384.0d0 * x(5) ** 3 * x(17) + 16.0d0 * x(5) ** 3 * x(18)

      else if ( ind .eq. 9 ) then
          c = - 6144.0d0 * x(5) ** 3 * x(15) - 384.0d0 * x(5) ** 3
     +      * x(16) + ( 6144.0d0 * x(5) ** 3 + 6144.0d0 * x(6) ** 3 )
     +      * x(17) + ( -384.0d0 * x(5) ** 3 + 384.0d0 * x(6) ** 3 )
     +      * x(18) - 6144.0d0 * x(6) ** 3 * x(19) + 384.0d0 * x(6) ** 3
     +      * x(20) + 0.125d0

      else if ( ind .eq. 10 ) then
          c = 384.0d0 * x(5) ** 3 * x(15) + 16.0d0 * x(5) ** 3 * x(16)
     +      + ( - 384.0d0 * x(5) ** 3 + 384.0d0 * x(6) ** 3 ) * x(17)
     +      + ( 32.0d0 * x(5) ** 3 + 32.0d0 * x(6) ** 3 ) * x(18)
     +      - 384.0d0 * x(6) ** 3 * x(19) + 16.0d0 * x(6) ** 3 * x(20)

      else if ( ind .eq. 11 ) then
          c = - 6144.0d0 * x(6) ** 3 * x(17) - 384.0d0 * x(6) ** 3
     +      * x(18) + ( 6144.0d0 * x(6) ** 3 + 6144.0d0 * x(7) ** 3 )
     +      * x(19) + ( - 384.0d0 * x(6) ** 3 + 384.0d0 * x(7) ** 3 )
     +      * x(20) - 6144.0d0 * x(7) ** 3 * x(21) + 384.0d0 * x(7) ** 3
     +      * x(22) + 0.125d0

      else if ( ind .eq. 12 ) then
          c = 384.0d0 * x(6) ** 3 * x(17) + 16.0d0 * x(6) ** 3 * x(18)
     +      + ( -384.0d0 * x(6) ** 3 + 384.0d0 * x(7) ** 3 ) * x(19)
     +      + ( 32.0d0 * x(6) ** 3 + 32.0d0 * x(7) ** 3 ) * x(20)
     +      - 384.0d0 * x(7) ** 3 * x(21) + 16.0d0 * x(7) ** 3 * x(22)

      else if ( ind .eq. 13 ) then
          c = - 6144.0d0 * x(7) ** 3 * x(19) - 384.0d0 * x(7) ** 3
     +      * x(20) + ( 6144.0d0 * x(7) ** 3 + 6144.0d0 * x(8) ** 3 )
     +      * x(21) + ( - 384.0d0 * x(7) ** 3 + 384.0d0 * x(8) ** 3 )
     +      * x(22) + 384.0d0 * x(8) ** 3 * x(23) + 0.125d0

      else if ( ind .eq. 14 ) then
          c = 384.0d0 * x(7) ** 3 * x(19) + 16.0d0 * x(7) ** 3 * x(20)
     +      + ( - 384.0d0 * x(7) ** 3 + 384.0d0 * x(8) ** 3 ) * x(21)
     +      + ( 32.0d0 * x(7) ** 3 + 32.0d0 * x(8) ** 3 ) * x(22)
     +      + 16.0d0 * x(8) ** 3 * x(23)

      else if ( ind .eq. 15 ) then
          c = 384.0d0 * x(8) ** 3 * x(21) + 16.0d0 * x(8) ** 3 * x(22)
     +      + 32.0d0 * x(8) ** 3 * x(23) - 0.0104167d0

      else if ( ind .eq. 16 ) then
          c = - 5.0d-02
          do j = 1,8
              c = c + 0.125d0 * x(j)
          end do

      else if ( 17 .le. ind .and. ind .le. 23 ) then
          c =   x(ind-15) - x(ind-16) - 0.05d0

      else if ( 24 .le. ind .and. ind .le. 30 ) then
          c = - x(ind-22) + x(ind-23) - 0.05d0
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

C     SCALAR ARGUMENTS
      integer j,k

      flag = 0

      if ( ind .eq. 1 ) then

          jcnnz = 6

          jcvar(1) = 1
          jcval(1) = x(1) ** 2 * ( 18432.0d0 * x(9) - 1152.0d0 * x(10) )

          jcvar(2) = 2
          jcval(2) = x(2) ** 2 * ( 18432.0d0 * x(9) + 1152.0d0 * x(10)
     +             - 18432.0d0 * x(11) + 1152.0d0 * x(12) )

          jcvar(3) = 9
          jcval(3) = 6144.0d0 * ( x(1) ** 3 + x(2) ** 3 )

          jcvar(4) = 10
          jcval(4) = 384.0d0 * ( - x(1) ** 3 + x(2) ** 3 )

          jcvar(5) = 11
          jcval(5) = - 6144.0d0 * x(2) ** 3

          jcvar(6) = 12
          jcval(6) = 384.0d0 * x(2) ** 3

      else if ( ind .eq. 2 ) then

          jcnnz = 6

          jcvar(1) = 1
          jcval(1) = x(1) ** 2 * ( - 1152.0d0 * x(9) + 96.0d0 * x(10) )

          jcvar(2) = 2
          jcval(2) = x(2) ** 2 * ( 1152.0d0 * x(9) + 96.0d0 * x(10)
     +             - 1152.0d0 * x(11) + 48.0d0 * x(12) )

          jcvar(3) = 9
          jcval(3) = 384.0d0 * ( - x(1) ** 3 + x(2) ** 3 )

          jcvar(4) = 10
          jcval(4) = 32.0d0 * ( x(1) ** 3 + x(2) ** 3 )

          jcvar(5) = 11
          jcval(5) = - 384.0d0 * x(2) ** 3

          jcvar(6) = 12
          jcval(6) = 16.0d0 * x(2) ** 3

      else if ( ind .eq. 3 .or. ind .eq. 5 .or. ind .eq. 7 .or.
     +          ind .eq. 9 .or. ind .eq. 11 ) then

          j = ( ind + 1 ) / 2
          k = ind + 6

          jcnnz = 8

          jcvar(1) = j
          jcval(1) = x(j) ** 2 * ( - 18432.0d0 * x(k) - 1152.0d0
     +             * x(k+1) + 18432.0d0 * x(k+2) -1152.0d0 * x(k+3) )

          jcvar(2) = j + 1
          jcval(2) = x(j+1) ** 2 * ( 18432.0d0 * x(k+2) + 1152.0d0
     +             * x(k+3) - 18432.0d0 * x(k+4) + 1152.0d0 * x(k+5) )

          jcvar(3) = k
          jcval(3) = - 6144.0d0 * x(j) ** 3

          jcvar(4) = k + 1
          jcval(4) = - 384.0d0 * x(j) ** 3

          jcvar(5) = k + 2
          jcval(5) = 6144.0d0 * ( x(j) ** 3 + x(j+1) ** 3 )

          jcvar(6) = k + 3
          jcval(6) = 384.0d0 * ( - x(j) ** 3 + x(j+1) ** 3 )

          jcvar(7) = k + 4
          jcval(7) = - 6144.0d0 * x(j+1) ** 3

          jcvar(8) = k + 5
          jcval(8) = 384.0d0 * x(j+1) ** 3

      else if ( ind .eq. 4 .or. ind .eq. 6 .or. ind .eq. 8 .or.
     +          ind .eq. 10 .or. ind .eq. 12 ) then

          j = ind / 2
          k = ind + 5

          jcnnz = 8

          jcvar(1) = j
          jcval(1) = x(j) ** 2 * ( 1152.0d0 * x(k) + 48.0d0 * x(k+1)
     +             - 1152.0d0 * x(k+2) + 96.0d0 * x(k+3) )

          jcvar(2) = j+1
          jcval(2) = x(j+1) ** 2 * ( 1152.0d0 * x(k+2) + 96.0d0 * x(k+3)
     +             - 1152.0d0 * x(k+4) + 48.0d0 * x(k+5) )

          jcvar(3) = k
          jcval(3) = 384.0d0 * x(j) ** 3

          jcvar(4) = k + 1
          jcval(4) = 16.0d0 * x(j) ** 3

          jcvar(5) = k + 2
          jcval(5) = 384.0d0 * ( - x(j) ** 3 + x(j+1) ** 3 )

          jcvar(6) = k + 3
          jcval(6) = 32.0d0 * ( x(j) ** 3 + x(j+1) ** 3 )

          jcvar(7) = k + 4
          jcval(7) = - 384.0d0 * x(j+1) ** 3

          jcvar(8) = k + 5
          jcval(8) = 16.0d0 * x(j+1) ** 3

      else if ( ind .eq. 13 ) then

          jcnnz = 7

          jcvar(1) = 7
          jcval(1) = x(7) ** 2 * ( - 18432.0d0 * x(19) - 1152.0d0
     +             * x(20) + 18432.0d0 * x(21) - 1152.0d0 * x(22) )

          jcvar(2) = 8
          jcval(2) = x(8) ** 2 * ( 18432.0d0 * x(21) + 1152.0d0
     +             * ( x(22) + x(23) ) )

          jcvar(3) = 19
          jcval(3) = - 6144.0d0 * x(7) ** 3

          jcvar(4) = 20
          jcval(4) = - 384.0d0 * x(7) ** 3

          jcvar(5) = 21
          jcval(5) = 6144.0d0 * ( x(7) ** 3 + x(8) ** 3 )

          jcvar(6) = 22
          jcval(6) = 384.0d0 * ( - x(7) ** 3 + x(8) ** 3 )

          jcvar(7) = 23
          jcval(7) = 384.0d0 * x(8) ** 3

      else if ( ind .eq. 14 ) then

          jcnnz = 7

          jcvar(1) = 7
          jcval(1) = x(7) ** 2 * ( 1152.0d0 * x(19) + 48.0d0 * x(20)
     +             - 1152.0d0 * x(21) + 96.0d0 * x(22) )

          jcvar(2) = 8
          jcval(2) = x(8) ** 2 * ( 1152.0d0 * x(21) + 96.0d0 * x(22)
     +             + 48.0d0 * x(23) )

          jcvar(3) = 19
          jcval(3) = 384.0d0 * x(7) ** 3

          jcvar(4) = 20
          jcval(4) = 16.0d0 * x(7) ** 3

          jcvar(5) = 21
          jcval(5) = 384.0d0 * ( - x(7) ** 3 + x(8) ** 3 )

          jcvar(6) = 22
          jcval(6) = 32.0d0 * ( x(7) ** 3 + x(8) ** 3 )

          jcvar(7) = 23
          jcval(7) = 16.0d0 * x(8) ** 3

      else if ( ind .eq. 15 ) then

          jcnnz = 4

          jcvar(1) = 8
          jcval(1) = x(8) ** 2 * ( 1152.0d0 * x(21) + 48.0d0 * x(22)
     +             + 96.0d0 * x(23) )

          jcvar(2) = 21
          jcval(2) = 384.0d0 * x(8) ** 3

          jcvar(3) = 22
          jcval(3) = 16.0d0 * x(8) ** 3

          jcvar(4) = 23
          jcval(4) = 32.0d0 * x(8) ** 3

      else if ( ind .eq. 16 ) then

          jcnnz = 8

          do j = 1,8
              jcvar(j) = j
              jcval(j) = 0.125d0
          end do

      else if ( 17 .le. ind .and. ind .le. 23 ) then

          jcnnz = 2

          jcvar(1) = ind - 15
          jcval(1) = 1.0d0

          jcvar(2) = ind - 16
          jcval(2) = - 1.0d0

      else if ( 24 .le. ind .and. ind .le. 30 ) then

          jcnnz = 2

          jcvar(1) = ind - 22
          jcval(1) = - 1.0d0

          jcvar(2) = ind - 23
          jcval(2) = 1.0d0
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

      subroutine endp(n,x,l,u,m,lambda,rho,equatn,linear)

      implicit none

C     SCALAR ARGUMENTS
      integer m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision l(n),lambda(m),rho(m),u(n),x(n)

      end
