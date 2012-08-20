C     *****************************************************************
C     *****************************************************************

      subroutine checkd(n,l,u,m,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n

C     ARRAY ARGUMENTS
      double precision l(n),u(n)

C     This subrotutine checks the user supplied first and second
C     derivatives subroutines (evalg, evalh, evaljac and evalhc) for
C     computing the objective function gradient and Hessian and the
C     constraints gradients and Hessians, respectively.

#include "fintrf.h"

#include "dim.par"
#include "outtyp.com"
#include "algparam.com"

C     LOCAL SCALARS
      character answer
      integer i,j,stat
      double precision drand,seed,smalll,smallu

C     LOCAL ARRAYS
      mwPointer plhs(1),prhs(2)
      double precision x(nmax)
      character*159 warn

C     EXTERNAL FUNCTIONS
      external drand
      integer mxGetString
      mwPointer mxCreateCharMatrixFromStrings
      mwSize integerToMwSize

      warn(1:45)    = '\n Warning! If you are in graphical mode, chec'
      warn(46:90)   = "kder output does not appear.\n Please type 'a'"
      warn(91:135)  = ' to finish running ALGENCAN and run Matlab wi'
      warn(136:159) = 'th -nodesktop\n option.\n\n'

      prhs(1) = mxCreateCharMatrixFromStrings(integerToMwSize(1), warn)

      call mexCallMATLAB(0,plhs,1,prhs,'fprintf')

      call mxDestroyArray(prhs(1))

      prhs(1) = mxCreateCharMatrixFromStrings(integerToMwSize(1),'')
      prhs(2) = mxCreateCharMatrixFromStrings(integerToMwSize(1),'s')

      call vunsetp()

C     SET A RANDOM POINT

      seed = 123456.0d0
      do i = 1,n
          smalll = max( l(i), - 10.0d0 )
          smallu = min( u(i),   10.0d0 )
          if ( .not. smalll .lt. smallu ) then
              smalll = l(i)
              smallu = u(i)
          end if
          x(i) = smalll + ( smallu - smalll ) * drand(seed)
      end do

      if ( iprintout .ge. 1 ) then
          write(* ,100)
          write(10,100)

          do i = 1,n
              write(* ,110) i,x(i)
              write(10,110) i,x(i)
          end do
      end if

C     CHECK OBJECTIVE FUNCTION GRADIENT

      if ( .not. gcoded ) then
          write(* ,160) 'evalg'
          write(10,160) 'evalg'

          go to 1000
      end if

      write(* ,120)
      write(10,120)

      call mexCallMATLAB(1,plhs,2,prhs,'input')
      stat = mxGetString(plhs(1),answer,integerToMwSize(1))

      if ( answer .eq. 'A' .or. answer .eq. 'a' ) then
          return

      else if ( answer .eq. 'Y' .or. answer .eq. 'y' ) then
          call checkg(n,x,inform)
          if ( inform .lt. 0 ) return
      end if

C     CHECK CONSTRAINTS GRADIENTS

 1000 continue

      if ( .not. jaccoded ) then
          write(* ,160) 'evaljac'
          write(10,160) 'evaljac'

          go to 1020
      end if

      j = 1

 1010 if ( j .le. m ) then

          write(* ,130) j
          write(10,130) j

          call mexCallMATLAB(1,plhs,2,prhs,'input')
          stat = mxGetString(plhs(1),answer,integerToMwSize(1))

          if ( answer .eq. 'A' .or. answer .eq. 'a' ) then
              return

          else if ( answer .eq. 'S' .or. answer .eq. 's' ) then
              go to 1020

          else if ( answer .eq. 'Y' .or. answer .eq. 'y' ) then
              call checkjac(n,x,j,inform)
              if ( inform .lt. 0 ) return
          end if

          j = j + 1

          go to 1010

      end if

C     CHECK HESSIAN OF THE OBJECTIVE FUNCTION

 1020 continue

      if ( .not. hcoded ) then
          write(* ,160) 'evalh'
          write(10,160) 'evalh'

          go to 1030
      end if

      write(* ,140)
      write(10,140)

      call mexCallMATLAB(1,plhs,2,prhs,'input')
      stat = mxGetString(plhs(1),answer,integerToMwSize(1))

      if ( answer .eq. 'A' .or. answer .eq. 'a' ) then
          return

      else if ( answer .eq. 'Y' .or. answer .eq. 'y' ) then
          call checkh(n,x,inform)
          if ( inform .lt. 0 ) return
      end if

C     CHECK HESSIANS OF THE CONSTRAINTS

 1030  continue

      if ( .not. hccoded ) then
          write(* ,160) 'evalhc'
          write(10,160) 'evalhc'

          go to 1050
      end if

      j = 1

 1040 if ( j .le. m ) then

          write(* ,150) j
          write(10,150) j

          call mexCallMATLAB(1,plhs,2,prhs,'input')
          stat = mxGetString(plhs(1),answer,integerToMwSize(1))

          if ( answer .eq. 'A' .or. answer .eq. 'a' ) then
              return

          else if ( answer .eq. 'S' .or. answer .eq. 's' ) then
              return

          else if ( answer .eq. 'Y' .or. answer .eq. 'y' ) then
              call checkhc(n,x,j,inform)
              if ( inform .lt. 0 ) return
          end if

          j = j + 1

          go to 1040

      end if

C     CHECK HESSIAN OF THE LAGRANGIAN

 1050 continue

      if ( .not. hlcoded ) then
          write(* ,160) 'evalhl'
          write(10,160) 'evalhl'

          go to 1060
      end if

      write(* ,*) 'Test of evalhl not implemented yet!'
      write(10,*) 'Test of evalhl not implemented yet!'

C     CHECK HESSIAN OF THE LAGRANGIAN TIMES A VECTOR

 1060 continue

      if ( .not. hlpcoded ) then
          write(* ,160) 'evalhlp'
          write(10,160) 'evalhlp'

          return
      end if

      write(* ,*) 'Test of evalhlp not implemented yet!'
      write(10,*) 'Test of evalhlp not implemented yet!'


      call mxDestroyArray(prhs(1))
      call mxDestroyArray(prhs(2))

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'Derivatives will be tested at the random point: ')
 110  format(  1X,'x(',I6,') = ',1P,D15.8)
 120  format(/,1X,'Check the gradient of the objective function?',
     +       /,1X,'Type Y(es), N(o) or A(bort checking): ')
 130  format(/,1X,'Check the gradient of constraint ',I5,'?',
     +       /,1X,'Type Y(es), N(o), A(bort checking) or ',
     +            'S(kip constraints gradients): ')
 140  format(/,1X,'Check the Hessian matrix of the objective function?',
     +       /,1X,'Type Y(es), N(o) or A(bort checking): ')
 150  format(/,1X,'Check the Hessian of constraint ',I5,'?',
     +       /,1X,'Type Y(es), N(o), A(bort checking) or ',
     +            'S(kip constraints gradients): ')
 160  format(/,1X,'Skipping test of uncoded ',A7,' subroutine.')

      end

C     *****************************************************************
C     *****************************************************************

      subroutine checkg(n,x,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n

C     ARRAY ARGUMENTS
      double precision x(n)

C     This subrotutine checks the user supplied subroutine evalg for
C     computing the gradient of the objective function using central
C     finite differences with two different discretization steps.

#include "dim.par"
#include "machconst.com"
#include "outtyp.com"

C     LOCAL SCALARS
      integer i
      double precision fminus,fplus,gdiff1,gdiff2,maxerr,step1,step2,tmp

C     LOCAL ARRAYS
      double precision g(nmax)

      call vevalg(n,x,g,inform)
      if ( inform .lt. 0 ) return

      if ( iprintout .ge. 1 ) then
          write(* ,100)
          write(10,100)
      end if

      maxerr = 0.0d0

      do i = 1,n
          tmp  = x(i)

          step1 = macheps13 * max( abs( tmp ), 1.0d0 )

          x(i) = tmp + step1
          call vevalf(n,x,fplus,inform)
          if ( inform .lt. 0 ) return

          x(i) = tmp - step1
          call vevalf(n,x,fminus,inform)
          if ( inform .lt. 0 ) return

          gdiff1 = ( fplus - fminus ) / ( 2.0d0 * step1 )

          step2 = macheps13 * max( abs( tmp ), 1.0d-03 )

          x(i) = tmp + step2
          call vevalf(n,x,fplus,inform)
          if ( inform .lt. 0 ) return

          x(i) = tmp - step2
          call vevalf(n,x,fminus,inform)
          if ( inform .lt. 0 ) return

          x(i) = tmp

          gdiff2 = ( fplus - fminus ) / ( 2.0d0 * step2 )

          tmp = min( abs( g(i) - gdiff1 ), abs( g(i) - gdiff2 ) )

          if ( iprintout .ge. 1 ) then
              write(* ,110) i,g(i),gdiff1,gdiff2,tmp
              write(10,110) i,g(i),gdiff1,gdiff2,tmp
          end if

          maxerr = max( maxerr, tmp )

      end do

      if ( iprintout .ge. 1 ) then
          write(* ,120) maxerr
          write(10,120) maxerr
      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'Gradient vector of the objective function.',
     +       /,1X,'Index',13X,'evalg',2X,'Central diff (two different ',
     +            'steps)',4X,'Absolute error')
 110  format(  1X,I5,4(3X,1P,D15.8))
 120  format(  1X,'Maximum absolute error = ',1P,D15.8)

      end

C     *****************************************************************
C     *****************************************************************

      subroutine checkh(n,x,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n

C     ARRAY ARGUMENTS
      double precision x(n)

C     This subrotutine checks the user supplied subroutine evalh for
C     computing the Hessian of the objective function using central
C     finite differences with two different discretization steps.

#include "dim.par"
#include "machconst.com"
#include "outtyp.com"

C     LOCAL SCALARS
      logical nullcol
      integer i,j,hnnz
      double precision elem,hdiff1,hdiff2,maxerr,step1,step2,tmp

C     LOCAL ARRAYS
      integer hlin(nsmax**2),hcol(nsmax**2)
      double precision g(nsmax),gplus1(nsmax),gplus2(nsmax),
     +        H(nsmax,nsmax),hval(nsmax**2),maxcoe(nsmax)

C     Check viability of the test

      if ( n .gt. nsmax ) then
          if ( iprintout .ge. 1 ) then
              write(*, 100) nsmax,nsmax
              write(10,100) nsmax,nsmax
          end if

          return
      end if

C     Compute the gradient of the objective function at x

      call vevalg(n,x,g,inform)
      if ( inform .lt. 0 ) return

C     Compute the Hessian of the objective function at x and save in a
C     dense matrix

      call vevalh(n,x,hlin,hcol,hval,hnnz,inform)
      if ( inform .lt. 0 ) return

      do j = 1,n
          do i = 1,n
              H(i,j) = 0.0d0
          end do
      end do

      do i = 1,hnnz
          H(hlin(i),hcol(i)) = H(hlin(i),hcol(i)) + hval(i)
      end do

C     Test column by column

      if ( iprintout .ge. 1 ) then
          write(* ,200)
          write(10,200)
      end if

      maxerr = 0.0d0

      do j = 1,n

          tmp  = x(j)

          step1 = macheps12 * max( abs( tmp ), 1.0d0 )

          x(j) = tmp + step1
          call vevalg(n,x,gplus1,inform)
          if ( inform .lt. 0 ) return

          step2 = macheps12 * max( abs( tmp ), 1.0d-03 )

          x(j) = tmp + step2
          call vevalg(n,x,gplus2,inform)
          if ( inform .lt. 0 ) return

          x(j) = tmp

          if ( iprintout .ge. 1 ) then
              write(* ,210) j
              write(10,210) j
          end if

          maxcoe(j) = 0.0d0

          nullcol = .true.

          do i = 1,n
              if ( i .ge. j ) then
                  elem = H(i,j)
              else
                  elem = H(j,i)
              end if
              hdiff1 = ( gplus1(i) - g(i) ) / step1
              hdiff2 = ( gplus2(i) - g(i) ) / step2
              tmp = min( abs( elem - hdiff1 ), abs( elem - hdiff2 ) )
              if ( elem   .ne. 0.0d0 .or.
     +             hdiff1 .ne. 0.0d0 .or.
     +             hdiff2 .ne. 0.0d0 ) then
                  if ( nullcol ) then
                      nullcol = .false.
                      if ( iprintout .ge. 1 ) then
                          write(* ,220)
                          write(10,220)
                      end if
                  end if
                  if ( iprintout .ge. 1 ) then
                      write(* ,230) i,elem,hdiff1,hdiff2,tmp
                      write(10,230) i,elem,hdiff1,hdiff2,tmp
                  end if
              end if
              maxcoe(j) = max( maxcoe(j), tmp )
          end do

          maxerr = max( maxerr, maxcoe(j) )

          if ( iprintout .ge. 1 ) then
              if ( nullcol ) then
                  write(* ,240)
                  write(10,240)
              else
                  write(* ,250) maxcoe(j)
                  write(10,250) maxcoe(j)
              end if
          end if

      end do

      if ( iprintout .ge. 1 ) then
          write(* ,*)
          write(10,*)

          do j = 1,n
              write(* ,260) j,maxcoe(j)
              write(10,260) j,maxcoe(j)
          end do

          write(* ,270) maxerr
          write(10,270) maxerr
      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'Subroutine CHECKH uses dense matrices up to ',
     +            'dimension ',I6,' times ',I6,'. The Hessian ',
     +            'checking will be skipped.')

 200  format(/,1X,'Hessian matrix of the objective function column by ',
     +            'column.')
 210  format(/,1X,'Column:  ',I6)
 220  format(/,1X,'Index',13X,'evalh',3X,'Incr. Quoc. (two different ',
     +            'steps)',4X,'Absolute error')
 230  format(  1X,I5,4(3X,1P,D15.8))
 240  format(  1X,'All the elements of this column are null.')
 250  format(  1X,'Maximum absolute error = ',1P,D15.8)
 260  format(  1X,'Column ',I6,' Maximum absolute error = ',1P,D15.8)
 270  format(/,1X,'Overall maximum absolute error = ',1P,D15.8)

      end

C     *****************************************************************
C     *****************************************************************

      subroutine checkjac(n,x,ind,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer ind,inform,n

C     ARRAY ARGUMENTS
      double precision x(n)

C     This subrotutine checks the user supplied subroutine evaljac for
C     computing the gradients of the constraints using central finite
C     differences with two different discretization steps.

#include "dim.par"
#include "machconst.com"
#include "outtyp.com"

C     LOCAL SCALARS
      logical nullcol
      integer i,jcnnz
      double precision cminus,cplus,jacdiff1,jacdiff2,maxerr,step1,
     +        step2,tmp

C     LOCAL ARRAYS
      integer jcvar(nmax)
      double precision g(nmax),jcval(nmax)

C     COMPUTE THE GRADIENT OF THE CONSTRAINT AND SAVE IT INTO A DENSE
C     VECTOR

      call vevaljac(n,x,ind,jcvar,jcval,jcnnz,inform)
      if ( inform .lt. 0 ) return

      do i = 1,n
          g(i) = 0.0d0
      end do

      do i = 1,jcnnz
          g(jcvar(i)) = g(jcvar(i)) + jcval(i)
      end do

C     COMPARE WITH CENTRAL FINITE DIFFERENCES

      if ( iprintout .ge. 1 ) then
          write(* ,100) ind
          write(10,100) ind
      end if

      maxerr = 0.0d0

      nullcol = .true.

      do i = 1,n
          tmp  = x(i)

          step1 = macheps13 * max( abs( tmp ), 1.0d0 )

          x(i) = tmp + step1
          call vevalc(n,x,ind,cplus,inform)
          if ( inform .lt. 0 ) return

          x(i) = tmp - step1
          call vevalc(n,x,ind,cminus,inform)
          if ( inform .lt. 0 ) return

          jacdiff1 = ( cplus - cminus ) / ( 2.0d0 * step1 )

          step2 = macheps13 * max( abs( tmp ), 1.0d-03 )

          x(i) = tmp + step2
          call vevalc(n,x,ind,cplus,inform)
          if ( inform .lt. 0 ) return

          x(i) = tmp - step2
          call vevalc(n,x,ind,cminus,inform)
          if ( inform .lt. 0 ) return

          x(i) = tmp

          jacdiff2 = ( cplus - cminus ) / ( 2.0d0 * step2 )

          tmp = min( abs( g(i) - jacdiff1 ), abs( g(i) - jacdiff2 ) )

          if ( g(i)     .ne. 0.0d0 .or.
     +         jacdiff1 .ne. 0.0d0 .or.
     +         jacdiff2 .ne. 0.0d0 ) then
              if ( nullcol ) then
                  nullcol = .false.
                  if ( iprintout .ge. 1 ) then
                      write(* ,110)
                      write(10,110)
                  end if
              end if
              if ( iprintout .ge. 1 ) then
                  write(* ,120) i,g(i),jacdiff1,jacdiff2,tmp
                  write(10,120) i,g(i),jacdiff1,jacdiff2,tmp
              end if
          end if

          maxerr = max( maxerr, tmp )
      end do

      if ( iprintout .ge. 1 ) then
          if ( nullcol ) then
              write(* ,130)
              write(10,130)
          else
              write(* ,140) maxerr
              write(10,140) maxerr
          end if
      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'Gradient vector of constraints ',I5,'.')
 110  format(/,1X,'Index',11X,'evaljac',2X,'Central diff (two ',
     +            'different steps)',4X,'Absolute error')
 120  format(  1X,I5,4(3X,1P,D15.8))
 130  format(  1X,'All the elements of this gradient are null.')
 140  format(  1X,'Maximum absolute error = ',1P,D15.8)

      end

C     *****************************************************************
C     *****************************************************************
      subroutine checkhc(n,x,ind,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer ind,inform,n

C     ARRAY ARGUMENTS
      double precision x(n)

C     This subrotutine checks the user supplied subroutine evalhc for
C     computing the Hessians of the constraints using finite
C     differences.

#include "dim.par"
#include "machconst.com"
#include "outtyp.com"

C     LOCAL SCALARS
      logical nullcol
      integer i,j,hnnz,jcnnz
      double precision elem,hdiff1,hdiff2,maxerr,step1,step2,tmp

C     LOCAL ARRAYS
      integer hlin(nsmax**2),hcol(nsmax**2),jcvar(nsmax)
      double precision g(nsmax),gplus1(nsmax),gplus2(nsmax),
     +        H(nsmax,nsmax),hval(nsmax**2),jcval(nsmax),maxcoe(nsmax)

C     Check viability of the test

      if ( n .gt. nsmax ) then
          if ( iprintout .ge. 1 ) then
              write(*, 100) nsmax,nsmax
              write(10,100) nsmax,nsmax
          end if

          return
      end if

C     Compute the gradient of constraint ind at x and save it in a
C     dense vector

      call vevaljac(n,x,ind,jcvar,jcval,jcnnz,inform)
      if ( inform .lt. 0 ) return

      do i = 1,n
          g(i) = 0.0d0
      end do

      do i = 1,jcnnz
          g(jcvar(i)) = g(jcvar(i)) + jcval(i)
      end do

C     Compute the Hessian of constraint ind at x and save it in a
C     dense matrix

      call vevalhc(n,x,ind,hlin,hcol,hval,hnnz,inform)
      if ( inform .lt. 0 ) return

      do j = 1,n
          do i = 1,n
              H(i,j) = 0.0d0
          end do
      end do

      do i = 1,hnnz
          H(hlin(i),hcol(i)) = H(hlin(i),hcol(i)) + hval(i)
      end do

      if ( iprintout .ge. 1 ) then
          write(* ,200) ind
          write(10,200) ind
      end if

      maxerr = 0.0d0

      do j = 1,n

          tmp  = x(j)

C         Compute the gradient of constraint ind at xplus1 and
C         save in a dense vector

          step1 = macheps12 * max( abs( tmp ), 1.0d0 )

          x(j) = tmp + step1
          call vevaljac(n,x,ind,jcvar,jcval,jcnnz,inform)
          if ( inform .lt. 0 ) return

          do i = 1,n
              gplus1(i) = 0.0d0
          end do

          do i = 1,jcnnz
              gplus1(jcvar(i)) = jcval(i)
          end do

C         Compute the gradient of constraint ind at xplus2 and
C         save in a dense vector

          step2 = macheps12 * max( abs( tmp ), 1.0d-03 )

          x(j) = tmp + step2
          call vevaljac(n,x,ind,jcvar,jcval,jcnnz,inform)
          if ( inform .lt. 0 ) return

          do i = 1,n
              gplus2(i) = 0.0d0
          end do

          do i = 1,jcnnz
              gplus2(jcvar(i)) = jcval(i)
          end do

          x(j) = tmp

          if ( iprintout .ge. 1 ) then
              write(* ,210) j
              write(10,210) j
          end if

          maxcoe(j) = 0.0d0

          nullcol = .true.

          do i = 1,n
              if ( i .ge. j ) then
                  elem = H(i,j)
              else
                  elem = H(j,i)
              end if
              hdiff1 = ( gplus1(i) - g(i) ) / step1
              hdiff2 = ( gplus2(i) - g(i) ) / step2
              tmp = min( abs( elem - hdiff1 ), abs( elem - hdiff2 ) )

              if ( elem   .ne. 0.0d0 .or.
     +             hdiff1 .ne. 0.0d0 .or.
     +             hdiff2 .ne. 0.0d0 ) then
                  if ( nullcol ) then
                      nullcol = .false.
                      if ( iprintout .ge. 1 ) then
                          write(* ,220)
                          write(10,220)
                      end if
                  end if
                  if ( iprintout .ge. 1 ) then
                      write(* ,230) i,elem,hdiff1,hdiff2,tmp
                      write(10,230) i,elem,hdiff1,hdiff2,tmp
                  end if
              end if

              maxcoe(j) = max( maxcoe(j), tmp )
          end do

          maxerr = max( maxerr, maxcoe(j) )

          if ( iprintout .ge. 1 ) then
              if ( nullcol ) then
                  write(* ,240)
                  write(10,240)
              else
                  write(* ,250) maxcoe(j)
                  write(10,250) maxcoe(j)
              end if
          end if

      end do

      if ( iprintout .ge. 1 ) then
          write(* ,*)
          write(10,*)

          do j = 1,n
              write(* ,260) j,maxcoe(j)
              write(10,260) j,maxcoe(j)
          end do

          write(* ,270) maxerr
          write(10,270) maxerr
      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'Subroutine CHECKHC uses dense matrices up to ',
     +            'dimension ',I6,' times ',I6,'. The Hessian ',
     +            'checking will be skipped.')

 200  format(/,1X,'Hessian matrix of constraint ',I5,' column by ',
     +            'column.')
 210  format(/,1X,'Column:  ',I6)
 220  format(/,1X,'Index',12X,'evalhc',3X,'Incr. Quoc. (two different ',
     +            'steps)',4X,'Absolute error')
 230  format(  1X,I5,4(3X,1P,D15.8))
 240  format(  1X,'All the elements of this column are null.')
 250  format(  1X,'Maximum absolute error = ',1P,D15.8)
 260  format(  1X,'Column ',I6,' Maximum absolute error = ',1P,D15.8)
 270  format(/,1X,'Overall maximum absolute error = ',1P,D15.8)

      end
