C     =================================================================
C     File: matlabwrapper64.f
C     =================================================================
C
C     =================================================================
C     Module: Matlab Wrapper
C     =================================================================

C     Last update: June 4th, 2008

C     ******************************************************************
C     ******************************************************************

#include "fintrf.h"

      subroutine mexFunction(nlhs,plhs,nrhs,prhs)

C     This subroutine is the gateway between MATLAB and ALGENCAN solver.
C     The function is defined in Matlab as:
C
C     [x,lambda,f,cnorm,snorm,nlpsupn,inform] =                   ...
C         algencan(evalf,evalg,evalh,evalc,evaljac,evalhc,evalfc, ...
C         evalgjac,evalhl,evalhlp,epsfeas,epsopt,iprint,ncomp,n,  ...
C         x,l,u,m,lambda,equatn,linear,coded,checkder);
C
C     but is called using 'algencanma' script, also defined in Matlab
C     (see how to call the solver, by typing 'help algencan' in Matlab
C     console). This subroutine expects that the functions inip and
C     param have already been called.

#include "dim.par"

C     PARAMETERS
      integer ARGS_NUMBER_IN,ARGS_NUMBER_OUT
      parameter( ARGS_NUMBER_IN  = 24 )
      parameter( ARGS_NUMBER_OUT =  7 )

C     SCALAR ARGUMENTS
      integer nlhs,nrhs

C     ARRAY ARGUMENTS
      mwPointer plhs(*),prhs(*)

C     LOCAL SCALARS
      integer countargs,status,i

      logical checkder
C      logical rhoauto
C      character*6 precond
      integer m,n,iprint,ncomp,inform
C     integer gtype,hptype,intype,maxtotit,maxoutit,maxtotfc,
      double precision epsfeas,epsopt,f,cnorm,snorm,nlpsupn
C     double precision rhofrac,rhomult,outiter,rhotype,totiter,
C    +        totfcnt,totgcnt,totcgcnt,
C      real time

C     LOCAL ARRAYS
      double precision tmparray(mmax)

C      integer wi1(nmax)
      logical coded(10),equatn(mmax),linear(mmax)
      double precision l(nmax),lambda(mmax),u(nmax),x(nmax)
C     double precision rho(mmax),wd1(mmax),
C     +        wd2(mmax),wd3(nmax),wd4(mmax),wd5(mmax),wd6(nmax),
C     +        wd7(nmax),wd8(nmax),wd9(nmax),wd10(nmax),wd11(nmax),
C     +        wd12(nmax),wd13(nmax),wd14(nmax),wd15(nmax),wd16(nmax),
C     +        wd17(nmax),wd18(nmax),wd19(nmax)

C     COMMON ARRAYS
      character*50 fnames(10),error(3)

C     COMMON BLOCKS
      common /matlab/ fnames

C     EXTERNAL FUNCTIONS
      integer mxGetString
      double precision mxGetScalar
      mwSize integerToMwSize
      mwPointer mxGetPr,mxCreateDoubleMatrix,mxCreateDoubleScalar


C     Error messages

      if ( nrhs .ne. ARGS_NUMBER_IN ) then
          call mexErrMsgTxt(error(1))
          return
      end if

      if ( nlhs .ne. ARGS_NUMBER_OUT ) then
          call mexErrMsgTxt(error(2))
          return
      end if

      countargs = 1


C     GET INFORMATION FROM MATLAB

C     User-defined function names

C      1 - evalf
C      2 - evalg
C      3 - evalh
C      4 - evalc
C      5 - evaljac
C      6 - evalhc
C      7 - evalfc
C      8 - evalgjac
C      9 - evalhl
C     10 - evalhlp


      do i = 1,10
          status = mxGetString(prhs(countargs),fnames(i),
     +    integerToMwSize(50))

          if ( status .ne. 0 ) then
              call mexErrMsgTxt(error(3))
              return
          end if

          countargs = countargs + 1

      end do

C     Other arguments

C     Feasibility epsilon

      epsfeas = mxGetScalar(prhs(countargs))
      countargs = countargs + 1

C     Optimality epsilon

      epsopt = mxGetScalar(prhs(countargs))
      countargs = countargs + 1


C     Output information detail

      iprint = int(mxGetScalar(prhs(countargs)))

      countargs = countargs + 1

C     ncomp

      ncomp = int(mxGetScalar(prhs(countargs)))

      countargs = countargs + 1

C     Number of variables

      n = int(mxGetScalar(prhs(countargs)))

      countargs = countargs + 1

C     Starting point

      call mxCopyPtrToReal8(mxGetPr(prhs(countargs)),x,
     +     integerToMwSize(n))

      countargs = countargs + 1

C     Lower and upper bounds

      call mxCopyPtrToReal8(mxGetPr(prhs(countargs)),l,
     +     integerToMwSize(n))

      countargs = countargs + 1

      call mxCopyPtrToReal8(mxGetPr(prhs(countargs)),u,
     +     integerToMwSize(n))

      countargs = countargs + 1

C     Number of constraints

      m = int(mxGetScalar(prhs(countargs)))

      countargs = countargs + 1

C     Lambda
      call mxCopyPtrToReal8(mxGetPr(prhs(countargs)),lambda,
     +     integerToMwSize(m))

      countargs = countargs + 1

C     Creating ALGENCAN's vector 'equatn'
C     (1 if i-th is an equality constraint 0 otherwise)

      call mxCopyPtrToReal8(mxGetPr(prhs(countargs)),tmparray,
     +     integerToMwSize(m))

      do i = 1,m
         if ( int(tmparray(i)) .eq. 1 ) then
            equatn(i) = .true.
         else
            equatn(i) = .false.
         end if
      end do

      countargs = countargs + 1

C     Creating ALGENCAN's vector 'linear'
C     (1 if i-th is a linear constraint 0 otherwise)

      call mxCopyPtrToReal8(mxGetPr(prhs(countargs)),tmparray,
     +     integerToMwSize(m))

      do i = 1,m
         if ( int(tmparray(i)) .eq. 1 ) then
            linear(i) = .true.
         else
            linear(i) = .false.
         end if
      end do

      countargs = countargs + 1

C     Creating ALGENCAN's vector 'coded'
C     (1 if i-th function was coded 0 otherwise)

      call mxCopyPtrToReal8(mxGetPr(prhs(countargs)),tmparray,
     +     integerToMwSize(10))

      do i = 1,10
         if ( int(tmparray(i)) .eq. 1 ) then
            coded(i) = .true.
         else
            coded(i) = .false.
         end if
      end do

      countargs = countargs + 1

C     Check derivatives = 0 or 1

      if ( int(mxGetScalar(prhs(countargs))) .eq. 1 ) then
         checkder = .true.
      else
         checkder = .false.
      end if


C     CALL OPTIMIZATION SOLVER

      call algencan(epsfeas,epsopt,iprint,ncomp,n,x,l,u,m,lambda,equatn,
     +linear,coded,checkder,f,cnorm,snorm,nlpsupn,inform)

C     RETURN INFORMATION TO MATLAB


      countargs = 1

C     Final estimation of the solution

      plhs(countargs) = mxCreateDoubleMatrix(integerToMwSize(n),
     +integerToMwSize(1),0)
      call mxCopyReal8ToPtr(x,mxGetPr(plhs(countargs)),
     +     integerToMwSize(n))

      countargs = countargs + 1

C     Final estimation of the Lagrange multipliers

      plhs(countargs) = mxCreateDoubleMatrix(integerToMwSize(n),
     +integerToMwSize(1),0)
      call mxCopyReal8ToPtr(lambda,mxGetPr(plhs(countargs)),
     +     integerToMwSize(n))

      countargs = countargs + 1

C     Objective functional value at the solution

      plhs(countargs) = mxCreateDoubleScalar(f)
      countargs = countargs + 1

C     Sup-norm of the ?

      plhs(countargs) = mxCreateDoubleScalar(cnorm)
      countargs = countargs + 1

C     Sup-norm of the constraints

      plhs(countargs) = mxCreateDoubleScalar(snorm)
      countargs = countargs + 1

C     Sup-norm of Augmented Lagrangian

      plhs(countargs) = mxCreateDoubleScalar(nlpsupn)
      countargs = countargs + 1

C     Status of ALGENCAN

      plhs(countargs) = mxCreateDoubleScalar(1.0d0 * inform)

      end


C     ******************************************************************
C     ******************************************************************

      subroutine evalf(n,x,f,flag)

      implicit none

C     This subroutine computes the objective function.

C     SCALAR ARGUMENTS
      integer flag,n
      double precision f

C     ARRAY ARGUMENTS
      double precision x(n)

C     LOCAL ARRAYS
      mwPointer plhs(2), prhs(2)

C     COMMON ARRAYS
      character*50 fnames(10)

C     COMMON BLOCKS
      common /matlab/ fnames

C     EXTERNAL FUNCTIONS
      mwPointer mxCreateDoubleMatrix,mxCreateDoubleScalar,mxGetPr
      mwSize integerToMwSize
      double precision mxGetScalar

C     INPUT VALUES

      prhs(1) = mxCreateDoubleScalar(1.0d0 * n)

      prhs(2) = mxCreateDoubleMatrix(integerToMwSize(n),
     +integerToMwSize(1),0)
      call mxCopyReal8ToPtr(x,mxGetPr(prhs(2)),integerToMwSize(n))

C     CALL MATLAB FUNCTION

      call mexCallMATLAB(2,plhs,2,prhs,'evalf')

C     OUTPUT VALUES

C     Objective function value.
      f = mxGetScalar(plhs(1))

C     Flag.
      flag = int(mxGetScalar(plhs(2)))

C     CLEAR UNUSED STRUCTURES

      call mxDestroyArray(prhs(1))
      call mxDestroyArray(prhs(2))

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalg(n,x,g,flag)

      implicit none

C     This subroutine computes the gradient vector of the objective
C     function.

C     SCALAR ARGUMENTS
      integer flag,n

C     ARRAY ARGUMENTS
      double precision g(n),x(n)

C     LOCAL ARRAYS
      mwPointer plhs(2), prhs(2)

C     COMMON ARRAYS
      character*50 fnames(10)

C     COMMON BLOCKS
      common /matlab/ fnames

C     EXTERNAL FUNCTIONS
      mwPointer mxCreateDoubleMatrix,mxCreateDoubleScalar,mxGetPr
      mwSize integerToMwSize
      double precision mxGetScalar

C     INPUT VALUES

      prhs(1) = mxCreateDoubleScalar(1.0d0 * n)

      prhs(2) = mxCreateDoubleMatrix(integerToMwSize(n),
     +integerToMwSize(1),0)
      call mxCopyReal8ToPtr(x,mxGetPr(prhs(2)),integerToMwSize(n))

C     CALL MATLAB FUNCTION

      call mexCallMATLAB(2,plhs,2,prhs,fnames(2))

C     OUTPUT VALUES

C     The gradient vector.
      call mxCopyPtrToReal8(mxGetPr(plhs(1)),g,integerToMwSize(n))

C     Flag.
      flag = int(mxGetScalar(plhs(2)))

C     CLEAR UNUSED STRUCTURES

      call mxDestroyArray(prhs(1))
      call mxDestroyArray(prhs(2))

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalh(n,x,hlin,hcol,hval,hnnz,flag)

      implicit none

C     This subroutine computes the Hessian matrix of the objective
C     function.

#include "dim.par"

C     SCALAR ARGUMENTS
      integer flag,n,hnnz

C     ARRAY ARGUMENTS
      integer hcol(*),hlin(*)
      double precision hval(*),x(n)

C     LOCAL SCALARS
      integer i,j,k
      mwPointer sparsematrix

C     LOCAL ARRAYS
      integer*8 ir(nmax),jc(nmax)
      mwPointer plhs(2),prhs(2)

C     COMMON ARRAYS
      character*50 fnames(10)

C     COMMON BLOCKS
      common /matlab/ fnames

C     EXTERNAL FUNCTIONS
      mwIndex mxGetIr,mxGetJc
      mwSize integerToMwSize
      mwPointer mxCreateDoubleMatrix,mxCreateDoubleScalar,mxGetPr
      integer mxIsSparse
      double precision mxGetScalar

C     INPUT VALUES

      prhs(1) = mxCreateDoubleScalar(1.0d0 * n)

      prhs(2) = mxCreateDoubleMatrix(integerToMwSize(n),
     +integerToMwSize(1),0)
      call mxCopyReal8ToPtr(x,mxGetPr(prhs(2)),integerToMwSize(n))

C     CALL MATLAB FUNCTION

      call mexCallMATLAB(2,plhs,2,prhs,fnames(3))

C     OUTPUT VALUES

      if ( mxIsSparse(plhs(1)) .eq. 0 ) then
C         Creating sparse matrix (if it is not sparse)
          call mexCallMATLAB(1,sparsematrix,1,plhs(1),'sparse')
      else

          sparsematrix = plhs(1)

      endif

C     Sparse matrix structure.

      call mxCopyPtrToInteger8(mxGetJc(sparsematrix),jc,
     +     integerToMwSize(n + 1))

      hnnz = jc(n + 1)

C     Conversion from MatLab sparsity mode to algencan's mode.

      call mxCopyPtrToReal8(mxGetPr(sparsematrix),hval,
     +     integerToMwSize(hnnz))

C     TAVA DANDO WARNING POR ISSO TROQUEI!!!!
      call mxCopyPtrToInteger8(mxGetIr(sparsematrix),ir,
     +     integerToMwSize(hnnz))

      k = 0
      do i = 1,n
          do j = 1,jc(i+1) - jc(i)
              k = k + 1
              hlin(k) = ir(k) + 1
              hcol(k) = i
          end do
      end do

C      call mxCopyPtrToInteger8(mxGetIr(sparsematrix),hlin,
C     +     integerToMwSize(hnnz))
C
C      k = 0
C      do i = 1,n
C          do j = 1,jc(i+1) - jc(i)
C              k = k + 1
C              hlin(k) = hlin(k) + 1
C              hcol(k) = i
C          end do
C      end do

C     Flag.
      flag = int(mxGetScalar(plhs(2)))

C     CLEAR UNUSED STRUCTURES

      call mxDestroyArray(prhs(1))
      call mxDestroyArray(prhs(2))

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalc(n,x,ind,c,flag)

      implicit none

C     This subroutine computes the ind-th constraint of your
C     problem.

C     SCALAR ARGUMENTS
      integer ind,flag,n
      double precision c

C     ARRAY ARGUMENTS
      double precision x(n)

C     LOCAL ARRAYS
      mwPointer plhs(2), prhs(3)

C     COMMON ARRAYS
      character*50 fnames(10)

C     COMMON BLOCKS
      common /matlab/ fnames

C     EXTERNAL FUNCTIONS
      mwPointer mxCreateDoubleMatrix,mxCreateDoubleScalar,mxGetPr
      mwSize integerToMwSize
      double precision mxGetScalar

C     INPUT VALUES

      prhs(1) = mxCreateDoubleScalar(1.0d0 * n)

      prhs(2) = mxCreateDoubleMatrix(integerToMwSize(n),
     +integerToMwSize(1),0)
      call mxCopyReal8ToPtr(x,mxGetPr(prhs(2)),integerToMwSize(n))

      prhs(3) = mxCreateDoubleScalar(1.0d0 * ind)

C     CALL MATLAB FUNCTION

      call mexCallMATLAB(2,plhs,3,prhs,fnames(4))

C     OUTPUT VALUES

C     ind-th contraint value.
      c = mxGetScalar(plhs(1))

C     Flag.
      flag = int(mxGetScalar(plhs(2)))

C     CLEAR UNUSED STRUCTURES

      call mxDestroyArray(prhs(1))
      call mxDestroyArray(prhs(2))
      call mxDestroyArray(prhs(3))

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evaljac(n,x,ind,jcvar,jcval,jcnnz,flag)

      implicit none

#include "dim.par"

C     This subroutine computes the gradient of the ind-th constraint.

C     SCALAR ARGUMENTS
      integer flag,ind,n,jcnnz

C     ARRAY ARGUMENTS
      integer jcvar(*)
      double precision x(n),jcval(*)

C     LOCAL SCALARS
      integer i
      mwPointer sparsematrix

C     LOCAL ARRAYS
      mwPointer plhs(2),prhs(3)
      integer*8 jccol(2),jclin(nmax)

C     COMMON ARRAYS
      character*50 fnames(10)

C     COMMON BLOCKS
      common /matlab/ fnames

C     EXTERNAL FUNCTIONS
      mwIndex mxGetIr,mxGetJc
      mwPointer mxCreateDoubleScalar,mxCreateDoubleMatrix,mxGetPr
      mwSize integerToMwSize
      integer mxIsSparse
      double precision mxGetScalar

C     INPUT VALUES

      prhs(1) = mxCreateDoubleScalar(1.0d0 * n)

      prhs(2) = mxCreateDoubleMatrix(integerToMwSize(n),
     +integerToMwSize(1),0)
      call mxCopyReal8ToPtr(x,mxGetPr(prhs(2)),integerToMwSize(n))

      prhs(3) = mxCreateDoubleScalar(1.0d0 * ind)

C     CALL MATLAB FUNCTION

      call mexCallMATLAB(2,plhs,3,prhs,fnames(5))

      if ( mxIsSparse(plhs(1)) .eq. 0 ) then
C         Creating sparse matrix (if it is not sparse)
          call mexCallMATLAB(1,sparsematrix,1,plhs(1),'sparse')
      else
          sparsematrix = plhs(1)
      endif

C     OUTPUT VALUES

C     Sparse matrix structure.

C     (Here, we suppose that the returned Jacobian is a column vector.)

      call mxCopyPtrToInteger8(mxGetJc(sparsematrix),jccol,
     +     integerToMwSize(2))

      jcnnz = jccol(2)

      call mxCopyPtrToInteger8(mxGetIr(sparsematrix),jclin,
     +     integerToMwSize(jcnnz))
      call mxCopyPtrToReal8(mxGetPr(sparsematrix),jcval,
     +     integerToMwSize(jcnnz))

      do i = 1,jcnnz
          jcvar(i) = jclin(i) + 1
      end do

C     Flag.
      flag = int(mxGetScalar(plhs(2)))

C     CLEAR UNUSED STRUCTURES

      call mxDestroyArray(prhs(1))
      call mxDestroyArray(prhs(2))
      call mxDestroyArray(prhs(3))

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalhc(n,x,ind,hclin,hccol,hcval,hcnnz,flag)

      implicit none

C     This subroutine computes the Hessian matrix of the ind-th
C     constraint.

#include "dim.par"

C     SCALAR ARGUMENTS
      integer flag,ind,n,hcnnz

C     ARRAY ARGUMENTS
      integer hccol(*),hclin(*)
      double precision hcval(*),x(n)

C     LOCAL SCALARS
      integer i,j,k
      mwPointer sparsematrix

C     LOCAL ARRAYS
      integer*8 jc(nmax),ir(nmax)
      mwPointer plhs(2),prhs(3)

C     COMMON ARRAYS
      character*50 fnames(10)

C     COMMON BLOCKS
      common /matlab/ fnames

C     EXTERNAL FUNCTIONS
      mwIndex mxGetIr,mxGetJc
      mwPointer mxCreateDoubleScalar,mxCreateDoubleMatrix,mxGetPr
      mwSize integerToMwSize
      integer mxIsSparse
      double precision mxGetScalar

C     INPUT VALUES

      prhs(1) = mxCreateDoubleScalar(1.0d0 * n)

      prhs(2) = mxCreateDoubleMatrix(integerToMwSize(n),
     +integerToMwSize(1),0)
      call mxCopyReal8ToPtr(x,mxGetPr(prhs(2)),integerToMwSize(n))

      prhs(3) = mxCreateDoubleScalar(1.0d0 * ind)

C     CALL MATLAB FUNCTION

      call mexCallMATLAB(2,plhs,3,prhs,fnames(6))

      if ( mxIsSparse(plhs(1)) .eq. 0 ) then
C         Creating sparse matrix (if it is not sparse)
          call mexCallMATLAB(1,sparsematrix,1,plhs(1),'sparse')
      else
          sparsematrix = plhs(1)
      endif

C     OUTPUT VALUES

C     Sparse matrix structure.

      call mxCopyPtrToInteger8(mxGetJc(sparsematrix),jc,
     +     integerToMwSize(n + 1))

      hcnnz = jc(n + 1)

      call mxCopyPtrToReal8(mxGetPr(sparsematrix),hcval,
     +     integerToMwSize(hcnnz))
      call mxCopyPtrToInteger8(mxGetIr(sparsematrix),ir,
     +     integerToMwSize(hcnnz))

      k = 0
      do i = 1,n
          do j = 1,jc(i+1) - jc(i)
              k = k + 1
              hclin(k) = ir(k) + 1
              hccol(k) = i
          end do
      end do

C     Flag.
      flag = int(mxGetScalar(plhs(2)))

C     CLEAR UNUSED STRUCTURES

      call mxDestroyArray(prhs(1))
      call mxDestroyArray(prhs(2))
      call mxDestroyArray(prhs(3))

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

C     LOCAL ARRAYS
      mwPointer plhs(3),prhs(3)

C     COMMON ARRAYS
      character*50 fnames(10)

C     COMMON BLOCKS
      common /matlab/ fnames

C     EXTERNAL FUNCTIONS
      mwIndex mxGetIr,mxGetJc
      mwPointer mxCreateDoubleScalar,mxCreateDoubleMatrix,mxGetPr
      mwSize integerToMwSize
      double precision mxGetScalar

C     INPUT VALUES

      prhs(1) = mxCreateDoubleScalar(1.0d0 * n)

      prhs(2) = mxCreateDoubleMatrix(integerToMwSize(n),
     +     integerToMwSize(1),0)
      call mxCopyReal8ToPtr(x,mxGetPr(prhs(2)),integerToMwSize(n))

      prhs(3) = mxCreateDoubleScalar(1.0d0 * m)

C     CALL MATLAB FUNCTION

      call mexCallMATLAB(3,plhs,3,prhs,fnames(7))

C     OUTPUT VALUES

C     Objective function value.
      f = mxGetScalar(plhs(1))

C     The vector of contraints.
      call mxCopyPtrToReal8(mxGetPr(plhs(2)),c,integerToMwSize(n))

C     Flag.
      flag = int(mxGetScalar(plhs(3)))

C     CLEAR UNUSED STRUCTURES

      call mxDestroyArray(prhs(1))
      call mxDestroyArray(prhs(2))
      call mxDestroyArray(prhs(3))

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,flag)

      implicit none

#include "dim.par"

C     SCALAR ARGUMENTS
      integer flag,jcnnz,m,n

C     ARRAY ARGUMENTS
      integer jcfun(*),jcvar(*)
      double precision g(n),jcval(*),x(n)

C     LOCAL SCALARS
      integer i,j,k
      mwPointer sparsematrix

C     LOCAL ARRAYS
      mwPointer plhs(3),prhs(3)
      integer*8 jc(nmax),jclin(nmax)

C     COMMON ARRAYS
      character*50 fnames(10)

C     COMMON BLOCKS
      common /matlab/ fnames

C     EXTERNAL FUNCTIONS
      mwIndex mxGetIr,mxGetJc
      mwPointer mxCreateDoubleScalar,mxCreateDoubleMatrix,mxGetPr
      mwSize integerToMwSize
      integer mxIsSparse
      double precision mxGetScalar


C     INPUT VALUES

      prhs(1) = mxCreateDoubleScalar(1.0d0 * n)

      prhs(2) = mxCreateDoubleMatrix(integerToMwSize(n),
     +     integerToMwSize(1),0)
      call mxCopyReal8ToPtr(x,mxGetPr(prhs(2)),integerToMwSize(n))

      prhs(3) = mxCreateDoubleScalar(1.0d0 * m)

C     CALL MATLAB FUNCTION

      call mexCallMATLAB(3,plhs,3,prhs,fnames(8))

      if ( mxIsSparse(plhs(2)) .eq. 0 ) then
C         Creating sparse matrix (if it is not sparse)
          call mexCallMATLAB(1,sparsematrix,1,plhs(2),'sparse')
      else
          sparsematrix = plhs(2)
      endif

C     OUTPUT VALUES

C     The gradient vector.
      call mxCopyPtrToReal8(mxGetPr(plhs(1)),g,integerToMwSize(n))

C     Sparse matrix structure.

      call mxCopyPtrToInteger8(mxGetJc(sparsematrix),jc,
     +     integerToMwSize(n+1))

      jcnnz = jc(n+1)

      call mxCopyPtrToReal8(mxGetPr(sparsematrix),jcval,
     +     integerToMwSize(jcnnz))
      call mxCopyPtrToInteger8(mxGetIr(sparsematrix),jclin,
     +     integerToMwSize(jcnnz))

      k = 0
      do i = 1,n
          do j = 1,jc(i+1) - jc(i)
              k = k + 1
              jcfun(k) = jclin(k) + 1
              jcvar(k) = i
          end do
      end do

C     Flag.
      flag = int(mxGetScalar(plhs(3)))

C     CLEAR UNUSED STRUCTURES

      call mxDestroyArray(prhs(1))
      call mxDestroyArray(prhs(2))
      call mxDestroyArray(prhs(3))

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalhl(n,x,m,lambda,scalef,scalec,hllin,hlcol,hlval,
     +hlnnz,flag)

      implicit none

#include "dim.par"

C     SCALAR ARGUMENTS
      integer flag,hlnnz,m,n
      double precision scalef

C     ARRAY ARGUMENTS
      integer hlcol(*),hllin(*)
      double precision hlval(*),lambda(m),scalec(m),x(n)

C     LOCAL SCALARS
      mwPointer sparsematrix
      integer i,j,k

C     LOCAL ARRAYS
      mwPointer plhs(2),prhs(6)
      integer*8 jc(nmax),ir(nmax)

C     COMMON ARRAYS
      character*50 fnames(10)

C     COMMON BLOCKS
      common /matlab/ fnames

C     EXTERNAL FUNCTIONS
      mwIndex mxGetIr,mxGetJc
      mwPointer mxCreateDoubleScalar,mxCreateDoubleMatrix,mxGetPr
      mwSize integerToMwSize
      integer mxIsSparse
      double precision mxGetScalar


C     INPUT VALUES

      prhs(1) = mxCreateDoubleScalar(1.0d0 * n)

      prhs(2) = mxCreateDoubleMatrix(integerToMwSize(n),
     +     integerToMwSize(1),0)
      call mxCopyReal8ToPtr(x,mxGetPr(prhs(2)),integerToMwSize(n))

      prhs(3) = mxCreateDoubleScalar(1.0d0 * m)

      prhs(4) = mxCreateDoubleMatrix(integerToMwSize(n),
     +     integerToMwSize(1),0)
      call mxCopyReal8ToPtr(lambda,mxGetPr(prhs(4)),
     +     integerToMwSize(m))

      prhs(5) = mxCreateDoubleScalar(scalef)

      prhs(6) = mxCreateDoubleMatrix(integerToMwSize(n),
     +     integerToMwSize(1),0)
      call mxCopyReal8ToPtr(scalec,mxGetPr(prhs(6)),
     +     integerToMwSize(m))

C     CALL MATLAB FUNCTION

      call mexCallMATLAB(2,plhs,6,prhs,fnames(9))

      if ( mxIsSparse(plhs(1)) .eq. 0 ) then
C         Creating sparse matrix (if it is not sparse)
          call mexCallMATLAB(1,sparsematrix,1,plhs(1),'sparse')
      else
          sparsematrix = plhs(1)
      endif

C     OUTPUT VALUES

C     Sparse matrix structure.

      call mxCopyPtrToInteger8(mxGetJc(sparsematrix),jc,
     +     integerToMwSize(n + 1))

      hlnnz = jc(n + 1)

      call mxCopyPtrToReal8(mxGetPr(sparsematrix),hlval,
     +     integerToMwSize(hlnnz))
      call mxCopyPtrToInteger8(mxGetIr(sparsematrix),ir,
     +     integerToMwSize(hlnnz))

      k = 0
      do i = 1,n
          do j = 1,jc(i+1) - jc(i)
              k = k + 1
              hllin(k) = ir(k) + 1
              hlcol(k) = i
          end do
      end do

C     Flag.
      flag = int(mxGetScalar(plhs(2)))

C     CLEAR UNUSED STRUCTURES

      do i = 1,6
          call mxDestroyArray(prhs(i))
      end do

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

C     LOCAL SCALARS
      integer i

C     LOCAL ARRAYS
      mwPointer plhs(3),prhs(8)

C     COMMON ARRAYS
      character*50 fnames(10)

C     COMMON BLOCKS
      common /matlab/ fnames

C     EXTERNAL FUNCTIONS
      mwPointer mxCreateDoubleScalar,mxCreateDoubleMatrix,mxGetPr
      mwSize integerToMwSize
      double precision mxGetScalar

C     INPUT VALUES

      prhs(1) = mxCreateDoubleScalar(1.0d0 * n)

      prhs(2) = mxCreateDoubleMatrix(integerToMwSize(n),
     +integerToMwSize(1),0)
      call mxCopyReal8ToPtr(x,mxGetPr(prhs(2)),integerToMwSize(n))

      prhs(3) = mxCreateDoubleScalar(1.0d0 * m)

      prhs(4) = mxCreateDoubleMatrix(integerToMwSize(m),
     +integerToMwSize(1),0)
      call mxCopyReal8ToPtr(lambda,mxGetPr(prhs(4)),integerToMwSize(m))

      prhs(5) = mxCreateDoubleScalar(sf)

      prhs(6) = mxCreateDoubleMatrix(integerToMwSize(m),
     +integerToMwSize(1),0)
      call mxCopyReal8ToPtr(sc,mxGetPr(prhs(6)),integerToMwSize(m))

      prhs(7) = mxCreateDoubleMatrix(integerToMwSize(n),
     +integerToMwSize(1),0)
      call mxCopyReal8ToPtr(p,mxGetPr(prhs(7)),integerToMwSize(n))

      if ( goth ) then
          prhs(8) = mxCreateDoubleScalar(1.0d0)
      else
          prhs(8) = mxCreateDoubleScalar(0.0d0)
      end if

C     CALL MATLAB FUNCTION

      call mexCallMATLAB(3,plhs,8,prhs,fnames(10))

C     OUTPUT VALUES

C     Hessian of the Lagrangian times vector p
      call mxCopyPtrToReal8(mxGetPr(plhs(1)),hp,integerToMwSize(n))

C     Goth
      if ( int(mxGetScalar(plhs(2))) .ne. 0 ) then
          goth = .true.
      else
          goth = .false.
      end if

C     Flag
      flag = int(mxGetScalar(plhs(3)))

C     CLEAR UNUSED STRUCTURES

      do i = 1,8
          call mxDestroyArray(prhs(i))
      end do

      end

C     ******************************************************************
C     ******************************************************************

      mwSize function integerToMwSize(number)

C     SCALAR ARGUMENTS
      integer number

      integerToMwSize = number

      end
