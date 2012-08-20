/* =================================================================
   File: cwrapper.c
   =================================================================

   =================================================================
   Module: C Wrapper
   =================================================================

   Last update of any of the component of this module:

   March 5, 2008.

   Users are encouraged to download periodically updated versions of
   this code at the TANGO home page:

   www.ime.usp.br/~egbirgin/tango/

   *****************************************************************

   TANGO LICENSE:
   --------------

   TANGO is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation. Non-free versions of TANGO are
   available under terms different from those of the General Public
   License. Professors J. M. Martínez (martinez@ime.unicamp.br,
   martinezimecc@gmail.com) or E. G. Birgin (egbirgin@ime.usp.br,
   egbirgin@gmail.com) should be contacted for more information
   related to such a license, future developments and/or technical
   support.

   Every published work that uses ALGENCAN should cite:

   R. Andreani, E. G. Birgin, J. M. Martínez and M. L. Schuverdt,
   "On Augmented Lagrangian methods with general lower-level
   constraints", to appear in SIAM Journal on Optimization.

   and

   R. Andreani, E. G. Birgin, J. M. Martínez and M. L. Schuverdt,
   "Augmented Lagrangian methods under the Constant Positive Linear
   Dependence constraint qualification", Mathematical
   Programming, 111, pp. 5-32, 2008.

   Every published work that uses GENCAN should cite:

   E. G. Birgin and J. M. Martínez, "Large-scale active-set
   box-constrained optimization method with spectral projected
   gradients", Computational Optimization and Applications 23, pp.
   101-125, 2002.

   (See other related works at the TANGO home page.)*/

/* ******************************************************************
   ****************************************************************** */

#include <stdio.h>
#include <string.h>
#include "cfortran.h"
#include "cwrapper.h"

/* ******************************************************************
   ****************************************************************** */

void Evalh(int n,double *x,int *hlin,int *hcol,double *hval,int *hnnz,
int *flag)
{

   int i;

   evalh(n,x,hlin,hcol,hval,hnnz,flag);

   for ( i = 0; i < *hnnz; i++ ) {

      hlin[i]++;
      hcol[i]++;
   }

}

/* ******************************************************************
   ****************************************************************** */

void Evaljac(int n,double *x,int ind,int *jcvar,double *jcval,int *jcnnz,
int *flag)
{

   int i;

   evaljac(n,x,ind,jcvar,jcval,jcnnz,flag);

   for ( i = 0; i < *jcnnz; i++ )
      jcvar[i]++;

}

/* ******************************************************************
   ****************************************************************** */

void Evalhc(int n,double *x,int ind,int *hclin,int *hccol,double *hcval,
int *hcnnz,int *flag)
{

   int i;

   evalhc(n,x,ind,hclin,hccol,hcval,hcnnz,flag);

   for ( i = 0; i < *hcnnz; i++ ) {

      hclin[i]++;
      hccol[i]++;
   }

}

/* ******************************************************************
   ****************************************************************** */

void Evalgjac(int n,double *x,double *g,int m,int *jcfun,int *jcvar,
double *jcval,int *jcnnz,int *flag)
{

   int i;

   evalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,flag);

   for ( i = 0; i < *jcnnz; i++ )
      jcvar[i]++;

}

/* ******************************************************************
   ****************************************************************** */

void Evalhl(int n,double *x,int m,double *lambda,double scalef,
double *scalec,int *hllin,int *hlcol,double *hlval,int *hlnnz,
int *flag)
{

   int i;

   evalhl(n,x,m,lambda,scalef,scalec,hllin,hlcol,hlval,hlnnz,flag);

   for ( i = 0; i < *hlnnz; i++ ) {

      hllin[i]++;
      hlcol[i]++;
   }

}

/* ******************************************************************
   ****************************************************************** */

void algencan(double epsfeas,double epsopt,int iprint,int ncomp,int n,
double *x,double *l,double *u,int m,double *lambda,int *equatn,int *linear,
int *coded,int checkder,double *f,double *cnorm,double *snorm,double *nlpsupn,
int *inform)
{

/* TRANSFORMS THE TYPE INT OF C IN THE TYPE LOGICAL OF FORTRAN
 * FOR THE ARRAYS equatn AND linear, AND VARIABLE checkder */
   C2FLOGICALV(equatn,m);
   C2FLOGICALV(linear,m);
   C2FLOGICALV(coded,10);

/* CALL OPTIMIZATION SOLVER */

   Algencan(epsfeas,epsopt,iprint,ncomp,n,x,l,u,m,lambda,equatn,
   linear,coded,checkder,*f,*cnorm,*snorm,*nlpsupn,*inform);

}
