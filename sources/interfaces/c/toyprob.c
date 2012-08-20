/* =================================================================
   File: toyprob.c
   =================================================================

   =================================================================
   Module: Subroutines that define the problem
   =================================================================

   Last update of any of the component of this module:

   March 5, 2008.

   Users are encouraged to download periodically updated versions of
   this code at the TANGO home page:

   www.ime.usp.br/~egbirgin/tango/

   ******************************************************************
   ****************************************************************** */

#include <stdio.h>
#include <stdlib.h>

void inip(int *n,double **x,double **l,double **u,int *m,double **lambda,
          int **equatn,int **linear,int *coded,int *checkder)
{

   /* This subroutine set some problem data. */

   int i;

   /* Number of variables */

   *n = 2;

   /* Number of constraints (equalities plus inequalities) */

   *m = 2;

   /* Memory allocation */

   *x      = (double *) malloc(*n * sizeof(double));
   *l      = (double *) malloc(*n * sizeof(double));
   *u      = (double *) malloc(*n * sizeof(double));
   *lambda = (double *) malloc(*m * sizeof(double));
   *equatn = (int    *) malloc(*m * sizeof(int   ));
   *linear = (int    *) malloc(*m * sizeof(int   ));

   if (      *x == NULL ||      *l == NULL ||      *u == NULL ||
        *lambda == NULL || *equatn == NULL || *linear == NULL ) {

      printf( "\nC INTERFACE ERROR: in inip, malloc returned NULL.\n" );
      exit( 0 );

   }

   /* Initial point */

   for(i = 0; i < *n; i++)
      (*x)[i] = 0.0;

   /* Lower and upper bounds */

   (*l)[0] = - 10.0;
   (*l)[1] = - 1.0e20;

   (*u)[0] = 10.0;
   (*u)[1] = 1.0e20;

   /* For each constraint i, set equatn[i] = 1. if it is an equality
      constraint of the form c_i(x) = 0, and set equatn[i] = 0 if
      it is an inequality constraint of the form c_i(x) <= 0. */

   (*equatn)[0] = 0;
   (*equatn)[1] = 0;

   /* For each constraint i, set linear[i] = 1 if it is a linear
      constraint, otherwise set linear[i] = 0 */

   (*linear)[0] = 0;
   (*linear)[1] = 1;

   /* Lagrange multipliers approximation. */

   for( i = 0; i < *m; i++ )
      (*lambda)[i] = 0.0;

   /* In this C interface evalf, evalg, evalh, evalc, evaljac and
      evalhc are present. evalfc, evalgjac, evalhl and evalhlp are not. */

   coded[0] = 1; /* evalf    */
   coded[1] = 1; /* evalg    */
   coded[2] = 1; /* evalh    */
   coded[3] = 1; /* evalc    */
   coded[4] = 1; /* evaljac  */
   coded[5] = 1; /* evalhc   */
   coded[6] = 0; /* evalfc   */
   coded[7] = 0; /* evalgjac */
   coded[8] = 0; /* evalhl   */
   coded[9] = 0; /* evalhlp  */

   *checkder = 1;
}

/* ******************************************************************
   ****************************************************************** */

void evalf(int n,double *x,double *f,int *flag) {

/* This subroutine must compute the objective function. For achieving
   this objective YOU MUST MODIFY it according to your problem. See
   below where your modifications must be inserted.

   Parameters of the subroutine:

   On Entry:

   n        integer,
            number of variables,

   x        double precision x(n),
            current point,

   On Return

   f        double precision,
            objective function value at x,

   flag     integer,
            You must set it to any number different of 0 (zero) if
            some error ocurred during the evaluation of the objective
            function. (For example, trying to compute the square root
            of a negative number, dividing by zero or a very small
            number, etc.) If everything was o.k. you must set it
            equal to zero. */

   *flag = 0;

   *f = x[n - 1];
}

/* ******************************************************************
   ****************************************************************** */

void evalg(int n,double *x,double *g,int *flag) {

/* This subroutine must compute the gradient vector of the objective
   function. For achieving these objective YOU MUST MODIFY it in the
   way specified below. However, if you decide to use numerical
   derivatives (we dont encourage this option at all!) you dont need
   to modify evalg.

   Parameters of the subroutine:

   On Entry:

   n        integer,
            number of variables,

   x        double precision x(n),
            current point,

   On Return

   g        double precision g(n),
            gradient vector of the objective function evaluated at x,

   flag     integer,
            You must set it to any number different of 0 (zero) if
            some error ocurred during the evaluation of any component
            of the gradient vector. (For example, trying to compute
            the square root of a negative number, dividing by zero or
            a very small number, etc.) If everything was o.k. you
            must set it equal to zero. */

   int i;

   *flag = 0;

   for (i = 0; i < n-1; i++)
       g[i] = 0.0;

   g[n-1] = 1.0;
}

/* ******************************************************************
   ****************************************************************** */

void evalh(int n,double *x,int *hlin,int *hcol,double *hval,int *hnnz,
int *flag) {

/* This subroutine might compute the Hessian matrix of the objective
   function. For achieving this objective YOU MAY MODIFY it according
   to your problem. To modify this subroutine IS NOT MANDATORY. See
   below where your modifications must be inserted.

   Parameters of the subroutine:

   On Entry:

   n        integer,
            number of variables,

   x        double precision x(n),
            current point,

   On Return

   hnnz     integer,
            number of perhaps-non-null elements of the computed
            Hessian,

   hlin     integer hlin(hnnz),
            see below,

   hcol     integer hcol(hnnz),
            see below,

   hval     double precision hval(hnnz),
            the non-null value of the (hlin(k),hcol(k)) position
            of the Hessian matrix of the objective function must
            be saved at hval(k). Just the lower triangular part of
            Hessian matrix must be computed,

   flag     integer,
            You must set it to any number different of 0 (zero) if
            some error ocurred during the evaluation of the Hessian
            matrix of the objective funtion. (For example, trying
            to compute the square root of a negative number,
            dividing by zero or a very small number, etc.) If
            everything was o.k. you must set it equal to zero. */

   *flag = 0;

   *hnnz = 0;
}

/* ******************************************************************
   ****************************************************************** */

void evalc(int n,double *x,int ind,double *c,int *flag) {

/* This subroutine must compute the ind-th constraint of your problem.
   For achieving this objective YOU MUST MOFIFY it according to your
   problem. See below the places where your modifications must be
   inserted.

   Parameters of the subroutine:

   On Entry:

   n        integer,
            number of variables,

   x        double precision x(n),
            current point,

   ind      integer,
            index of the constraint to be computed,

   On Return

   c        double precision,
            ind-th constraint evaluated at x,

   flag     integer
            You must set it to any number different of 0 (zero) if
            some error ocurred during the evaluation of the
            constraint. (For example, trying to compute the square
            root of a negative number, dividing by zero or a very
            small number, etc.) If everything was o.k. you must set
            it equal to zero. */

   *flag = 0;

   if ( ind == 1 )
     *c = x[0] * x[0] + 1.0 - x[n-1];

   else if ( ind == 2 )
     *c = 2.0 - x[0] - x[n-1];

   else
     *flag = -1;
}

/* ******************************************************************
   ****************************************************************** */

void evaljac(int n,double *x,int ind,int *jcvar,double *jcval,
int *jcnnz,int *flag) {

/* This subroutine must compute the gradient of the ind-th constraint.
   For achieving these objective YOU MUST MODIFY it in the way specified
   below.

   Parameters of the subroutine:

   On Entry:

   n        integer,
            number of variables,

   x        double precision x(n),
            current point,

   ind      integer,
            index of the constraint whose gradient will be computed,

   On Return

   jcnnz   integer,
            number of perhaps-non-null elements of the computed
            gradient,

   jcvar   integer jcvar(jcnnz),
            see below,

   jcval   double precision jcval(jcnnz),
            the non-null value of the partial derivative of the ind-th
            constraint with respect to the jcvar(k)-th variable must
            be saved at jcval(k).

   flag     integer
            You must set it to any number different of 0 (zero) if
            some error ocurred during the evaluation of the
            constraint. (For example, trying to compute the square
            root of a negative number, dividing by zero or a very
            small number, etc.) If everything was o.k. you must set
            it equal to zero. */

   *flag = 0;

   if ( ind == 1 ) {
     *jcnnz = 2;

     jcvar[0] = 0;
     jcval[0] = 2.0 * x[0];

     jcvar[1] = n - 1;
     jcval[1] = - 1.0;
   }

   else if ( ind == 2 ) {
     *jcnnz = 2;

     jcvar[0] =   0;
     jcval[0] = - 1.0;

     jcvar[1] = n - 1;
     jcval[1] = - 1.0;
   }

   else
     *flag = -1;
}

/* ******************************************************************
   ****************************************************************** */

void evalhc(int n,double *x,int ind,int *hclin,int *hccol,double *hcval,
int *hcnnz,int *flag) {

/* This subroutine might compute the Hessian matrix of the ind-th
   constraint. For achieving this objective YOU MAY MODIFY it
   according to your problem. To modify this subroutine IS NOT
   MANDATORY. See below where your modifications must be inserted.

   Parameters of the subroutine:

   On Entry:

   n        integer,
            number of variables,

   x        double precision x(n),
            current point,

   ind      integer,
            index of the constraint whose Hessian will be computed,

   On Return

   hcnnz    integer,
            number of perhaps-non-null elements of the computed
            Hessian,

   hclin    integer hclin(hcnnz),
            see below,

   hccol    integer hccol(hcnnz),
            see below,

   hcval    double precision hcval(hcnnz),
            the non-null value of the (hclin(k),hccol(k)) position
            of the Hessian matrix of the ind-th constraint must
            be saved at hcval(k). Just the lower triangular part of
            Hessian matrix must be computed,

   flag     integer,
            You must set it to any number different of 0 (zero) if
            some error ocurred during the evaluation of the Hessian
            matrix of the ind-th constraint. (For example, trying
            to compute the square root of a negative number,
            dividing by zero or a very small number, etc.) If
            everything was o.k. you must set it equal to zero. */

   *flag = 0;

   if ( ind == 1 ) {
     *hcnnz = 1;

     hclin[0] = 0;
     hccol[0] = 0;
     hcval[0] = 2.0;
   }

   else if ( ind == 2 )
     *hcnnz = 0;

   else
     *flag = -1;
}

/* *****************************************************************
   ***************************************************************** */

void evalfc(int n,double *x,double *f,int m,double *c,int *flag) {

   *flag = -1;
}

/* *****************************************************************
   ***************************************************************** */

void evalgjac(int n,double *x,double *g,int m,int *jcfun,int *jcvar,
double *jcval,int *jcnnz,int *flag) {

   *flag = -1;
}

/* *****************************************************************
   ***************************************************************** */

void evalhl(int n,double *x,int m,double *lambda,double scalef,
double *scalec,int *hllin,int *hlcol,double *hlval,int *hlnnz,
int *flag) {

   *flag = -1;
}

/* *****************************************************************
   ***************************************************************** */

void evalhlp(int n,double *x,int m,double *lambda,double scalef,
double *scalec,double *p,double *hp,int *goth,int *flag) {

   *flag = -1;
}

/* ******************************************************************
   ****************************************************************** */

void endp(int n,double *x,double *l,double *u,int m,double *lambda,
int *equatn,int *linear)
{
   free(equatn);
   free(linear);
   free(l     );
   free(lambda);
   free(u     );
   free(x     );
}
