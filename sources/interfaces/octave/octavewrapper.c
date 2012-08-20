/* ==================================================================
   Module: Interface between Octave and ALGENCAN
   ==================================================================

   Coded by: Ricardo Andrade

   Last update of any of the component of this module:

   September 11, 2009.

   Users are encouraged to download periodically updated versions of
   this code at the TANGO home page:

   www.ime.usp.br/~egbirgin/tango/

   ******************************************************************
   ******************************************************************

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
   Programming 111, pp. 5-32, 2008.

   Every published work that uses GENCAN should cite:

   E. G. Birgin and J. M. Martínez, "Large-scale active-set
   box-constrained optimization method with spectral projected
   gradients", Computational Optimization and Applications 23, pp.
   101-125, 2002.

   (See other related works at the TANGO home page.)

   ****************************************************************** */


/* ******************************************************************
   ****************************************************************** */

#include "../c/cfortran.h"
#include "octavewrapper.h"
#include <string.h>
#include <iostream>
#include <oct.h>
#include <parse.h>

/* ******************************************************************
   ****************************************************************** */

   void inip (int* n,double** x,double** l, double** u,int* m,
	      double** lambda,int** equatn,int** linear, int* coded,
              int* checkder) {

   int i;
   NDArray x_array, l_array, u_array;
   NDArray lambda_array, equatn_array, linear_array, coded_array;
   octave_value_list input_arguments, output_values;

   output_values = feval("inip",input_arguments,0);

   *n        = output_values(0).int_value();
   *m        = output_values(4).int_value();
   *coded    = output_values(8).int_value();
   *checkder = output_values(9).int_value();

   *x      = (double *) malloc(*n * sizeof(double));
   *l      = (double *) malloc(*n * sizeof(double));
   *u      = (double *) malloc(*n * sizeof(double));
   *lambda = (double *) malloc(*m * sizeof(double));
   *equatn = (int    *) malloc(*m * sizeof(int   ));
   *linear = (int    *) malloc(*m * sizeof(int   ));


   x_array = output_values(1).array_value();
   l_array = output_values(2).array_value();
   u_array = output_values(3).array_value();
   for(i = 0; i < *n; i++) {
     (*x)[i] = x_array.elem(i);
     (*l)[i] = l_array.elem(i);
     (*u)[i] = u_array.elem(i);
   }

   lambda_array = output_values(5).array_value();
   equatn_array = output_values(6).array_value();
   linear_array = output_values(7).array_value();
   for(i = 0; i < *m; i++) {
     (*lambda)[i] = lambda_array.elem(i);
     (*equatn)[i] = (int) equatn_array.elem(i);
     (*linear)[i] = (int) linear_array.elem(i);
   }

   coded_array = output_values(8).array_value();
   for(i = 0; i < 10; i++) {
     coded[i] = (int) coded_array.elem(i);
   }

   }

/* ******************************************************************
   ****************************************************************** */

   void evalf(int n,double* x,double* f,int* flag) {

   int i;
   octave_value_list input_arguments, output_values;
   NDArray x_array;

   x_array.resize_no_fill(n);
   for (i=0; i < n; ++i){
     x_array(i) = x[i];
   }

   input_arguments(0) = n;
   input_arguments(1) = x_array;

   output_values = feval("evalf",input_arguments,2);

   *flag = output_values(1).int_value();

   if (*flag == 0) {
     *f = output_values(0).double_value();
   }


   }

/* ******************************************************************
   ****************************************************************** */

   void evalg(int n,double* x,double* g,int* flag) {

   int i;
   octave_value_list input_arguments, output_values;
   NDArray x_array, g_array;

   x_array.resize_no_fill(n);
   for (i=0; i < n; ++i){
     x_array(i) = x[i];
   }

   input_arguments(0) = n;
   input_arguments(1) = x_array;

   output_values = feval("evalg",input_arguments,2);

   *flag = output_values(1).int_value();

   if (*flag == 0) {
     g_array = output_values(0).array_value();
     for (i = 0; i < n; i++)
       g[i] = g_array.elem(i);
   }

   }

/* ******************************************************************
   ****************************************************************** */

   void evalh(int n,double* x,int* hlin,int* hcol,double* hval,int* hnnz,
   int* flag) {

   int i;
   NDArray x_array;
   NDArray hlin_array, hcol_array, hval_array;
   octave_value_list input_arguments, output_values;

   x_array.resize_no_fill(n);
   for (i=0; i < n; ++i){
     x_array(i) = x[i];
   }

   input_arguments(0) = n;
   input_arguments(1) = x_array;

   output_values = feval("evalh",input_arguments,2);

   *flag = output_values(4).int_value();

   if (*flag == 0){
     *hnnz = output_values(3).int_value();
     if (*hnnz > 0 ){
       hlin_array = output_values(0).array_value();
       hcol_array = output_values(1).array_value();
       hval_array = output_values(2).array_value();
       for (i = 0; i < *hnnz; i++) {
	 hlin[i] = (int) hlin_array.elem(i);
	 hcol[i] = (int) hcol_array.elem(i);
	 hval[i] = hval_array.elem(i);
       }
     }
   }


   }

/* ******************************************************************
   ****************************************************************** */

   void evalc(int n,double* x,int ind,double* cind,int* flag) {

   int i;
   octave_value_list input_arguments, output_values;
   NDArray x_array;

   x_array.resize_no_fill(n);
   for (i=0; i < n; ++i){
     x_array(i) = x[i];
   }

   input_arguments(0) = n;
   input_arguments(1) = x_array;
   input_arguments(2) = ind;

   output_values = feval("evalc",input_arguments,3);

   *flag = output_values(1).int_value();

   if (*flag == 0){
     *cind = output_values(0).double_value();
   }


   }

/* ******************************************************************
   ****************************************************************** */

   void evaljac(int n,double* x,int ind,int* jcvar,double* jcval,int* jcnnz,
   int* flag) {

   int i;
   octave_value_list input_arguments, output_values;
   NDArray jcvar_array, jcval_array;
   NDArray x_array;

   x_array.resize_no_fill(n);
   for (i=0; i < n; ++i){
     x_array(i) = x[i];
   }

   input_arguments(0) = n;
   input_arguments(1) = x_array;
   input_arguments(2) = ind;

   output_values = feval("evaljac",input_arguments,3);

   *flag = output_values(3).int_value();

   if (*flag == 0){
     *jcnnz = output_values(2).int_value();
     if (*jcnnz > 0) {
       jcvar_array = output_values(0).array_value();
       jcval_array = output_values(1).array_value();

       for (i = 0; i < *jcnnz; i++) {
	 jcvar[i] = (int) jcvar_array.elem(i);
	 jcval[i] = jcval_array.elem(i);
       }
     }
   }

   }

/* ******************************************************************
   ****************************************************************** */

   void evalhc(int n,double* x,int ind,int* hclin,int* hccol,double* hcval,
   int* hcnnz,int* flag) {

   int i;

   octave_value_list input_arguments, output_values;
   NDArray hclin_array, hccol_array, hcval_array;
   NDArray x_array;

   x_array.resize_no_fill(n);
   for (i=0; i < n; ++i){
     x_array(i) = x[i];
   }

   input_arguments(0) = n;
   input_arguments(1) = x_array;
   input_arguments(2) = ind;

   output_values = feval("evalhc",input_arguments,3);

   *flag = output_values(4).int_value();

   if (*flag == 0 ){
     *hcnnz = output_values(3).int_value();
     if (*hcnnz > 0 ){
       hclin_array = output_values(0).array_value();
       hccol_array = output_values(1).array_value();
       hcval_array = output_values(2).array_value();
       for (i = 0; i < *hcnnz; i++) {
	 hclin[i] = (int) hclin_array.elem(i);
	 hccol[i] = (int) hccol_array.elem(i);
	 hcval[i] = hcval_array.elem(i);
       }
     }
   }

   }

/* ******************************************************************
   ****************************************************************** */

   void evalfc(int n,double* x,double* f,int m,double* constr,
   int* flag) {

   int i;
   octave_value_list input_arguments, output_values;
   NDArray x_array, constr_array;

   x_array.resize_no_fill(n);
   for (i=0; i < n; ++i){
     x_array(i) = x[i];
   }

   input_arguments(0) = n;
   input_arguments(1) = x_array;
   input_arguments(2) = m;

   output_values = feval("evalfc",input_arguments,3);

   *flag = output_values(2).int_value();

   if(*flag == 0){
     *f  = output_values(0).double_value();
     constr_array = output_values(1).array_value();
     for (i = 0; i < m; i++) {
       constr[i] = constr_array.elem(i);
     }
   }

   }

/* ******************************************************************
   ****************************************************************** */

   void evalgjac(int n,double *x,double *g,int m,int *jcfun,int *jcvar,
   double *jcval,int *jcnnz,int *flag) {

   int i;
   octave_value_list input_arguments, output_values;
   NDArray jcfun_array, jcvar_array, jcval_array;
   NDArray x_array, g_array;

   x_array.resize_no_fill(n);
   for (i=0; i < n; ++i){
     x_array(i) = x[i];
   }

   input_arguments(0) = n;
   input_arguments(1) = x_array;
   input_arguments(2) = m;

   output_values = feval("evalgjac",input_arguments,3);

   *flag = output_values(5).int_value();

   if (*flag == 0 ) {
     g_array = output_values(0).array_value();
     for (i = 0; i < n; i++)
       g[i] = g_array.elem(i);

     *jcnnz = output_values(4).int_value();

     if (*jcnnz > 0 ){
       jcfun_array = output_values(1).array_value();
       jcvar_array = output_values(2).array_value();
       jcval_array = output_values(3).array_value();

       for (i = 0; i < *jcnnz; i++) {
	 jcfun[i] = (int) jcfun_array.elem(i);
	 jcvar[i] = (int) jcvar_array.elem(i);
	 jcval[i] = jcval_array.elem(i);
       }
     }
   }

   }

/* ******************************************************************
   ****************************************************************** */

   void evalhl(int n,double *x,int m,double *lambda,double sf,
   double *sc,int *hllin,int *hlcol,double *hlval,int *hlnnz,
   int *flag) {

   int i;
   octave_value_list input_arguments, output_values;
   NDArray  x_array, lambda_array, sc_array;
   NDArray hllin_array, hlcol_array, hlval_array;

   x_array.resize_no_fill(n);
   for (i=0; i < n; ++i){
     x_array(i) = x[i];
   }

   lambda_array.resize_no_fill(m);
   sc_array.resize_no_fill(m);
   for (i=0; i < m; ++i){
     lambda_array(i) = lambda[i];
     sc_array(i) = sc[i];
   }

   input_arguments(0) = n;
   input_arguments(1) = x_array;
   input_arguments(2) = m;
   input_arguments(3) = lambda_array;
   input_arguments(4) = sf;
   input_arguments(5) = sc_array;

   output_values = feval("evalhl",input_arguments,6);

   *flag = output_values(4).int_value();

   if (*flag == 0 ){
     *hlnnz = output_values(3).int_value();
     if (*hlnnz > 0 ){
       hllin_array = output_values(0).array_value();
       hlcol_array = output_values(1).array_value();
       hlval_array = output_values(2).array_value();
       for (i = 0; i < *hlnnz; i++) {
	 hllin[i] = (int) hllin_array.elem(i);
	 hlcol[i] = (int) hlcol_array.elem(i);
	 hlval[i] = hlval_array.elem(i);
       }
     }
   }

   }

/* ******************************************************************
   ****************************************************************** */

   void evalhlp(int n,double *x,int m,double *lambda,double sf,
   double *sc, double *p, double *hp,int *gothl,int *flag) {

   int i;
   octave_value_list input_arguments, output_values;
   NDArray  x_array, lambda_array, sc_array, p_array, hp_array;

   x_array.resize_no_fill(n);
   p_array.resize_no_fill(n);
   hp_array.resize_no_fill(n);
   for (i=0; i < n; ++i){
     x_array(i) = x[i];
     p_array(i) = p[i];
     hp_array(i) = hp[i];
   }

   lambda_array.resize_no_fill(m);
   sc_array.resize_no_fill(m);
   for (i=0; i < m; ++i){
     lambda_array(i) = lambda[i];
     sc_array(i) = sc[i];
   }

   input_arguments(0) = n;
   input_arguments(1) = x_array;
   input_arguments(2) = m;
   input_arguments(3) = lambda_array;
   input_arguments(4) = sf;
   input_arguments(5) = sc_array;
   input_arguments(6) = p_array;
   input_arguments(7) = *gothl;

   output_values = feval("evalhlp",input_arguments,8);

   *flag   = output_values(2).int_value();

   if(*flag == 0){
     hp_array = output_values(0).array_value();
     for (i = 0; i < n; i++)
       hp[i] = hp_array.elem(i);

     *gothl  = output_values(1).int_value();
   }

   }

/* ******************************************************************
****************************************************************** */

   void endp (int n,double* x,double* l,double* u,int m,double* lambda,
   int* equatn, int* linear) {

   int i;
   octave_value_list input_arguments;
   NDArray x_array, l_array, u_array;
   NDArray lambda_array, equatn_array, linear_array;

   x_array.resize_no_fill(n);
   l_array.resize_no_fill(n);
   u_array.resize_no_fill(n);
   for (i=0; i < n; ++i){
     x_array(i) = x[i];
     l_array(i) = l[i];
     u_array(i) = u[i];
   }

   lambda_array.resize_no_fill(m);
   equatn_array.resize_no_fill(m);
   linear_array.resize_no_fill(m);
   for (i=0; i < m; ++i){
     lambda_array(i) = lambda[i];
     equatn_array(i) = equatn[i];
     linear_array(i) = linear[i];
   }

   input_arguments(0) =  n;
   input_arguments(1) =  x_array;
   input_arguments(2) =  l_array;
   input_arguments(3) =  u_array;
   input_arguments(4) =  m;
   input_arguments(5) =  lambda_array;
   input_arguments(6) =  equatn_array;
   input_arguments(7) =  linear_array;

   feval("endp",input_arguments,8);

   }

/* ******************************************************************
   ****************************************************************** */
   void param(double *epsfeas,double *epsopt,int *iprint,int *ncomp) {

   int i;
   NDArray rho_array;
   octave_value_list input_arguments, output_values;

   output_values = feval("param",input_arguments,0);
   printf("param\n");
   *epsfeas  = output_values(0).double_value();
   *epsopt   = output_values(1).double_value();
   *iprint   = output_values(2).int_value();
   *ncomp    = output_values(3).int_value();

   }

/* ******************************************************************
   ****************************************************************** */


void solver(int n,double *x,double *l,double *u,int m,double *lambda,
int *equatn,int *linear,int* coded,int checkder,double epsfeas,double epsopt,
int iprint,int ncomp,double *f,double* cnorm,double *snorm,double *nlpsupn,
int *inform){

/* TRANSFORMS THE TYPE INT OF C IN THE TYPE LOGICAL OF FORTRAN */
   C2FLOGICALV(equatn,m);
   C2FLOGICALV(linear,m);
   C2FLOGICALV(coded,10);

/* CALL OPTIMIZATION SOLVER */
   Algencan(epsfeas,epsopt,iprint,ncomp,n,x,l,u,m,lambda,equatn,
   linear,coded,checkder,*f,*cnorm,*snorm,*nlpsupn,*inform);

}

/* ******************************************************************
   ****************************************************************** */
DEFUN_DLD (algencan,arguments,nargout," This calls the ALGENCAN subroutine.\n \
You have to define the following octave functions before call it:\n \
\tevalf\n \tevalg\n \tevalh\n \tevalc\n \tevaljac\n \tevalhc\n \tevalhlp\n \
\tinip\n \tendp\n \tparam\n \
For more information go to http://www.ime.usp.br/~egbirgin/tango/\n")
{

   int checkder,inform,iprint,m,n,ncomp;
   double cnorm,f,nlpsupn,epsfeas,epsopt,snorm;

   int i;
   int coded[10];
   int *equatn,*linear;
   double *l,*lambda,*u,*x;

   NDArray l_array,lambda_array,coded_array,u_array,x_array, equatn_array,
           linear_array;
   octave_value_list return_value;

   n = arguments(0).int_value();
   m = arguments(4).int_value();

   x      = (double *) malloc(n * sizeof(double));
   l      = (double *) malloc(n * sizeof(double));
   u      = (double *) malloc(n * sizeof(double));
   lambda = (double *) malloc(m * sizeof(double));
   equatn = (int    *) malloc(m * sizeof(int   ));
   linear = (int    *) malloc(m * sizeof(int   ));

   x_array = arguments(1).array_value();
   l_array = arguments(2).array_value();
   u_array = arguments(3).array_value();
   for(i = 0; i < n; i++) {
     x[i] = x_array.elem(i);
     l[i] = l_array.elem(i);
     u[i] = u_array.elem(i);
   }

   lambda_array = arguments(5).array_value();
   equatn_array = arguments(6).array_value();
   linear_array = arguments(7).array_value();
   coded_array  = arguments(8).array_value();
   for(i = 0; i < m; i++) {
     lambda[i] = lambda_array.elem(i);
     equatn[i] = (int) equatn_array.elem(i);
     linear[i] = (int) linear_array.elem(i);
   }

   for(i = 0; i < 10; i++)
     coded[i] = (int) coded_array.elem(i);

   checkder  = arguments(9).int_value();
   epsfeas   = arguments(10).double_value();
   epsopt    = arguments(11).double_value();
   iprint    = arguments(12).int_value();
   ncomp     = arguments(13).int_value();

/* CALL OPTIMIZATION SOLVER */
   solver(n,x,l,u,m,lambda,equatn,linear,coded,checkder,epsfeas,epsopt,
   iprint,ncomp,&f,&cnorm,&snorm,&nlpsupn,&inform);

   x_array.resize_no_fill(n);
   for (i=0; i < n; ++i){
     x_array(i) = x[i];
   }

   lambda_array.resize_no_fill(m);
   for (i=0; i < m; ++i){
     lambda_array(i) = lambda[i];
   }

   return_value(0) = x_array;
   return_value(1) = lambda_array;
   return_value(2) = f;
   return_value(3) = cnorm;
   return_value(4) = snorm;
   return_value(5) = nlpsupn;
   return_value(6) = inform;

   return return_value;

}
