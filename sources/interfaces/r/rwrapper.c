/* ==================================================================
   Module: Interface between R and ALGENCAN
   ==================================================================

   Last update of any of the component of this module:

   September 10, 2008.

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

   R. Andreani, E. G. Birgin, J. M. Martínez and M. L. Schuverdt, "On
   Augmented Lagrangian methods with general lower-level constraints",
   SIAM Journal on Optimization 18, pp. 1286-1309, 2007.

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

   (See other related works at the TANGO home page.)

   ****************************************************************** */

/* ******************************************************************
   ****************************************************************** */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <R.h>
#include <Rdefines.h>
#include "../c/cfortran.h"
#include "rwrapper.h"

/* ******************************************************************
   ****************************************************************** */

   SEXP createRRealScalar(double x) {

   SEXP ans;

   PROTECT(ans = allocVector(REALSXP, 1));
   REAL(ans)[0] = x;
   UNPROTECT(1);

   return ans;

   }

/* ******************************************************************
   ****************************************************************** */

   SEXP createRRealVector(int size,double* x) {

   SEXP ans;
   int i;

   if(x == NULL){
     PROTECT(ans = allocVector(REALSXP, 1));
     REAL(ans)[0] = 0.0;
     UNPROTECT(1);
   }
   else{
     PROTECT(ans = allocVector(REALSXP, size));
     for (i = 0; i < size; i++)
       REAL(ans)[i] = x[i];
     UNPROTECT(1);
   }

   return ans;

   }

/* ******************************************************************
   ****************************************************************** */

   SEXP createRIntScalar(int x) {

   SEXP ans;

   PROTECT(ans = allocVector(INTSXP, 1));
   INTEGER(ans)[0] = x;
   UNPROTECT(1);

   return ans;

   }

/* ******************************************************************
   ****************************************************************** */

   SEXP createRIntVector(int size,int* x) {

   SEXP ans;
   int i;

   if(x == NULL){
     PROTECT(ans = allocVector(INTSXP, 1));
     INTEGER(ans)[0] = 0;
     UNPROTECT(1);
   }
   else{
     PROTECT(ans = allocVector(INTSXP, size));
     for (i = 0; i < size; i++)
       INTEGER(ans)[i] = x[i];
     UNPROTECT(1);
   }

   return ans;

   }

/* ******************************************************************
   ****************************************************************** */
   void inip (int* n,double** x,double** l, double** u,int* m,
	      double** lambda,int** equatn,int** linear,int* coded,
              int* checkder) {

   int i;

   SEXP n_r,m_r,x_r,l_r,u_r,lambda_r,rho_r,equatn_r,linear_r,coded_r,
        checkder_r;

   *n = 0;
   *m = 0;

   defineVar(install("x")       ,createRRealVector(*n,NULL) ,environment_r);
   defineVar(install("l")       ,createRRealVector(*n,NULL) ,environment_r);
   defineVar(install("u")       ,createRRealVector(*n,NULL) ,environment_r);
   defineVar(install("lambda")  ,createRRealVector(*m,NULL) ,environment_r);
   defineVar(install("equatn")  ,createRIntVector(*m,NULL)  ,environment_r);
   defineVar(install("linear")  ,createRIntVector(*m,NULL)  ,environment_r);
   defineVar(install("coded")   ,createRIntVector(10,NULL)  ,environment_r);
   defineVar(install("checkder"),createRIntScalar(*checkder),environment_r);

   EVAL(inip_r);

   n_r        = findVar(install("n")       ,environment_r);
   x_r        = findVar(install("x")       ,environment_r);
   l_r        = findVar(install("l")       ,environment_r);
   u_r        = findVar(install("u")       ,environment_r);
   m_r        = findVar(install("m")       ,environment_r);
   lambda_r   = findVar(install("lambda")  ,environment_r);
   equatn_r   = findVar(install("equatn")  ,environment_r);
   linear_r   = findVar(install("linear")  ,environment_r);
   coded_r    = findVar(install("coded")   ,environment_r);
   checkder_r = findVar(install("checkder"),environment_r);

   *n        = (INTEGER(AS_INTEGER(EVAL(n_r))))[0];
   *m        = (INTEGER(AS_INTEGER(EVAL(m_r))))[0];
   *checkder = (INTEGER(AS_INTEGER(EVAL(checkder_r))))[0];

   *x      = (double *) malloc(*n * sizeof(double));
   *l      = (double *) malloc(*n * sizeof(double));
   *u      = (double *) malloc(*n * sizeof(double));
   *lambda = (double *) malloc(*m * sizeof(double));
   *equatn = (int    *) malloc(*m * sizeof(int   ));
   *linear = (int    *) malloc(*m * sizeof(int   ));

   for(i = 0; i < *n; i++) {
     (*x)[i] = (REAL(EVAL(x_r)))[i];
     (*l)[i] = (REAL(EVAL(l_r)))[i];
     (*u)[i] = (REAL(EVAL(u_r)))[i];
   }

   for(i = 0; i < *m; i++) {
     (*lambda)[i] = (REAL(EVAL(lambda_r)))[i];
     (*equatn)[i] = (INTEGER(AS_INTEGER(EVAL(equatn_r))))[i];
     (*linear)[i] = (INTEGER(AS_INTEGER(EVAL(linear_r))))[i];
   }

   for(i = 0; i < 10; i++) {
     coded[i] = (INTEGER(AS_INTEGER(EVAL(coded_r))))[i];
   }

   }

/* ******************************************************************
   ****************************************************************** */

   void evalf(int n,double* x,double* f,int* flag) {

   SEXP flag_r,f_r;

   defineVar(install("n"),createRIntScalar(n)   ,environment_r);
   defineVar(install("x"),createRRealVector(n,x),environment_r);

   EVAL(evalf_r);

   f_r    = findVar(install("f"),   environment_r);
   flag_r = findVar(install("flag"),environment_r);

   *f    = (REAL(EVAL(f_r)))[0];
   *flag = (INTEGER(AS_INTEGER(EVAL(flag_r))))[0];

   }

/* ******************************************************************
   ****************************************************************** */

   void evalg(int n,double* x,double* g,int* flag) {

   int i;
   SEXP flag_r,g_r;

   defineVar(install("n"),createRIntScalar(n)      ,environment_r);
   defineVar(install("x"),createRRealVector(n,x)   ,environment_r);
   defineVar(install("g"),createRRealVector(n,NULL),environment_r);

   EVAL(evalg_r);

   g_r    = findVar(install("g"),   environment_r);
   flag_r = findVar(install("flag"),environment_r);

   for (i = 0; i < n; i++)
     g[i] = (REAL(EVAL(g_r)))[i];

   *flag = (INTEGER(AS_INTEGER(EVAL(flag_r))))[0];

   }

/* ******************************************************************
   ****************************************************************** */

   void evalh(int n,double* x,int* hlin,int* hcol,double* hval,int* hnnz,
   int* flag) {

   int i;
   double zero[1] = {0};
   SEXP hlin_r,hcol_r,hval_r,hnnz_r,flag_r;

   defineVar(install("n")   ,createRIntScalar(n)      ,environment_r);
   defineVar(install("x")   ,createRRealVector(n,x)   ,environment_r);
   defineVar(install("hlin"),createRIntVector(1,NULL) ,environment_r);
   defineVar(install("hcol"),createRIntVector(1,NULL) ,environment_r);
   defineVar(install("hval"),createRRealVector(1,NULL),environment_r);

   EVAL(evalh_r);

   hnnz_r = findVar(install("hnnz"),environment_r);
   flag_r = findVar(install("flag"),environment_r);
   hlin_r = findVar(install("hlin"),environment_r);
   hcol_r = findVar(install("hcol"),environment_r);
   hval_r = findVar(install("hval"),environment_r);

   *hnnz = (INTEGER(AS_INTEGER(EVAL(hnnz_r))))[0];

   for (i = 0; i < *hnnz; i++) {
     hlin[i] = (INTEGER(AS_INTEGER(EVAL(hlin_r))))[i];
     hcol[i] = (INTEGER(AS_INTEGER(EVAL(hcol_r))))[i];
     hval[i] = (REAL(EVAL(hval_r)))[i];
   }

   *flag = (INTEGER(AS_INTEGER(EVAL(flag_r))))[0];

   }

/* ******************************************************************
   ****************************************************************** */

   void evalc(int n,double* x,int ind,double* cind,int* flag) {

   SEXP flag_r,cind_r;

   defineVar(install("n")  ,createRIntScalar(n)   ,environment_r);
   defineVar(install("x")  ,createRRealVector(n,x),environment_r);
   defineVar(install("ind"),createRIntScalar(ind) ,environment_r);

   EVAL(evalc_r);

   cind_r = findVar(install("cind"),environment_r);
   flag_r = findVar(install("flag"),environment_r);

   *cind = (REAL(EVAL(cind_r)))[0];
   *flag = (INTEGER(AS_INTEGER(EVAL(flag_r))))[0];

   }

/* ******************************************************************
   ****************************************************************** */

   void evaljac(int n,double* x,int ind,int* jcvar,double* jcval,int* jcnnz,
   int* flag) {

   int i;
   double zero[1] = {0};
   SEXP flag_r,jcvar_r,jcval_r,jcnnz_r;

   defineVar(install("n")     ,createRIntScalar(n)      ,environment_r);
   defineVar(install("x")     ,createRRealVector(n,x)   ,environment_r);
   defineVar(install("ind")   ,createRIntScalar(ind)    ,environment_r);
   defineVar(install("jcvar") ,createRIntVector(1,NULL) ,environment_r);
   defineVar(install("jcval") ,createRRealVector(1,NULL),environment_r);

   EVAL(evaljac_r);

   jcnnz_r = findVar(install("jcnnz"),environment_r);
   jcvar_r = findVar(install("jcvar"),environment_r);
   jcval_r = findVar(install("jcval"),environment_r);
   flag_r  = findVar(install("flag")  ,environment_r);

   *jcnnz = (INTEGER(AS_INTEGER(EVAL(jcnnz_r))))[0];

   for (i = 0; i < *jcnnz; i++) {
     jcvar[i] = (INTEGER(AS_INTEGER(EVAL(jcvar_r))))[i];
     jcval[i] = (REAL(EVAL(jcval_r)))[i];
   }

   *flag = (INTEGER(AS_INTEGER(EVAL(flag_r))))[0];

   }

/* ******************************************************************
   ****************************************************************** */

   void evalhc(int n,double* x,int ind,int* hclin,int* hccol,double* hcval,
   int* hcnnz,int* flag) {

   int i;
   double zero[1] = {0};

   SEXP hclin_r,hccol_r,hcval_r,hcnnz_r,flag_r;

   defineVar(install("n")    ,createRIntScalar(n)      ,environment_r);
   defineVar(install("x")    ,createRRealVector(n,x)   ,environment_r);
   defineVar(install("ind")  ,createRIntScalar(ind)    ,environment_r);
   defineVar(install("hclin"),createRIntVector(1,NULL) ,environment_r);
   defineVar(install("hccol"),createRIntVector(1,NULL) ,environment_r);
   defineVar(install("hcval"),createRRealVector(1,NULL),environment_r);

   EVAL(evalhc_r);

   hcnnz_r = findVar(install("hcnnz"),environment_r);
   hclin_r = findVar(install("hclin"),environment_r);
   hccol_r = findVar(install("hccol"),environment_r);
   hcval_r = findVar(install("hcval"),environment_r);
   flag_r  = findVar(install("flag") ,environment_r);

   *hcnnz = (INTEGER(AS_INTEGER(EVAL(hcnnz_r))))[0];

   for (i = 0; i < *hcnnz; i++) {
     hclin[i] = (INTEGER(AS_INTEGER(EVAL(hclin_r))))[i];
     hccol[i] = (INTEGER(AS_INTEGER(EVAL(hccol_r))))[i];
     hcval[i] = (REAL(EVAL(hcval_r)))[i];
   }

   *flag = (INTEGER(AS_INTEGER(EVAL(flag_r))))[0];

   }

/* ******************************************************************
   ****************************************************************** */

   void evalfc(int n,double* x,double* f,int m,double* constr,
   int* flag) {

   int i;
   SEXP flag_r,constr_r,f_r;

   defineVar(install("n"),createRIntScalar(n)           ,environment_r);
   defineVar(install("x"),createRRealVector(n,x)        ,environment_r);
   defineVar(install("m"),createRIntScalar(m)           ,environment_r);
   defineVar(install("constr"),createRRealVector(1,NULL),environment_r);

   EVAL(evalfc_r);

   f_r      = findVar(install("f")   ,environment_r);
   constr_r = findVar(install("constr"),environment_r);
   flag_r   = findVar(install("flag"),environment_r);

   *f = (REAL(EVAL(f_r)))[0];
   for (i = 0; i < m; i++) {
     constr[i] = (REAL(EVAL(constr_r)))[i];
   }

   *flag = (INTEGER(AS_INTEGER(EVAL(flag_r))))[0];

   }

/* ******************************************************************
   ****************************************************************** */

   void evalgjac(int n,double *x,double *g,int m,int *jcfun,int *jcvar,
   double *jcval,int *jcnnz,int *flag) {

   int i;

   SEXP g_r,jcfun_r,jcvar_r,jcval_r,jcnnz_r,flag_r;

   defineVar(install("n")    ,createRIntScalar(n)      ,environment_r);
   defineVar(install("x")    ,createRRealVector(n,x)   ,environment_r);
   defineVar(install("g")    ,createRRealVector(n,x)   ,environment_r);
   defineVar(install("m")    ,createRIntScalar(m)      ,environment_r);
   defineVar(install("jcfun"),createRIntVector(1,NULL) ,environment_r);
   defineVar(install("jcvar"),createRIntVector(1,NULL) ,environment_r);
   defineVar(install("jcval"),createRRealVector(1,NULL),environment_r);
   defineVar(install("jcnnz"),createRIntScalar(0)      ,environment_r);

   EVAL(evalgjac_r);

   g_r     = findVar(install("g")    ,environment_r);
   jcnnz_r = findVar(install("jcnnz"),environment_r);
   jcfun_r = findVar(install("jcfun"),environment_r);
   jcvar_r = findVar(install("jcvar"),environment_r);
   jcval_r = findVar(install("jcval"),environment_r);
   flag_r  = findVar(install("flag") ,environment_r);

   *jcnnz = (INTEGER(AS_INTEGER(EVAL(jcnnz_r))))[0];

   for (i = 0; i < n; i++)
     g[i] = (REAL(EVAL(g_r)))[i];

   for (i = 0; i < *jcnnz; i++) {
     jcfun[i] = (INTEGER(AS_INTEGER(EVAL(jcfun_r))))[i];
     jcvar[i] = (INTEGER(AS_INTEGER(EVAL(jcvar_r))))[i];;
     jcval[i] = (REAL(EVAL(jcval_r)))[i];;
   }

   *flag = (INTEGER(AS_INTEGER(EVAL(flag_r))))[0];


   }

/* ******************************************************************
   ****************************************************************** */

   void evalhl(int n,double *x,int m,double *lambda,double sf,
   double *sc,int *hllin,int *hlcol,double *hlval,int *hlnnz,
   int *flag) {

   int i;
   SEXP hlnnz_r,hlcol_r,hllin_r,hlval_r,flag_r;

   defineVar(install("n")     ,createRIntScalar(n)        ,environment_r);
   defineVar(install("x")     ,createRRealVector(n,x)     ,environment_r);
   defineVar(install("m")     ,createRIntScalar(m)        ,environment_r);
   defineVar(install("lambda"),createRRealVector(m,lambda),environment_r);
   defineVar(install("sf")    ,createRRealScalar(sf)      ,environment_r);
   defineVar(install("sc")    ,createRRealVector(m,sc)    ,environment_r);
   defineVar(install("hllin") ,createRIntVector(1,NULL)   ,environment_r);
   defineVar(install("hlcol") ,createRIntVector(1,NULL)   ,environment_r);
   defineVar(install("hlval") ,createRRealVector(1,NULL)  ,environment_r);

   EVAL(evalhl_r);

   hlnnz_r = findVar(install("hlnnz"),environment_r);
   hllin_r = findVar(install("hllin"),environment_r);
   hlcol_r = findVar(install("hlcol"),environment_r);
   hlval_r = findVar(install("hlval"),environment_r);
   flag_r  = findVar(install("flag") ,environment_r);

   *hlnnz = (INTEGER(AS_INTEGER(EVAL(hlnnz_r))))[0];

   for (i = 0; i < *hlnnz; i++){
     hllin[i] = (INTEGER(AS_INTEGER(EVAL(hllin_r))))[i];
     hlcol[i] = (INTEGER(AS_INTEGER(EVAL(hlcol_r))))[i];
     hlval[i] = (REAL(EVAL(hlval_r)))[i];
   }

   *flag = (INTEGER(AS_INTEGER(EVAL(flag_r))))[0];

   }

/* ******************************************************************
   ****************************************************************** */

   void evalhlp(int n,double *x,int m,double *lambda,double sf,
   double *sc,double *p,double *hp,int *gothl,int *flag) {

   int i;
   SEXP hp_r,gothl_r,flag_r;

   defineVar(install("n")     ,createRIntScalar(n)        ,environment_r);
   defineVar(install("x")     ,createRRealVector(n,x)     ,environment_r);
   defineVar(install("m")     ,createRIntScalar(m)        ,environment_r);
   defineVar(install("lambda"),createRRealVector(m,lambda),environment_r);
   defineVar(install("sf")    ,createRRealScalar(sf)      ,environment_r);
   defineVar(install("sc")    ,createRRealVector(m,sc)    ,environment_r);
   defineVar(install("p")     ,createRRealVector(n,p)     ,environment_r);
   defineVar(install("hp")    ,createRRealVector(n,hp)    ,environment_r);
   defineVar(install("gothl") ,createRIntScalar(*gothl)   ,environment_r);

   EVAL(evalhlp_r);

   hp_r    = findVar(install("hp")    ,environment_r);
   gothl_r = findVar(install("gothl") ,environment_r);
   flag_r  = findVar(install("flag")  ,environment_r);

   for (i = 0; i < n; i++)
     hp[i] = (REAL(EVAL(hp_r)))[i];

   *gothl = (INTEGER(AS_INTEGER(EVAL(gothl_r))))[0];
   *flag  = (INTEGER(AS_INTEGER(EVAL(flag_r))))[0];

   }

/* ******************************************************************
   ****************************************************************** */

   void endp (int n,double* x,double* l,double* u,int m,double* lambda,
   int* equatn, int* linear) {

   defineVar(install("n")     ,createRIntScalar(n)        ,environment_r);
   defineVar(install("x")     ,createRRealVector(n,x)     ,environment_r);
   defineVar(install("l")     ,createRRealVector(n,l)     ,environment_r);
   defineVar(install("u")     ,createRRealVector(n,u)     ,environment_r);
   defineVar(install("m")     ,createRIntScalar(m)        ,environment_r);
   defineVar(install("lambda"),createRRealVector(m,lambda),environment_r);
   defineVar(install("equatn"),createRIntVector(m,equatn) ,environment_r);
   defineVar(install("linear"),createRIntVector(m,linear) ,environment_r);

   EVAL(endp_r);

   free(x);
   free(l);
   free(u);
   free(lambda);
   free(equatn);
   free(linear);

   x      = NULL;
   l      = NULL;
   u      = NULL;
   lambda = NULL;
   equatn = NULL;
   linear = NULL;

   }

/* ******************************************************************
   ****************************************************************** */

   void param(double *epsfeas,double *epsopt,int *iprint,int *ncomp) {
   int i;

   SEXP epsfeas_r,epsopt_r,iprint_r,ncomp_r;

   EVAL(param_r);

   epsfeas_r = findVar(install("epsfeas") ,environment_r);
   epsopt_r  = findVar(install("epsopt")  ,environment_r);
   iprint_r  = findVar(install("iprint")  ,environment_r);
   ncomp_r   = findVar(install("ncomp")   ,environment_r);

   *epsfeas  = (REAL(EVAL(epsfeas_r)))[0];
   *epsopt   = (REAL(EVAL(epsopt_r)))[0];
   *iprint   = (INTEGER(AS_INTEGER(EVAL(iprint_r))))[0];
   *ncomp    = (INTEGER(AS_INTEGER(EVAL(ncomp_r))))[0];

   }

/* ******************************************************************
   ****************************************************************** */

   SEXP ralgencan(SEXP evalf_ptr,SEXP evalg_ptr,SEXP evalh_ptr,
   SEXP evalc_ptr,SEXP evaljac_ptr,SEXP evalhc_ptr,SEXP evalfc_ptr,
   SEXP evalgjac_ptr, SEXP evalhl_ptr, SEXP evalhlp_ptr,
   SEXP inip_ptr, SEXP endp_ptr,SEXP param_ptr,
   SEXP environment_ptr) {

   int checkder,inform,iprint,m,n,ncomp;
   double cnorm,f,nlpsupn,epsfeas,epsopt,snorm;

   int coded[10];
   int *equatn,*linear;
   double *l,*lambda,*u,*x;
   SEXP return_value;

   evalf_r       = evalf_ptr;
   evalg_r       = evalg_ptr;
   evalh_r       = evalh_ptr;
   evalc_r       = evalc_ptr;
   evaljac_r     = evaljac_ptr;
   evalhc_r      = evalhc_ptr;
   evalfc_r      = evalfc_ptr;
   evalgjac_r    = evalgjac_ptr;
   evalhl_r      = evalhl_ptr;
   evalhlp_r     = evalhlp_ptr;
   inip_r        = inip_ptr;
   endp_r        = endp_ptr;
   param_r       = param_ptr;
   environment_r = environment_ptr;

/* SET SOME SOLVER ARGUMENTS */
   param(&epsfeas,&epsopt,&iprint,&ncomp);

/* SET UP PROBLEM DATA */
   inip(&n,&x,&l,&u,&m,&lambda,&equatn,&linear,coded,&checkder);

   C2FLOGICALV(equatn,m);
   C2FLOGICALV(linear,m);
   C2FLOGICALV(coded,10);

   Algencan(epsfeas,epsopt,iprint,ncomp,n,x,l,u,m,lambda,equatn,
   linear,coded,checkder,f,cnorm,snorm,nlpsupn,inform);

/* WRITE ADDITIONAL OUTPUT INFORMATION CODED BY THE USER */
   endp(n,x,l,u,m,lambda,equatn,linear);


   defineVar(install("AlgencanReturnValue"),createRIntScalar(0),
   environment_r);
   defineVar(install("f"),createRRealScalar(f),environment_r);
   defineVar(install("cnorm"),createRRealScalar(cnorm),environment_r);
   defineVar(install("snorm"),createRRealScalar(snorm),environment_r);
   defineVar(install("nlpsupn"),createRRealScalar(nlpsupn),environment_r);
   defineVar(install("inform"),createRIntScalar(inform),environment_r);

   return_value = findVar(install("AlgencanReturnValue"),environment_r);

   return return_value;

   }
