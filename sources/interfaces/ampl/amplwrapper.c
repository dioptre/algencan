/* =================================================================
   File: amplwrapper.c
   =================================================================

   =================================================================
   Module: AMPL Wrapper
   =================================================================

   Last update of any of the component of this module:

   March 4, 2011. (Previous one: December 17th, 2007.)

   Users are encouraged to download periodically updated versions of
   this code at the TANGO home page:

   http://www.ime.usp.br/~egbirgin/tango/

   *****************************************************************
   *****************************************************************

   TANGO Project

   License:

   All the TANGO Project components are free software; you can
   redistribute it and/or modify it under the terms of the GNU General
   Public License as published by the Free Software Foundation; either
   version 2 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA. You can also find the GPL on the GNU web site.

   Non-free versions of TANGO are available under terms different from
   those of the General Public License. Professors J. M. Martínez
   (martinez@ime.unicamp.br, martinezimecc@gmail.com) or E. G. Birgin
   (egbirgin@ime.usp.br, egbirgin@gmail.com) should be contacted for
   more information related to such a license, future developments
   and/or technical support.

   In addition, we kindly ask you to acknowledge the TANGO Project and
   its authors in any program or publication in which you use of the
   TANGO Project components. For published works that use ALGENCAN we
   suggest referencing:

   R. Andreani, E. G. Birgin, J. M. Martínez and M. L. Schuverdt, "On
   Augmented Lagrangian methods with general lower-level constraints",
   SIAM Journal on Optimization 18, pp. 1286-1309, 2007.

   and

   R. Andreani, E. G. Birgin, J. M. Martínez and M. L. Schuverdt,
   "Augmented Lagrangian methods under the Constant Positive Linear
   Dependence constraint qualification", Mathematical Programming 111,
   pp. 5-32, 2008

   For published works that use GENCAN we suggest referencing:

   E. G. Birgin and J. M. Martínez, "Large-scale active-set
   box-constrained optimization method with spectral projected
   gradients", Computational Optimization and Applications 23,
   pp. 101-125, 2002.

   For published works that use SPG we suggest referencing:

   E. G. Birgin, J. M. Martínez and M. Raydan, "Nonmonotone spectral
   projected gradient methods on convex sets", SIAM Journal on
   Optimization 10, pp. 1196-1211, 2000,

   and

   E. G. Birgin, J. M. Martínez and M. Raydan, "Algorithm 813: SPG -
   software for convex-constrained optimization", ACM Transactions on
   Mathematical Software 27, pp. 340-349, 2001.

   (See also other related works in the TANGO Project home page:

   http://www.ime.usp.br/~egbirgin/tango/)

   ***************************************************************** */

/* *****************************************************************
   ***************************************************************** */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "../c/cfortran.h"
#include "amplwrapper.h"
#include "asl.h"
#include "getstub.h"

#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))

ASL *asl;
int *cmap,*slaind,useslacks;
double *ca,*cb,objsign;
char *stub,message[200];

/* *****************************************************************
   ***************************************************************** */

void inip(int *n,double **x,double **l,double **u,int *m,double **lambda,
int **equatn,int **linear,int **coded,int *checkder) {

/* This subroutine set some problem data. */

   int flag,i,k,nnzh,nconstr,nslacks;
   double cl,cu;
   FILE *nl;

   strcpy( message, "" );

   /* Sparse jacobian option. */

   asl->i.congrd_mode = 1;

   /* Name of file that contains the description of the problem. */

   if ( stub == NULL ) {

     write_sol( "\nThe 'stub' file of the problem is missing.\n",
                NULL, NULL, NULL );
     exit( 0 );
   }

   nl = jac0dim(stub, (fint)strlen(stub));

   /* Initial point (if provided by ampl). */

   X0     = (real *) malloc( n_var * sizeof (real) );
   havex0 = (char *) malloc( n_var * sizeof (char) );

   if ( X0 == NULL || havex0 == NULL ) {

     write_sol( "\nAMPL INTERFACE ERROR: in inip, malloc returned NULL.\n",
                NULL, NULL, NULL );
     exit( 0 );
   }

   pfgh_read(nl,ASL_findgroups);

   if ( n_obj == 0 )
     nnzh = sphsetup(-1,0,1,1);
   else
     nnzh = sphsetup(-1,1,1,1);

   if ( n_obj > 0 && objtype[0] )
     objsign = - 1.0;
   else
     objsign =   1.0;

   /* Verifies if the problem has integer variables. */

   if ( nbv > 0 || niv > 0 || nlvci > 0 || nlvoi > 0 ) {

     write_sol( "\nALGENCAN can not handle integer variables.\n",
                NULL, NULL, NULL );
     exit( 0 );
   }

   /* Number of variables. */

   useslacks = 1;

   //printf("\nnranges= %d",nranges);
   //exit( 0 );

   *n = n_var + ( useslacks ? nranges : 0 );

   /* Number of constraints (equalities plus inequalities). */

   *m = n_con + ( useslacks ? 0 : nranges );

   /* Memory allocation. */

   *x      = (double *) malloc(((*n)+(*m)) * sizeof(double));
   *l      = (double *) malloc(((*n)+(*m)) * sizeof(double));
   *u      = (double *) malloc(((*n)+(*m)) * sizeof(double));
   *lambda = (double *) malloc((*m) * sizeof(double));
   *equatn = (int *)    malloc((*m) * sizeof(int));
   *linear = (int *)    malloc((*m) * sizeof(int));
   *coded  = (int *)    malloc( 10  * sizeof(int));

   cmap    = (int *)    malloc((*m) * sizeof(int));
   ca      = (double *) malloc((*m) * sizeof(double));
   cb      = (double *) malloc((*m) * sizeof(double));
   slaind  = (int *)    malloc((*m) * sizeof(int));

   if (      *x == NULL ||      *l == NULL ||      *u == NULL ||
        *lambda == NULL || *equatn == NULL || *linear == NULL ||
         *coded == NULL ||  slaind == NULL ||    cmap == NULL ||
             ca == NULL ||      cb == NULL ) {

     write_sol( "\nAMPL INTERFACE ERROR: in inip, malloc returned NULL.\n",
                NULL, NULL, NULL );
     exit( 0 );

   }

   /* Initial point. */

   for(i = 0; i < n_var; i++)
     if ( havex0[i] )
       (*x)[i] = X0[i];
     else
       (*x)[i] = 0.0;

   /* Lower and upper bounds. */

   for(i = 0; i < n_var; i++) {
     (*l)[i] = LUv[2 * i];
     (*u)[i] = LUv[2 * i + 1];
   }

  /* For each constraint i, set equatn(i) = 1. if it is an equality
     constraint of the form c_i(x) = 0, and set equatn(i) = 0 if it is
     an inequality constraint of the form c_i(x) <= 0. Moreover, set
     linear(i) = 1 if it is a linear constraint, otherwise set
     linear(i) = 0. */

   nconstr = 0;
   nslacks = 0;

   for(i = 0; i < n_con; i++) {

      cl = LUrhs[2 * i];
      cu = LUrhs[2 * i + 1];

      /* Equality constraint. */

      if ( cl == cu ) {
	 (*equatn)[nconstr] = 1;
	 (*linear)[nconstr] = (i < nlc ? 0 : 1);

	 slaind[nconstr]    = - 1;
         cmap[nconstr]      = i;
	 ca[nconstr]        = 1.0;
	 cb[nconstr]        = - cl;

	 nconstr++;
      }

      /* Ranged inequality constraint: add slack or split it. */

      else if ( cl > negInfinity && cu < Infinity ) {

	 /* Replace by c(x) - s = 0 and l <= s <= u. */

	 if ( useslacks ) {
	    (*equatn)[nconstr] = 1;
	    (*linear)[nconstr] = (i < nlc ? 0 : 1);

	    k = n_var + nslacks;

	    slaind[nconstr] = k;
	    cmap[nconstr]   = i;
	    ca[nconstr]     = 1.0;
	    cb[nconstr]     = 0.0;

	    (*l)[k] = cl;
	    (*u)[k] = cu;

	    (*x)[k] = conival(i, (real *)(*x), (fint *)&flag );

	    if ( flag != 0 ) {
	       sprintf(message,"AMPL INTERFACE ERROR: conival returned %d",(int)flag);
	       write_sol(message,NULL, NULL, NULL );
	       exit( 0 );
	    }

	    (*x)[k] = max(cl, min( (*x)[k], cu));

	    nconstr++;
	    nslacks++;
	 }

	 /* Split into c(x) - u <= 0 and l - c(x) <= 0. */

	 else {
	    (*equatn)[nconstr] = 0;
	    (*linear)[nconstr] = (i < nlc ? 0 : 1);

	    slaind[nconstr]    = - 1;
	    cmap[nconstr]      = i;
	    ca[nconstr]        = 1.0;
	    cb[nconstr]        = - cu;

	    nconstr++;

	    (*equatn)[nconstr] = 0;
	    (*linear)[nconstr] = (i < nlc ? 0 : 1);

	    slaind[nconstr]    = - 1;
	    cmap[nconstr]      = i;
	    ca[nconstr]        = - 1.0;
	    cb[nconstr]        = cl;

	    nconstr++;
	 }
      }

      /* Inequality constraint of type c(x) <= u. */

      else if ( cu < Infinity ) {
	 (*equatn)[nconstr] = 0;
	 (*linear)[nconstr] = (i < nlc ? 0 : 1);

	 slaind[nconstr]    = - 1;
	 cmap[nconstr]      = i;
	 ca[nconstr]        = 1.0;
	 cb[nconstr]        = - cu;

	 nconstr++;
      }
      /* Inequality constraint of type l <= c(x). */

      else if ( cl > negInfinity ) {
	 (*equatn)[nconstr] = 0;
	 (*linear)[nconstr] = (i < nlc ? 0 : 1);

	 slaind[nconstr]    = - 1;
	 cmap[nconstr]      = i;
	 ca[nconstr]        = - 1.0;
	 cb[nconstr]        = cl;

	 nconstr++;
      }
   }

   /* Lagrange multipliers approximation. */

   for( i = 0; i < *m; i++ )
       (*lambda)[i] = 0.0;

   /* In this AMPL interface evalf, evalg, evalc, evaljac and evalhl are
      present. evalh, evalhc, evalhlp, evalfc and evalgjac are not. */

   (*coded)[0] = 1; /* evalf    */
   (*coded)[1] = 1; /* evalg    */
   (*coded)[2] = 0; /* evalh    */
   (*coded)[3] = 1; /* evalc    */
   (*coded)[4] = 1; /* evaljac  */
   (*coded)[5] = 0; /* evalhc   */
   (*coded)[6] = 0; /* evalfc   */
   (*coded)[7] = 0; /* evalgjac */
   (*coded)[8] = 1; /* evalhl   */
   (*coded)[9] = 0; /* evalhlp  */

   *checkder = 0;
}

/* *****************************************************************
   ***************************************************************** */

void evalf(int n,double *x,double *f,int *flag) {

/* This subroutine computes the objective function. */

   *flag = 0;

   if ( n_obj == 0 )
     *f = 0;

   else {
     *f = objval(0, (real *)x, (fint *)flag);

     if ( *flag != 0 ) {
       sprintf( message, "AMPL INTERFACE ERROR: objval returned %d",*flag );
       //return;
     }

     *f *= objsign;
   }

}

/* *****************************************************************
   ***************************************************************** */

void evalg(int n,double *x,double *g,int *flag) {

/* This subroutine computes the gradient vector of the objective
   function. */

   int i;

   *flag = 0;

   if ( n_obj == 0 ) {
     for ( i = 0; i < n; i++ )
       g[i] = 0;
   }

   else {
     objgrd(0, (real *)x, (real *)g, (fint *)flag);

     if ( *flag != 0 ) {
       sprintf( message, "AMPL INTERFACE ERROR: objgrd returned %d",*flag);
       //return;
     }

     for ( i = 0; i < n_var; i++ )
	g[i] *= objsign;

     for ( i = n_var; i < n; i++ )
	g[i] = 0.0;
   }
}

/* *****************************************************************
   ***************************************************************** */

void evalh(int n,double *x,int *hlin,int *hcol,double *hval,
int *hnnz,int *flag) {

/* This subroutine might compute the Hessian matrix of the objective
   function. */

   *flag = -1;
}

/* ******************************************************************
   ****************************************************************** */

void evalc(int n,double *x,int ind,double *c,int *flag) {

/* This subroutine computes the ind-th constraint. */

   *flag = 0;

   *c = conival( cmap[ind - 1], (real *)x, (fint *)flag );

   if ( *flag != 0 ) {
     sprintf( message, "AMPL INTERFACE ERROR: conival returned %d",*flag );
     //return;
   }

   if ( slaind[ind - 1] == -1 )
      *c = *c * ca[ind - 1] + cb[ind - 1];

   else
      *c = *c - x[ slaind[ind - 1] ];
}

/* *****************************************************************
   ***************************************************************** */

void evaljac(int n,double *x,int ind,int *jcvar,double *jcval,
int *jcnnz,int *flag) {

/* This subroutine computes the sparse gradient of the ind-th
   constraint. */

   struct cgrad *og;

   *flag = 0;

   congrd( cmap[ind - 1], (real *)x, (real *)jcval, (fint *)flag );

   if ( *flag != 0 ) {
      sprintf( message, "AMPL INTERFACE ERROR: congrd returned %d",*flag );
      //return;
   }

   *jcnnz = 0;
   for(og = Cgrad[ cmap[ind - 1] ]; og != NULL; og = og->next) {
     jcvar[*jcnnz]  = og->varno + 1;
     jcval[*jcnnz] *= ca[ind - 1];
     (*jcnnz)++;
   }

   if ( slaind[ind - 1] != -1 ) {
      jcvar[*jcnnz] = slaind[ind - 1] + 1;
      jcval[*jcnnz] = - 1.0;
      (*jcnnz)++;
   }
}

/* *****************************************************************
   ***************************************************************** */

void evalhc(int n,double *x,int ind,int *hclin,int *hccol,
double *hcval,int *hcnnz,int *flag) {

   *flag = -1;
}

/* *****************************************************************
   ***************************************************************** */

void evalhl(int n,double *x,int m,double *lambda,double scalef,
double *scalec,int *hllin,int *hlcol,double *hlval,int *hlnnz,
int *flag) {

/* This subroutine computes the sparse Hessian of the Lagrangian. */

   int i,j;
   double *coeff,*OW;

   *flag = 0;

   OW    = (real *) calloc( n_obj, sizeof (real) );
   coeff = (real *) calloc( n_con, sizeof (real) );

   if ( OW == NULL || coeff == NULL ) {

     write_sol( "\nAMPL INTERFACE ERROR: in evalhl, calloc returned NULL.\n",
                NULL, NULL, NULL );
     exit( 0 );

   }

   OW[0] = objsign * scalef;

   for ( i = 0; i < m; i++ )
     coeff[cmap[i]] += ca[i] * lambda[i] * scalec[i];

   if ( n_obj == 0 )
     sphes(hlval,-1,NULL,coeff);
   else
     sphes(hlval,-1,OW,coeff);

   for ( j = 0; j < n_var; j++ ) {
     for ( i = sputinfo->hcolstarts[j]; i < sputinfo->hcolstarts[j+1]; i++ ) {

       /* The AMPL subroutine returns the upper triangle (as lists of columns
          with elements from the first row to the diagonal). Invert indices to
          return the lower triangle. */

       hllin[i] = j + 1;
       hlcol[i] = sputinfo->hrownos[i] + 1;
     }
   }

   *hlnnz = sputinfo->hcolstarts[n_var];

   free( OW );
   free( coeff );
}

/* *****************************************************************
   ***************************************************************** */

void evalhlp(int n,double *x,int m,double *lambda,double scalef,
double *scalec,double *p,double *hp,int *goth,int *flag) {

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

void endp(int n,double *x,double *l,double *u,int m,double *lambda,
int *equatn,int *linear) {

/* This subroutine can be used to do some extra job after the solver
   has found the solution. */

   int i;
   double * ampl_lambda;

   ampl_lambda = (real *) calloc( n_con, sizeof (real) );

   if ( ampl_lambda == NULL) {
     write_sol( "\nAMPL INTERFACE ERROR: in endp, calloc returned NULL.\n",
                NULL, NULL, NULL );
     exit( 0 );
   }

   for ( i = 0; i < m; i++ )
     ampl_lambda[cmap[i]] += ca[i] * lambda[i];

   write_sol( message, x, ampl_lambda, NULL );

   free(ampl_lambda);
}

/* *****************************************************************
   ***************************************************************** */

int main(int argc,char *argv[]) {

   int checkder,inform,m,n,iprint,ncomp;
   double f,nlpsupn,epsfeas,epsopt,cnorm,snorm;

   int *coded,*equatn,*linear;
   double *l,*lambda,*u,*x;

   asl = ASL_alloc(ASL_read_pfgh);

/* SET SOME SOLVER ARGUMENTS */
   param(&epsfeas,&epsopt,&iprint,&ncomp);

/* GET THE OPTIONS AND STUB'S FILE NAME. */
   getOptions(argv,&epsfeas,&epsopt,&iprint,&ncomp);

/* SET UP PROBLEM DATA */
   inip(&n,&x,&l,&u,&m,&lambda,&equatn,&linear,&coded,&checkder);

/* TRANSFORMS THE TYPE INT OF C IN THE TYPE LOGICAL OF FORTRAN
 * FOR THE ARRAYS equatn AND linear, AND VARIABLE checkder */
   C2FLOGICALV(equatn,m);
   C2FLOGICALV(linear,m);
   C2FLOGICALV(coded,10);

/* CALL ALGENCAN */
   algencan(epsfeas,epsopt,iprint,ncomp,n,x,l,u,m,lambda,equatn,
   linear,coded,checkder,f,cnorm,snorm,nlpsupn,inform);

/* WRITE ADDITIONAL OUTPUT INFORMATION CODED BY THE USER */
   endp(n,x,l,u,m,lambda,equatn,linear);

/* FREEING MEMORY */
   memFree(&coded,&equatn,&linear,&l,&lambda,&u,&x);

   return 0;
}

/* *****************************************************************
   ***************************************************************** */

void param(double *epsfeas,double *epsopt,int *iprint,int *ncomp) {

/* This subroutine sets some algencan arguments related to
   stopping criteria and output. */

   *epsfeas  = 1.0e-08;
   *epsopt   = 1.0e-08;
   *iprint   = 10;
   *ncomp    = 6;
}

/* *****************************************************************
   ***************************************************************** */

void setp(int n,double * x) {

   xknown((real *)x);
}

/* *****************************************************************
   ***************************************************************** */

void unsetp() {

   xunknown();
}

/* *****************************************************************
   ***************************************************************** */

void memFree(int **coded,int **equatn,int **linear,double **l,
double **lambda,double **u,double **x) {

   free(*coded );
   free(*equatn);
   free(*linear);
   free(*l     );
   free(*lambda);
   free(*u     );
   free(*x     );
   free(cmap   );
   free(ca     );
   free(cb     );
   free(slaind );
}

/* *****************************************************************
   ***************************************************************** */

void getOptions(char *argv[],double *epsfeas,double *epsopt,int *iprint,
int *ncomp) {

  char* epsfeasDescription  = "DOUBLE PRECISION\n\n"
"Feasibility tolerance for the sup-norm of the constraints. (Ignored in the     \n"
"unconstrained and bound-constrained cases.)                                    \n";

  char* epsoptDescription = "DOUBLE PRECISION\n\n"
"Optimality tolerance for the sup-norm of the projected gradient of the         \n"
"Lagrangian in the constrained case and the sup-norm of the projected           \n"
"gradient of the objective function in the unconstrained and the bound-         \n"
"constrained cases.                                                             \n";

  char* iprintDescription = "INTEGER\n\n"
"Controls the ammount of information of the output.                             \n";

  char* ncompDescription = "INTEGER\n\n"
"Every time a vector is printed, just its first ncomp component will be         \n"
"displayed.                                                                     \n";

  /* Must be in alphabetical order. */
  keyword keywds[] = {
    KW("epsfeas",D_val,epsfeas,epsfeasDescription),

    KW("epsopt",D_val,epsopt,epsoptDescription),

    KW("iprint",I_val,iprint,iprintDescription),

    KW("ncomp",I_val,ncomp,ncompDescription)
  };

  Option_Info Oinfo = {"algencan","ALGENCAN","algencan_options",
                       keywds,nkeywds,1,"October, 2007",NULL,NULL,NULL,
                       NULL,20071022};

  /* File name that contains the description of the problem. */
  stub = getstops(argv,&Oinfo);

}
