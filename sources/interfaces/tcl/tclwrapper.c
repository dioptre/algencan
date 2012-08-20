/* =================================================================
   File: tclwrapper.c
   =================================================================

   =================================================================
   Module: TCL Wrapper
   =================================================================

   Last update of any of the component of this module:

   April 28th, 2010.

   Users are encouraged to download periodically updated versions of
   this code at the TANGO home page:

   http://www.ime.usp.br/~egbirgin/tango/

   =================================================================
   This interface was developed by

   Juan José Miranda Bront,
   Pablo Matías Factorovich,
   Daniel Negrotto,
   Agustín Pecorari, and
   Daniel Esteban Severin,

   and modified by

   Francisco N. C. Sobral.
   =================================================================

   *****************************************************************
   *****************************************************************

   TANGO Project
   
   License:
   
   All the TANGO Project components are free software; you can redistribute
   it and/or modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or (at your option) any later version.
   
   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
   Public License for more details.
   
   You should have received a copy of the GNU General Public License along
   with this program; if not, write to the Free Software Foundation, Inc.,
   59 Temple Place, Suite 330, Boston, MA 02111-1307 USA. You can also find
   the GPL on the GNU web site.
   
   Non-free versions of TANGO are available under terms different from
   those of the General Public License. Professors J. M. Martínez
   (martinez@ime.unicamp.br, martinezimecc@gmail.com) or E. G. Birgin
   (egbirgin@ime.usp.br, egbirgin@gmail.com) should be contacted for more
   information related to such a license, future developments and/or
   technical support.
   
   In addition, we kindly ask you to acknowledge the TANGO Project and its
   authors in any program or publication in which you use of the TANGO
   Project components. For published works that use ALGENCAN we suggest
   referencing:
   
   R. Andreani, E. G. Birgin, J. M. Martínez and M. L. Schuverdt, "On
   Augmented Lagrangian methods with general lower-level constraints",
   SIAM Journal on Optimization 18, pp. 1286-1309, 2007
   
   and
   
   R. Andreani, E. G. Birgin, J. M. Martínez and M. L. Schuverdt,
   "Augmented Lagrangian methods under the Constant Positive Linear
   Dependence constraint qualification", Mathematical Programming 111,
   pp. 5-32, 2008
   
   For published works that use GENCAN we suggest referencing:
   
   E. G. Birgin and J. M. Martínez, "Large-scale active-set
   box-constrained optimization method with spectral projected gradients",
   Computational Optimization and Applications 23, pp. 101-125, 2002.
   
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
   
   *****************************************************************/


#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include <string.h>
#include "../c/cfortran.h"
#include "tclwrapper.h"

#include <tcl.h>

int Tclwrapper_Init(Tcl_Interp *interp);
static int wrapper(ClientData clientData, Tcl_Interp *interp, int objc,
Tcl_Obj *CONST objv[]);

static Tcl_Interp *tcl_interp;  

char *fnames[10];
char chind[50],command[200];

/***********************************************************************
 ***********************************************************************/

int Tclwrapper_Init(Tcl_Interp *interp) {

  if (Tcl_Init(interp) == TCL_ERROR) return TCL_ERROR;

  tcl_interp = interp;
  Tcl_CreateObjCommand(interp, "algencan", wrapper,
                       (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);

  return TCL_OK;
}

/***********************************************************************
 ***********************************************************************/

static int wrapper(ClientData clientData, Tcl_Interp *interp, int objc,
Tcl_Obj *CONST objv[]) {
	
  Tcl_Obj *tcl_result,*lambdai,*vobj,*xi;

  double cnorm,epsfeas,epsopt,f,snorm,nlpsupn;
  int checkder,inform,iprint,m,n,ncomp;
  double *l,*lambda,*u,*x;
  int *equatn,*linear,*coded;

  char *lambdaname,*svector,*xname; 
  int err,i,tempint,countargs = 1;
  

  tcl_result = Tcl_GetObjResult(interp);

  /* Verify the number of arguments */

  if ( objc != 30 ) {
    Tcl_SetStringObj(tcl_result,
                     "Wrong number of input arguments to the solver. Usage:\n"
                     " algencan  evalfname evalgname evalhname "
                     "evalcname evaljacname evalhcname evalfcname evalgjacname"
                     " evalhlname evalhlpname epsfeas epsopt iprint ncomp n x "
                     "l u m lambda equatn linear coded checkder f "
                     "cnorm snorm nlpsupn inform\n", -1);
    return TCL_ERROR;
  }


  /* 
     Get information from TCL
  */

  /* User-defined function names

     0 - evalf
     1 - evalg
     2 - evalh
     3 - evalc
     4 - evaljac
     5 - evalhc
     6 - evalfc
     7 - evalgjac
     8 - evalhl
     9 - evalhlp
  */

  for (i = 0; i < 10; i++) {
    fnames[i] = Tcl_GetStringFromObj(objv[countargs], &tempint);
    if ( fnames[i] == NULL ) return TCL_ERROR;
    
    countargs++;
  }

  /* Feasibility epsilon */

  err = Tcl_GetDoubleFromObj(interp, objv[countargs], &epsfeas);
  if ( err == TCL_ERROR ) return TCL_ERROR;

  countargs++;

  /* Optimality epsilon */

  err = Tcl_GetDoubleFromObj(interp, objv[countargs], &epsopt);
  if ( err == TCL_ERROR ) return TCL_ERROR;

  countargs++;

  /* Output information detail */

  err = Tcl_GetIntFromObj(interp, objv[countargs], &iprint); 
  if ( err == TCL_ERROR ) return TCL_ERROR;

  countargs++;

  /* Number of components */

  err = Tcl_GetIntFromObj(interp, objv[countargs], &ncomp);
  if ( err == TCL_ERROR) return TCL_ERROR;

  countargs++;

  /* Number of variables */

  err = Tcl_GetIntFromObj(interp, objv[countargs], &n);
  if ( err == TCL_ERROR ) return TCL_ERROR;

  x = (double *) malloc(sizeof(double) * n);
  l = (double *) malloc(sizeof(double) * n);
  u = (double *) malloc(sizeof(double) * n);

  if ( x == NULL || l == NULL || u == NULL ) {
    Tcl_SetStringObj(tcl_result,"\nTCL INTERFACE ERROR: "
                     "malloc returned NULL.\n", -1);
    return TCL_ERROR;
  }

  countargs++;

  /* Initial point */

  svector = Tcl_GetString(objv[countargs]);

  xname = (char *) malloc(strlen(svector) * sizeof(char));

  if ( xname == NULL ) {
    Tcl_SetStringObj(tcl_result,"\nTCL INTERFACE ERROR: "
                     "malloc returned NULL.\n", -1);
    return TCL_ERROR;
  }

  strcpy(xname,svector);

  for (i = 0; i < n; i++) {
    sprintf(chind,"%d",i);

    vobj = Tcl_GetVar2Ex(tcl_interp,svector,chind,TCL_LEAVE_ERR_MSG);
    if ( vobj == NULL ) return TCL_ERROR;
      
    err = Tcl_GetDoubleFromObj(tcl_interp,vobj,&(x[i]));
    if ( err == TCL_ERROR) return TCL_ERROR;
  }

  countargs++;

  /* Lower bounds */

  svector = Tcl_GetString(objv[countargs]);
  
  for (i = 0; i < n; i++) {
    sprintf(chind,"%d",i);

    vobj = Tcl_GetVar2Ex(tcl_interp,svector,chind,TCL_LEAVE_ERR_MSG);
    if ( vobj == NULL ) return TCL_ERROR;

    err = Tcl_GetDoubleFromObj(tcl_interp,vobj,&(l[i]));
    if ( err == TCL_ERROR ) return TCL_ERROR;
  }

  countargs++;

  /* Upper bounds */

  svector = Tcl_GetString(objv[countargs]);
  
  for (i = 0; i < n; i++) {
    sprintf(chind,"%d",i);

    vobj = Tcl_GetVar2Ex(tcl_interp,svector,chind,TCL_LEAVE_ERR_MSG);
    if ( vobj == NULL ) return TCL_ERROR;
    
    err = Tcl_GetDoubleFromObj(tcl_interp,vobj,&(u[i]));
    if ( err == TCL_ERROR ) return TCL_ERROR;
  }

  countargs++;

  /* Number of constraints */

  err = Tcl_GetIntFromObj(interp, objv[countargs], &m);
  if ( err == TCL_ERROR ) return TCL_ERROR;

  lambda = (double *) malloc(sizeof(double) * m);
  equatn = (int    *) malloc(sizeof(int)    * m);
  linear = (int    *) malloc(sizeof(int)    * m);
  coded  = (int    *) malloc(sizeof(int)   * 10);            

  if ( lambda == NULL || equatn == NULL || linear == NULL ||
       coded  == NULL ) {
    Tcl_SetStringObj(tcl_result,"\nTCL INTERFACE ERROR: "
                     "malloc returned NULL.\n", -1);
    return TCL_ERROR;
  }

  countargs++;

  /* Lambda */

  svector = Tcl_GetString(objv[countargs]);

  lambdaname = (char *) malloc(strlen(svector) * sizeof(char));

  if ( lambdaname == NULL ) {
    Tcl_SetStringObj(tcl_result,"\nTCL INTERFACE ERROR: "
                     "malloc returned NULL.\n", -1);
    return TCL_ERROR;
  }

  strcpy(lambdaname,svector);
  
  for (i = 0; i < m; i++) {
    sprintf(chind,"%d",i);
    vobj = Tcl_GetVar2Ex(tcl_interp,svector,chind,TCL_LEAVE_ERR_MSG);
    
    if ( vobj == NULL ) return TCL_ERROR;
    
    err = Tcl_GetDoubleFromObj(tcl_interp,vobj,&(lambda[i]));
    if ( err == TCL_ERROR ) return TCL_ERROR;
  }

  countargs++;

  /* Equatn */

  svector = Tcl_GetString(objv[countargs]);
  
  for (i = 0; i < m; i++) {
    sprintf(chind,"%d",i);

    vobj = Tcl_GetVar2Ex(tcl_interp,svector,chind,TCL_LEAVE_ERR_MSG);
    if ( vobj == NULL ) return TCL_ERROR;
    
    err = Tcl_GetIntFromObj(tcl_interp,vobj,&(equatn[i]));
    if ( err == TCL_ERROR ) return TCL_ERROR;
  }

  countargs++;

  /* Linear */

  svector = Tcl_GetString(objv[countargs]);
  
  for (i = 0; i < m; i++) {
    sprintf(chind,"%d",i);

    vobj = Tcl_GetVar2Ex(tcl_interp,svector,chind,TCL_LEAVE_ERR_MSG);
    if ( vobj == NULL ) return TCL_ERROR;
    
    err = Tcl_GetIntFromObj(tcl_interp,vobj,&(linear[i]));
    if ( err == TCL_ERROR ) return TCL_ERROR;
  }

  countargs++;

  /* Coded */

  svector = Tcl_GetString(objv[countargs]);
  
  for (i = 0; i < 10; i++) {
    sprintf(chind,"%d",i);

    vobj = Tcl_GetVar2Ex(tcl_interp,svector,chind,TCL_LEAVE_ERR_MSG);
    if ( vobj == NULL ) return TCL_ERROR;
    
    err = Tcl_GetIntFromObj(tcl_interp,vobj,&(coded[i]));
    if ( err == TCL_ERROR ) return TCL_ERROR;
  }

  countargs++;

  /* Check derivatives? */

  err = Tcl_GetIntFromObj(interp, objv[countargs], &checkder);
  if ( err == TCL_ERROR ) return TCL_ERROR;

  
  /* TRANSFORMS THE TYPE INT OF C IN THE TYPE LOGICAL OF FORTRAN
   * FOR THE ARRAYS equatn AND linear, AND VARIABLE checkder */
  C2FLOGICALV(equatn,m);
  C2FLOGICALV(linear,m);
  C2FLOGICALV(coded,10);
  
  /*
    CALL OPTIMIZATION SOLVER
  */
  
  Algencan(epsfeas,epsopt,iprint,ncomp,n,x,l,u,m,lambda,equatn,linear,
  coded,checkder,f,cnorm,snorm,nlpsupn,inform);
  

  /* 
     RETRIEVE INFORMATION
  */

  vobj = Tcl_NewDoubleObj(f);
  if ( Tcl_SetVar2Ex(tcl_interp,      "f",NULL,vobj,
                     TCL_LEAVE_ERR_MSG) == NULL ) { return TCL_ERROR; }

  vobj = Tcl_NewDoubleObj(cnorm);
  if ( Tcl_SetVar2Ex(tcl_interp,  "cnorm",NULL,vobj,
                     TCL_LEAVE_ERR_MSG) == NULL ) { return TCL_ERROR; }

  vobj = Tcl_NewDoubleObj(snorm);
  if ( Tcl_SetVar2Ex(tcl_interp,  "snorm",NULL,vobj,
                     TCL_LEAVE_ERR_MSG) == NULL ) { return TCL_ERROR; }

  vobj = Tcl_NewDoubleObj(nlpsupn);
  if ( Tcl_SetVar2Ex(tcl_interp,"nlpsupn",NULL,vobj,
                     TCL_LEAVE_ERR_MSG) == NULL ) { return TCL_ERROR; }

  vobj = Tcl_NewIntObj(inform);
  if ( Tcl_SetVar2Ex(tcl_interp, "inform",NULL,vobj,
                     TCL_LEAVE_ERR_MSG) == NULL ) { return TCL_ERROR; }

  for (i = 0; i < n; i++) {
    sprintf(chind,"%d",i);

    xi = Tcl_NewDoubleObj(x[i]);
    if ( xi == NULL ) return TCL_ERROR;
    
    if ( Tcl_SetVar2Ex(tcl_interp,xname,chind,xi,
                       TCL_LEAVE_ERR_MSG) == NULL ) { return TCL_ERROR; }
  }

  for (i = 0; i < m; i++) {
    sprintf(chind,"%d",i);

    lambdai = Tcl_NewDoubleObj(lambda[i]);
    if ( lambdai == NULL ) return TCL_ERROR;
    
    if ( Tcl_SetVar2Ex(tcl_interp,lambdaname,chind,lambdai,
                       TCL_LEAVE_ERR_MSG) == NULL ) { return TCL_ERROR; }
  }

  /*
    FREE MEMORY.
  */

  free(coded     );
  free(equatn    );
  free(linear    );
  free(l         );
  free(lambda    );
  free(lambdaname);
  free(u         );
  free(x         );
  free(xname     );

  return TCL_OK;
}

/***********************************************************************
 ***********************************************************************/

void evalf(int n,double *x,double *f,int *flag) {

  Tcl_Obj *tf,*tflag,*xi;
  int i,err;

  sprintf(command,"%s %d x f flag",fnames[0],n);

  for ( i = 0; i < n; i++) {
    sprintf(chind,"%d",i);
    xi = Tcl_NewDoubleObj(x[i]);

    if ( xi == NULL ||
         Tcl_SetVar2Ex(tcl_interp,"x",chind,xi,
                       TCL_LEAVE_ERR_MSG) == NULL ) {
      *flag = - 1;
      return;
    }
  }

  err = Tcl_EvalEx(tcl_interp,command,-1,TCL_EVAL_DIRECT);

  if ( err == TCL_ERROR) {
    fprintf(stderr,"\nTCL INTERFACE ERROR: %s\n",
            Tcl_GetStringResult(tcl_interp));
    *flag = - 1;
    return;
  }

  tf = Tcl_GetVar2Ex(tcl_interp,"f",NULL,TCL_LEAVE_ERR_MSG);

  if ( tf == NULL ||
       Tcl_GetDoubleFromObj(tcl_interp,tf,f) == TCL_ERROR ) {
    *flag = - 1;
    return;
  }

  tflag = Tcl_GetVar2Ex(tcl_interp,"flag",NULL,TCL_LEAVE_ERR_MSG);

  if ( tflag == NULL ||
       Tcl_GetIntFromObj(tcl_interp,tflag,flag) == TCL_ERROR ) {
    *flag = - 1;
    return;
  }

}


/***********************************************************************
 ***********************************************************************/

void evalg(int n,double *x,double *g,int *flag) {

  Tcl_Obj *tflag,*tg,*xi;
  int i,err;

  sprintf(command,"%s %d x g flag",fnames[1],n);

  for ( i = 0; i < n; i++) {
    sprintf(chind,"%d",i);
    xi = Tcl_NewDoubleObj(x[i]);

    if ( xi == NULL ||
         Tcl_SetVar2Ex(tcl_interp,"x",chind,xi,
                       TCL_LEAVE_ERR_MSG) == NULL ) {
      *flag = - 1;
      return;
    }
  }

  err = Tcl_EvalEx(tcl_interp,command,-1,TCL_EVAL_DIRECT);

  if ( err == TCL_ERROR) {
    fprintf(stderr,"\nTCL INTERFACE ERROR: %s\n",
            Tcl_GetStringResult(tcl_interp));
    *flag = - 1;
    return;
  }

  for (i = 0; i < n; i++) {
    sprintf(chind,"%d",i);
    tg = Tcl_GetVar2Ex(tcl_interp,"g",chind,TCL_LEAVE_ERR_MSG);

    if ( tg == NULL ||
         Tcl_GetDoubleFromObj(tcl_interp,tg,&(g[i])) == TCL_ERROR ) {
      *flag = - 1;
      return;
    }
  }

  tflag = Tcl_GetVar2Ex(tcl_interp,"flag",NULL,TCL_LEAVE_ERR_MSG);

  if ( tflag == NULL ||
       Tcl_GetIntFromObj(tcl_interp,tflag,flag) == TCL_ERROR ) {
    *flag = - 1;
    return;
  }

}

/***********************************************************************
 ***********************************************************************/

void evalh(int n,double *x,int *hlin,int *hcol,double *hval,int *hnnz,
int *flag) {

  Tcl_Obj *tflag,*thlin,*thcol,*thnnz,*thval,*xi;
  int i,err;

  sprintf(command,"%s %d x hlin hcol hval hnnz flag",fnames[2],n);

  for ( i = 0; i < n; i++) {
    sprintf(chind,"%d",i);
    xi = Tcl_NewDoubleObj(x[i]);

    if ( xi == NULL ||
         Tcl_SetVar2Ex(tcl_interp,"x",chind,xi,
                       TCL_LEAVE_ERR_MSG) == NULL ) {
      *flag = - 1;
      return;
    }
  }

  err = Tcl_EvalEx(tcl_interp,command,-1,TCL_EVAL_DIRECT);

  if ( err == TCL_ERROR) {
    fprintf(stderr,"\nTCL INTERFACE ERROR: %s\n",
            Tcl_GetStringResult(tcl_interp));
    *flag = - 1;
    return;
  }

  thnnz = Tcl_GetVar2Ex(tcl_interp,"hnnz",NULL,TCL_LEAVE_ERR_MSG);

  if ( thnnz == NULL ||
       Tcl_GetIntFromObj(tcl_interp,thnnz,hnnz) == TCL_ERROR ) {
    *flag = - 1;
    return;
  }

  for (i = 0; i < *hnnz; i++) {
    sprintf(chind,"%d",i);
    thlin = Tcl_GetVar2Ex(tcl_interp,"hlin",chind,TCL_LEAVE_ERR_MSG);
    thcol = Tcl_GetVar2Ex(tcl_interp,"hcol",chind,TCL_LEAVE_ERR_MSG);
    thval = Tcl_GetVar2Ex(tcl_interp,"hval",chind,TCL_LEAVE_ERR_MSG);

    if ( thlin == NULL || thcol == NULL || thval == NULL ) {
      *flag = - 1;
      return;
    }

    if ( Tcl_GetIntFromObj(   tcl_interp,thlin,&(hlin[i])) == TCL_ERROR ||
         Tcl_GetIntFromObj(   tcl_interp,thcol,&(hcol[i])) == TCL_ERROR ||
         Tcl_GetDoubleFromObj(tcl_interp,thval,&(hval[i])) == TCL_ERROR ) {
      *flag = - 1;
      return;
    }
  }

  tflag = Tcl_GetVar2Ex(tcl_interp,"flag",NULL,TCL_LEAVE_ERR_MSG);

  if ( tflag == NULL ||
       Tcl_GetIntFromObj(tcl_interp,tflag,flag) == TCL_ERROR) {
    *flag = - 1;
    return;
  }

}

/***********************************************************************
 ***********************************************************************/

void evalc(int n,double *x,int ind,double *cind,int *flag) {

  Tcl_Obj *tcind,*tind,*tflag,*xi;
  int i,err;

  sprintf(command,"%s %d x %d cind flag",fnames[3],n,ind);

  for ( i = 0; i < n; i++) {
    sprintf(chind,"%d",i);
    xi = Tcl_NewDoubleObj(x[i]);

    if ( xi == NULL ||
         Tcl_SetVar2Ex(tcl_interp,"x",chind,xi,
                       TCL_LEAVE_ERR_MSG) == NULL ) {
      *flag = - 1;
      return;
    }
  }

  err = Tcl_EvalEx(tcl_interp,command,-1,TCL_EVAL_DIRECT);

  if ( err == TCL_ERROR) {
    fprintf(stderr,"\nTCL INTERFACE ERROR: %s\n",
            Tcl_GetStringResult(tcl_interp));
    *flag = - 1;
    return;
  }

  tcind = Tcl_GetVar2Ex(tcl_interp,"cind",NULL,TCL_LEAVE_ERR_MSG);

  if ( tcind == NULL ||
       Tcl_GetDoubleFromObj(tcl_interp,tcind,cind) == TCL_ERROR ) {
    *flag = - 1;
    return;
  }

  tflag = Tcl_GetVar2Ex(tcl_interp,"flag",NULL,TCL_LEAVE_ERR_MSG);

  if ( tflag == NULL ||
       Tcl_GetIntFromObj(tcl_interp,tflag,flag) == TCL_ERROR ) {
    *flag = - 1;
    return;
  }

}

/***********************************************************************
 ***********************************************************************/

void evaljac(int n,double *x,int ind,int *jcvar,double *jcval,int *jcnnz, 
int *flag) {

  Tcl_Obj *tjcnnz,*tjcval,*tjcvar,*tflag,*xi;
  int i,err;

  sprintf(command,"%s %d x %d jcvar jcval jcnnz flag",fnames[4],n,ind);

  for ( i = 0; i < n; i++) {
    sprintf(chind,"%d",i);
    xi = Tcl_NewDoubleObj(x[i]);

    if ( xi == NULL ||
         Tcl_SetVar2Ex(tcl_interp,"x",chind,xi,
                       TCL_LEAVE_ERR_MSG) == NULL ) {
      *flag = - 1;
      return;
    }
  }

  err = Tcl_EvalEx(tcl_interp,command,-1,TCL_EVAL_DIRECT);

  if ( err == TCL_ERROR) {
    fprintf(stderr,"\nTCL INTERFACE ERROR: %s\n",
            Tcl_GetStringResult(tcl_interp));
    *flag = - 1;
    return;
  }

  tjcnnz = Tcl_GetVar2Ex(tcl_interp,"jcnnz",NULL,TCL_LEAVE_ERR_MSG);

  if ( tjcnnz == NULL ||
       Tcl_GetIntFromObj(tcl_interp,tjcnnz,jcnnz) == TCL_ERROR ) {
    *flag = - 1;
    return;
  }

  for (i = 0; i < *jcnnz; i++) {
    sprintf(chind,"%d",i);
    tjcvar = Tcl_GetVar2Ex(tcl_interp,"jcvar",chind,TCL_LEAVE_ERR_MSG);
    tjcval = Tcl_GetVar2Ex(tcl_interp,"jcval",chind,TCL_LEAVE_ERR_MSG);

    if ( tjcvar == NULL || tjcval == NULL ) {
      *flag = - 1;
      return;
    }
    
    if ( Tcl_GetIntFromObj(   tcl_interp,tjcvar,&(jcvar[i])) == TCL_ERROR ||
         Tcl_GetDoubleFromObj(tcl_interp,tjcval,&(jcval[i])) == TCL_ERROR ) {
      *flag = - 1;
      return;
    }
  }

  tflag = Tcl_GetVar2Ex(tcl_interp,"flag",NULL,TCL_LEAVE_ERR_MSG);

  if ( tflag == NULL ||
       Tcl_GetIntFromObj(tcl_interp,tflag,flag) == TCL_ERROR ) {
    *flag = - 1;
    return;
  }

}

/***********************************************************************
 ***********************************************************************/

void evalhc(int n,double *x,int ind,int *hclin,int *hccol,double *hcval,
int *hcnnz,int *flag) {

  Tcl_Obj *thccol,*thclin,*thcnnz,*thcval,*tflag,*xi;
  int i,err;

  sprintf(command,"%s %d x %d hclin hccol hcval hcnnz flag",
          fnames[5],n,ind);

  for ( i = 0; i < n; i++) {
    sprintf(chind,"%d",i);
    xi = Tcl_NewDoubleObj(x[i]);

    if ( xi == NULL ||
         Tcl_SetVar2Ex(tcl_interp,"x",chind,xi,
                       TCL_LEAVE_ERR_MSG) == NULL ) {
      *flag = - 1;
      return;
    }
  }

  err = Tcl_EvalEx(tcl_interp,command,-1,TCL_EVAL_DIRECT);

  if ( err == TCL_ERROR) {
    fprintf(stderr,"\nTCL INTERFACE ERROR: %s\n",
            Tcl_GetStringResult(tcl_interp));
    *flag = - 1;
    return;
  }

  thcnnz = Tcl_GetVar2Ex(tcl_interp,"hcnnz",NULL,TCL_LEAVE_ERR_MSG);

  if ( thcnnz == NULL ||
       Tcl_GetIntFromObj(tcl_interp,thcnnz,hcnnz) == TCL_ERROR ) {
    *flag = - 1;
    return;
  }

  for (i = 0; i < *hcnnz; i++) {
    sprintf(chind,"%d",i);
    thclin = Tcl_GetVar2Ex(tcl_interp,"hclin",chind,TCL_LEAVE_ERR_MSG);
    thccol = Tcl_GetVar2Ex(tcl_interp,"hccol",chind,TCL_LEAVE_ERR_MSG);
    thcval = Tcl_GetVar2Ex(tcl_interp,"hcval",chind,TCL_LEAVE_ERR_MSG);

    if ( thclin == NULL || thccol == NULL || thcval == NULL ) {
      *flag = - 1;
      return;
    }

    if ( Tcl_GetIntFromObj(   tcl_interp,thclin,&(hclin[i])) == TCL_ERROR ||
         Tcl_GetIntFromObj(   tcl_interp,thccol,&(hccol[i])) == TCL_ERROR ||
         Tcl_GetDoubleFromObj(tcl_interp,thcval,&(hcval[i])) == TCL_ERROR ) {
      *flag = - 1;
      return;
    }
  }

  tflag = Tcl_GetVar2Ex(tcl_interp,"flag",NULL,TCL_LEAVE_ERR_MSG);

  if ( tflag == NULL ||
       Tcl_GetIntFromObj(tcl_interp,tflag,flag) == TCL_ERROR ) {
    *flag = - 1;
    return;
  }

}

/***********************************************************************
 ***********************************************************************/

void evalfc(int n,double *x,double *f,int m,double *c,int *flag) {

  Tcl_Obj *tc,*tf,*tflag,*xi;
  int i,err;

  sprintf(command,"%s %d x f %d c flag",fnames[6],n,m);

  for ( i = 0; i < n; i++) {
    sprintf(chind,"%d",i);
    xi = Tcl_NewDoubleObj(x[i]);

    if ( xi == NULL ||
         Tcl_SetVar2Ex(tcl_interp,"x",chind,xi,
                       TCL_LEAVE_ERR_MSG) == NULL ) {
      *flag = - 1;
      return;
    }
  }

  err = Tcl_EvalEx(tcl_interp,command,-1,TCL_EVAL_DIRECT);

  if ( err == TCL_ERROR) {
    fprintf(stderr,"\nTCL INTERFACE ERROR: %s\n",
            Tcl_GetStringResult(tcl_interp));
    *flag = - 1;
    return;
  }

  tf = Tcl_GetVar2Ex(tcl_interp,"f",NULL,TCL_LEAVE_ERR_MSG);

  if ( tf == NULL ||
       Tcl_GetDoubleFromObj(tcl_interp,tf,f) == TCL_ERROR ) {
    *flag = - 1;
    return;
  }
      
  for (i = 0; i < m; i++) {
    sprintf(chind,"%d",i);
    tc = Tcl_GetVar2Ex(tcl_interp,"c",chind,TCL_LEAVE_ERR_MSG);

    if ( tc == NULL ||
         Tcl_GetDoubleFromObj(tcl_interp,tc,&(c[i])) == TCL_ERROR ) {
      *flag = - 1;
      return;
    }
  }

  tflag = Tcl_GetVar2Ex(tcl_interp,"flag",NULL,TCL_LEAVE_ERR_MSG);

  if ( tflag == NULL ||
       Tcl_GetIntFromObj(tcl_interp,tflag,flag) == TCL_ERROR ) {
    *flag = - 1;
    return;
  }

}

/***********************************************************************
 ***********************************************************************/

void evalgjac(int n,double *x,double *g,int m,int *jcfun,int *jcvar,
double *jcval,int *jcnnz,int *flag) {

  Tcl_Obj *tflag,*tg,*tjcfun,*tjcnnz,*tjcvar,*tjcval,*xi;
  int i;

  sprintf(command,"%s %d x g %d jcfun jcvar jcval jcnnz flag",
          fnames[7],n,m);

  for ( i = 0; i < n; i++) {
    sprintf(chind,"%d",i);
    xi = Tcl_NewDoubleObj(x[i]);

    if ( xi == NULL ||
         Tcl_SetVar2Ex(tcl_interp,"x",chind,xi,
                       TCL_LEAVE_ERR_MSG) == NULL ) {
      *flag = - 1;
      return;
    }
  }

  if ( Tcl_EvalEx(tcl_interp,command,-1,TCL_EVAL_DIRECT) == TCL_ERROR) {
    fprintf(stderr,"\nTCL INTERFACE ERROR: %s\n",
            Tcl_GetStringResult(tcl_interp));
    *flag = - 1;
    return;
  }

  for (i = 0; i < n; i++) {
    sprintf(chind,"%d",i);
    tg = Tcl_GetVar2Ex(tcl_interp,"g",chind,TCL_LEAVE_ERR_MSG);

    if ( tg == NULL ||
         Tcl_GetDoubleFromObj(tcl_interp,tg,&(g[i])) == TCL_ERROR ) {
      *flag = - 1;
      return;
    }
  }

  tjcnnz = Tcl_GetVar2Ex(tcl_interp,"jcnnz",NULL,TCL_LEAVE_ERR_MSG);

  if ( tjcnnz == NULL ||
       Tcl_GetIntFromObj(tcl_interp,tjcnnz,jcnnz) == TCL_ERROR ) {
    *flag = - 1;
    return;
  }

  for (i = 0; i < *jcnnz; i++) {
    sprintf(chind,"%d",i);
    tjcfun = Tcl_GetVar2Ex(tcl_interp,"jcfun",chind,TCL_LEAVE_ERR_MSG);
    tjcvar = Tcl_GetVar2Ex(tcl_interp,"jcvar",chind,TCL_LEAVE_ERR_MSG);
    tjcval = Tcl_GetVar2Ex(tcl_interp,"jcval",chind,TCL_LEAVE_ERR_MSG);

    if ( tjcfun == NULL || tjcvar == NULL || tjcval == NULL ) {
      *flag = - 1;
      return;
    }

    if ( Tcl_GetIntFromObj(   tcl_interp,tjcfun,&(jcfun[i])) == TCL_ERROR ||
         Tcl_GetIntFromObj(   tcl_interp,tjcvar,&(jcvar[i])) == TCL_ERROR ||
         Tcl_GetDoubleFromObj(tcl_interp,tjcval,&(jcval[i])) == TCL_ERROR ) {
      *flag = - 1;
      return;
    }
  }

  tflag = Tcl_GetVar2Ex(tcl_interp,"flag",NULL,TCL_LEAVE_ERR_MSG);

  if ( tflag == NULL ||
       Tcl_GetIntFromObj(tcl_interp,tflag,flag) == TCL_ERROR ) {
    *flag = - 1;
    return;
  }

}


/***********************************************************************
 ***********************************************************************/

void evalhl(int n,double *x,int m,double *lambda,double sf,double *sc,
int *hllin,int *hlcol,double *hlval,int *hlnnz,int *flag) {

  Tcl_Obj *tflag,*thlcol,*thllin,*thlnnz,*tlambdai,*tsci,*tsf,*thlval,*xi;
  int i;

  sprintf(command,"%s %d x %d lambda $sf sc hllin hlcol hlval hlnnz flag",
          fnames[8],n,m);

  for ( i = 0; i < n; i++) {
    sprintf(chind,"%d",i);
    xi = Tcl_NewDoubleObj(x[i]);

    if ( xi == NULL ||
         Tcl_SetVar2Ex(tcl_interp,"x",chind,xi,
                       TCL_LEAVE_ERR_MSG) == NULL ) {
      *flag = - 1;
      return;
    }
  }

  tsf = Tcl_NewDoubleObj(sf);

  if ( tsf == NULL ||
       Tcl_SetVar2Ex(tcl_interp,"sf",NULL,tsf,
                     TCL_LEAVE_ERR_MSG) == NULL ) {
    *flag = - 1;
    return;
  }
  
  for (i = 0; i < m; i++) {
    sprintf(chind,"%d",i);

    tsci = Tcl_NewDoubleObj(sc[i]);

    if ( tsci == NULL ||
         Tcl_SetVar2Ex(tcl_interp,"sc",chind,tsci,
                       TCL_LEAVE_ERR_MSG) == NULL ) {
      *flag = - 1;
      return;
    }

    tlambdai = Tcl_NewDoubleObj(lambda[i]);

    if ( tlambdai == NULL ||
         Tcl_SetVar2Ex(tcl_interp,"lambda",chind,tlambdai,
                       TCL_LEAVE_ERR_MSG) == NULL ) {
      *flag = - 1;
      return;
    }
  }

  if ( Tcl_EvalEx(tcl_interp,command,-1,TCL_EVAL_DIRECT) == TCL_ERROR) {
    fprintf(stderr,"\nTCL INTERFACE ERROR: %s\n",
            Tcl_GetStringResult(tcl_interp));
    *flag = - 1;
    return;
  }

  thlnnz = Tcl_GetVar2Ex(tcl_interp,"hlnnz",NULL,TCL_LEAVE_ERR_MSG);

  if ( thlnnz == NULL ||
       Tcl_GetIntFromObj(tcl_interp,thlnnz,hlnnz) == TCL_ERROR ) {
    *flag = - 1;
    return;
  }

  for (i = 0; i < *hlnnz; i++) {
    sprintf(chind,"%d",i);
    thllin = Tcl_GetVar2Ex(tcl_interp,"hllin",chind,TCL_LEAVE_ERR_MSG);
    thlcol = Tcl_GetVar2Ex(tcl_interp,"hlcol",chind,TCL_LEAVE_ERR_MSG);
    thlval = Tcl_GetVar2Ex(tcl_interp,"hlval",chind,TCL_LEAVE_ERR_MSG);

    if ( thllin == NULL || thlcol == NULL || thlval == NULL ) {
      *flag = - 1;
      return;
    }

    if ( Tcl_GetIntFromObj(   tcl_interp,thllin,&(hllin[i])) == TCL_ERROR ||
         Tcl_GetIntFromObj(   tcl_interp,thlcol,&(hlcol[i])) == TCL_ERROR ||
         Tcl_GetDoubleFromObj(tcl_interp,thlval,&(hlval[i])) == TCL_ERROR ) {
      *flag = - 1;
      return;
    }
  }

  tflag = Tcl_GetVar2Ex(tcl_interp,"flag",NULL,TCL_LEAVE_ERR_MSG);

  if ( tflag == NULL ||
       Tcl_GetIntFromObj(tcl_interp,tflag,flag) == TCL_ERROR ) {
    *flag = - 1;
    return;
  }

}

/***********************************************************************
 ***********************************************************************/

void evalhlp(int n,double *x,int m,double *lambda,double sf,double *sc,
double *p,double *hp,int gothl,int *flag) {

  Tcl_Obj *tflag,*tlambdai,*tpi,*thp,*tsci,*tsf,*xi;
  int i;

  sprintf(command,"%s %d x %d lambda $sf sc p hp %d flag",
          fnames[9],n,m,gothl);

  for ( i = 0; i < n; i++) {
    sprintf(chind,"%d",i);

    xi = Tcl_NewDoubleObj(x[i]);

    if ( xi == NULL ||
         Tcl_SetVar2Ex(tcl_interp,"x",chind,xi,
                       TCL_LEAVE_ERR_MSG) == NULL ) {
      *flag = - 1;
      return;
    }

    tpi = Tcl_NewDoubleObj(p[i]);

    if ( tpi == NULL ||
         Tcl_SetVar2Ex(tcl_interp,"p",chind,tpi,
                       TCL_LEAVE_ERR_MSG) == NULL ) {
      *flag = - 1;
      return;
    }
  }

  tsf = Tcl_NewDoubleObj(sf);
  Tcl_SetVar2Ex(tcl_interp,"sf",NULL,tsf,TCL_LEAVE_ERR_MSG);

  for (i = 0; i < m; i++) {
    sprintf(chind,"%d",i);

    tsci = Tcl_NewDoubleObj(sc[i]);

    if ( tsci == NULL ||
         Tcl_SetVar2Ex(tcl_interp,"sc",chind,tsci,
                       TCL_LEAVE_ERR_MSG) == NULL ) {
      *flag = - 1;
      return;
    }

    tlambdai = Tcl_NewDoubleObj(lambda[i]);

    if ( tlambdai == NULL ||
         Tcl_SetVar2Ex(tcl_interp,"lambda",chind,tlambdai,
                       TCL_LEAVE_ERR_MSG) == NULL ) {
      *flag = - 1;
      return;
    }
  }

  if ( Tcl_EvalEx(tcl_interp,command,-1,TCL_EVAL_DIRECT) == TCL_ERROR) {
    fprintf(stderr,"\nTCL INTERFACE ERROR: %s\n",
            Tcl_GetStringResult(tcl_interp));
    *flag = - 1;
    return;
  }

  for (i = 0; i < n; i++) {
    sprintf(chind,"%d",i);
    thp = Tcl_GetVar2Ex(tcl_interp,"hp",chind,TCL_LEAVE_ERR_MSG);

    if ( thp == NULL ||
         Tcl_GetDoubleFromObj(tcl_interp,thp,&(hp[i])) == TCL_ERROR ) {
      *flag = - 1;
      return;
    }
  }

  tflag = Tcl_GetVar2Ex(tcl_interp,"flag",NULL,TCL_LEAVE_ERR_MSG);

  if ( tflag == NULL ||
       Tcl_GetIntFromObj(tcl_interp,tflag,flag) == TCL_ERROR ) {
    *flag = - 1;
    return;
  }

}
