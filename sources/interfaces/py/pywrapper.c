/* ==================================================================
   Module: Interface between Python and ALGENCAN
   ==================================================================

   Last update of any of the components of this module:

   March 25, 2008.

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

#include <Python.h>
#include <numpy/arrayobject.h>
#include "../c/cfortran.h"
#include "pywrapper.h"

static PyMethodDef pywrapper_methods[] = {
  {"solver",pywrapper_solver,METH_VARARGS,"Call optimization solver."},
  {NULL,NULL,0,NULL}
};

/* ******************************************************************
   ****************************************************************** */

PyMODINIT_FUNC initpywrapper(void) {

  (void) Py_InitModule("pywrapper",pywrapper_methods);
  import_array();

}

/* ******************************************************************
   ****************************************************************** */

static PyObject *pywrapper_solver(PyObject *self,PyObject *args) {

  int checkder,inform,iprint,m,n,ncomp;
  double cnorm,f,nlpsupn,epsfeas,epsopt,snorm;

  int coded[10];
  int *equatn,*linear;
  double *l,*lambda,*u,*x;

  return_value = Py_None;
  setbuf(stderr,(char *) malloc(BUFSIZ));

  if (!PyArg_ParseTuple(args,"O!O!O!O!O!O!O!O!O!O!O!O!O!:solver",
                        &PyFunction_Type,&evalf_py,
                        &PyFunction_Type,&evalg_py,
                        &PyFunction_Type,&evalh_py,
                        &PyFunction_Type,&evalc_py,
                        &PyFunction_Type,&evaljac_py,
                        &PyFunction_Type,&evalhc_py,
                        &PyFunction_Type,&evalfc_py,
                        &PyFunction_Type,&evalgjac_py,
                        &PyFunction_Type,&evalhl_py,
                        &PyFunction_Type,&evalhlp_py,
                        &PyFunction_Type,&inip_py,
                        &PyFunction_Type,&endp_py,
                        &PyDict_Type,&param_py)) goto cleanup;

  Py_INCREF(evalf_py   );
  Py_INCREF(evalg_py   );
  Py_INCREF(evalh_py   );
  Py_INCREF(evalc_py   );
  Py_INCREF(evaljac_py );
  Py_INCREF(evalhc_py  );
  Py_INCREF(evalfc_py  );
  Py_INCREF(evalgjac_py);
  Py_INCREF(evalhl_py  );
  Py_INCREF(evalhlp_py );
  Py_INCREF(inip_py    );
  Py_INCREF(endp_py    );
  Py_INCREF(param_py   );

  param(&epsfeas,&epsopt,&iprint,&ncomp);
  if (return_value == NULL) goto cleanup;

  inip(&n,&x,&l,&u,&m,&lambda,&equatn,&linear,coded,&checkder);
  if (return_value == NULL) goto cleanup;

  C2FLOGICALV(equatn,m);
  C2FLOGICALV(linear,m);
  C2FLOGICALV(coded,10);

  Algencan(epsfeas,epsopt,iprint,ncomp,n,x,l,u,m,lambda,equatn,
  linear,coded,checkder,f,cnorm,snorm,nlpsupn,inform);
  if (return_value == NULL) goto cleanup;

  endp(n,x,l,u,m,lambda,equatn,linear);
  if (return_value == NULL) goto cleanup;

 cleanup:
  Py_XDECREF(evalf_py   );
  Py_XDECREF(evalg_py   );
  Py_XDECREF(evalh_py   );
  Py_XDECREF(evalc_py   );
  Py_XDECREF(evaljac_py );
  Py_XDECREF(evalhc_py  );
  Py_XDECREF(evalfc_py  );
  Py_XDECREF(evalgjac_py);
  Py_XDECREF(evalhl_py  );
  Py_XDECREF(evalhlp_py );
  Py_XDECREF(inip_py    );
  Py_XDECREF(endp_py    );
  Py_XDECREF(param_py   );

  fflush(stdout);
  fflush(stderr);

  Py_XINCREF(return_value);
  return return_value;

}

/* ******************************************************************
   ****************************************************************** */

void inip(int *n,double **x,double **l,double **u,int *m,double **lambda,
          int **equatn,int **linear,int *coded,int *checkder) {

  PyObject *result = NULL,*x_py,*l_py,*u_py,*lambda_py,*equatn_py,*linear_py,
           *coded_py;

  if ((result = PyEval_CallFunction(inip_py,"()")) == NULL) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, PyEval_CallFunction "
                   "returned NULL\n\n",__func__);
    return_value = NULL;
    goto cleanup;
  }

  if (!PyArg_ParseTuple(result,"iOOOiOOOOi:inip",n,&x_py,&l_py,&u_py,m,
                        &lambda_py,&equatn_py,&linear_py,&coded_py,checkder)) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, PyArg_ParseTuple "
                   "returned false\n\n",__func__);
    return_value = NULL;
    goto cleanup;
  }

  if ((*x      = (double *) malloc(*n * sizeof(double))) == NULL ||
      (*l      = (double *) malloc(*n * sizeof(double))) == NULL ||
      (*u      = (double *) malloc(*n * sizeof(double))) == NULL ||
      (*lambda = (double *) malloc(*m * sizeof(double))) == NULL ||
      (*equatn = (int    *) malloc(*m * sizeof(int   ))) == NULL ||
      (*linear = (int    *) malloc(*m * sizeof(int   ))) == NULL) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, malloc "
                   "returned NULL\n\n",__func__);
    return_value = NULL;
    goto cleanup;
  }

  if (BuildRealArray(*n,x_py,*x          ) == -1 ||
      BuildRealArray(*n,l_py,*l          ) == -1 ||
      BuildRealArray(*n,u_py,*u          ) == -1 ||
      BuildRealArray(*m,lambda_py,*lambda) == -1 ||
      BuildIntArray( *m,equatn_py,*equatn) == -1 ||
      BuildIntArray( *m,linear_py,*linear) == -1 ||
      BuildIntArray( 10, coded_py,  coded) == -1) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, Build(Int|Real)Array "
                   "returned -1\n\n",__func__);
    return_value = NULL;
    goto cleanup;
  }

 cleanup:
  Py_XDECREF(result);

}

/* ******************************************************************
   ****************************************************************** */

void evalf(int n,double *x,double *f,int *flag) {

  PyObject *result = NULL,*x_py = NULL;

  *flag = -1;

  if (BuildRealPyArray(n,x,&x_py) == -1) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, BuildRealPyArray "
                   "returned -1\n\n",__func__);
    return_value = NULL;
    goto cleanup;
  }

  if ((result = PyEval_CallFunction(evalf_py,"(O)",x_py)) == NULL) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, PyEval_CallFunction "
                   "returned NULL\n\n",__func__);
    return_value = NULL;
    goto cleanup;
  }

  if (!PyArg_ParseTuple(result,"di:evalf",f,flag)) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, PyArg_ParseTuple "
                   "returned false\n\n",__func__);
    return_value = NULL; *flag = -1;
    goto cleanup;
  }

 cleanup:
  Py_XDECREF(x_py  );
  Py_XDECREF(result);

}

/* ******************************************************************
   ****************************************************************** */

void evalg(int n,double *x,double *g,int *flag) {

  PyObject *result = NULL,*x_py = NULL,*g_py;

  *flag = -1;

  if (BuildRealPyArray(n,x,&x_py) == -1) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, BuildRealPyArray "
                   "returned -1\n\n",__func__);
    return_value = NULL;
    goto cleanup;
  }

  if ((result = PyEval_CallFunction(evalg_py,"(O)",x_py)) == NULL) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, PyEval_CallFunction "
                   "returned NULL\n\n",__func__);
    return_value = NULL;
    goto cleanup;
  }

  if (!PyArg_ParseTuple(result,"Oi:evalg",&g_py,flag)) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, PyArg_ParseTuple "
                   "returned false\n\n",__func__);
    return_value = NULL; *flag = -1;
    goto cleanup;
  }

  if (BuildRealArray(n,g_py,g) == -1) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, BuildRealArray "
                   "returned -1\n\n",__func__);
    return_value = NULL; *flag = -1;
    goto cleanup;
  }

 cleanup:
  Py_XDECREF(x_py  );
  Py_XDECREF(result);

}

/* ******************************************************************
   ****************************************************************** */

void evalh(int n,double *x,int *hlin,int *hcol,double *hval,int *hnnz,
           int *flag) {

  int i;
  PyObject *result = NULL,*x_py = NULL,*hlin_py,*hcol_py,*hval_py;

  *flag = -1;

  if (BuildRealPyArray(n,x,&x_py) == -1) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, BuildRealPyArray "
                   "returned -1\n\n",__func__);
    return_value = NULL;
    goto cleanup;
  }

  if ((result = PyEval_CallFunction(evalh_py,"(O)",x_py)) == NULL) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, PyEval_CallFunction "
                   "returned NULL\n\n",__func__);
    return_value = NULL;
    goto cleanup;
  }

  if (!PyArg_ParseTuple(result,"OOOii:evalh",&hlin_py,&hcol_py,&hval_py,hnnz,
                        flag)) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, PyArg_ParseTuple "
                   "returned false\n\n",__func__);
    return_value = NULL; *flag = -1;
    goto cleanup;
  }

  if (BuildIntArray( *hnnz,hlin_py,hlin) == -1 ||
      BuildIntArray( *hnnz,hcol_py,hcol) == -1 ||
      BuildRealArray(*hnnz,hval_py,hval) == -1) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, Build(Int|Real)Array "
                   "returned -1\n\n",__func__);
    return_value = NULL; *flag = -1;
    goto cleanup;
  }

  for (i = 0; i < *hnnz; i++) {
     hlin[i]++;
     hcol[i]++;
  }

 cleanup:
  Py_XDECREF(x_py  );
  Py_XDECREF(result);

}

/* ******************************************************************
   ****************************************************************** */

void evalc(int n,double *x,int ind,double *c,int *flag) {

  PyObject *result = NULL,*x_py = NULL;

  *flag = -1;

  if (BuildRealPyArray(n,x,&x_py) == -1) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, BuildRealPyArray "
                   "returned -1\n\n",__func__);
    return_value = NULL;
    goto cleanup;
  }

  if ((result = PyEval_CallFunction(evalc_py,"(Oi)",x_py,ind)) == NULL) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, PyEval_CallFunction "
                   "returned NULL\n\n",__func__);
    return_value = NULL;
    goto cleanup;
  }

  if (!PyArg_ParseTuple(result,"di:evalc",c,flag)) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, PyArg_ParseTuple "
                   "returned false\n\n",__func__);
    return_value = NULL; *flag = -1;
    goto cleanup;
  }

 cleanup:
  Py_XDECREF(x_py  );
  Py_XDECREF(result);

}

/* ******************************************************************
   ****************************************************************** */

void evaljac(int n,double *x,int ind,int *jcvar,double *jcval,int *jcnnz,
             int *flag) {

  int i;
  PyObject *result = NULL,*x_py = NULL,*jcvar_py,*jcval_py;

  *flag = -1;

  if (BuildRealPyArray(n,x,&x_py) == -1) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, BuildRealPyArray "
                   "returned -1\n\n",__func__);
    return_value = NULL;
    goto cleanup;
  }

  if ((result = PyEval_CallFunction(evaljac_py,"(Oi)",x_py,ind)) == NULL) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, PyEval_CallFunction "
                   "returned NULL\n\n",__func__);
    return_value = NULL;
    goto cleanup;
  }

  if (!PyArg_ParseTuple(result,"OOii:evaljac",&jcvar_py,&jcval_py,jcnnz,
                        flag)) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, PyArg_ParseTuple "
                   "returned false\n\n",__func__);
    return_value = NULL; *flag = -1;
    goto cleanup;
  }

  if (BuildIntArray( *jcnnz,jcvar_py,jcvar) == -1 ||
      BuildRealArray(*jcnnz,jcval_py,jcval) == -1) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, Build(Int|Real)Array "
                   "returned -1\n\n",__func__);
    return_value = NULL; *flag = -1;
    goto cleanup;
  }

  for (i = 0; i < *jcnnz; i++)
     jcvar[i]++;

 cleanup:
  Py_XDECREF(x_py  );
  Py_XDECREF(result);

}

/* ******************************************************************
   ****************************************************************** */

void evalhc(int n,double *x,int ind,int *hclin,int *hccol,double *hcval,
            int *hcnnz,int *flag) {

  int i;
  PyObject *result = NULL,*x_py = NULL,*hclin_py,*hccol_py,*hcval_py;

  *flag = -1;

  if (BuildRealPyArray(n,x,&x_py) == -1) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, BuildRealPyArray "
                   "returned -1\n\n",__func__);
    return_value = NULL;
    goto cleanup;
  }

  if ((result = PyEval_CallFunction(evalhc_py,"(Oi)",x_py,ind)) == NULL) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, PyEval_CallFunction "
                   "returned NULL\n\n",__func__);
    return_value = NULL;
    goto cleanup;
  }

  if (!PyArg_ParseTuple(result,"OOOii:evalhc",&hclin_py,&hccol_py,&hcval_py,
                        hcnnz,flag)) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, PyArg_ParseTuple "
                   "returned false\n\n",__func__);
    return_value = NULL; *flag = -1;
    goto cleanup;
  }

  if (BuildIntArray( *hcnnz,hclin_py,hclin) == -1 ||
      BuildIntArray( *hcnnz,hccol_py,hccol) == -1 ||
      BuildRealArray(*hcnnz,hcval_py,hcval) == -1) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, Build(Int|Real)Array "
                   "returned -1\n\n",__func__);
    return_value = NULL; *flag = -1;
    goto cleanup;
  }

  for (i = 0; i < *hcnnz; i++) {
     hclin[i]++;
     hccol[i]++;
  }

 cleanup:
  Py_XDECREF(x_py  );
  Py_XDECREF(result);

}

/* *****************************************************************
   ***************************************************************** */

void evalfc(int n,double *x,double *f,int m,double *c,int *flag) {

  PyObject *result = NULL,*x_py = NULL,*c_py;

  *flag = -1;

  if (BuildRealPyArray(n,x,&x_py) == -1) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, BuildRealPyArray "
                   "returned -1\n\n",__func__);
    return_value = NULL;
    goto cleanup;
  }

  if ((result = PyEval_CallFunction(evalfc_py,"(Oi)",x_py,m)) == NULL) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, PyEval_CallFunction "
                   "returned NULL\n\n",__func__);
    return_value = NULL;
    goto cleanup;
  }

  if (!PyArg_ParseTuple(result,"dOi:evalfc",f,&c_py,flag)) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, PyArg_ParseTuple "
                   "returned false\n\n",__func__);
    return_value = NULL; *flag = -1;
    goto cleanup;
  }

  if (BuildRealArray(m,c_py,c) == -1) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, BuildRealArray "
                   "returned -1\n\n",__func__);
    return_value = NULL; *flag = -1;
    goto cleanup;
  }

 cleanup:
  Py_XDECREF(x_py     );
  Py_XDECREF(result   );

}

/* *****************************************************************
   ***************************************************************** */

void evalgjac(int n,double *x,double *g,int m,int *jcfun,int *jcvar,
              double *jcval,int *jcnnz,int *flag) {

  int i;
  PyObject *result = NULL,*x_py = NULL,*g_py,*jcfun_py,*jcvar_py,*jcval_py;

  *flag = -1;

  if (BuildRealPyArray(n,x,&x_py) == -1) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, BuildRealPyArray "
                   "returned -1\n\n",__func__);
    return_value = NULL;
    goto cleanup;
  }

  if ((result = PyEval_CallFunction(evalgjac_py,"(Oi)",x_py,m)) == NULL) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, PyEval_CallFunction "
                   "returned NULL\n\n",__func__);
    return_value = NULL;
    goto cleanup;
  }

  if (!PyArg_ParseTuple(result,"OOOOii:evalgjac",&g_py,&jcfun_py,&jcvar_py,
                        &jcval_py,jcnnz,flag)) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, PyArg_ParseTuple "
                   "returned false\n\n",__func__);
    return_value = NULL; *flag = -1;
    goto cleanup;
  }

  if (BuildRealArray(n,g_py,g             ) == -1 ||
      BuildIntArray( *jcnnz,jcfun_py,jcfun) == -1 ||
      BuildIntArray( *jcnnz,jcvar_py,jcvar) == -1 ||
      BuildRealArray(*jcnnz,jcval_py,jcval) == -1) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, Build(Int|Real)Array "
                   "returned -1\n\n",__func__);
    return_value = NULL; *flag = -1;
    goto cleanup;
  }

  for (i = 0; i < *jcnnz; i++)
     jcvar[i]++;

 cleanup:
  Py_XDECREF(x_py  );
  Py_XDECREF(result);

}

/* *****************************************************************
   ***************************************************************** */

void evalhl(int n,double *x,int m,double *lambda,double scalef,
            double *scalec,int *hllin,int *hlcol,double *hlval,int *hlnnz,
            int *flag) {

  int i;
  PyObject *result = NULL,*x_py = NULL,*lambda_py = NULL,*scalec_py = NULL,
           *hllin_py,*hlcol_py,*hlval_py;

  *flag = -1;

  if (BuildRealPyArray(n,x,&x_py          ) == -1 ||
      BuildRealPyArray(m,lambda,&lambda_py) == -1 ||
      BuildRealPyArray(m,scalec,&scalec_py) == -1) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, BuildRealPyArray "
                   "returned -1\n\n",__func__);
    return_value = NULL;
    goto cleanup;
  }

  if ((result = PyEval_CallFunction(evalhl_py,"(OiOdO)",x_py,m,lambda_py,scalef,
                                    scalec_py)) == NULL) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, PyEval_CallFunction "
                   "returned NULL\n\n",__func__);
    return_value = NULL;
    goto cleanup;
  }

  if (!PyArg_ParseTuple(result,"OOOii:evalhl",&hllin_py,&hlcol_py,&hlval_py,
                        hlnnz,flag)) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, PyArg_ParseTuple "
                   "returned false\n\n",__func__);
    return_value = NULL; *flag = -1;
    goto cleanup;
  }

  if (BuildIntArray( *hlnnz,hllin_py,hllin) == -1 ||
      BuildIntArray( *hlnnz,hlcol_py,hlcol) == -1 ||
      BuildRealArray(*hlnnz,hlval_py,hlval) == -1) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, Build(Int|Real)Array "
                   "returned -1\n\n",__func__);
    return_value = NULL; *flag = -1;
    goto cleanup;
  }

  for (i = 0; i < *hlnnz; i++) {
     hllin[i]++;
     hlcol[i]++;
  }

 cleanup:
  Py_XDECREF(x_py     );
  Py_XDECREF(lambda_py);
  Py_XDECREF(scalec_py);
  Py_XDECREF(result   );

}

/* ******************************************************************
   ****************************************************************** */

void evalhlp(int n,double *x,int m,double *lambda,double scalef,
             double *scalec,double *p,double *hp,int *goth,int *flag) {

  PyObject *result = NULL,*x_py = NULL,*lambda_py = NULL,*scalec_py = NULL,
           *p_py = NULL,*hp_py;

  *flag = -1;

  if (BuildRealPyArray(n,x,&x_py          ) == -1 ||
      BuildRealPyArray(m,lambda,&lambda_py) == -1 ||
      BuildRealPyArray(m,scalec,&scalec_py) == -1 ||
      BuildRealPyArray(n,p,&p_py          ) == -1) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, BuildRealPyArray "
                   "returned -1\n\n",__func__);
    return_value = NULL;
    goto cleanup;
  }

  if ((result = PyEval_CallFunction(evalhlp_py,"(OiOdOi)",x_py,m,lambda_py,
                                    scalef,scalec_py,p_py,*goth)) == NULL) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, PyEval_CallFunction "
                   "returned NULL\n\n",__func__);
    return_value = NULL;
    goto cleanup;
  }

  if (!PyArg_ParseTuple(result,"Oii:evalhlp",&hp_py,goth,flag)) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, PyArg_ParseTuple "
                   "returned false\n\n",__func__);
    return_value = NULL; *flag = -1;
    goto cleanup;
  }

  if (BuildRealArray(n,hp_py,hp) == -1) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, BuildRealArray "
                   "returned -1\n\n",__func__);
    return_value = NULL; *flag = -1;
    goto cleanup;
  }

 cleanup:
  Py_XDECREF(x_py     );
  Py_XDECREF(lambda_py);
  Py_XDECREF(scalec_py);
  Py_XDECREF(p_py     );
  Py_XDECREF(result   );

}

/* ******************************************************************
   ****************************************************************** */

void endp(int n,double *x,double *l,double *u,int m,double *lambda,
          int *equatn,int *linear) {

  PyObject *result = NULL,*x_py = NULL,*l_py = NULL,*u_py = NULL,
           *lambda_py = NULL,*equatn_py = NULL,*linear_py = NULL;

  if (BuildRealPyArray(n,x,&x_py          ) == -1 ||
      BuildRealPyArray(n,l,&l_py          ) == -1 ||
      BuildRealPyArray(n,u,&u_py          ) == -1 ||
      BuildRealPyArray(m,lambda,&lambda_py) == -1 ||
      BuildIntPyArray( m,equatn,&equatn_py) == -1 ||
      BuildIntPyArray( m,linear,&linear_py) == -1) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, Build(Int|Real)PyArray "
                   "returned -1\n\n",__func__);
    return_value = NULL;
    goto cleanup;
  }

  if ((result = PyEval_CallFunction(endp_py,"(OOOiOOO)",x_py,l_py,u_py,m,
                                    lambda_py,equatn_py,linear_py)) == NULL) {
    fprintf(stderr,"\nPYTHON INTERFACE ERROR: in %s, PyEval_CallFunction "
                   "returned NULL\n\n",__func__);
    return_value = NULL;
    goto cleanup;
  }

 cleanup:
  free(x     );
  free(l     );
  free(u     );
  free(lambda);
  free(equatn);
  free(linear);

  Py_XDECREF(x_py     );
  Py_XDECREF(l_py     );
  Py_XDECREF(u_py     );
  Py_XDECREF(lambda_py);
  Py_XDECREF(equatn_py);
  Py_XDECREF(linear_py);
  Py_XDECREF(result   );

}

/* ******************************************************************
   ****************************************************************** */

void param(double *epsfeas,double *epsopt,int *iprint,int *ncomp) {

  if (!PyDict_Check(param_py)) {
    PyErr_SetString(PyExc_TypeError,"param must be a dictionary");
    return_value = NULL;
    return;
  }

  *epsfeas = PyFloat_AsDouble(PyDict_GetItemString(param_py,"epsfeas"));
  if (PyErr_Occurred() != NULL) {
    PyErr_SetString(PyExc_ValueError,"value of key \'epsfeas\' is invalid");
    return_value = NULL;
    return;
  }

  *epsopt = PyFloat_AsDouble(PyDict_GetItemString(param_py,"epsopt"));
  if (PyErr_Occurred() != NULL) {
    PyErr_SetString(PyExc_ValueError,"value of key \'epsopt\' is invalid");
    return_value = NULL;
    return;
  }

  *iprint = PyInt_AsLong(PyDict_GetItemString(param_py,"iprint"));
  if (PyErr_Occurred() != NULL) {
    PyErr_SetString(PyExc_ValueError,"value of key \'iprint\' is invalid");
    return_value = NULL;
    return;
  }

  *ncomp = PyInt_AsLong(PyDict_GetItemString(param_py,"ncomp"));
  if (PyErr_Occurred() != NULL) {
    PyErr_SetString(PyExc_ValueError,"value of key \'ncomp\' is invalid");
    return_value = NULL;
    return;
  }

}

/* ******************************************************************
   ****************************************************************** */

int BuildIntArray(int size,PyObject *input,int *array) {

  int i,result = -1;
  PyArrayObject *array_py = NULL;

  Py_INCREF(input);

  if ((array_py = (PyArrayObject *)
       PyArray_ContiguousFromAny(input,PyArray_LONG,1,1)) == NULL)
    goto cleanup;

  if (PyArray_DIM(array_py,0) < size) {
    PyErr_SetString(PyExc_ValueError,"array shorter than expected");
    goto cleanup;
  }

  for (i = 0; i < size; i++)
    array[i] = ((long *) PyArray_DATA(array_py))[i];

  result = 0;

 cleanup:
  Py_XDECREF(input   );
  Py_XDECREF(array_py);

  return result;

}

/* ******************************************************************
   ****************************************************************** */

int BuildRealArray(int size,PyObject *input,double *array) {

  int i,result = -1;
  PyArrayObject *array_py = NULL;

  Py_INCREF(input);

  if ((array_py = (PyArrayObject *)
       PyArray_ContiguousFromAny(input,PyArray_DOUBLE,1,1)) == NULL)
    goto cleanup;

  if (PyArray_DIM(array_py,0) < size) {
    PyErr_SetString(PyExc_ValueError,"array shorter than expected");
    goto cleanup;
  }

  for (i = 0; i < size; i++)
    array[i] = ((double *) PyArray_DATA(array_py))[i];

  result = 0;

 cleanup:
  Py_XDECREF(input   );
  Py_XDECREF(array_py);

  return result;

}

/* ******************************************************************
   ****************************************************************** */

int BuildIntPyArray(int size,int *input,PyObject **array) {

  npy_intp dimensions[1];

  dimensions[0] = size;

  /* ENSURE THAT input IS NOT FREED BEFORE THE ARRAY OBJECT IS DESTROYED */
  *array = PyArray_SimpleNewFromData(1,dimensions,PyArray_INT,(char *) input);
  if (*array == NULL)
    return -1;

  return 0;

}

/* ******************************************************************
   ****************************************************************** */

int BuildRealPyArray(int size,double *input,PyObject **array) {

  npy_intp dimensions[1];

  dimensions[0] = size;

  /* ENSURE THAT input IS NOT FREED BEFORE THE ARRAY OBJECT IS DESTROYED */
  *array = PyArray_SimpleNewFromData(1,dimensions,PyArray_DOUBLE,
                                     (char *) input);
  if (*array == NULL)
    return -1;

  return 0;

}
