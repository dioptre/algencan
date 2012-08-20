static PyObject *evalf_py    = NULL;
static PyObject *evalg_py    = NULL;
static PyObject *evalh_py    = NULL;
static PyObject *evalc_py    = NULL;
static PyObject *evaljac_py  = NULL;
static PyObject *evalhc_py   = NULL;
static PyObject *evalfc_py   = NULL;
static PyObject *evalgjac_py = NULL;
static PyObject *evalhl_py   = NULL;
static PyObject *evalhlp_py  = NULL;
static PyObject *inip_py     = NULL;
static PyObject *endp_py     = NULL;
static PyObject *param_py    = NULL;

static PyObject *return_value;

static PyObject *pywrapper_solver(PyObject *self,PyObject *args);

void inip(int *n,double **x,double **l,double **u,int *m,double **lambda,
          int **equatn,int **linear,int *coded,int *checkder);

void evalf(int n,double *x,double *f,int *flag);

void evalg(int n,double *x,double *g,int *flag);

void evalh(int n,double *x,int *hlin,int *hcol,double *hval,int *hnnz,
           int *flag);

void evalc(int n,double *x,int ind,double *c,int *flag);

void evaljac(int n,double *x,int ind,int *jcvar,double *jcval,int *jcnnz,
             int *flag);

void evalhc(int n,double *x,int ind,int *hclin,int *hccol,double *hcval,
            int *hcnnz,int *flag);

void evalfc(int n,double *x,double *f,int m,double *c,int *flag);

void evalgjac(int n,double *x,double *g,int m,int *jcfun,int *jcvar,
              double *jcval,int *jcnnz,int *flag);

void evalhl(int n,double *x,int m,double *lambda,double scalef,
            double *scalec,int *hllin,int *hlcol,double *hlval,int *hlnnz,
            int *flag);

void evalhlp(int n,double *x,int m,double *lambda,double scalef,
             double *scalec,double *p,double *hp,int *goth,int *flag);

void endp(int n,double *x,double *l,double *u,int m,double *lambda,
int *equatn,int *linear);

void param(double *epsfeas,double *epsopt,int *iprint,int *ncomp);

int BuildIntArray(int size,PyObject *input,int *array);

int BuildRealArray(int size,PyObject *input,double *array);

int BuildIntPyArray(int size,int *input,PyObject **array);

int BuildRealPyArray(int size,double *input,PyObject **array);

FCALLSCSUB4(evalf,EVALF,evalf,INT,DOUBLEV,PDOUBLE,PINT)

FCALLSCSUB4(evalg,EVALG,evalg,INT,DOUBLEV,DOUBLEV,PINT)

FCALLSCSUB7(evalh,EVALH,evalh,INT,DOUBLEV,INTV,INTV,DOUBLEV,PINT,PINT)

FCALLSCSUB5(evalc,EVALC,evalc,INT,DOUBLEV,INT,PDOUBLE,PINT)

FCALLSCSUB7(evaljac,EVALJAC,evaljac,INT,DOUBLEV,INT,INTV,DOUBLEV,PINT,PINT)

FCALLSCSUB8(evalhc,EVALHC,evalhc,INT,DOUBLEV,INT,INTV,INTV,DOUBLEV,PINT,PINT)

FCALLSCSUB6(evalfc,EVALFC,evalfc,INT,DOUBLEV,PDOUBLE,INT,DOUBLEV,PINT)

FCALLSCSUB9(evalgjac,EVALGJAC,evalgjac,INT,DOUBLEV,DOUBLEV,INT,INTV,INTV,\
DOUBLEV,PINT,PINT)

FCALLSCSUB11(evalhl,EVALHL,evalhl,INT,DOUBLEV,INT,DOUBLEV,DOUBLE,DOUBLEV,\
INTV,INTV,DOUBLEV,PINT,PINT)

FCALLSCSUB10(evalhlp,EVALHLP,evalhlp,INT,DOUBLEV,INT,DOUBLEV,DOUBLE,DOUBLEV,\
DOUBLEV,DOUBLEV,PLOGICAL,PINT)

PROTOCCALLSFSUB19(ALGENCAN,algencan,DOUBLE,DOUBLE,INT,INT,INT,DOUBLEV,\
DOUBLEV,DOUBLEV,INT,DOUBLEV,LOGICALV,LOGICALV,LOGICALV,LOGICAL,PDOUBLE,\
PDOUBLE,PDOUBLE,PDOUBLE,PINT)

#define Algencan(epsfeas,epsopt,iprint,ncomp,n,x,l,u,m,lambda,\
equatn,linear,coded,checkder,f,cnorm,snorm,nlpsupn,inform) \
CCALLSFSUB19(ALGENCAN,algencan,DOUBLE,DOUBLE,INT,INT,INT,DOUBLEV,DOUBLEV,\
DOUBLEV,INT,DOUBLEV,LOGICALV,LOGICALV,LOGICALV,LOGICAL,PDOUBLE,\
PDOUBLE,PDOUBLE,PDOUBLE,PINT,\
epsfeas,epsopt,iprint,ncomp,n,x,l,u,m,lambda,equatn,linear,coded,checkder,\
f,cnorm,snorm,nlpsupn,inform)
