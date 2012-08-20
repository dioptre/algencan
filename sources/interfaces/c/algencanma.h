void algencan(double epsfeas,double epsopt,int iprint,int ncomp,int n,
              double *x,double *l,double *u,int m,double *lambda,int *equatn,
              int *linear,int *coded,int checkder,double *f,double *cnorm,
              double *snorm,double *nlpsupn,int *inform);

void param(double *epsfeas,double *epsopt,int *iprint,int *ncomp);

void inip(int *n,double **x,double **l,double **u,int *m,double **lambda,
	  int **equatn,int **linear,int *coded,int *checkder);

void endp(int n,double *x,double *l,double *u,int m,double *lambda,
	  int *equatn,int *linear);

