#include <jni.h>
#include "../c/cfortran.h"
#include "javawrapper.h"

JNIEnv * ENV;

jdoubleArray jlambda,jsc,jp,jx;
jmethodID    jevalf,jevalg,jevalh,jevalc,jevaljac,jevalhc,jevalfc,
             jevalgjac,jevalhl,jevalhlp;
jobject algProb;


/***********************************************************************
 ***********************************************************************/

/*
 * Class:     tango_algencan_base_ALGENCANProblem
 * Method:    Algencan
 * Signature: (DDIII[D[D[DI[D[Z[Z[ZZ)Ltango/algencan/base/Status;
 */

JNIEXPORT jobject JNICALL Java_tango_algencan_base_ALGENCANProblem_Algencan
(JNIEnv * env,jobject jobj,jdouble epsfeas,jdouble epsopt,jint iprint,
 jint ncomp,jint n,jdoubleArray x,jdoubleArray l,jdoubleArray u,jint m,
 jdoubleArray lambda,jbooleanArray equatn,jbooleanArray linear,
 jbooleanArray coded,jboolean checkder) {

  double cnorm,f,nlpsupn,snorm;
  int i,inform, *jie, *jil, *jic;
  jdouble *jdx, *jdl, *jdu, *jdlambda;
  jboolean *jbe, *jbl, *jbc;
  jclass clazz;
  jmethodID constructor;
  jobject status;

  ENV     = env;
  algProb = jobj;

  /*
   * Load user-supplied functions into global variables.
   */

  clazz  = (*ENV)->GetObjectClass(ENV, algProb);

  jevalf = (*ENV)->GetMethodID(ENV, clazz, "evalf",            "([D)D");
  jevalg = (*ENV)->GetMethodID(ENV, clazz, "evalg",           "([D)[D");
  jevalh = (*ENV)->GetMethodID(ENV, clazz, "evalh", 
                                   "([D)Ltango/algencan/base/Hessian;");

  jevalc   = (*ENV)->GetMethodID(ENV, clazz,   "evalc",       "([DI)D");
  jevaljac = (*ENV)->GetMethodID(ENV, clazz, "evaljac",
                                 "([DI)Ltango/algencan/base/Jacobian;");
  jevalhc  = (*ENV)->GetMethodID(ENV, clazz,  "evalhc",
                                  "([DI)Ltango/algencan/base/Hessian;");

  jevalfc   = (*ENV)->GetMethodID(ENV, clazz,   "evalfc",
                   "([D)Ltango/algencan/base/ObjectiveAndConstraints;");
  jevalgjac = (*ENV)->GetMethodID(ENV, clazz, "evalgjac",
                      "([D)Ltango/algencan/base/GradientAndJacobian;" );
  jevalhl   = (*ENV)->GetMethodID(ENV, clazz,   "evalhl",
                             "([D[DD[D)Ltango/algencan/base/Hessian;" );
  jevalhlp  = (*ENV)->GetMethodID(ENV, clazz,  "evalhlp",
                    "([D[DD[D[DZ)Ltango/algencan/base/HLTimesVector;" );

  /*
   * Alocate data structure (global and local).
   */

  jie = (int*) malloc(m * sizeof(int));
  jil = (int*) malloc(m * sizeof(int));
  jic = (int*) malloc(10 * sizeof(int));

  jlambda = (*ENV)->NewGlobalRef(ENV, (*ENV)->NewDoubleArray(ENV, m));
  jx      = (*ENV)->NewGlobalRef(ENV, (*ENV)->NewDoubleArray(ENV, n));
  jsc     = (*ENV)->NewGlobalRef(ENV, (*ENV)->NewDoubleArray(ENV, m));
  jp      = (*ENV)->NewGlobalRef(ENV, (*ENV)->NewDoubleArray(ENV, n));

  jdx = (*ENV)->GetDoubleArrayElements(ENV, x, (jboolean*) JNI_FALSE);
  jdl = (*ENV)->GetDoubleArrayElements(ENV, l, (jboolean*) JNI_FALSE);
  jdu = (*ENV)->GetDoubleArrayElements(ENV, u, (jboolean*) JNI_FALSE);
  jdlambda = (*ENV)->GetDoubleArrayElements(ENV, lambda, (jboolean*) JNI_FALSE);

  jbe = (*ENV)->GetBooleanArrayElements(ENV, equatn, (jboolean*) JNI_FALSE);
  jbl = (*ENV)->GetBooleanArrayElements(ENV, linear, (jboolean*) JNI_FALSE);
  jbc = (*ENV)->GetBooleanArrayElements(ENV,  coded, (jboolean*) JNI_FALSE);


  /* TRANSFORMS THE TYPE INT OF C IN THE TYPE LOGICAL OF FORTRAN 
   * FOR THE ARRAYS equatn AND linear, AND VARIABLE checkder */

  for ( i = 0; i < m; i++ ) {
    jie[i] = jbe[i];
    jil[i] = jbl[i];
  }

  for ( i = 0; i < 10; i++ ) {
    jic[i] = jbc[i];
  }

  C2FLOGICALV(jie,m);
  C2FLOGICALV(jil,m);
  C2FLOGICALV(jic,10);

  /*
   * Call ALGENCAN
   */

  algencan((double) epsfeas,(double) epsopt,(int) iprint,(int) ncomp,
           (int) n,(double*) jdx,(double*) jdl,(double *) jdu,(int) m,
           (double*) jdlambda,jie,jil,jic,(int) checkder,f,cnorm,snorm,
           nlpsupn,inform);

  /*
   * Retrieve information and create Status object.
   */
  
  (*ENV)->SetDoubleArrayRegion(ENV, jx, 0, n, jdx);
  (*ENV)->SetDoubleArrayRegion(ENV, jlambda, 0, m, jdlambda);
  
  clazz       = (*ENV)->FindClass(  ENV,  "Ltango/algencan/base/Status;");
  constructor = (*ENV)->GetMethodID(ENV, clazz, "<init>", "(I[DI[DDDDDI)V");
  status      = (*ENV)->NewObject(ENV, clazz, constructor, (jint) n,jx,
                                  (jint) m,jlambda,(jdouble) f,
                                  (jdouble) cnorm,(jdouble) snorm,
                                  (jdouble) nlpsupn,(jint) inform);
  
  /*
   * Release data structure.
   */

  (*ENV)->ReleaseDoubleArrayElements(ENV, x, jdx, 0);
  (*ENV)->ReleaseDoubleArrayElements(ENV, l, jdl, 0);
  (*ENV)->ReleaseDoubleArrayElements(ENV, u, jdu, 0);
  (*ENV)->ReleaseDoubleArrayElements(ENV, u, jdlambda, 0);

  (*ENV)->ReleaseBooleanArrayElements(ENV, equatn, jbe, 0);
  (*ENV)->ReleaseBooleanArrayElements(ENV, linear, jbl, 0);
  (*ENV)->ReleaseBooleanArrayElements(ENV,  coded, jbc, 0);

  free(jie);
  free(jil);
  free(jic);

  (*ENV)->DeleteGlobalRef(ENV, jlambda);
  (*ENV)->DeleteGlobalRef(ENV, jsc    );
  (*ENV)->DeleteGlobalRef(ENV,      jx);
  (*ENV)->DeleteGlobalRef(ENV,      jp);

  (*ENV)->DeleteLocalRef(ENV, algProb);

  return status;
}

/***********************************************************************
 ***********************************************************************/

void evalf(int n,double *x,double *f,int *flag) {

   *flag = 0;

   (*ENV)->SetDoubleArrayRegion(ENV, jx, 0, n, (jdouble*) x);
   *f = (*ENV)->CallDoubleMethod(ENV, algProb, jevalf, jx);

   if ( (*ENV)->ExceptionOccurred(ENV) != NULL ) {
      *flag = - 1;
      return;
   }

}


/***********************************************************************
 ***********************************************************************/

void evalg(int n,double *x,double *g,int *flag) {

   jmethodID method;
   jclass clazz;
   jdoubleArray jg;

   *flag = 0;

   (*ENV)->SetDoubleArrayRegion(ENV, jx, 0, n, (jdouble*) x);

   if ( (*ENV)->ExceptionOccurred(ENV) != NULL ) {
      *flag = - 1;
      return;
   }

   jg = (jdoubleArray) (*ENV)->CallObjectMethod(ENV, algProb, jevalg, jx);

   if ( (*ENV)->ExceptionOccurred(ENV) != NULL ) {
      *flag = - 1;
      return;
   }

   (*ENV)->GetDoubleArrayRegion(ENV, jg, 0, n, (jdouble*) g);
   (*ENV)->DeleteLocalRef(ENV, jg);

}

/***********************************************************************
 ***********************************************************************/

void evalh(int n,double *x,int *hlin,int *hcol,double *hval,int *hnnz,
int *flag) {

   int i;
   jclass clazz;
   jmethodID method;
   jobject hessian;
   jintArray jhcol,jhlin;
   jdoubleArray jhval;

   *flag = 0;

   (*ENV)->SetDoubleArrayRegion(ENV, jx, 0, n, (jdouble*) x);
   hessian = (*ENV)->CallObjectMethod(ENV, algProb, jevalh, jx);

   if ( (*ENV)->ExceptionOccurred(ENV) != NULL ) {
      *flag = - 1;
      return;
   }

   clazz = (*ENV)->GetObjectClass(ENV, hessian);

   method = (*ENV)->GetMethodID(ENV, clazz, "getHnnz", "()I");
   *hnnz  = (*ENV)->CallIntMethod(ENV, hessian, method);

   method = (*ENV)->GetMethodID(ENV, clazz, "getHlin", "()[I");
   jhlin  = (jintArray) (*ENV)->CallObjectMethod(ENV, hessian, method);

   method = (*ENV)->GetMethodID(ENV, clazz, "getHcol", "()[I");
   jhcol  = (jintArray) (*ENV)->CallObjectMethod(ENV, hessian, method);

   method = (*ENV)->GetMethodID(ENV, clazz, "getHval", "()[D");
   jhval  = (jdoubleArray) (*ENV)->CallObjectMethod(ENV, hessian, method);

   (*ENV)->GetDoubleArrayRegion(ENV, jhval, 0, *hnnz, (jdouble*) hval);
   (*ENV)->GetIntArrayRegion(ENV, jhlin, 0, *hnnz, (jint*) hlin);
   (*ENV)->GetIntArrayRegion(ENV, jhcol, 0, *hnnz, (jint*) hcol);

   for ( i = 0; i < *hnnz; i++ ) {
      hlin[i]++;
      hcol[i]++;
   }

   (*ENV)->DeleteLocalRef(ENV, hessian);
   (*ENV)->DeleteLocalRef(ENV,   jhval);
   (*ENV)->DeleteLocalRef(ENV,   jhlin);
   (*ENV)->DeleteLocalRef(ENV,   jhcol);

}

/***********************************************************************
 ***********************************************************************/

void evalc(int n,double *x,int ind,double *c,int *flag) {

   *flag = 0;

   (*ENV)->SetDoubleArrayRegion(ENV, jx, 0, n, (jdouble*) x);
   *c = (*ENV)->CallDoubleMethod(ENV, algProb, jevalc, jx, (jint) ind);

   if ( (*ENV)->ExceptionOccurred(ENV) != NULL ) {
      *flag = - 1;
      return;
   }

}

/***********************************************************************
 ***********************************************************************/

void evaljac(int n,double *x,int ind,int *jcvar,double *jcval,
int *jcnnz,int *flag) {

   int i;
   jobject jacobian;
   jclass clazz;
   jmethodID method;
   jintArray jjcvar;
   jdoubleArray jjcval;

   *flag = 0;

   (*ENV)->SetDoubleArrayRegion(ENV, jx, 0, n, (jdouble*) x);
   jacobian = (*ENV)->CallObjectMethod(ENV, algProb, jevaljac, jx, (jint) ind);

   if ( (*ENV)->ExceptionOccurred(ENV) != NULL ) {
      *flag = - 1;
      return;
   }

   clazz  = (*ENV)->GetObjectClass(ENV, jacobian);

   method = (*ENV)->GetMethodID(ENV, clazz, "getJCnnz", "()I");
   *jcnnz = (*ENV)->CallIntMethod(ENV, jacobian, method);

   method = (*ENV)->GetMethodID(ENV, clazz, "getJCvar", "()[I");
   jjcvar = (jintArray) (*ENV)->CallObjectMethod(ENV, jacobian, method);

   method = (*ENV)->GetMethodID(ENV, clazz, "getJCval", "()[D");
   jjcval = (jdoubleArray) (*ENV)->CallObjectMethod(ENV, jacobian, method);
              
   (*ENV)->GetIntArrayRegion(ENV, jjcvar, 0, *jcnnz, (jint*) jcvar);
   (*ENV)->GetDoubleArrayRegion(ENV, jjcval, 0, *jcnnz, (jdouble*) jcval);

    for (i = 0; i < *jcnnz; i++)
       jcvar[i]++;

    (*ENV)->DeleteLocalRef(ENV, jacobian);
    (*ENV)->DeleteLocalRef(ENV,   jjcvar);
    (*ENV)->DeleteLocalRef(ENV,   jjcval);

}

/***********************************************************************
 ***********************************************************************/

void evalhc(int n,double *x,int ind,int *hclin,int *hccol,double *hcval,
int *hcnnz,int *flag) {

   int i;
   jobject hessian;
   jclass clazz;
   jmethodID method;
   jintArray jhccol,jhclin;
   jdoubleArray jhcval;

   *flag = 0;

   (*ENV)->SetDoubleArrayRegion(ENV, jx, 0, n, (jdouble*) x);
   hessian = (*ENV)->CallObjectMethod(ENV, algProb, jevalhc, jx, (jint) ind);

   if ( (*ENV)->ExceptionOccurred(ENV) != NULL ) {
      *flag = - 1;
      return;
   }

   clazz = (*ENV)->GetObjectClass(ENV, hessian);

   method = (*ENV)->GetMethodID(ENV, clazz, "getHnnz", "()I");
   *hcnnz = (*ENV)->CallIntMethod(ENV, hessian, method);

   method = (*ENV)->GetMethodID(ENV, clazz, "getHlin", "()[I");
   jhclin = (jintArray) (*ENV)->CallObjectMethod(ENV, hessian, method);

   method = (*ENV)->GetMethodID(ENV, clazz, "getHcol", "()[I");
   jhccol = (jintArray) (*ENV)->CallObjectMethod(ENV, hessian, method);

   method = (*ENV)->GetMethodID(ENV, clazz, "getHval", "()[D");
   jhcval = (jdoubleArray) (*ENV)->CallObjectMethod(ENV, hessian, method);

   (*ENV)->GetDoubleArrayRegion(ENV, jhcval, 0, *hcnnz, (jdouble*) hcval);
   (*ENV)->GetIntArrayRegion(ENV, jhclin, 0, *hcnnz, (jint*) hclin);
   (*ENV)->GetIntArrayRegion(ENV, jhccol, 0, *hcnnz, (jint*) hccol);

   for ( i = 0; i < *hcnnz; i++ ) {
      hclin[i]++;
      hccol[i]++;
   }

   (*ENV)->DeleteLocalRef(ENV, hessian);
   (*ENV)->DeleteLocalRef(ENV,  jhcval);
   (*ENV)->DeleteLocalRef(ENV,  jhclin);
   (*ENV)->DeleteLocalRef(ENV,  jhccol);

}

/***********************************************************************
 ***********************************************************************/

void evalfc(int n,double *x,double *f,int m,double *c,int *flag) {

   jobject fc;
   jclass clazz;
   jmethodID method;
   jdoubleArray jc;

   *flag = 0;

   (*ENV)->SetDoubleArrayRegion(ENV, jx, 0, n, (jdouble*) x);
   fc = (*ENV)->CallObjectMethod(ENV, algProb, jevalfc, jx);

   if ( (*ENV)->ExceptionOccurred(ENV) != NULL ) {
      *flag = - 1;
      return;
   }

   clazz = (*ENV)->GetObjectClass(ENV, fc);
    
   method = (*ENV)->GetMethodID(ENV, clazz, "getF", "()D");
   *f     = (*ENV)->CallDoubleMethod(ENV, fc, method);

   method = (*ENV)->GetMethodID(ENV, clazz, "getC", "()[D");
   jc     = (jdoubleArray) (*ENV)->CallObjectMethod(ENV, fc, method);

   (*ENV)->GetDoubleArrayRegion(ENV, jc, 0, m, (jdouble*) c);
    
   (*ENV)->DeleteLocalRef(ENV, fc);
   (*ENV)->DeleteLocalRef(ENV, jc);

}

/***********************************************************************
 ***********************************************************************/

void evalgjac(int n,double *x,double *g,int m,int *jcfun,int *jcvar,
double *jcval,int *jcnnz,int *flag) {

   int i;
   jobject gj;
   jclass clazz;
   jmethodID method;
   jintArray jjcfun,jjcvar;
   jdoubleArray jg,jjcval;

   *flag = 0;

   (*ENV)->SetDoubleArrayRegion(ENV, jx, 0, n, (jdouble*) x);
   gj = (*ENV)->CallObjectMethod(ENV, algProb, jevalgjac, jx);

   if ( (*ENV)->ExceptionOccurred(ENV) != NULL ) {
      *flag = - 1;
      return;
   }

   clazz = (*ENV)->GetObjectClass(ENV, gj);
    
   method = (*ENV)->GetMethodID(ENV, clazz, "getG", "()[D");
   jg     = (jdoubleArray) (*ENV)->CallObjectMethod(ENV, gj, method);
    
   method = (*ENV)->GetMethodID(ENV, clazz, "getJCnnz", "()I");
   *jcnnz = (*ENV)->CallIntMethod(ENV, gj, method);

   method = (*ENV)->GetMethodID(ENV, clazz, "getJCfun", "()[I");
   jjcfun = (jintArray) (*ENV)->CallObjectMethod(ENV, gj, method);

   method = (*ENV)->GetMethodID(ENV, clazz, "getJCvar", "()[I");
   jjcvar = (jintArray) (*ENV)->CallObjectMethod(ENV, gj, method);

   method = (*ENV)->GetMethodID(ENV, clazz, "getJCval", "()[D");
   jjcval = (jdoubleArray) (*ENV)->CallObjectMethod(ENV, gj, method);
    
   (*ENV)->GetDoubleArrayRegion(ENV, jg, 0, n, (jdouble*) g);
   (*ENV)->GetIntArrayRegion(ENV, jjcfun, 0, *jcnnz, (jint*) jcfun);
   (*ENV)->GetIntArrayRegion(ENV, jjcvar, 0, *jcnnz, (jint*) jcvar);
   (*ENV)->GetDoubleArrayRegion(ENV, jjcval, 0, *jcnnz, (jdouble*) jcval);

   for (i = 0; i < *jcnnz; i++) {
      jcfun[i]++;
      jcvar[i]++;
   }

   (*ENV)->DeleteLocalRef(ENV,     gj);
   (*ENV)->DeleteLocalRef(ENV,     jg);
   (*ENV)->DeleteLocalRef(ENV, jjcfun);
   (*ENV)->DeleteLocalRef(ENV, jjcvar);
   (*ENV)->DeleteLocalRef(ENV, jjcval);

}

/***********************************************************************
 ***********************************************************************/

void evalhl(int n,double *x,int m,double *lambda,double sf,double *sc,
int *hllin,int *hlcol,double *hlval,int *hlnnz,int *flag) {

   int i;
   jobject hessianL;
   jclass clazz;
   jmethodID method;
   jintArray jhllin,jhlcol;
   jdoubleArray jhlval;

   *flag = 0;

   (*ENV)->SetDoubleArrayRegion(ENV, jx, 0, n, (jdouble*) x);
   (*ENV)->SetDoubleArrayRegion(ENV, jlambda, 0, m, (jdouble*) lambda);
   (*ENV)->SetDoubleArrayRegion(ENV, jsc, 0, m, (jdouble*) sc);

   hessianL = (*ENV)->CallObjectMethod(ENV, algProb, jevalhl, jx, jlambda, (jdouble) sf, jsc);

   if ( (*ENV)->ExceptionOccurred(ENV) != NULL ) {
      *flag = - 1;
      return;
   }

   clazz = (*ENV)->GetObjectClass(ENV, hessianL);

   method = (*ENV)->GetMethodID(ENV, clazz, "getHnnz", "()I");
   *hlnnz = (*ENV)->CallIntMethod(ENV, hessianL, method);

   method = (*ENV)->GetMethodID(ENV, clazz, "getHlin", "()[I");
   jhllin = (jintArray) (*ENV)->CallObjectMethod(ENV, hessianL, method);

   method = (*ENV)->GetMethodID(ENV, clazz, "getHcol", "()[I");
   jhlcol = (jintArray) (*ENV)->CallObjectMethod(ENV, hessianL, method);

   method = (*ENV)->GetMethodID(ENV, clazz, "getHval", "()[D");
   jhlval = (jdoubleArray) (*ENV)->CallObjectMethod(ENV, hessianL, method);

   (*ENV)->GetDoubleArrayRegion(ENV, jhlval, 0, *hlnnz, (jdouble*) hlval);
   (*ENV)->GetIntArrayRegion(ENV, jhllin, 0, *hlnnz, (jint*) hllin);
   (*ENV)->GetIntArrayRegion(ENV, jhlcol, 0, *hlnnz, (jint*) hlcol);

   for ( i = 0; i < *hlnnz; i++ ) {
      hllin[i]++;
      hlcol[i]++;      
   }

   (*ENV)->DeleteLocalRef(ENV, hessianL);
   (*ENV)->DeleteLocalRef(ENV,   jhlval);
   (*ENV)->DeleteLocalRef(ENV,   jhllin);
   (*ENV)->DeleteLocalRef(ENV,   jhlcol);

}

/***********************************************************************
 ***********************************************************************/

void evalhlp(int n,double *x,int m,double *lambda,double sf,double *sc,
double *p,double *hp,int *goth,int *flag) {

   jobject hlp;
   jclass clazz;
   jmethodID method;
   jdoubleArray jhp;

   *flag = 0;

   (*ENV)->SetDoubleArrayRegion(ENV, jx, 0, n, (jdouble*) x);
   (*ENV)->SetDoubleArrayRegion(ENV, jp, 0, n, (jdouble*) p); 
   (*ENV)->SetDoubleArrayRegion(ENV, jlambda, 0, m, (jdouble*) lambda);
   (*ENV)->SetDoubleArrayRegion(ENV, jsc, 0, m, (jdouble*) sc);

   hlp = (*ENV)->CallObjectMethod(ENV, algProb, jevalhlp, jx, jlambda, (jdouble) sf, jsc, jp,*goth);

   if ( (*ENV)->ExceptionOccurred(ENV) != NULL ) {
      *flag = - 1;
      return;
   }

   clazz = (*ENV)->GetObjectClass(ENV, hlp);

   method = (*ENV)->GetMethodID(ENV, clazz, "getHP", "()[D");
   jhp    = (jdoubleArray) (*ENV)->CallObjectMethod(ENV, hlp, method);
    
   method = (*ENV)->GetMethodID(ENV, clazz, "getGoth", "()Z");
   *goth  = (*ENV)->CallBooleanMethod(ENV, hlp, method);

   (*ENV)->GetDoubleArrayRegion(ENV, jhp, 0, n, (jdouble*) hp);

   (*ENV)->DeleteLocalRef(ENV, hlp);
   (*ENV)->DeleteLocalRef(ENV, jhp);

}
