#ifndef __EXPO_H_
#define __EXPO_H_

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>

/************************************************************/ 

double dnrm2_(int* n, double* x, int* incx);
double dscal_(int* n, double* alpha, double* x, int* incx);

/************************************************************/ 

void matFuncExpmvShifted(char trans, int n1, int n2, 
                  double alpha, double* pin, double* pout);
double matFuncOneNormShifted();
double matFuncTrace();
void transEvent(double* pin, double* pout);
void sampleEvent(double* pin, double* pout);

/************************************************************/ 

typedef struct _expo {
  int N;
  int ki;
  double beta;
  double mu;
  double psi;
  double rho;
  double gamma;
  double shift;

  int vflag;
  int rescale;
  double tol;
  double cutoff;

  int SImodel;
  double (*lambda)(int);

  int offset;
  int useLog;
  int includeZero;

  void (*matvec)(char,int,int,double,double*,double*);
  double (*norm)();
  double (*trace)();

  void (*fs)(double*,double*);
  void (*ft)(double*,double*);

  int p_max;
  int m_max;
} expo_type;

static expo_type *expo;

expo_type* init_expo(expo_type* e) {
  e->N = 0;
  e->ki = 0; 
  e->beta = 0.0;
  e->mu = 0.0;
  e->psi = 0.0;
  e->rho = 0.0;
  e->gamma = 0.0;
  e->shift = 0.0;
  e->vflag = 0;
  e->rescale = 1;
  e->tol = 0.0;
  e->cutoff = 0.0;
  e->SImodel = 1;
  e->lambda = NULL;
  e->offset = 0;
  e->useLog = 0;
  e->includeZero = 0;
  e->matvec = &matFuncExpmvShifted;
  e->norm = &matFuncOneNormShifted;
  e->trace = &matFuncTrace;
  e->ft = &transEvent;
  e->fs = &sampleEvent;
  e->p_max = 8;  // +1 than used in the loop calculation
  e->m_max = 55; // +1 than used in the loop calculation
  return e;
}

/************************************************************/ 
/* Lambda functions                                         */
/************************************************************/ 

double lambdaSI(int I) { 
  if (I > 0 && I <= expo->N) 
    return expo->beta/expo->N*(expo->N-I); 
  else 
    return 0.0;
}

/************************************************************/ 

double lambdaInf(int I) { 
  if (I > 0 && I <= expo->N) 
    return expo->beta;
  else 
    return 0.0;
}

/************************************************************/ 
/* Column and row sums                                      */
/************************************************************/ 

double matRowSum(int m, int SI) {
  double l = expo->lambda(m);
  double rs = 0.0;
  if (m > 0) rs = m*(l+expo->mu+expo->psi) 
    + (m+expo->ki)*l + (m-expo->ki)*expo->mu;
  return rs;
}

/************************************************************/ 

double matColSum(int m, int SI) {
  double cs = 0.0;
  if (m >= expo->ki) {
    cs = m*(expo->lambda(m)+expo->mu+expo->psi);
    if (m < expo->N) cs += ((m+1)+expo->ki)*expo->lambda(m+1);
    if (m > 0) cs += ((m-1)-expo->ki)*expo->mu;
  }
  return cs;
}

/************************************************************/ 

double matColSumShifted(int m, int SI) {
  double cs = 0.0;
  if (m >= expo->ki) {
    cs = fabs(-m*(expo->lambda(m)+expo->mu+expo->psi)-expo->shift);
    if (m < expo->N) cs += ((m+1)+expo->ki)*expo->lambda(m+1);
    if (m > 0) cs += ((m-1)-expo->ki)*expo->mu;
  }
  return cs;
}

/************************************************************/ 
/* Matrix-vector product required for EXPMV                 */
/* - matrix-matrix product                                  */
/* - shifted matrix-matrix product                          */
/* - shifted matrix norm                                    */
/* - matrix trace                                           */
/************************************************************/ 

void matFuncExpmv(char trans, int n1, int n2, 
                  double alpha, double* pin, double* pout) 
{
  int m;
  double a, b, c;
  for (m = 0; m <= expo->N; ++m) {
    if (m < expo->ki) pout[m] = 0.0;
    else {
      a = -m*(expo->lambda(m)+expo->psi+expo->mu)*pin[m];
      b = (m < expo->N) ? (m+expo->ki)*expo->lambda(m)*pin[m+1] : 0.0;
      c = (m > 0) ? (m-expo->ki)*expo->mu*pin[m-1] : (m-expo->ki)*expo->mu;
      pout[m] = alpha*(a+b+c);
      if (isnan(pout[m])) pout[m] = 0.0;
    }
  }
}

/************************************************************/ 

double matFuncTrace() {
  int m = 0;
  double trace = 0.0;
  for (m = expo->ki; m <= expo->N; ++m) {
    trace += -m*(expo->lambda(m)+expo->psi+expo->mu);
  }
  return trace;
}

/************************************************************/ 

void matFuncExpmvShifted(char trans, int n1, int n2, 
                  double alpha, double* pin, double* pout) 
{
  int m;
  double a, b, c, l;
  for (m = 0; m <= expo->N; ++m) {
    if (m < expo->ki) pout[m] = 0.0;
    else {
      l = expo->lambda(m);
      a = (-m*(l+expo->psi+expo->mu)-expo->shift)*pin[m];
      b = (m < expo->N) ? (m+expo->ki)*l*pin[m+1] : 0.0;
      c = (m > 0) ? (m-expo->ki)*expo->mu*pin[m-1] : (m-expo->ki)*expo->mu;
      pout[m] = alpha*(a+b+c);
      if (isnan(pout[m])) pout[m] = 0.0;
    }
  }
}


/************************************************************/ 

double matFuncOneNorm() {
  double maxM = .5*expo->N*(1.+(expo->mu+.5*expo->psi)/expo->beta)
                - .5 - .25*expo->ki;
  int maxMi = rint(maxM);
  if (maxMi < expo->ki) maxMi = expo->ki;
  if (maxMi > expo->N) maxMi = expo->N;
  return matColSum(maxMi,expo->SImodel);
}

/************************************************************/ 

double matFuncOneNormShifted() {
  double maxM = .5*expo->N*(1.+(expo->mu+.5*expo->psi)/expo->beta)
                - .5 - .25*expo->ki;
  int maxMi = rint(maxM);
  if (maxMi < expo->ki) maxMi = expo->ki;
  if (maxMi > expo->N) maxMi = expo->N;
  return matColSumShifted(maxMi,expo->SImodel);
}

/************************************************************/ 

double matFuncInfNorm() {
  double maxRow = .5*expo->N + .25*expo->ki 
                  + .25*expo->N/expo->beta*(expo->mu+expo->psi);
  int maxRi = rint(maxRow);
  if (maxRi < expo->ki) maxRi = expo->ki;
  if (maxRi > expo->N) maxRi = expo->N;
  return matRowSum(maxRi,expo->SImodel);
}

/************************************************************/ 

void transEvent(double* pin, double* pout) {
  int m;
  double l;
  for (m = 1; m < expo->N; ++m) {
    l = expo->lambda(m);
    pout[m] = 2*l*pin[m+1];
  }
  pout[expo->N] = 0.0;
}

/************************************************************/ 

void sampleEvent(double* pin, double* pout) {
  int m;
  pout[1] = 0.0;
  for (m = 2; m <= expo->N; ++m) pout[m] = expo->psi*pin[m-1];
}

#endif // __EXPO_H_
