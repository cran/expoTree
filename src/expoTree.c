#include "expoTree.h"

SEXP expoTreeEval(SEXP parameters, SEXP times, SEXP ttypes) {
  PROTECT(parameters = AS_NUMERIC(parameters));
  PROTECT(times = AS_NUMERIC(times));
  PROTECT(ttypes = AS_INTEGER(ttypes));

  double* pars = NUMERIC_POINTER(parameters);

  int N = (int) pars[0];
  double beta = pars[1];
  double mu = pars[2];
  double psi = pars[3];
  double rho = pars[4];

  int SI = 1;
  int vf = 0;
  int rs = 1;

  double* ptimes = NUMERIC_POINTER(times);
  int* pttypes = INTEGER_POINTER(ttypes);
  int nt = LENGTH(times);
  int extant = 0;
  int maxExtant = 0;

  int i = 0;
  int ki = 0;

  SEXP p;
  PROTECT(p = NEW_NUMERIC(N+1));
  double* p0 = NUMERIC_POINTER(p);

  int root = (ptimes[nt-1] == ptimes[nt-2]) ? 1 : 0;

  double t0 = 0.0;

  for (i = nt-1; i >= 0; --i) {
    if (pttypes[i] == 0) --extant;
    else ++extant;
    if (maxExtant < extant) maxExtant = extant;
  }

  if (N < maxExtant || beta <= 0.0 || mu < 0.0 || 
      psi < 0.0 || rho < 0.0 || rho > 1.0) {
    for (i = 0; i < N; ++i) p0[i] = R_NegInf;
  } else {
    // set initial value of p
    if (extant == 0) {
      p0[0] = 0.0;
      for (i = 1; i <= N; ++i) p0[i] = psi;
      ki = 1;
      t0 = ptimes[0];
      ptimes = ptimes+1;
      pttypes = pttypes+1;
      --nt;
    } else {
      ki = extant;
      p0[0] = 0.0;
      for (i = 1; i <= N; ++i) {
        if (i < extant) p0[i] = 0.0;
        else p0[i] = pow(rho,1.*extant)*pow(1.-rho,i-extant);
      }
    }
    rExpoTree(&N,&ki,&beta,&mu,&psi,&nt,ptimes,pttypes,p0,&t0,&SI,&vf,&rs);
    if (root) {
      for (i = 0; i < N; ++i) p0[i] -= M_LN2 + log(beta) + log(1.-1./N);
    }
  }

  UNPROTECT(4);
  return p;
}

/****************************************************************************/

SEXP expoTreeSurvival(SEXP parameters, SEXP torig) {
  double fx = R_NegInf;

  PROTECT(parameters = AS_NUMERIC(parameters));
  PROTECT(torig = AS_NUMERIC(torig));

  double* pars = NUMERIC_POINTER(parameters);

  int N = (int) pars[0];
  double beta = pars[1];
  double mu = pars[2];
  double psi = pars[3];
  double rho = pars[4];

  int SI = 1;
  int vf = 0;
  int rs = 1;

  double* ptorig = NUMERIC_POINTER(torig);
  int extant = 0;
  int maxExtant = 0;

  int ttype = 1;
  double t0 = 0.0;
  SEXP p;
  PROTECT(p = NEW_NUMERIC(N+1));
  double* p0 = NUMERIC_POINTER(p);
  int ki = 0;
  int nt = 1;
  int i = 0;

  if (beta <= 0.0 || mu < 0.0 || psi < 0.0 || rho < 0.0 || rho > 1.0) {
    for (i = 0; i <= N; ++i) p0[i] = R_NegInf;
  } else {
    for (i = 0; i <= N; ++i) p0[i] = 1.0;
    rExpoTree(&N,&ki,&beta,&mu,&psi,&nt,ptorig,&ttype,p0,&t0,&SI,&vf,&rs);
    for (i = 0; i <= N; ++i) p0[i] = (p0[i] < 0) ? log(1.-exp(p0[i])) : R_NegInf;
  }

  UNPROTECT(3);
  return p;
}

/****************************************************************************/

double bdss_q(double t, double c1, double c2) {
  double q = (1-c2)*exp(-t*c1/2.0) + ((1+c2)*exp(t*c1/2));
  return 0.25*q*q;
}

/****************************************************************************/

SEXP infTreeEval(SEXP parameters, SEXP times, SEXP ttypes, SEXP survival) {
  PROTECT(parameters = AS_NUMERIC(parameters));
  PROTECT(times = AS_NUMERIC(times));
  PROTECT(ttypes = AS_INTEGER(ttypes));
  PROTECT(survival = AS_INTEGER(survival));

  double* pars = NUMERIC_POINTER(parameters);
  double beta = pars[0];
  double mu   = pars[1];
  double psi  = pars[2];
  double rho  = pars[3];

  double* ptimes = NUMERIC_POINTER(times);
  int* pttypes = INTEGER_POINTER(ttypes);
  int nt = LENGTH(times);
  int extant = 0;
  int maxExtant = 0;

  SEXP fx;
  PROTECT(fx = NEW_NUMERIC(1));
  double* lik = NUMERIC_POINTER(fx);

  int root = (ptimes[nt-1] == ptimes[nt-2]) ? 1 : 0;

  double t0 = 0.0;
  int i;

  for (i = nt-1; i >= 0; --i) {
    if (pttypes[i] == 0) --extant;
    else ++extant;
    if (maxExtant < extant) maxExtant = extant;
  }

  if (beta <= 0.0 || mu < 0.0 || psi < 0.0 || rho < 0.0 || rho > 1.0) {
    *lik = R_NegInf;
  } else {
    double c1;
    double c2;
    double p0;
    double mp = mu+psi;

    c1 = beta-mp;
    c1 = sqrt(c1*c1 + 4*beta*psi);
    c2 = -(beta-mp-2*beta*rho)/c1;
    *lik = -log(2*beta);

    if (survival) {
      p0 = exp(-c1*ptimes[nt-1])*(1.-c2);
      p0 = beta+mp+c1*(p0-(1.+c2))/(p0+1.+c2);
      p0 = p0/(2.*beta);
      *lik -= log(1.-p0);
    }

    if (extant > 0) *lik += extant*log(rho);
    for (int i = 0; i < nt; ++i) {
      if (pttypes[i]) *lik += log(2*beta/bdss_q(ptimes[i],c1,c2));
      else *lik += log(psi*bdss_q(ptimes[i],c1,c2));
    }
  }

  UNPROTECT(5);
  return fx;
}


