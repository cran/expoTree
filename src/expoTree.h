#ifndef __EXPOTREE_H__
#define __EXPOTREE_H__

#include <math.h>
#include <R.h>
#include <Rdefines.h>

void rExpoTree(int* RN, int* Rki, double* Rbeta, double* Rmu,
    double* Rpsi, int* Rn, double* times, int* ttypes, 
    double* p, double* t0, int* RSImodel, int* Rvflag, int* Rrescale);

void expoTree(int n, double* times, int* ttypes, 
    double* p, int wrklen, double* wrk, int iwrklen, int* iwrk);

SEXP expoTreeEval(SEXP parameters, SEXP times, SEXP ttypes);
SEXP expoTreeSurvival(SEXP parameters, SEXP torig);

/* N = Inf reimplementation */
/* ======================== */

double bdss_q(double t, double c1, double c2);

SEXP infTreeEval(SEXP parameters, SEXP times, SEXP ttypes, SEXP survival);


#endif
