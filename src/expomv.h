#ifndef __EXPOMV_H__
#define __EXPOMV_H__

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <R.h>

#include "expo.h"
#include "expmv.h"

void rExpoTree(int* RN, int* Rki, double* Rbeta, double* Rmu,
    double* Rpsi, int* Rn, double* times, int* ttypes, double* p,
    double* t0, int* RSImodel, int* Rvflag, int* Rrescale);

void expoTree(int n, double* times, int* ttypes, 
    double* p, int wrklen, double* wrk, int iwrklen, int* iwrk);

#endif
