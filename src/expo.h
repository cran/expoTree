#ifndef __EXPO_H_
#define __EXPO_H_

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>

#include "expo_type.h"

/************************************************************/ 

double dnrm2_(int* n, double* x, int* incx);
double dscal_(int* n, double* alpha, double* x, int* incx);

/************************************************************/ 

void initMat();

double lambdaSI(int I);
double lambdaInf(int I);

double matRowSum(int m);
double matColSum(int m);
double matColSumShifted(int m, int SI);

void matFuncExpmv(char trans, int n1, int n2, 
                  double alpha, double* pin, double* pout);
void matFuncExpmvShifted(char trans, int n1, int n2, 
                  double alpha, double* pin, double* pout);

double matFuncOneNorm();
double matFuncOneNormShifted();
double matFuncInfNorm();
double matFuncTrace();

void transEvent(double* pin, double* pout);
void sampleEvent(double* pin, double* pout);
void sampleEventNoShift(double* pin, double* pout);

/************************************************************/ 

extern expo_type* expo;
expo_type* init_expo(expo_type* e);

#endif // __EXPO_H_
