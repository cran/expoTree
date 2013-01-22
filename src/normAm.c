/*
  function [c,mv] = normAm(A,m)
  %NORMAM   Estimate of 1-norm of power of matrix.
  %   NORMAM(A,m) estimates norm(A^m,1).
  %   If A has nonnegative elements the estimate is exact.
  %   [C,MV] = NORMAM(A,m) returns the estimate C and the number MV of
  %   matrix-vector products computed involving A or A^*.

  %   Reference: A. H. Al-Mohy and N. J. Higham, A New Scaling and Squaring
  %   Algorithm for the Matrix Exponential, SIAM J. Matrix Anal. Appl. 31(3):
  %   970-989, 2009.

  %   Awad H. Al-Mohy and Nicholas J. Higham, September 7, 2010.
 *
 *
 * Calculate power of matrix-vector product
 *   trans : 'N' no, 'T' transpose
 *   fmv   : matrix-vector product
 *   n     : dimension of square matrix
 *   m     : power
 *   x     : input vector (dim = n)
 *   wrk   : workspace (dim = 2*n)
 *           on output the first n element of wrk are the result of y = (A^m)*x
 */

#include <stdlib.h>
#include <math.h>
#include "expmv.h"

void afun_power(char trans, matMat fmv, double t, int n1, int n2,
                int m, double* x, double* wrk) 
{
  int i;
  double* y = wrk;
  double* yout = wrk+n1*n2;
  double* ytmp;
  memcpy(y,x,n1*n2*sizeof(double));
  for (i = 0; i < m; ++i) {
    fmv(trans,n1,n2,t,y,yout);
    // memcpy(y,yout,n1*n2*sizeof(double));
    ytmp = y; y = yout; yout = ytmp;
  }
  if (y != wrk) memcpy(wrk,y,n1*n2*sizeof(double));
}

void normAm(int n, matMat fmv, double t, int m, double* c, int* mv,
            int dwrklen, double* dwrk, int iwrklen, int* iwrk)
{
  int tcol = 1; // Number of columns used by DLACN1
  char trans = 'N';

  // int memlen = 3*n + (4*n+1)*tcol;
  // double* wrk = (double*) malloc(memlen*sizeof(double));
  // int imemlen = 2*n+4;
  // int* iwrk = (int*) malloc(imemlen*sizeof(int));

  double* mvwrk = dwrk;
  double* v     = mvwrk + 2*n*tcol;
  double* x     = v + n;
  double* xold  = x + n*tcol;
  double* wrk   = xold + n*tcol;
  double* h     = wrk + tcol;

  int* ind   = iwrk;
  int* indh  = ind + n;
  int* iseed = indh + n;

  int kase = 0;
  int info = 0;

  iseed[0] = 153;
  iseed[1] = 1673;
  iseed[2] = 2;
  iseed[3] = 3567;

  mv = 0;

  *c = 0.0;
  dlacn1_(&n,&tcol,v,x,&n,xold,&n,wrk,h,ind,indh,c,&kase,iseed,&info);

  while (kase != 0) {
    if (kase == 1) trans = 'N';
    else if (kase == 2) trans = 'T';
    // call matrix-matrix product
    afun_power(trans,fmv,t,n,tcol,m,x,mvwrk);
    memcpy(x,mvwrk,n*tcol*sizeof(double));
    dlacn1_(&n,&tcol,v,x,&n,xold,&n,wrk,h,ind,indh,c,&kase,iseed,&info);
    mv += m*tcol;
  }

  // free(iwrk);
  // free(wrk);
}
