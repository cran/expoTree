/*
 * Select degree of Taylor approximation. 
 * ======================================
 * 
 * Adapted from MATLAB code of 
 *    Awad H. Al-Mohy and Nicholas J. Higham, October 26, 2010.
 *  Reference: A. H. Al-Mohy and N. J. Higham, Computing the action of
 *  the matrix exponential, with an application to exponential
 *  integrators. MIMS EPrint 2010.30, The University of Manchester, 2010.
 *
 * Gabriel E Leventhal, November 13, 2011.
 *
 * PARAMETERS
 *
 * n (input)    : dimension of matrix
 * fmv (intput) : matrix-matrix product function
 * b (intput)   : matrix b in x = A*b
 * m (intput)   : number of columns of b
 * m_max (input) : default = 8
 * p_max (input) : default = 55
 * prec (input)  : desired precision ('d' = double, 's' = single, 'h' = half)
 * alpha (output) : array of dimension p_max-1
 * eta (output)   : array of dimension p_max
 * M (output)     : approximation matrix (dimension (m_max,p_max))
 * shift (input) : set to true if a shift of the matrix is desired
 *                 (not implemented; shifting on matrix-vector product?)
 * bal (input)   : set to true is balancing of the matrix is desired
 *                 (not implemented; balancing on matrix-vector product?)
 * force_estm (intput) : default = false
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "expmv.h"

#ifdef THETA_HALF
#include "theta_taylor_half.h"
#elif defined THETA_SINGLE
#include "theta_taylor_single.h"
#else
#include "theta_taylor.h"
#endif

void get_theta(char prec, double* theta) {
  FILE* pFile;
  char mystring[100];
  int i;
  switch (prec) {
    case 'h': pFile = fopen("theta_taylor_half.txt","r"); break;
    case 's': pFile = fopen("theta_taylor_single.txt","r"); break;
    case 'd':
    default: pFile = fopen("theta_taylor.txt","r"); break;
  }
  i = 0;
  if (pFile == NULL) perror ("Error opening file");
  else {
    while (! feof(pFile)) {
      if (fgets(mystring, 100, pFile) != NULL) {
        theta[i] = atof(mystring);
        ++i;
      }
    }
    fclose (pFile);
  }
}

int select_taylor_degree(int n, matMat fmv, double t, double* b, int m, 
    normFunc nf, int m_max, int p_max, char prec, double* alpha, double* eta, 
    double* tm, int* mv, int* unA, int shift, char bal, int force_estm, 
    int wrklen, double* wrk, int iwrklen, int* iwrk)
{
  if (p_max < 2 || m_max > 60 || m_max + 1 < p_max*(p_max-1)) {
    return -1;
  }

  int k = 0;
  double normA = 0.0;
  double c = 0.0;
  int i;
  int j;

  double normLim = 4*theta[m_max-1]*p_max*(p_max+3)/(m_max*m);

  mv = 0;

  if (!force_estm) normA = t*nf();

  if (!force_estm && normA <= normLim) {
    *unA = 1;
    c = normA;
    for (i = 0; i < p_max-1; ++i) alpha[i] = c;
  } else {
    *unA = 0;
    for (i = 0; i < p_max; ++i) {
      normAm(n,fmv,t,i+1,&c,&k,wrklen,wrk,iwrklen,iwrk);
      c = pow(c,1./(i+1.));
      mv = mv + k;
      eta[i] = c;
    }

    for (i = 0; i < p_max-1; ++i) {
      alpha[i] = (eta[i] > eta[i+1]) ? eta[i] : eta[i+1];
    }
  }

  for (i = 0; i < m_max; ++i) {
    for (j = 0; j < p_max-1; ++j) {
      tm[j*m_max+i] = 0.0;
    }
  }

  for (j = 1; j < p_max; ++j) {
    for (i = (j+1)*j-2; i < m_max; ++i) {
      tm[(j-1)*m_max+i] = alpha[j-1]/theta[i];
    }
  }

  return 0;
}

