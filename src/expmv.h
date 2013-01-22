/*
 * input
 * -----
 *
 * t         : time
 * A         : matrix of derivatives
 * b         : starting vector
 * M         : taylor degree
 * prec      : required accuracy (double, single, half)
 * shift     : shift matrix (default true)
 * bal       : balance matrix (default false)
 *               N = none
 *               P = permute only
 *               S = scale only
 *               B = scale + permute
 * full_term : default false
 * prnt      : default false
 *
 *
 *
 * output
 * ------
 *  mv  : number of matrix-vector products
 *  mvd : (mvd < mv) number of MV products that were used for norm
 *        estimation
 *
 * ORIGINAL FUNCTION DESCRIPTION (MATLAB CODE)
 * ===========================================
 *
 * function [f,s,m,mv,mvd,unA] = ...
 *  %EXPMV   Matrix exponential times vector or matrix.
 *  %   [F,S,M,MV,MVD] = EXPMV(t,A,B,[],PREC) computes EXPM(t*A)*B without
 *  %   explicitly forming EXPM(t*A). PREC is the required accuracy, 'double',
 *  %   'single' or 'half', and defaults to CLASS(A).
 *  %   A total of MV products with A or A^* are used, of which MVD are
 *  %   for norm estimation.
 *  %   The full syntax is
 *  %     [f,s,m,mv,mvd,unA] = expmv(t,A,b,M,prec,shift,bal,full_term,prnt).
 *  %   unA = 1 if the alpha_p were used instead of norm(A).
 *  %   If repeated invocation of EXPMV is required for several values of t
 *  %   or B, it is recommended to provide M as an external parameter as
 *  %   M = SELECT_TAYLOR_DEGREE(A,m_max,p_max,prec,shift,bal,true).
 *  %   This also allows choosing different m_max and p_max.
 *
 *  %   Reference: A. H. Al-Mohy and N. J. Higham, Computing the action of
 *  %   the matrix exponential, with an application to exponential
 *  %   integrators. MIMS EPrint 2010.30, The University of Manchester, 2010.
 *
 *  %   Awad H. Al-Mohy and Nicholas J. Higham, October 26, 2010.
 *
 *
 *
 */
/*
 * PARAMETRS
 *
 * t (input) : scalar in y = exp(tA)*b
 * fmv (intput) : matrix-matrix product function
 * nf (input)   : function to calculate matrix 1-norm
 * tf (input)   : function to calculate trace
 * b (input)    : matrix b in y = exp(tA)*b
 * M (input)    : length of Taylor truncation 
 *                (will be calculated if M <= 0 or M > n)
 * prec (input) : desired precision
 * shift (input) : should the matrix be shifted?
 * bal (intput)  : should the matrix be balanced?
 *
 */

#ifndef __EXPMV_H_
#define __EXPMV_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <R.h>

typedef void (*matMat)(char trans, int n1, int n2, double alpha, double* x, double* y);
typedef double (*normFunc)();
typedef double (*traceFunc)();

double inf_norm(int n1, int n2, double* A);

void expmv(double t, int n, matMat fmv, normFunc nf, traceFunc tf,
    double* b, int ncol, int tmrow, int tmcol, double* tm, int recalcm,
    char prec, double shift, char bal, int full_term, int prnt, int* info,
    int wrklen, double* wrk, int iwrklen, int* iwrk);

int select_taylor_degree(int n, matMat fmv, double t, double* b, int m, 
    normFunc nf, int m_max, int p_max, char prec, double* alpha, double* eta, 
    double* M, int* mv, int* unA, int shift, char bal, int force_estm,
    int wrklen, double* wrk, int iwrklen, int* iwrk);

void normAm(int n, matMat fmv, double t, int m, double* c, int* mv,
            int wrklen, double* wrk, int iwrklen, int* iwrk);

void afun_power(char trans, matMat fmv, double t, int n1, int n2,
                int m, double* x, double* wrk);

void dlacon_(int* n, double* v, double* x, int* isgn, double* est, int* kase);

void dlacn1_(int* n, int* t, double* v, double* x, int* ldx, double* xold, 
    int* ldxold, double* wrk, double* h, int* ind, int* indh, double* est, 
    int* kase, int* iseed, int* info);

double dnrm2_(int*,double*,int*);

#endif
