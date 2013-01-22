#include "expmv.h"

void expmv(double t, int n, matMat fmv, normFunc nf, traceFunc tf,
    double* b, int ncol, int m_max, int p_max, double* tm, int recalcm,
    char prec, double shift, char bal, int full_term, int prnt, int* info,
    int wrklen, double* wrk, int iwrklen, int* iwrk)
{
  int mv;
  int mvd;
  int unA;

  int i;
  int j;
  int k;

  double tt;
  // double trace;
  double tol;

  /* size of Taylor approximation */
  int tcol;

  int wlen;
  int nalen;

  double* talpha;
  double* teta;
  double* C;
  double* b1;
  double* b2;
  double* nawrk;

  wlen = 2*p_max + 2*n*ncol + m_max*(p_max-1);
  nalen = 3*n + (4*n+1)*ncol;
  if (wrklen < wlen+nalen || iwrklen < 2*n + 4) {
    if (prnt) Rprintf("Not enough workspace supplied!\n");
    *info = -2;
    return;
  }

  talpha = wrk;
  teta   = talpha + p_max;
  C      = teta + p_max;
  b1     = C + m_max*(p_max-1);
  b2     = b1 + n*ncol;
  nawrk  = b2 + n*ncol;

  switch (prec) {
    case 's': 
      tol = pow(2.,-24.);
      break;
    case 'h': 
      tol = pow(2.,-10.);
      break;
    case 'd': 
    default:
      tol = pow(2.,-53.);
      break;
  }

  // get required Taylor truncation (if not set)
  if (recalcm) {
    if (prnt > 2) Rprintf("Calculating required Taylor truncation...");
    tt = 1.;
    select_taylor_degree(n,fmv,t,b,ncol,nf,m_max,p_max,
                         prec,talpha,teta,tm,&mvd,&unA,shift,bal,0,
                         wrklen-wlen,nawrk,iwrklen,iwrk);
    if (prnt > 3) {
      for (i = 0; i < m_max; ++i) {
        for (j = 0; j < p_max-1; ++j) {
          Rprintf("%16.8e ",tm[j*m_max+i]);
        }
        Rprintf("\n");
      }
    }
    mv = mvd;
    if (prnt > 2) Rprintf("done.\n");
  } else {
    tt = t; 
    mv = 0; 
    mvd = 0;
  }

  double s = 1.0;     // cost per column
  double cost = 0.0;
  double ccost = 0.0;

  if (t == 0.0) { 
    tcol = 0; 
  } else {
    k = 0;
    for (i = 0; i < m_max; ++i) {
      for (j = 0; j < p_max-1; ++j) {
        C[k] = ceil(fabs(tt)*tm[j*(m_max)+i])*(i+1.);
        if (C[k] == 0.0) C[k] = 1./0.;
        ++k;
      }
    }

    cost = INFINITY;
    tcol = 0;
    for (i = 0; i < m_max; ++i) {
      ccost = INFINITY;
      for (j = 0; j < p_max-1; ++j) {
        if (C[i*(p_max-1)+j] < ccost) {
          ccost = C[i*(p_max-1)+j];
        }
      }
      if (ccost < cost) {
        cost = ccost;
        tcol = i+1;
      }
    }

    s = (cost/tcol > 1) ? cost/tcol : 1.0;
  }

  if (tcol == 0) {
    if (prnt) {
      Rprintf("Cannot calculate matrix exponential (under-/overflow?).\n");
      Rprintf("Returned results may be gibberish.\n");
    }
    *info = -1;
    return;
  }

  double eta = 1.0;
  if (shift != 0.0) eta = exp(t*shift/s);

  double* btmp = NULL;

  memcpy(b1,b,n*ncol*sizeof(double));

  if (prnt > 2) Rprintf("m = %2d, s = %g\n", tcol, s);

  double c1, c2;
  double bnorm;
  int ss;
  for (ss = 1; ss <= s; ++ss) {
    c1 = inf_norm(n,ncol,b1);
    if (prnt > 3) Rprintf("s = %d, ",ss);
    if (prnt > 4) Rprintf("\n");
    for (k = 1; k <= tcol; ++k) {
      fmv('n',n,ncol,t/(s*k),b1,b2);
      btmp = b1; b1 = b2; b2 = btmp; btmp = NULL;
      ++mv;
      for (i = 0; i < n; ++i) {
        for (j = 0; j < ncol; ++j) {
          b[j*n+i] += b1[j*n+i];
        }
      }
      c2 = inf_norm(n,ncol,b1);
      if (prnt > 4) Rprintf("k=%3d: %9.2e %9.2e %9.2e",k,bnorm,c1,c2);
      if (!full_term) {
        bnorm = inf_norm(n,ncol,b);
        if (prnt > 4) Rprintf(" %9.2e, \n",(c1+c2)/bnorm);
        if (c1+c2 <= tol*bnorm) {
          if (prnt == 4) {
            Rprintf("k=%3d: %9.2e %9.2e %9.2e",k,bnorm,c1,c2);
            Rprintf(" %9.2e, ",(c1+c2)/bnorm);
          }
          if (prnt > 3) Rprintf("m_actual = %2d\n", k);
          break;
        }
        c1 = c2;
      }
    }
    for (i = 0; i < n*ncol; ++i) b[i] *= eta;

    memcpy(b1,b,n*ncol*sizeof(double));
  }

  iwrk[0] = tcol;
  iwrk[1] = k;
}

double inf_norm(int n1, int n2, double* A) {
  int i, j;
  double rowsum;
  double c = 0.0;
  for (i = 0; i < n1; ++i) {
    rowsum = 0.0;
    for (j = 0; j < n2; ++j) rowsum += fabs(A[j*n1+i]);
    if (rowsum > c) c = rowsum;
  }
  return c;
}


