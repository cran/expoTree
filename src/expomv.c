#include "expomv.h"

/************************************************************/ 

void rExpoTree(int* RN, int* Rki, double* Rbeta, double* Rmu,
    double* Rpsi, int* Rn, double* times, int* ttypes, double* p,
    double* t0, int* RSImodel, int* Rvflag, int* Rrescale) 
{
  expo = (expo_type*) malloc(sizeof(expo_type));
  init_expo(expo);

  expo->N       = *RN;
  expo->ki      = *Rki;
  expo->beta    = *Rbeta;
  expo->mu      = *Rmu;
  expo->psi     = *Rpsi;
  expo->SImodel = *RSImodel;
  expo->rescale = *Rrescale;
  expo->vflag   = *Rvflag;
  expo->cutoff  = 0.0;

  int n = *Rn;
  int N = *RN+1;

  int ncol = 1;
  int p_max = expo->p_max;
  int m_max = expo->m_max;

  int memlen = (2*N+(p_max-1)*m_max);
  int wlen   = 2*p_max + (6*ncol+3)*N + ncol + m_max*(p_max-1);
  int iwlen  = 2*N + 4;

  double* wrk = (double*) malloc((memlen+wlen)*sizeof(double));
  int* iwrk = (int*) malloc(iwlen*sizeof(int));

  wrk[0] = *t0;

  expoTree(n,times,ttypes,p,memlen+wlen,wrk,iwlen,iwrk);

  free(iwrk);
  free(wrk);
  free(expo);
}

/************************************************************/ 

void expoTree(int n, double* times, int* ttypes, 
    double* p, int wrklen, double* wrk, int iwrklen, int* iwrk) 
{
  /***************** SI PARAMETERS *****************/
  double t, dt;
  int i, j;

  double nrm;
  double scale = 0.0;
  
  int info = 0;
  int one = 1;

  int ncol = 1;
  int p_max = expo->p_max;
  int m_max = expo->m_max;

  if (expo->SImodel == 0) {
    expo->lambda = &lambdaInf;
  } else {
    expo->lambda = &lambdaSI;
  }

  /********************* EXPMV *********************/

  int N = expo->N+1;

  int wrk_steps = 1;
  double wrk_dt = 0.0;
  int wrk_info = 0;
  double wrk_scale = 0.0;
  int max_wrk_steps = 1000;

  int memlen = (2*N+(p_max-1)*m_max);
  int wlen   = 2*p_max + (6*ncol+3)*N + ncol + m_max*(p_max-1);
  int iwlen  = 2*N + 4;

  if (wrklen < memlen+wlen || iwrklen < iwlen) {
    Rprintf("Not enough memory supplied!");
    return;
  }

  double* p0 = wrk;
  double* pT = wrk + N;
  double* tm = wrk + 2*N;
  double* expowrk = wrk + memlen;

  if (expo->vflag > 1) {
    Rprintf("Running expoTree with parameters:\n");
    Rprintf(" N    = %d\n",expo->N);
    Rprintf(" ki   = %d\n",expo->ki);
    Rprintf(" beta = %g\n",expo->beta);
    Rprintf(" mu   = %g\n",expo->mu);
    Rprintf(" psi  = %g\n",expo->psi);
    Rprintf(" n    = %d\n",n);
    Rprintf(" |v|  = %f\n",dnrm2_(&N,p,&one));
    Rprintf("\n");
  }

  /******** INTEGRATE DIFFERENTIAL EQUATION ********/

  t = wrk[0];
  scale = 0.0;
  nrm = 1.0;

  memcpy(p0,p,N*sizeof(double));
  nrm = dnrm2_(&N,p0,&one);
  if (expo->rescale && nrm != 1.0) {
    scale += log(nrm);
    nrm = 1./nrm;
    dscal_(&N,&nrm,p0,&one);
    nrm = dnrm2_(&N,p0,&one);
  }

  for (i = 0; i < n; ++i) {
    dt = times[i] - t;

    if (dt < 0.0) {
      Rprintf("Negative dt (%f) at time step %d! Aborting.\n",dt,i);
      Rprintf("Numer of time points: n = %d\n",n);
      Rprintf("t(0)   = %8.4e\n",times[0]);
      Rprintf("t(i)   = %8.4e\n",t);
      Rprintf("t(i+1) = %8.4e\n",times[i]);
      Rprintf("dt     = %8.4e\n",dt);
      p[0] = R_NegInf;
      return;
    }

    if (dt > 0.0) {
      expo->offset = 0;
      expo->shift = expo->trace()/N;
      wrk_steps = 1;
      wrk_info = 0;

      // save a backup of the vector
      memcpy(pT,p0,N*sizeof(double));

      while (wrk_info == 0) {
        if (wrk_steps > max_wrk_steps) {
          Rprintf("Maximum time step intervals reached.\n");
          for (j = 0; j < N; ++j) p[j] = R_NegInf;
          return;
        }

        wrk_dt = (wrk_steps > 1) ? dt/wrk_steps : dt;
        wrk_scale = 0.0;

        int k;
        for (k = 0; k < wrk_steps; ++k) {
          info = 0;
          expmv(wrk_dt,N,expo->matvec,expo->norm,expo->trace,p0+expo->offset,
                1,m_max,p_max,tm,1,'d',expo->shift,0,0,expo->vflag,&info,
                wrklen-memlen,expowrk,iwrklen,iwrk);

          // Error during calculation. Return INF.
          if (info < 0) {
            for (j = 0; j < N; ++j) p[j] = R_NegInf;
            return;
          }

          // check for negative values
          for (j = 0; j < N; ++j) if (p0[j] < 0.0) p0[j] = 0.0;

          // calculate norm of vector
          nrm = dnrm2_(&N,p0,&one);

          // validate 2-norm of the vector
          if (nrm < expo->cutoff) {
            if (expo->vflag > 1) Rprintf("Vector norm is zero. Aborting.\n");
            for (i = 0; i < N; ++i) p[i] = R_NegInf;
            return;
          }

          if (expo->rescale) {
            wrk_scale += log(nrm);
            nrm = 1./nrm;
            dscal_(&N,&nrm,p0,&one);
          }

          if (expo->vflag > 1) {
            if (wrk_steps == 1) {
              Rprintf("%6d/%1d %8f %5d %8.4e % 8.4e | %2d %2d\n",
                  i,ttypes[i],dt,expo->ki,nrm,scale+wrk_scale,iwrk[0],iwrk[1]);
            } else {
              Rprintf("%8d %8f %5d %8.4e % 8.4e %2d %2d\n",
                  k,wrk_dt,expo->ki,nrm,scale+wrk_scale,iwrk[0],iwrk[1]);
            }
          }

          if ((iwrk[0] > 40 && iwrk[0] < iwrk[1]) && wrk_steps < max_wrk_steps) {
            memcpy(p0,pT,N*sizeof(double));
            ++wrk_steps;
            wrk_info = 0;
            if (expo->vflag > 2) {
              Rprintf("Decreasing time step (%d).\n",wrk_steps);
            }
            break;
          } else {
            wrk_info = 1;
          }
        }
      }

      if (expo->rescale) scale += wrk_scale;
    }

    // update initial condition with event information
    if (i < n-1) {
      if (ttypes[i] == 1) {
        (*expo->ft)(p0,pT);
        --expo->ki;
      } else {
        (*expo->fs)(p0,pT);
        ++expo->ki;
      }
      memcpy(p0,pT,N*sizeof(double));
    }

    nrm = dnrm2_(&N,p0,&one);

    // rescale likelihoods for numerical reasons
    if (expo->rescale) {
      if (nrm > 1e20) {
        Rprintf("Problem with 2-norm of vector in rescaling!\n");
        Rprintf("Most likely ki >= N.\n");
        Rprintf("||pT|| = 0.0\n");
        Rprintf("||p0|| = %14.8e\n",dnrm2_(&N,p0,&one));
        if (expo->vflag > 2) {
          for (i = 0; i < N; ++i) {
            Rprintf("  pT(%3d) = %12.6e\n",i,pT[i]);
          }
        }
        p[0] = R_NegInf;
        return;
      }
      scale += log(nrm);
      nrm = 1./nrm;
      dscal_(&N,&nrm,p0,&one);
    }

    t = times[i];
  }

  for (i = 0; i < N; ++i) p[i] = log(p0[i])+scale;
}
