#include "expomv.h"

expo_type* expo = NULL;

/************************************************************/ 

void rExpoTree(double* RN, int* Rki, double* Rbeta, double* Rmu,
    double* Rpsi, int* Rn, int* parVecLen, double* times, 
    int* ttypes, double* p, double* t0, int* RSImodel, 
    int* Rvflag, int* Rrescale) 
{
  // allocate parameter structure and initialize
  expo = (expo_type*) malloc(sizeof(expo_type));
  init_expo(expo);

  expo->N       = (int) ceil(RN[0]); /* total population size */
  expo->K       = RN[0];             /* carrying capacity */
  expo->ki      = *Rki;              /* current number of lineages */
  expo->beta    = Rbeta[0];          /* branching rate */
  expo->mu      = Rmu[0];            /* extinction rate */
  expo->psi     = Rpsi[0];           /* sampling rate */
  expo->SImodel = *RSImodel;         /* Model: DD (SI=1), INF (SI=0) */
  expo->rescale = *Rrescale;         /* rescale probability vector */
  expo->vflag   = *Rvflag;           /* verbosity level */
  expo->cutoff  = 0.0;               /* precision for zero */

  expo->parVecLen = *parVecLen;      /* number of parameter sets */
  expo->curPar  = 0;                 /* current parameter set */
  expo->NVec    = RN;                /* N parameters */
  expo->betaVec = Rbeta;             /* beta parameters */
  expo->muVec   = Rmu;               /* mu parameters */
  expo->psiVec  = Rpsi;              /* psi parameters */

  expo->offset = 0;                  /* don't calculate for zero entries */

  int n = *Rn;                       /* number of events in the tree */

  /* get maximum carrying capacity */
  int j;
  int N = (int) ceil(RN[0]);
  for (j = 0; j < *parVecLen; ++j) {
    if (N < (int) ceil(RN[j])) N = (int) ceil(RN[j]);
  }
  expo->N_max = N;
  N += 1;

  expo->mat = (double*) malloc(3*N*sizeof(double));

  int ncol = 1;                  /* EXPMV parameters */
  int p_max = expo->p_max;
  int m_max = expo->m_max;

  /* calculate required memory and allocate */
  int memlen = (2*N+(p_max-1)*m_max);
  int wlen   = 2*p_max + (6*ncol+3)*N + ncol + m_max*(p_max-1);
  int iwlen  = 2*N + 4;
  double* wrk = (double*) malloc((memlen+wlen)*sizeof(double));
  int* iwrk = (int*) malloc(iwlen*sizeof(int));

  /* set present time */
  wrk[0] = *t0;

  /* call algorithm */
  expoTree(n,times,ttypes,p,memlen+wlen,wrk,iwlen,iwrk);

  free(expo->mat);
  expo->mat = NULL;

  /* clean up */
  free(iwrk);
  free(wrk);
  free(expo);

  expo = NULL;
}

/************************************************************/ 

void expoTree(int n, double* times, int* ttypes, 
    double* p, int wrklen, double* wrk, int iwrklen, int* iwrk) 
{
  /***************** SI PARAMETERS *****************/
  double t = 0.0;
  double dt = 0.0;

  double nrm = 0.0;         /* vector norm */
  double scale = 0.0;       /* log scaling factor */
  
  int info = 0;             /* info variable */
  int one = 1;              /* integer 1 (FORTRAN corpse) */

  int ncol = 1;             /* EXPMV parameters */
  int p_max = expo->p_max;
  int m_max = expo->m_max;

  /* choose lambda function */

  if (expo->SImodel == 0) {
    expo->lambda = &lambdaInf;
  } else {
    expo->lambda = &lambdaSI;
  }

  /********************* EXPMV *********************/

  int N = expo->N_max + 1;  /* maximum dimension */
  int Ncur = expo->N + 1;   /* current vector length */

  /* split up calculation when precision is bad */
  int    wrk_steps = 1;
  double wrk_dt = 0.0;
  int    wrk_info = 0;
  double wrk_scale = 0.0;
  int    max_wrk_steps = 1000;

  /* calculate required memory and check */
  int memlen = (2*N+(p_max-1)*m_max);
  int wlen   = 2*p_max + (6*ncol+3)*N + ncol + m_max*(p_max-1);
  int iwlen  = 2*N + 4;
  if (wrklen < memlen+wlen || iwrklen < iwlen) {
    Rprintf("Not enough memory supplied!");
    return;
  }

  /* set pointers for easier access to variables */
  double* p0 = wrk;
  double* pT = wrk + N;
  double* tm = wrk + 2*N;
  double* expowrk = wrk + memlen;

  /* verbose output */
  if (expo->vflag > 1) {
    Rprintf("Running expoTree with parameters:\n");
    Rprintf(" N    = %d\n",Ncur);
    Rprintf(" ki   = %d\n",expo->ki);
    Rprintf(" beta = %g\n",expo->beta);
    Rprintf(" mu   = %g\n",expo->mu);
    Rprintf(" psi  = %g\n",expo->psi);
    Rprintf(" n    = %d\n",n);
    Rprintf(" |v|  = %f\n",dnrm2_(&Ncur,p,&one));
    Rprintf("\n");
  }

  /******** INTEGRATE DIFFERENTIAL EQUATION ********/

  /* initialize variables */
  t = wrk[0];
  scale = 0.0;
  nrm = 1.0;

  /* copy probability vector and calculate norm */
  memcpy(p0,p,N*sizeof(double));
  nrm = dnrm2_(&Ncur,p0,&one);

  /* rescale probability vector if requested */
  if (expo->rescale && nrm != 1.0) {
    scale += log(nrm);
    nrm = 1./nrm;
    dscal_(&Ncur,&nrm,p0,&one);
    nrm = dnrm2_(&Ncur,p0,&one);
  }

  int Nold;
  /* loop over tree events */
  int i, j;
  for (i = 0; i < n; ++i) {
    /* get time interval between events */
    dt = times[i] - t;

    if (dt < 0.0) {
      Rprintf("Negative dt (%f) at time step %d! Aborting.\n",dt,i);
      Rprintf("Numer of time points: n = %d\n",n);
      Rprintf("t(0)   = %8.4e\n",times[0]);
      Rprintf("t(i)   = %8.4e\n",t);
      Rprintf("t(i+1) = %8.4e\n",times[i]);
      Rprintf("dt     = %8.4e\n",dt);
      p[1] = R_NegInf;
      return;
    }

    /* don't do anything for dt = 0 */
    if (dt > 0.0) {
      expo->shift = expo->trace()/Ncur;  /* shift matrix */
      wrk_steps = 1;
      wrk_info = 0;

      initMat(); /* setup matrix */

      /* save a backup of the vector */
      memcpy(pT,p0,Ncur*sizeof(double));

      while (wrk_info == 0) {
        if (wrk_steps > max_wrk_steps) {
          Rprintf("Maximum time step intervals reached.\n");
          for (j = 0; j < Ncur; ++j) p[j] = R_NegInf;
          return;
        }

        wrk_dt = (wrk_steps > 1) ? dt/wrk_steps : dt;
        wrk_scale = 0.0;

        int k;
        for (k = 0; k < wrk_steps; ++k) {
          info = 0;
          expmv(wrk_dt,Ncur,expo->matvec,expo->norm,expo->trace,p0+expo->offset,
                1,m_max,p_max,tm,1,'d',expo->shift,0,0,expo->vflag,&info,
                wrklen-memlen,expowrk,iwrklen,iwrk);

          // Error during calculation. Return INF.
          if (info < 0) {
            for (j = 0; j < Ncur; ++j) p[j] = R_NegInf;
            return;
          }

          // check for negative values
          for (j = 0; j < Ncur; ++j) {
            if (p0[j] < 0.0) {
              if (expo->vflag > 1) {
                Rprintf("Negative probabilities! p(%d) = %g\n",j,p0[j]);
              }
              p0[j] = 0.0;
            }
          }

          // calculate norm of vector
          nrm = dnrm2_(&Ncur,p0,&one);

          // validate 2-norm of the vector
          if (nrm < expo->cutoff) {
            if (expo->vflag > 1) Rprintf("Vector norm is zero. Aborting.\n");
            for (j = 0; j < Ncur; ++j) p[j] = R_NegInf;
            return;
          }

          if (expo->rescale) {
            wrk_scale += log(nrm);
            nrm = 1./nrm;
            dscal_(&Ncur,&nrm,p0,&one);
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
            memcpy(p0,pT,Ncur*sizeof(double));
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
      int m;
      switch (ttypes[i]) {
        case 1:
          // branching event
          (*expo->ft)(p0,pT);
          --expo->ki;
          break;

        case 0:
          // sampling-removal event
          (*expo->fs)(p0,pT);
          ++expo->ki;
          break;

        case 2:
          // sampling-only event
          ++expo->ki;
          (*expo->fs2)(p0,pT);
          break;

        case 3:
          // rate shift
          if (++expo->curPar < expo->parVecLen) {
            if (expo->vflag > 1) {
              Rprintf("Rate shift!\n");
              Rprintf("  OLD: N = %d, beta = %f, mu = %f, psi = %f\n",expo->N,expo->beta,expo->mu,expo->psi);
            }
            Nold = expo->N;
            expo->N    = expo->NVec[expo->curPar];
            expo->beta = expo->betaVec[expo->curPar];
            expo->mu   = expo->muVec[expo->curPar];
            expo->psi  = expo->psiVec[expo->curPar];
            Ncur = expo->N + 1;
            if (Nold > expo->N) {
              /* carrying capacity decreased */
              for (j = expo->N+1; j <= Nold; ++j) p0[j] = 0.0;
            } else if (Nold < expo->N) {
              for (j = Nold+1; j <= expo->N; ++j) p0[j] = 0.0;
            }
            if (expo->vflag > 1) {
              Rprintf("  NEW: N = %d, beta = %f, mu = %f, psi = %f\n",expo->N,expo->beta,expo->mu,expo->psi);
            }
          } else {
            if (expo->vflag > 1) {
              Rprintf("Not enough parameters supplied to shift rates. Ignoring rate shift!\n");
            }
          }

        case 4:
          ++expo->ki;
          memcpy(pT,p0,N*sizeof(double));
          for (j = 0; j < expo->ki; ++j) pT[j] = 0.0;
          break;

        case 99:
        default:
          memcpy(pT,p0,N*sizeof(double));
          break;
      }
      memcpy(p0,pT,N*sizeof(double));
    }

    nrm = dnrm2_(&Ncur,p0,&one);

    // rescale likelihoods for numerical reasons
    if (expo->rescale) {
      if (nrm > 1e20) {
        Rprintf("Problem with 2-norm of vector in rescaling!\n");
        Rprintf("Most likely ki >= N.\n");
        Rprintf("||pT|| = 0.0\n");
        Rprintf("||p0|| = %14.8e\n",dnrm2_(&Ncur,p0,&one));
        if (expo->vflag > 2) {
          for (j = 0; j < Ncur; ++j) {
            Rprintf("  pT(%3d) = %12.6e\n",i,pT[j]);
          }
        }
        p[0] = R_NegInf;
        return;
      }
      scale += log(nrm);
      nrm = 1./nrm;
      dscal_(&Ncur,&nrm,p0,&one);
    }

    t = times[i];
  }

  for (j = 0; j < Ncur; ++j) {
    p[j] = ((p0[j] < 1.0) ? log(p0[j]) : 0.0) + scale;
  }
  if (expo->vflag > 1) {
    Rprintf("p0(%d) = %g, scale = %g\n",1,p0[1],scale);
    Rprintf(" p(%d) = %g\n",1,p[1]);
  }
}
