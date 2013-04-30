#ifndef __EXPO_TYPE_H__
#define __EXPO_TYPE_H__

typedef struct _expo {
  int N;        /* maximum population size */
  int ki;       /* current number of lineages */

  double beta;  /* infection/speciation rate */
  double mu;    /* recovery/extinction rate */
  double psi;   /* per-lineage sampling rate */
  double rho;   /* proability of discovery at present */
  double K;     /* real-valued carrying capacity */

  int     parVecLen; /* number of parameter sets */
  int     curPar;    /* current parameter set */
  double* NVec;      /* carrying capacity parameters */
  double* betaVec;   /* infection rate parameters */
  double* muVec;     /* recovery rate parameters */
  double* psiVec;    /* sampling rate parameters */

  double gamma;
  double shift;

  int vflag;      /* verbosity flag */
  int rescale;    /* rescale probability vector at each iteration */
  double tol;     /* numerical tolerance */
  double cutoff;  /* numerical cutoff for zero */

  int SImodel;    /* use desity-dependent model */
  double (*lambda)(int); /* infection rate function */

  int offset;
  int useLog;
  int includeZero;

  double* mat;

  void (*matvec)(char,int,int,double,double*,double*);
  double (*norm)();
  double (*trace)();

  void (*ft)(double*,double*);
  void (*fs)(double*,double*);
  void (*fs2)(double*,double*);

  int p_max;
  int m_max;
  int N_max;
} expo_type;

#endif /* __EXPO_TYPE_H__ */
