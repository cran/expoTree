#include <iostream>
#include <vector>
#include <list>
#include <cmath>
#include <typeinfo>
#include <string>
#include <sstream>
#include <algorithm>
#include <getopt.h>
using namespace std;

#include "simulate_trees.h"

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

extern "C" {
  int RSimTrees(int* N, double* beta, double* mu, double* psi,
      int* max_samples, int* min_outbreak, double* max_time,
      int* maxlen, double* times, int* ttypes);

  SEXP RSimEpi(SEXP parameters, SEXP max_samples, SEXP min_outbreak, SEXP max_time);
}

// ===========================================================================

int RSimTrees(int* N, double* beta, double* mu, double* psi,
      int* max_samples, int* min_outbreak, double* max_time,
      int* maxlen, double* times, int* ttypes) {
  Pars pars;
  pars.N    = *N;
  pars.mu   = *mu;
  pars.psi  = *psi;
  pars.b() = *beta;

  if (*max_time < 0) *max_time = R_PosInf;

  vector<Individual> pop;
  list<int> inf;
  vector<int> samples;

  double t = sim_trees(pars,pop,inf,samples,*max_samples,*max_time,*min_outbreak);
  to_array(t,pop,samples,*maxlen,times,ttypes);

  return 1;
}

// ===========================================================================

SEXP RSimEpi(SEXP parameters, SEXP max_samples, SEXP min_outbreak, SEXP max_time)
{
  PROTECT(parameters = AS_NUMERIC(parameters));
  PROTECT(max_samples = AS_INTEGER(max_samples));
  PROTECT(min_outbreak = AS_INTEGER(min_outbreak));
  PROTECT(max_time = AS_NUMERIC(max_time));

  Pars pars;
  pars.N    = (int) (NUMERIC_POINTER(parameters)[0]);
  pars.b()  = NUMERIC_POINTER(parameters)[1];
  pars.mu   = NUMERIC_POINTER(parameters)[2];
  pars.psi  = NUMERIC_POINTER(parameters)[3];

  int max_samp = INTEGER_VALUE(max_samples);
  int min_outb = INTEGER_VALUE(min_outbreak);
  double max_t = NUMERIC_VALUE(max_time);
  
  if (max_t < 0) max_t = R_PosInf;

  vector<Individual> pop;
  list<int> inf;
  vector<int> samples;

  double t = sim_trees(pars,pop,inf,samples,max_samp,max_t,min_outb);

  int maxlen = 2*samples.size();
  int itlen  = pop.size();

  SEXP times = PROTECT(NEW_NUMERIC(maxlen));
  SEXP ttypes = PROTECT(NEW_INTEGER(maxlen));
  SEXP itimes = PROTECT(NEW_NUMERIC(itlen));
  SEXP dtimes = PROTECT(NEW_NUMERIC(itlen));
  SEXP dtypes = PROTECT(NEW_INTEGER(itlen));
  SEXP id = PROTECT(NEW_INTEGER(itlen));
  SEXP parent = PROTECT(NEW_INTEGER(itlen));

  double *ptimes = NUMERIC_POINTER(times);
  int    *pttypes = INTEGER_POINTER(ttypes);
  double *pitimes = NUMERIC_POINTER(itimes);
  double *pdtimes = NUMERIC_POINTER(dtimes);
  int    *pdtypes = INTEGER_POINTER(dtypes);
  int    *pid = INTEGER_POINTER(id);
  int    *pparent = INTEGER_POINTER(parent);

  for (int i(0); i < itlen; ++i) {
    pitimes[i] = 0.0;
    pdtimes[i] = 0.0;
    pdtypes[i] = 0;
    pid[i] = 0;
    pparent[i] = 0;
  }

  to_array(t,pop,samples,maxlen,ptimes,pttypes);
  inf_times(t,pop,itlen,pitimes,pdtimes,pdtypes,pid,pparent);

  const char* names[7] = { "times", "ttypes", "itimes", "dtimes", 
                           "dtypes", "id", "parent" };
  SEXP res_names = PROTECT(allocVector(STRSXP,7));
  for (int i(0); i < 7; ++i) SET_STRING_ELT(res_names,i,mkChar(names[i]));

  SEXP res = PROTECT(allocVector(VECSXP,7));
  SET_VECTOR_ELT(res,0,times);
  SET_VECTOR_ELT(res,1,ttypes);
  SET_VECTOR_ELT(res,2,itimes);
  SET_VECTOR_ELT(res,3,dtimes);
  SET_VECTOR_ELT(res,4,dtypes);
  SET_VECTOR_ELT(res,5,id);
  SET_VECTOR_ELT(res,6,parent);
  setAttrib(res, R_NamesSymbol, res_names);

  UNPROTECT(13);
  return res;
}


