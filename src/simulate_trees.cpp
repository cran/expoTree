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

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

// ===========================================================================

struct _pars {
  int N;
  double beta;
  double mu;
  double psi;
  double t;
} pars;

// ===========================================================================

double lambda(int I) { return pars.beta*(1.-1.*I/pars.N); }

// ===========================================================================

struct Individual {
  Individual(double it = 0.0, int p = 0) 
    : itime(it), dtime(-1.0), stime(-1.0), parent(p), id(-1)
  {}

  virtual ~Individual() {}

  double itime;  // infection time
  double dtime;  // death time
  double stime;  // sample time
  int parent;    // parent individual
  int id;        // my id
  vector<int> children;
};

// ===========================================================================

int getBranchingTimes(vector<double>& times, int subclade, 
    const vector<Individual>& pop) 
{
  int ret = 0;
  // check if the leaf is sampled
  if (pop[subclade].stime >= 0.0) ret = 1;
  // check if at least one of the children is sampled
  for (int i(pop[subclade].children.size()-1); i >= 0; --i) {
    if (getBranchingTimes(times,pop[subclade].children[i],pop)) {
      if (ret) {
        times.push_back(pop[pop[subclade].children[i]].itime);
      } else {
        ret = 1;
      }
    }
  }
  return ret;
}

// ===========================================================================

string newickString(int subclade, const vector<Individual>& pop, 
    bool sampled = false) 
{
  ostringstream str;
  double t(pars.t);
  double dtime(0.0);
  if (pop[subclade].stime >= 0.0) dtime = pop[subclade].stime;
  else if (pop[subclade].dtime >= 0.0) dtime = pop[subclade].dtime;
  else dtime = t;
  if (pop[subclade].children.size() == 0) {
    str << "'" << subclade;
    if (pop[subclade].stime >= 0.0) str << "+";
    str << "':";
    str << dtime-pop[subclade].itime;
  } else {
    t = dtime;
    for (size_t i(0); i < pop[subclade].children.size(); ++i) str << "(";
    str << "'" << subclade;
    if (pop[subclade].stime >= 0.0) str << "+";
    str << "'";
    for (int i(pop[subclade].children.size()-1); i >= 0; --i) {
      str << ":" << t-pop[pop[subclade].children[i]].itime;
      str << "," << newickString(pop[subclade].children[i],pop,sampled) << ")";
      // str << "'" << subclade << "->" << pop[subclade].children[i] << "'";
      t = pop[pop[subclade].children[i]].itime;
    }
    str << ":" << t-pop[subclade].itime;
  }
  return str.str();
}

// ===========================================================================

double sim_trees(vector<Individual>& pop, list<int>& inf, vector<int>& samples,
    int max_samples, double max_time, int min_outbreak) {
  GetRNGstate();

  double t(0.0);
  
  double totalRate(0.0);
  double nextEventTime(0.0);
  double r(0.0);
  double l(0.0);

  int S(pars.N-1); // susceptibles
  int I(1);        // infecteds
  int J(0);        // sampleds
  int node(0);

  list<int>::iterator infit;

  int max_tries = 1000;
  int tries = 0;
  while (tries++ < max_tries) {
    pop.clear();
    inf.clear();
    samples.clear();
    t = 0.0;
    nextEventTime = 0.0;
    r = 0.0;
    l = 0.0;
    S = pars.N-1;
    I = 1;
    J = 0;
    pop.push_back(Individual(0.0,-1));
    pop.back().id = 0;
    inf.push_back(0);
    while (I > 0 && samples.size() < max_samples && t < max_time) {
      if (pars.psi <= 0.0 && I >= max_samples) {
        for (infit = inf.begin(); infit != inf.end(); ++infit) {
          pop[*infit].stime = t;
          samples.push_back(*infit);
        }
        break;
      }
      if (pars.N > 0) l = lambda(I);
      else l = pars.beta;
      totalRate = l*I + (pars.mu+pars.psi)*I;
      nextEventTime = t-log(unif_rand())/totalRate;
      r = unif_rand();
      node = floor(unif_rand()*I);
      infit = inf.begin();
      advance(infit,node);
      if (r < l*I/totalRate) {
        // BRANCHING EVENT
        pop.push_back(Individual());
        pop.back().itime = nextEventTime;
        pop.back().parent = *infit;
        pop.back().id = pop.size()-1;
        inf.push_back(pop.size()-1);
        pop[*infit].children.push_back(pop.size()-1);
        ++I;
      } else {
        // REMOVAL EVENT
        if (r < (l*I+pars.mu*I)/totalRate) {
          pop[*infit].dtime = nextEventTime;
        } else {
          pop[*infit].stime = nextEventTime;
          samples.push_back(*infit);
          ++J;
        }
        inf.erase(infit);
        --I;
      }

      t = nextEventTime;
    }
    if (samples.size() >= min_outbreak) break;
  }

  PutRNGstate();

  return t;
}

// ===========================================================================

void to_array(double t, vector<Individual>& pop, vector<int>& samples, 
    int len, double* x, int* xtype) {
  vector<double> times;
  getBranchingTimes(times,0,pop);
  times.push_back(0.0);
  int j(0);
  if (len >= times.size() + samples.size()) {
    for (int i(0); i < times.size(); ++i) {
      x[j] = t-times[i];
      xtype[j] = 1;
      ++j;
    }
    for (int i(0); i < samples.size(); ++i) {
      x[j] = t-pop[samples[i]].stime;
      xtype[j] = 0;
      ++j;
    }
  }
}

// ===========================================================================

void inf_times(double t, vector<Individual>& pop, 
    int len, double* x0, double* x1, int* xtype, 
    int* id, int* parent) 
{
  vector<Individual>::iterator ind(pop.begin());
  int i(0);
  while (ind != pop.end() && i < len) {
    if (ind->itime >= 0) {
      x0[i] = t-ind->itime;
      id[i] = ind->id;
      parent[i] = ind->parent;
      if (ind->stime >= 0) {
        x1[i] = t-ind->stime;
        xtype[i] = 0;
      } else if (ind->dtime >= 0) {
        x1[i] = t-ind->dtime;
        xtype[i] = -1;
      } else {
        x1[i] = 0.0;
        xtype[i] = 1;
      }
      ++i;
      ++ind;
    }
  }
}

// ===========================================================================

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
  pars.N    = *N;
  pars.beta = *beta;
  pars.mu   = *mu;
  pars.psi  = *psi;

  if (*max_time < 0) *max_time = R_PosInf;

  vector<Individual> pop;
  list<int> inf;
  vector<int> samples;

  double t = sim_trees(pop,inf,samples,*max_samples,*max_time,*min_outbreak);
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

  pars.N    = (int) (NUMERIC_POINTER(parameters)[0]);
  pars.beta = NUMERIC_POINTER(parameters)[1];
  pars.mu   = NUMERIC_POINTER(parameters)[2];
  pars.psi  = NUMERIC_POINTER(parameters)[3];

  int max_samp = INTEGER_VALUE(max_samples);
  int min_outb = INTEGER_VALUE(min_outbreak);
  double max_t = NUMERIC_VALUE(max_time);
  
  if (max_t < 0) max_t = R_PosInf;

  vector<Individual> pop;
  list<int> inf;
  vector<int> samples;

  double t = sim_trees(pop,inf,samples,max_samp,max_t,min_outb);

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

  const char* names[7] = { "times", "ttypes", "itimes", "dtimes", "dtypes", "id", "parent" };
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


