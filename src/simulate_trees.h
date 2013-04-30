#ifndef __simulate_trees_h__
#define __simulate_trees_h__

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

#ifdef FAKER
#include <GSLRng.h>
#else
#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#endif

class Pars {
public:
  Pars() {}
  virtual ~Pars() {}

  double& b() { return beta; }
  double& u() { return mu; }
  double& s() { return psi; }

  int N;
  double beta;
  double mu;
  double psi;
  double t;
};

// ===========================================================================

inline double lambda(int I, const Pars& pars) { 
  return pars.beta*(1.-1.*I/pars.N);
}

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
    const vector<Individual>& pop);

string newickString(int subclade, const vector<Individual>& pop, 
    double t, bool sampled = false);

double sim_trees(const Pars& pars, vector<Individual>& pop, list<int>& inf, vector<int>& samples,
    int max_samples, double max_time, int min_outbreak);

void to_array(double t, vector<Individual>& pop, vector<int>& samples, 
    int len, double* x, int* xtype);

void inf_times(double t, vector<Individual>& pop, 
    int len, double* x0, double* x1, int* xtype, 
    int* id, int* parent);

#endif // __simulate_trees_h__

