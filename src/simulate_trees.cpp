#include "simulate_trees.h"

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
    double t, bool sampled) 
{
  ostringstream str;
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
      str << "," << newickString(pop[subclade].children[i],pop,t,sampled) << ")";
      // str << "'" << subclade << "->" << pop[subclade].children[i] << "'";
      t = pop[pop[subclade].children[i]].itime;
    }
    str << ":" << t-pop[subclade].itime;
  }
  return str.str();
}

// ===========================================================================

double sim_trees(const Pars& pars, vector<Individual>& pop, list<int>& inf, 
    vector<int>& samples, int max_samples, double max_time, int min_outbreak) {
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
      if (pars.N > 0) l = lambda(I,pars);
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
  // loop through all individuals of the population
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

