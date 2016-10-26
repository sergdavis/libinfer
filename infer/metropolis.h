/*
 *
 *
 *
 */

#ifndef __METROPOLIS_H__
#define __METROPOLIS_H__

#include <iostream>
#include <vector>
#include <cmath>

#include "random.h"
#include "block.h"
#include "average.h"
#include "function.h"

inline double CheckDrift(double * w, int n)
{
 int S = 0;
 for (int i=0;i<n-1;++i) 
 {
  if (fabs(w[i+1]-w[i]) < 1.0E-6) continue;
  S += ( (w[i+1]-w[i] < 0.0) ? -1 : 1 );
 }
 return double(S)/double(n-1);
}

inline double RejectionSigmoid(double R)
{
 return (1.0+0.2*(1.0/(1.0+exp(7.0*(R-0.7)))-0.5));
}

template <class T> class Metropolis
{
 public:
    Metropolis() { Reset(); ClearProperties(); delta = 1.0; }

    void Reset() { MCrej = 0; MCcount = 0; }

    double RejectionRate() { return 100.0*double(MCrej)/MCcount; }

    virtual bool OnBurnInStep(T & x, long int step) { return true; }
    
    virtual bool OnProductionStep(T & x, long int step) { return true; }

    void SetDelta(double d) { delta = d; }

    bool Move(const RealFunction<T> & logprob, T & s, double & logprob_s)
    {
     s.Mutate(delta);
     double logpnew = logprob(s);
     double Q = logpnew-logprob_s;
     if (log(Random()) < Q) { logprob_s = logpnew; return true; }
     else { s.UndoMutation(); return false; }
    }

    void Sweep(const RealFunction<T> & prob, T & s, double & probs)
    {
     for (int q=0;q<s.Size();++q)
     {
      if (!Move(prob, s, probs)) MCrej++;
      MCcount++;
     }
    }

    void ClearProperties() 
    {
     properties.clear();
     datasets.clear();
    }

    void AddProperty(const RealFunction<T> & property, Block<double> & data)
    {
     properties.push_back(&property);
     datasets.push_back(&data);
    }

    void Simulate(const RealFunction<T> & logmodel, T & x, long int steps, long int burnin=0, long int adapt=0)
    {
     double logprobx = logmodel(x);
     Reset();
     long int n = 0;
     double * buffer = new double[NBUFFER];
     for (int i=0;i<NBUFFER;++i) buffer[i] = 0.0;
     int bc = 0;
     SimpleAverage xav(100);
     adapt_steps = adapt;
     std::cerr << "DEBUG Started Burn-in stage\n";
     while (1)
     {
      Sweep(logmodel, x, logprobx);
      if ((adapt_steps > 0) && (n % adapt_steps == 0))
      {
       delta *= RejectionSigmoid(RejectionRate()/100.0);
      }
      if (!OnBurnInStep(x, n)) { delete [] buffer; return; }
      xav.Add(logprobx);
      if (xav.Full())
      { 
       if (bc == NBUFFER) bc = 0;
       buffer[bc] = xav.Average();
       bc++; 
       double drift = CheckDrift(buffer, bc); 
       if (((fabs(drift) > 0.0) && (fabs(drift) < fabs(buffer[bc-1])*CONVERGENCE_THRESHOLD))) break;
      }
      n++;
      if (n >= burnin) break;
     }
     std::cerr << "DEBUG Ended Burn-in stage, rejection rate: " << RejectionRate() << "\n";
     for (long int n=0;n<steps;++n)
     {
      auto it2 = datasets.begin();
      for (auto it = properties.begin();it != properties.end();it++)
      {
       RealFunction<T> & property = *(*it);
       Block<double> & data = *(*it2);
       data[n] = property(x);
       it2++;
      }
      Sweep(logmodel, x, logprobx);
      if (!OnProductionStep(x, n)) { return; }
     }
    }
  
   double delta;

 private:
   long int MCrej;
   long int MCcount;
   long int adapt_steps;
   std::vector< RealFunction<T> * > properties;
   std::vector< Block<double> * > datasets;

   const int NBUFFER = 10000;
   const double CONVERGENCE_THRESHOLD = 0.001;
};

#endif

