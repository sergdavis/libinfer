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
    Metropolis() { Reset(); }

    void Reset() { MCrej = 0; MCcount = 0; delta = 2.0; }

    double RejectionRate() { return 100.0*double(MCrej)/MCcount; }

    virtual bool OnBurnInStep(T & x, long int step) { return true; }
    
    virtual bool OnProductionStep(T & x, long int step) { return true; }

    bool Move(const RealFunction<T> & prob, T & s, double & probs)
    {
     double d = Gaussian(0.0, delta);
     s.Mutate(d);
     double g1 = -(d*d)/(2.0*delta*delta);
     double g2 = -(d*d)/(2.0*delta*delta);
     double pnew = prob(s);
     double Q = log(pnew)-log(probs)+g1-g2; // these cancel...???
     if (log(Random()) < Q) { probs = pnew; return true; }
     else { s.UndoMutation(); return false; }
    }

    void Sweep(const RealFunction<T> & prob, T & s, double & probs)
    {
     int d = s.Size();
     for (int q=0;q<d;++q)
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

    void Simulate(const RealFunction<T> & model, T & x, long int steps)
    {
     double probx = model(x);
     Reset();
     long int n = 0;
     double * buffer = new double[NBUFFER];
     for (int i=0;i<NBUFFER;++i) buffer[i] = 0.0;
     int bc = 0;
     SimpleAverage xav(100);
     std::cerr << "DEBUG Started Burn-in stage\n";
     while (1)
     {
      Sweep(model, x, probx);
      if (n % 50 == 0)
      {
       double f = RejectionSigmoid(RejectionRate()/100.0);
       delta *= f;
      }
      if (!OnBurnInStep(x, n)) { delete [] buffer; return; }
      xav.Add(log(probx));
      if (xav.Full()) 
      { 
       if (bc == NBUFFER) bc = 0;
       buffer[bc] = xav.Average();
       bc++; 
       double drift = CheckDrift(buffer, bc); 
       if (((fabs(drift) > 0.0) && (fabs(drift) < fabs(buffer[bc-1])*CONVERGENCE_THRESHOLD))) break;
      }
      n++;
      if (n == MAX_BURNIN) break;
     }
     std::cerr << "DEBUG Ended Burn-in stage\n";
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
      Sweep(model, x, probx);
      if (!OnProductionStep(x, n)) { return; }
     }
    }

 private: 
   double delta;
   long int MCrej;
   long int MCcount;
   std::vector< RealFunction<T> * > properties;
   std::vector< Block<double> * > datasets;

   const int NBUFFER = 10000;
   const double CONVERGENCE_THRESHOLD = 0.001;
   const int MAX_BURNIN = 200000;
};

#endif

