/*
 *
 *
 *
 */

#ifndef __MAXENT_H__
#define __MAXENT_H__

#include "function.h"
#include "state.h"
#include <list>

template <class T> class MaxEntModel: public RealFunction<T>
{
 public:

    MaxEntModel() { logprior = NULL; }

    virtual ~MaxEntModel() { }

    void SetPrior(const RealFunction<T> & lprior) { logprior = &lprior; }

    void AddConstraint(const RealFunction<T> & R) { f.push_back(&R); }

    void AddConstraint(const RealFunction<T> & R, double Rval)
    {
     AddConstraint(R);
     F.push_back(Rval);
    }

    void SetParams(const State & p) { params = p; }

    double operator()(const T & x) const override
    {
     double logp = 0.0;
     int i = 0;
     for (auto it=f.begin();it!=f.end();++it)
     {
      logp -= params[i++]*(*(*it))(x);
     }
     if (logprior == NULL) return logp;
     else return logp+(*logprior)(x);
    }

 private:
   const RealFunction<T> * logprior;
   std::list<const RealFunction<T> *> f;
   std::list<double> F;
   State params;
};

#endif

