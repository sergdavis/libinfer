/*
 *
 *
 *
 */

#ifndef __MAXENT_H__
#define __MAXENT_H__

#include "function.h"
#include "state.h"

template <class T> class FlatPrior: public RealFunction<T>
{
 public:
   double operator()(const T & x) const override { return 0.0; }
};

template <class T> class MaxEntModel: public RealFunction<T>
{
 public:

    MaxEntModel(int nparams)
    {
     f = new const RealFunction<T> *[nparams];
     param_index = 0;
     logprior = NULL;
    }

    virtual ~MaxEntModel() { delete [] f; }

    void SetPrior(const RealFunction<T> & lprior) { logprior = &lprior; }

    void AddConstraint(const RealFunction<T> & R) { f[param_index++] = &R; }

    void SetParams(const State & p) { params = p; }

    double operator()(const T & x) const override
    {
     double logp = 0.0;
     for (int i=0;i<params.Size();++i) logp -= params[i]*(*(f[i]))(x);
     if (logprior == NULL) return logp;
     else return logp+(*logprior)(x);
    }

 private:
   const RealFunction<T> * logprior;
   const RealFunction<T> ** f;
   State params;
   int param_index;
};

#endif

