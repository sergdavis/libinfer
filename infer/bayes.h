/*
 *
 *
 *
 */

#ifndef __BAYES_H__
#define __BAYES_H__

#include "state.h"

template <class T> class BayesModel: public RealFunction<T>
{
 public:
    BayesModel() { logprior = NULL; loglikelihood = NULL; }

    virtual ~BayesModel() { }

    void SetPrior(const RealFunction<T> & lprior) { logprior = &lprior; }

    void SetLikelihood(const RealFunction<T> & llikelihoodR) { loglikelihood = &llikelihood; }

    double operator()(const T & x) const override
    {
     if (logprior == NULL) return (*loglikelihood)(x);
     else return (*loglikelihood)(x)+(*logprior)(x);
    }

 private:
   const RealFunction<T> * logprior;
   const RealFunction<T> * loglikelihood;
};

#endif

