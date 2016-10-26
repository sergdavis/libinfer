/*
 *
 *
 *
 */

#ifndef __STATMECHANICS_H__
#define __STATMECHANICS_H__

#include "hamiltonian.h"

template <class T> class CanonicalModel: public RealFunction<T>
{
 public:
   CanonicalModel(const Hamiltonian<T> & H, double beta): H(H), beta(beta) { }

   double operator()(const T & s) const override { return -beta*H(s); }

 private:
   const Hamiltonian<T> & H;
   double beta;
};

#endif

