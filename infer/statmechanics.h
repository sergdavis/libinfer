/*
 *
 *
 *
 */

#ifndef __STATMECHANICS_H__
#define __STATMECHANICS_H__

#include "maxent.h"
#include "hamiltonian.h"

template <class T> class CanonicalModel: public MaxEntModel<T>
{
 public:
   CanonicalModel(const Hamiltonian<T> & H, double beta): MaxEntModel<T>(1)
   {
    MaxEntModel<T>::SetParams(State({beta}));
    MaxEntModel<T>::AddConstraint(H);
   }
};

#endif

