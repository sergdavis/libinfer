/*
 *
 *
 *
 */

#ifndef __PSO_H__
#define __PSO_H__

#include "objfunction.h"

//
//
//
namespace pso
{

 class Minimizer
 {
  public:
    Minimizer(int params, int particles);
    virtual ~Minimizer();

    const State & Minimize(const ObjectiveFunction & obj, const State & seed, double tolerance=1.0e-10);

    double tolerance, omega, c1, c2;
  
    virtual void AdvanceParticles(const ObjectiveFunction & obj, const State & globalmin);

    virtual void OnIteration(int step, double objval, const State & current) const { }

  private:
   int nparams, nparticles, nprocs;
   double * minibuffer, * buffer;
   State current;
   State * x, *v, *localmin;
   void UpdateGlobalMinimum(const ObjectiveFunction & obj);
 };

}

#endif

