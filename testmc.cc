/*
 *
 *
 *
 */

#include "infer/metropolis.h"
#include "infer/state.h"

class CanonicalModel: public RealFunction<State>
{
 public:
   double operator()(const State & s) const
   {
    return 0.0;
   }
};

int main()
{
 Metropolis<State> m;

 CanonicalModel P;

 State seed(100);
 
 m.Simulate(P, seed, 1000);

 return 0;
}

