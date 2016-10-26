/*
 *
 *
 *
 */

#ifndef __GAMMA_MODEL_H__
#define __GAMMA_MODEL_H__

#include "state.h"

class GammaModel: public RealFunction<PositiveRealState>
{
 public:
   GammaModel(double k, double theta): k(k), theta(theta) { }

   double operator()(const PositiveRealState & s) const override
   {
    const double & x = s[0];
    return (-x/theta + (k-1)*log(x));
   }

 private:
   double k, theta;
};

#endif

