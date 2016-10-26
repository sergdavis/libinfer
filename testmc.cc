/*
 *
 *
 *
 */

#include "infer/metropolis.h"
#include "infer/state.h"

class SingleState: public State
{
 public:
   SingleState(double x): State(1) { (*this)[0] = xold = x; }

   void Mutate(double delta) override 
   { 
    xold = (*this)[0];
    while (1)
    {
     (*this)[0] = xold + SignRandom()*delta;
     if ((*this)[0] >= 0.0) break;
    }
   }

   void UndoMutation() override { (*this)[0] = xold; }

 private:
  double xold;
};

class GammaModel: public RealFunction<SingleState>
{
 public:
   GammaModel(double k, double theta): k(k), theta(theta) { }

   double operator()(const SingleState & s) const override
   {
    const double & x = s[0];
    return (-x/theta + (k-1)*log(x));
   }

 private:
   double k, theta;
};

WindowAverage x_av(300000);
WindowAverage logx_av(300000);

class MyMetropolis: public Metropolis<SingleState>
{
 public:
    MyMetropolis() { }
    bool OnProductionStep(SingleState & s, long int step) override
    {  
     if (step % 500 == 0)
     {
      x_av.Add(s[0]);
      logx_av.Add(log(s[0]));
      std::cout << step << "  " << x_av.Average() << "  " << logx_av.Average() << "  " << delta << "\n";
     }
     return true; 
    }
};

int main()
{
 InitRandom();
 MyMetropolis m;

 GammaModel G(4.5, 3.0);
 SingleState seed(10.0);
 m.SetDelta(20.0);
 m.Simulate(G, seed, 5000000, 3000000, 100);
 std::cout << "Rejection: " << m.RejectionRate() << "\n";
 std::cout << "<x> = " << x_av.Average() << "\n";
 std::cout << "<ln x> = " << logx_av.Average() << "\n";

 return 0;
}

