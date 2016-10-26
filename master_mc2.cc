/*
 *
 *
 *
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <cassert>

double f1(double x) { return log(x); }

double f2(double x) { return x; }

double Random() { return (rand() % 1000000)/1000000.0; }

double SignRandom() { return 2.0*Random()-1.0; }

double LogProb(double x, double L1, double L2, double F1, double F2, double sigma)
{
 double u1 = (f1(x)-F1)/sigma;
 double u2 = (f2(x)-F2)/sigma;
 //return (-0.5*(u1*u1+u2*u2) -sigma*(L1/u1 + L2/u2) - log(fabs(u1)) - log(fabs(u2)));
 return (-0.5*(u1*u1+u2*u2) -sigma*(L1/u1 + L2/u2));
}

void DoMonteCarlo(double seed, double F1, double F2)
{
 long int rej_x = 0, count_x = 0;
 long int rej_L1 = 0, count_L1 = 0;
 long int rej_L2 = 0, count_L2 = 0;
 double x = seed, L1 = 0.0, L2 = 0.0;
 double sigma = 5.0;
 double delta_x = 0.001, delta_L1 = 0.001, delta_L2 = 0.001;
 long int n = 0;
 while (1)
 { 
  double logp_old = LogProb(x, L1, L2, F1, F2, sigma);
  double r = Random();
  if (r < (1.0/3.0))
  {
   double xnew;
   while (1)
   {
    xnew = x + SignRandom()*delta_x;
    if (xnew >= 0.0) break;
   }
   double logp_new = LogProb(xnew, L1, L2, F1, F2, sigma);
   if (log(Random()) < (logp_new-logp_old)) { x = xnew; }
   else { rej_x++; }
   count_x++;
  }
  else if (r < (2.0/3.0))
  {
   double L1new;
   while (1)
   {
    L1new = L1 + SignRandom()*delta_L1;
    if ((L1new < 1.0) && (L1new > -100.0)) break;
   }
   double logp_new = LogProb(x, L1new, L2, F1, F2, sigma);
   if (log(Random()) < (logp_new-logp_old)) { L1 = L1new; }
   else { rej_L1++; }
   count_L1++;
  }
  else
  {
   double L2new;
   while (1)
   {
    L2new = L2 + SignRandom()*delta_L2;
    if ((L2new >= 0.0) && (L2new < 100.0)) break;
   }
   double logp_new = LogProb(x, L1, L2new, F1, F2, sigma);
   if (log(Random()) < (logp_new-logp_old)) { L2 = L2new; }
   else { rej_L2++; }
   count_L2++;
  }
  if ((n > 0) && (n % 1000 == 0))
  {
   if ( (float(rej_x)/float(count_x)) > 0.75) delta_x *= 0.998;
   if ( (float(rej_x)/float(count_x)) < 0.65) delta_x *= 1.002;
   if ( (float(rej_L1)/float(count_L1)) > 0.75) delta_L1 *= 0.998;
   if ( (float(rej_L1)/float(count_L1)) < 0.65) delta_L1 *= 1.002;
   if ( (float(rej_L2)/float(count_L2)) > 0.75) delta_L2 *= 0.998;
   if ( (float(rej_L2)/float(count_L2)) < 0.65) delta_L2 *= 1.002;
   sigma *= (1.0-1.0E-05);
  }
  if (n % 5000 == 0)
  {
   std::cout << n << "  " << x << "  (" << (f1(x)-F1)/sigma << "  " << (f2(x)-F1)/sigma << ")  " << L1 << "  " << L2;
   std::cout << "  " << float(rej_x)/float(count_x) << "  " << float(rej_L1)/float(count_L1) << "  " << float(rej_L2)/float(count_L2);
   std::cout << "  " << sigma << std::endl;
  }
  n++;
 }
}

int main()
{
 srand(time(NULL));
 std::cout << std::setprecision(10);

 double F1 = -0.25;
 double F2 = 1.10;   

 double x = 4.0; // initial seed

 //
 DoMonteCarlo(x, F1, F2);
 return 0;
}

