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

#define MAXSTEPS 1000
#define XSCALE   1.0

double Random() { return (rand() % 1000000)/1000000.0; }

double SignRandom() { return 2.0*Random()-1.0; }

double DoMonteCarlo(double seed, double L1, double L2, double & F1, double & F2, double & var1, double & var2, double & covar)
{
 long int n = 1500000;
 long int rej = 0, count = 0;
 double f1av = 0.0, f2av = 0.0;
 double f1f1av = 0.0, f2f2av = 0.0, f1f2av = 0.0;
 double x = seed;
 double MCSTEP = 7.0;
 for (int nstep=0;nstep < n; ++nstep)
 { 
  double xnew = x + SignRandom()*MCSTEP*XSCALE;
  if (log(Random()) < (-L1*(f1(xnew)-f1(x))-L2*(f2(xnew)-f2(x))))
  {
   x = xnew;
  }
  else { rej++; }
  if ((nstep > 200000) && (nstep % 200 == 0))
  {
   f1av += f1(x);
   f2av += f2(x); 
   f1f1av += f1(x)*f2(x); 
   f2f2av += f2(x)*f2(x); 
   f1f2av += f1(x)*f2(x); 
   count++;
  }
  if ((nstep > 0) && (nstep % 2000 == 0))
  {
   if ( (float(rej)/float(nstep)) > 0.7) MCSTEP *= 0.998;
   if ( (float(rej)/float(nstep)) < 0.7) MCSTEP *= 1.002;
  }
 }
 f1av /= double(count);
 f2av /= double(count);
 F1 = f1av;
 F2 = f2av;
 f1f1av /= double(count);
 f2f2av /= double(count);
 f1f2av /= double(count);
 var1 = f1f1av - f1av*f1av;
 var2 = f2f2av - f2av*f2av;
 covar = f1f2av - f1av*f2av;
 return 100.0*float(rej)/float(n);
}

double DoReverseMonteCarlo(double seed, double F1, double F2, double & L1, double & L2)
{
 long int rej_x = 0, count_x = 0;
 long int rej_L1 = 0, count_L1 = 0;
 long int rej_L2 = 0, count_L2 = 0;
 double rejstat_x, rejstat_L1, rejstat_L2;
 double x = seed, sigma = 0.0001;
 double delta_x = 0.01, delta_L1 = 4.0, delta_L2 = 4.0;
 long int nstep = 0;
 while (1)
 { 
  if (Random() < 0.5)
  {
   double H = -L1*f1(x)-L2*f2(x) - 0.5*(1.0/(sigma*sigma))*(pow(F1-f1(x),2.0)+pow(F2-f2(x),2.0));
   double xnew = x + SignRandom()*delta_x;
   double Hnew = -L1*f1(xnew)-L2*f2(xnew) - 0.5*(1.0/(sigma*sigma))*(pow(F1-f1(xnew),2.0)+pow(F2-f2(xnew),2.0));
   if (log(Random()) < (Hnew-H)) { x = xnew; }
   else { rej_x++; }
   count_x ++;
  }
  else
  {
   if (Random() < 0.5)
   {
    double H = -L1*f1(x);
    double L1new = 1.01;
    while (L1new >= 1.0)
    {
     L1new = L1 + SignRandom()*delta_L1;
    }
    double Hnew = -L1new*f1(x);
    if (log(Random() < (Hnew-H))) { L1 = L1new; }
    else { rej_L1++; }
    count_L1++;
   }
   else
   {
    double H = -L2*f2(x);
    double L2new = -1.0; 
    while (L2new < 0.0)
    {
     L2new = L2 + SignRandom()*delta_L2;
    }
    double Hnew = -L2new*f2(x);
    if (log(Random() < (Hnew-H))) { L2 = L2new; }
    else { rej_L2++; }
    count_L2++;
   }
  }
  //if ((nstep > 2000000) && (nstep % 200 == 0))
  //{
  // f1av += f1(x);
  // f2av += f2(x); 
  // count++;
  //}
  if ((nstep > 0) && (nstep % 1000 == 0))
  {
   rejstat_x = float(rej_x)/float(count_x);
   rejstat_L1 = float(rej_L1)/float(count_L1);
   rejstat_L2 = float(rej_L2)/float(count_L2);
  /*
   if ( rejstat_x > 0.75) delta_x *= 0.998;
   if ( rejstat_x < 0.65) delta_x *= 1.002;
   if ( rejstat_L1 > 0.75) delta_L1 *= 0.998;
   if ( rejstat_L1 < 0.65) delta_L1 *= 1.002;
   if ( rejstat_L2 > 0.75) delta_L2 *= 0.998;
   if ( rejstat_L2 < 0.65) delta_L2 *= 1.002;
   //
  */
   sigma *= 0.998;
  }
  if (nstep % 200 == 0) 
     std::cout << x << "  " << L1 << "  " << L2 << "  " << sigma << "  |||||   " << rejstat_x << "  " << rejstat_L1 << "  " << rejstat_L2 << "\n";
  nstep++;
 }
}

int main()
{
 srand(time(NULL));
 std::cout << std::setprecision(10);
 // These correspond to L1 = -0.6 and theta = 1.42857

 double F1 = -0.23011711743689367;
 double F2 = 1.1206967632639822;   

 //double L1 = -0.6, L2 = 1.42857;
 double L1 = -0.8, L2 = 1.1;
 double eta = 0.01;
 double x = 4.0; // initial seed
 double f1av, f2av, var1, var2, covar;
 double rejstat = DoMonteCarlo(x, L1, L2, f1av, f2av, var1, var2, covar);
 double error0 = 0.5*((fabs(f1av-F1)/fabs(F1)) + (fabs(f2av-F2)/fabs(F2)));
 double beta = 0.5, delta_L1 = 0.01, delta_L2 = 0.01;
 int niter = 0;
 while (1)
 {
  std::cout << niter << "  " << rejstat << "  " << L1 << "  " << L2 << "  " << "  " << error0 << "  " << beta << "  " << eta << "\n";
  double zzz = eta*Random();
  double L1new = L1 - 2.0*zzz*((F1-f1av)*var1 + (F2-f2av)*covar);
  double L2new = L2 - 2.0*zzz*((F1-f1av)*covar + (F2-f2av)*var2);
  //double L1new = L1 + SignRandom()*delta_L1;
  //double L2new = L2 + SignRandom()*delta_L2;
  rejstat = DoMonteCarlo(x, L1, L2, f1av, f2av, var1, var2, covar);
  double error = 0.5*((fabs(f1av-F1)/fabs(F1)) + (fabs(f2av-F2)/fabs(F2)));
  if ((error < error0) || (log(Random()) < -beta*(error-error0)))
  {
   error0 = error;
   L1 = L1new;
   L2 = L2new;
  }
  else { }
  //L1 = 0.9*L1 - 2.0*eta*((F1-f1av)*var1 + (F2-f2av)*covar);
  //L2 = 0.9*L2 - 2.0*eta*((F1-f1av)*covar + (F2-f2av)*var2);
  //eta *= 0.99;
  beta *= 1.01;
  eta *= 0.999;
  niter++;
 }
 return 0;
}

