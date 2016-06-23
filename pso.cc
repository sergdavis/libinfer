/*
 *
 *
 *
 */

#include "infer/pso.h"
#include "infer/random.h"

#include <mpi.h>
#include <cassert>

using namespace pso;

Minimizer::Minimizer(int N, int Nparticles)
{
 nprocs = 0;
 MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
 int rank;
 MPI_Comm_size(MPI_COMM_WORLD, &rank);
 srand(time(0)+rank);
 nparticles = Nparticles/nprocs;
 x = new State[nparticles];
 v = new State[nparticles];
 localmin = new State[nparticles];
 nparams = N;
 minibuffer = new double[nparams+1];
 buffer = new double[nprocs*(nparams+1)];
 omega = 0.7;
 c1 = c2 = 1.4;
}

Minimizer::~Minimizer()
{
 MPI_Barrier(MPI_COMM_WORLD);
 delete [] buffer;
 delete [] minibuffer;
 delete [] localmin;
 delete [] x;
 delete [] v;
 MPI_Finalize();
}

void Minimizer::UpdateGlobalMinimum(const ObjectiveFunction & obj)
{
 int nprocs = 0, rank;
 MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
 MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 for (int i=0;i<nparticles;++i)
 {
  //std::cout << "DEBUG rank=" << rank << ": x[" << i << "]=" << x[i] << " with obj " << obj(x[i]) << "\n";
  if (obj(x[i]) < obj(localmin[i])) localmin[i] = x[i];
 }
 State globalmin_proc(localmin[0]);
 for (int i=0;i<nparticles;++i)
 {
  if (obj(localmin[i]) < obj(globalmin_proc)) globalmin_proc = localmin[i];
 }
 State globalmin(globalmin_proc);
 // Gather all globalmin_proc to rank 0
 for (int i=0;i<nparams;++i) minibuffer[i] = globalmin_proc[i];
 minibuffer[nparams] = obj(globalmin_proc);
 double globalmin_value = buffer[nparams];
 int i = 0, globalmin_pos = 0;
 MPI_Barrier(MPI_COMM_WORLD);
 MPI_Gather(minibuffer, nparams+1, MPI_DOUBLE, buffer, nparams+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
 // Choose globalmin among all globalmin_proc 
 while (i < (nprocs-1)*(nparams+1)) 
 {
  if (buffer[i+nparams] < globalmin_value)
  {
   globalmin_value = buffer[i+nparams];
   globalmin_pos = i;
  } 
  i += (nparams+1);
 }
 for (int i=0;i<nparams;++i) minibuffer[i] = buffer[globalmin_pos+i];
 // Broadcast globalmin from 0 to all ranks
 MPI_Barrier(MPI_COMM_WORLD);
 MPI_Bcast(minibuffer, nparams, MPI_DOUBLE, 0, MPI_COMM_WORLD);
 for (int i=0;i<nparams;++i) globalmin[i] = minibuffer[i]; 
 current = globalmin;
}

void Minimizer::AdvanceParticles(const ObjectiveFunction & obj, const State & globalmin)
{
 for (int i=0;i<nparticles;++i)
 {
  double r1 = Random();
  double r2 = Random();
  v[i] = v[i]*omega + (localmin[i]-x[i])*c1*r1 + (globalmin-x[i])*c2*r2;
  x[i] += v[i];
 }
}

const State & Minimizer::Minimize(const ObjectiveFunction & obj, const State & seed, double tolerance)
{
 nparams = seed.Size();
 int n = 0;
 double oldval = obj(seed);
 for (int i=0;i<nparticles;++i)
 {
  x[i].Initialize(nparams);
  v[i].Initialize(nparams);
  for (int q=0;q<nparams;++q)
  {
   x[i][q] = 10.0*SignRandom();
   v[i][q] = SignRandom();
  }
  localmin[i] = x[i];
 }
 x[0] = seed;
 while (n < 100000)
 {
  UpdateGlobalMinimum(obj);
  double objc = obj(current);
  assert (objc <= oldval);
  if (objc < tolerance) break;
  oldval = objc;
  AdvanceParticles(obj, current);
  n++;
 }
 MPI_Barrier(MPI_COMM_WORLD);
 return current;
}

