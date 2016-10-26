/*
 *
 *
 *
 */

#ifndef __ISING2D_H__
#define __ISING2D_H__

#include "state.h"
#include "binarymatrix.h"

class Ising2DState: public BinaryMatrix, public State
{
 public: 
   Ising2DState(int N): BinaryMatrix(N) { }

   void Mutate(double delta) override 
   { 
    i0 = rand() % Size();
    j0 = rand() % Size();
    int v0 = Get(i0, j0);
    Set(i0, j0, 1-v0);
   } 

   void UndoMutation() override 
   {
    int v0 = Get(i0, j0);
    Set(i0, j0, 1-v0);
   }

 private:
   int i0, j0;
};

#endif

