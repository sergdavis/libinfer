/*
 *
 *
 *
 */

#ifndef __STATE_H__
#define __STATE_H__

#include <iostream>
#include <cassert>
#include <initializer_list>

//
//
//
class State
{
 public:
   State() { buffer = NULL; N = 0; }

   State(const std::initializer_list<double> & L)
   {
    N = L.size();
    Initialize(N);
    int i = 0;
    for (auto it=L.begin();it!=L.end();++it) buffer[i++] = *it;
   }
 
   State(int N): N(N)
   {
    Initialize(N);
    for (int i=0;i<N;++i) buffer[i] = 0.0;
   }

   State(const State & s)
   {
    N = s.Size();
    Initialize(N);
    for (int i=0;i<N;++i) buffer[i] = s[i];
   }

   void Initialize(int N)
   {
    this->N = N;
    buffer = new double[N];
   }

   virtual ~State() { delete [] buffer; }

   const State & operator=(const State & s)
   {
    if (&s != this)
    {
     delete [] buffer;
     N = s.Size();
     Initialize(N);
     for (int i=0;i<N;++i) buffer[i] = s[i];
    }
    return s;
   }

   State operator*(double a) const
   {
    State s = (*this);
    for (int i=0;i<N;++i) s[i] *= a;
    return s;
   } 

   State operator+(const State & a) const
   {
    State s = (*this);
    for (int i=0;i<N;++i) s[i] += a[i];
    return s;
   } 

   const State & operator+=(const State & a)
   {
    for (int i=0;i<N;++i) buffer[i] += a[i];
    return (*this);
   } 

   State operator-(const State & a) const
   {
    State s = (*this);
    for (int i=0;i<N;++i) s[i] -= a[i];
    return s;
   } 

   int Size() const { return N; }

   double operator[](unsigned long int i) const { return buffer[i]; }

   double & operator[](unsigned long int i) { return buffer[i]; }

   void Mutate(double delta) { assert (false); }

   void UndoMutation() { assert (false); }

  private:
    int N;
    double * buffer;
 };

inline std::ostream & operator<<(std::ostream & os, const State & s)
{
 os << "<";
 for (int q=0;q<s.Size();++q) os << s[q] << ", ";
 os << ">";
 return os;
}

#endif

