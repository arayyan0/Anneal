///  @file     common.cpp
///  @author   Ahmed Rayyan
///  @date     December 2, 2019
///  @brief    useful functions for the simulated annealing algorithm
#include "common.hpp"


void BravaisIndicesToFlat(const uint& x, const uint&y, const uint & period1,
                                         const uint&z, const uint & period2,
                                         uint& flat_index)
{
  flat_index = x + y*period1 + z*period2;
}

long double Lorentzian(const long double& r, const long double& l){
  long double x = r - 0.5;
  long double lmax = 2*l;
  long double result = x<lmax ? 1 / ( 1.0 + pow(x/l,2) ) : 0;
  return result;
}

// long double Gaussian(const long double& s, const long double& x, const long double& l){
//   long double cutoff=1;
//   long double dist = (x-0.5)/l;
//   long double result = dist<cutoff ? s*exp(-pow(dist,2)/2.0) : 0;
//   return result;
// }

void SphericalAnglesToCubic(const long double& theta, const long double& phi,
                                 Vector3LD& some_spin)
{
  long double s = sin(theta);
  some_spin <<  s*cos(phi), s*sin(phi), cos(theta);
}

void PBCIndices(const uint& i, const uint& length, int& lower_i, int& higher_i){
  higher_i = (i+1)%length;
  lower_i = (i-1);
  if (lower_i >= 0){lower_i = lower_i%length;}
  else if (lower_i < 0){lower_i = length-1;}
}
