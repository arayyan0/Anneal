///  @file     common.cpp
///  @author   Ahmed Rayyan
///  @date     December 2, 2019
///  @brief    useful functions for the simulated annealing algorithm
#include "common.hpp"

void ThreeDBravaisIndicesToFlat(const uint& w, const uint&x, const uint & periodx,
                                               const uint&y, const uint & periody,
                                               const uint&z, const uint & periodz,
                                               uint& flat_index)
{
  flat_index = w + x*periodx + y*periody + z*periodz;
}

long double Lorentzian(const long double& r, const long double& l)
{
  long double x = r - 0.5;
  long double lmax = 2*l;
  long double result = x<lmax ? 1 / ( 1.0 + pow(x/l,2) ) : 0;
  return result;
}

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

bool EqualWithinEpsilon(const long double& a,
                        const long double& b,
                        const long double& eps)
{
  return abs(a - b) <= eps;
}

bool IsThisFloatAnIntegerWithinEpsilon(const long double& a,
                                       const long double& eps)
{
  return EqualWithinEpsilon(a, round(a), eps);
}
