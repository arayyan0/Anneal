///  @file     spins.cpp
///  @author   Ahmed Rayyan
///  @date     December 2, 2019
///  @brief    defining the spin degree of freedom and associated functions
#include "spin.hpp"

Spin::Spin(){}

Spin::Spin(const long double& theta, const long double& phi)
{
  long double s = sin(theta);
  Vector3LD v(s*cos(phi), s*sin(phi), cos(theta));
  VectorXYZ = v;
}

Spin::Spin(Vector3LD& v)
{
  v.normalize();
  VectorXYZ = v;
}

void SpherePointPicker(Spin& some_spin)
{
  long double u = MyRandom::unit_interval(MyRandom::RNG);
  long double v = MyRandom::unit_interval(MyRandom::RNG);
  long double theta = acos(2*u-1); long double phi = 2*pi*v;
  some_spin = Spin(theta, phi);
}
