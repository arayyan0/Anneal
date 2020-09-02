///  @file     spins.cpp
///  @author   Ahmed Rayyan
///  @date     December 2, 2019
///  @brief    defining the spin degree of freedom and associated functions
#include "spin.hpp"

Spin::Spin(){}

Spin::Spin(const double& theta, const double& phi)
{
  double s = sin(theta);
  Eigen::Vector3d v(s*cos(phi), s*sin(phi), cos(theta));
  VectorXYZ = v;
}

Spin::Spin(Eigen::Vector3d& v)
{
  v.normalize();
  VectorXYZ = v;
}

void SpherePointPicker(Spin& some_spin)
{
  double u = MyRandom::unit_interval(MyRandom::RNG); double v =  MyRandom::unit_interval(MyRandom::RNG);
  double theta = acos(2*u-1); double phi = 2*pi*v;
  some_spin = Spin(theta, phi);
}
