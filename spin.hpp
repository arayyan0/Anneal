///  @file     spins.cpp
///  @author   Ahmed Rayyan
///  @date     December 2, 2019
///  @brief    declaring the spin degree of freedom and associated functions
#ifndef SPIN_HPP
#define SPIN_HPP

#include "common.hpp"

class Spin
{
  public:
    Vector3LD VectorXYZ;
    //Really should add in VectorABC here...

    Spin();
    Spin(const long double& theta, const long double& phi);
    Spin(Vector3LD& vector);
};

Vector3LD SphericalAnglesToCubic(const long double& theta, const long double& phi);

void SpherePointPicker(Vector3LD& some_spin);

#endif // SPIN_HPP
