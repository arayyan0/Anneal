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
    Eigen::Vector3d VectorXYZ;
    //Really should add in VectorABC here...

    Spin();
    Spin(const double& theta, const double& phi);
    Spin(Eigen::Vector3d& vector);
};

void SpherePointPicker(Spin& some_spin);

#endif // SPIN_HPP
