///  @file     hamiltonian.cpp
///  @author   Ahmed Rayyan
///  @date     December 2, 2019
///  @brief    declaring the spin model
#ifndef HAMILTONIAN_HPP
#define HAMILTONIAN_HPP

#include "common.hpp"

class Hamiltonia
{
public:
  Matrix3LD Hx, Hy, Hz;
  Matrix3LD SpinBasis;
  const Vector3LD hField;
  std::string ParameterOutput;

  //kga honeycomb
  Hamiltonia(const long double& phi, const long double& g, const long double& a,
             const long double& h_magnitude, const Vector3LD& h_direction);
  //jtau triangular
  Hamiltonia(const long double& jtau, const long double& lambda,
             const long double& jquad, const long double& jocto,
             const long double& h_magnitude, const Vector3LD& h_direction);
  //jtau honeycomb
  Hamiltonia(const long double& phi, const long double& theta,
             const long double& h_magnitude, const Vector3LD& h_direction);
private:
  const long double hMagnitude;
  const Vector3LD hDirection;
  Matrix3LD ReturnJTauHamiltonian(const long double& angle);

  ABC_DIR;
};


#endif // HAMILTONIAN_HPP
