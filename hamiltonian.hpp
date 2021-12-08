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
  Vector3LD hField;
  std::string ParameterOutput;

  // kg anisotropy
  Hamiltonia(const long double& phi, const long double& g, const long double& a,
             const long double& h_magnitude, const Vector3LD& h_direction);
  // multipole hamiltonian
  Hamiltonia(Eigen::Array<long double, 5, 1>& params, Vector3LD& h_direction);
private:
  long double hMagnitude;
  Vector3LD hDirection;
  Matrix3LD BondDependentQuadQuad(const long double& angle);
  Matrix3LD BondDependentQuadOcto(const long double& angle);

  ABC_DIR;
};


#endif // HAMILTONIAN_HPP
