///  @file     hamiltonian.cpp
///  @author   Ahmed Rayyan
///  @date     December 2, 2019
///  @brief    declaring the spin model
#ifndef HAMILTONIAN_HPP
#define HAMILTONIAN_HPP

#include "common.hpp"
#include "spin.hpp"

struct Parameters
{
public:
  long double Kx, Ky, Kz, Gx, Gy, Gz;
  const long double Gp, J1, h, hTheta, hPhi;
  Vector3LD hDir;

  Parameters(const long double& kx, const long double& ky, const long double& kz,
             const long double& gx, const long double& gy, const long double& gz,
             const long double& gammap, const long double& heisenberg,
             const long double& hstrength, const long double& htheta, const long double& hphi);
  void Anisotropy(const long double& scale, const long double& g, const long double& a,
                    const bool& kitaev_or_gamma, const int& sign);
  void PrintParameters(std::ostream &out);

private:
  FIELD_DIR
};

class Bond
{
public:
  long double BondEnergy;
  Vector3LD MolecFieldContribution;

  Bond(const Spin& spin_i, const Spin& spin_j, const uint& bond_type, const Parameters& p);
private:
  Spin iSpin, jSpin;
  uint BondType;
  Parameters Pa;
  Matrix3LD BondHamiltonian;

  void SpecifyBondHamiltonian();
};

#endif // HAMILTONIAN_HPP
