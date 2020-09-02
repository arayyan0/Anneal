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
  double Kx, Ky, Kz, Gx, Gy, Gz;
  const double Gp, J1, h, hTheta, hPhi;
  Eigen::Vector3d hDir;

  Parameters(const double& kx, const double& ky, const double& kz,
             const double& gx, const double& gy, const double& gz,
             const double& gammap, const double& heisenberg,
             const double& hstrength, const double& htheta, const double& hphi);
  void Anisotropy(const double& scale, const double& g, const double& a, 
                    const bool& kitaev_or_gamma, const int& sign);
  void PrintParameters(std::ostream &out);

private:
  FIELD_DIR
};

class Bond
{
public:
  double BondEnergy;
  Eigen::Vector3d MolecFieldContribution;

  Bond(const Spin& spin_i, const Spin& spin_j, const uint& bond_type, const Parameters& p);
private:
  Spin iSpin, jSpin;
  uint BondType;
  Parameters Pa;
  Eigen::Matrix3d BondHamiltonian;

  void SpecifyBondHamiltonian();
};

#endif // HAMILTONIAN_HPP
