///  @file     hamiltonian.cpp
///  @author   Ahmed Rayyan
///  @date     December 2, 2019
///  @brief    defining the spin model
#include "hamiltonian.hpp"

Parameters::Parameters(const long double& kx, const long double& ky, const long double& kz,
                       const long double& gx, const long double& gy, const long double& gz,
                       const long double& gammap, const long double& heisenberg,
                       const long double& hstrength, const long double& htheta, const long double& hphi):
Kx(kx), Ky(ky), Kz(kz),
Gx(gx), Gy(gy), Gz(gz),
Gp(gammap), J1(heisenberg), h(hstrength), hTheta(htheta), hPhi(hphi)
{
  long double s = sin(hTheta/180*pi);
  hDir = s*cos(hPhi/180*pi)*A_Dir + s*sin(hPhi/180*pi)*B_Dir + cos(hTheta/180*pi)*C_Dir;
}

void Parameters::Anisotropy(const long double& scale, const long double& g, const long double& a,
                              const bool& kitaev_or_gamma, const int& sign)
// recall that 0 is false, 1 is true
// kitaev_or_gamma == 0 -> kitaev, kitaev_or_gamma == 1 -> gamma
{
  if (kitaev_or_gamma){
    Gx = sign*scale*(2*(1-a)*(1-g));
    Gy = sign*scale*(2*(1-a)*g);
    Gz = sign*scale*(1+2*a);
  } else{
    Kx = sign*scale*(2*(1-a)*(1-g));
    Ky = sign*scale*(2*(1-a)*g);
    Kz = sign*scale*(1+2*a);
  }
}

void Parameters::PrintParameters(std::ostream &out)
{
  out << "------------------------Hamiltonian Parameters------------------------\n";
  out << "Kx Ky Kz\n";
  out << Kx << " " << Ky << " " << Kz << "\n";
  out << "Gx Gy Gz\n";
  out << Gx << " " << Gy << " " << Gz << "\n";
  out << "Gp\n";
  out << Gp << "\n";
  out << "J1\n";
  out << J1 << "\n";
  out << "h\n";
  out << h << "\n";
  out << "h_theta\n";
  out << hTheta << "\n";
  out << "h_phi\n";
  out << hPhi << "\n";
}

Bond::Bond(const Spin& spin_i, const Spin& spin_j, const uint& bond_type, const Parameters& p):
iSpin(spin_i), jSpin(spin_j), BondType(bond_type), Pa(p)
{
  SpecifyBondHamiltonian();
  MolecFieldContribution = -BondHamiltonian*jSpin.VectorXYZ + Pa.h*Pa.hDir/3.0;
  BondEnergy = - iSpin.VectorXYZ.transpose().dot(MolecFieldContribution + Pa.h*Pa.hDir/3.0);
}

void Bond::SpecifyBondHamiltonian()
{
  Matrix3LD int_matrix;
  switch (BondType+1)
  {
    case 1:
      int_matrix << Pa.Kx + Pa.J1, Pa.Gp , Pa.Gp,
                           Pa.Gp , Pa.J1 , Pa.Gx,
                           Pa.Gp , Pa.Gx , Pa.J1;
      break;
    case 2:
      int_matrix << Pa.J1 , Pa.Gp        , Pa.Gy,
                     Pa.Gp, Pa.Ky + Pa.J1, Pa.Gp,
                     Pa.Gy,   Pa.Gp      , Pa.J1;
      break;
    case 3:
      int_matrix << Pa.J1 , Pa.Gz,         Pa.Gp,
                     Pa.Gz, Pa.J1,         Pa.Gp,
                     Pa.Gp, Pa.Gp, Pa.J1 + Pa.Kz;
      break;
  }
  BondHamiltonian = int_matrix;
}
