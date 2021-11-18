///  @file     hamiltonian.cpp
///  @author   Ahmed Rayyan
///  @date     December 2, 2019
///  @brief    defining the spin model
#include "hamiltonian.hpp"

Hamiltonia::Hamiltonia(const long double& phi, const long double& g, const long double& a,
                       const long double& h_magnitude, const Vector3LD& h_direction):
                       hMagnitude(h_magnitude),
                       hDirection(h_direction),
                       hField(h_magnitude*h_direction)
{
  const long double G0 = sin(phi*pi);
  const long double K0 = cos(phi*pi);

  const long double gx = +G0*(2*(1-a)*(1-g));
  const long double gy = +G0*(2*(1-a)*g);
  const long double gz = +G0*(1+2*a);

  const long double kx = -K0*(2*(1-a)*(1-g));
  const long double ky = -K0*(2*(1-a)*g);
  const long double kz = -K0*(1+2*a);

  const long double gp = 0;
  const long double j1 = 0;

  Matrix3LD matrix1, matrix2, matrix3;
  matrix1  <<  j1+kx,    gp,    gp,
                  gp,    j1,    gx,
                  gp,    gx,    j1;
  matrix2  <<     j1,    gp,    gy,
                  gp, j1+ky,    gp,
                  gy,    gp,    j1;
  matrix3  <<     j1,    gz,    gp,
                  gz,    j1,    gp,
                  gp,    gp, j1+kz;

  Hx = matrix1;
  Hy = matrix2;
  Hz = matrix3;
  SpinBasis = Matrix3LD::Identity();

  std::stringstream paramout;
  paramout << std::fixed << std::setprecision(14);
  paramout << "------------------------Hamiltonian Parameters------------------------\n";
  paramout << "Kx Ky Kz\n";
  paramout << kx << " " << ky << " " << kz << "\n";
  paramout << "Gx Gy Gz\n";
  paramout << gx << " " << gy << " " << gz << "\n";
  paramout << "Gp\n";
  paramout << gp << "\n";
  paramout << "J1\n";
  paramout << j1 << "\n";
  paramout << "hMagnitude\n";
  paramout << hMagnitude << "\n";
  paramout << "hDirection\n";
  paramout << hDirection.transpose() << "\n";
  ParameterOutput = paramout.str();
}

Hamiltonia::Hamiltonia(Vector3LD& params, const uint entry,
                       const long double& h_magnitude, const Vector3LD& h_direction):
                       hMagnitude(h_magnitude),
                       hDirection(h_direction),
                       hField(h_magnitude*h_direction)
{
  double jtau, jb, jquad, jocto;
  if (entry == 0){
    jtau  = cos(params(0)*pi);
    jb    = sin(params(0)*pi)*cos(params(1)*pi);
    jquad = sin(params(0)*pi)*sin(params(1)*pi)*cos(params(2)*pi);
    jocto = sin(params(0)*pi)*sin(params(1)*pi)*sin(params(2)*pi);
  } else if (entry == 1){
    jtau  = 1;
    jb    = params(0)*jtau;
    jquad = params(1)*jtau;
    jocto = params(2)*jtau;
  }

  Matrix3LD matrix1;
  matrix1  <<   jquad,     0,     0,
                    0, jocto,     0,
                    0,     0, jquad;

  long double pz = 0;
  long double px = 2*pi/3.0;
  long double py = 4*pi/3.0;
  Hz = matrix1 + jtau*BondDependentQuadQuad(pz) + jb*BondDependentQuadOcto(pz);
  Hx = matrix1 + jtau*BondDependentQuadQuad(px) + jb*BondDependentQuadOcto(px);
  Hy = matrix1 + jtau*BondDependentQuadQuad(py) + jb*BondDependentQuadOcto(py);

  Matrix3LD mat_1, mat_2;
  mat_1 << A_Dir, B_Dir, C_Dir;
  mat_2 <<   0, 0, -1,
            -1, 0,  0,
             0, 1,  0;
  SpinBasis = mat_1 * mat_2;

  std::stringstream paramout;
  paramout << std::fixed << std::setprecision(14);
  paramout << "------------------------Hamiltonian Parameters------------------------\n";
  paramout << "bond-dependent: JTau, JB\n";
  paramout << jtau << " " << jb << "\n";
  paramout << "bond-independent: JQuad, JOcto, lambda\n";
  paramout << jquad << " " << jocto << "\n";
  paramout << "hMagnitude\n";
  paramout << hMagnitude << "\n";
  paramout << "hDirection\n";
  paramout << hDirection.transpose() << "\n";
  ParameterOutput = paramout.str();
}

Matrix3LD Hamiltonia::BondDependentQuadQuad(const long double& angle)
{
  Matrix3LD matrix1;
  long double c = cos(angle);
  long double s = sin(angle);
  matrix1  << 1-c, 0,  -s,
                0, 0,   0,
               -s, 0, 1+c;
  return matrix1/2.0;
}

Matrix3LD Hamiltonia::BondDependentQuadOcto(const long double& angle)
{
  Matrix3LD matrix2;
  long double c = cos(angle);
  long double s = sin(angle);
  matrix2  <<  0, s, 0,
               s, 0, c,
               0, c, 0;
  return matrix2;
}
