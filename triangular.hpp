///  @file     triangular.hpp
///  @author   Ahmed Rayyan
///  @date     June 7, 2021
///  @brief    defining the triangular Lattice class
#ifndef TRIANGULAR_HPP
#define TRIANGULAR_HPP

#include "common.hpp"
#include "spin.hpp"
// #include "hamiltonian.hpp"

struct SiteInfo
{
  //nn_1, nn_2, bond-dep Hamiltonian
  vector<std::tuple<int, int, Matrix3LD>> NearestNeighbours;
};

class TriangularLattice
{
public:
  const uint L1, L2, NumSites, NumDefects;
  const long double JTau, Lambda, IsingY, Defect, HField;
  vector<vector<SiteInfo> > ClusterInfo;
  vector<vector<uint> > Defects;
  Matrix3LD Hx, Hy, Hz, Hdefect1, Hdefect2;
  Vector3LD HDirection;
  std::mt19937 RNG;
  std::uniform_real_distribution<double> unit_interval;
  std::uniform_int_distribution<uint> L1Dist, L2Dist;
  ArrayXXLD StripySignsX, StripySignsY, StripySignsZ;

  uint overrelaxMCratio;

  vector<vector<Vector3LD> > Cluster;                      //replica-dependent

  long double ClusterEnergy;                               //replica-dependent
  long double EBar, E2Bar, E3Bar, E4Bar;                   //replica-dependent

  Eigen::Matrix<long double, 3, 2> ClusterStripyOPMatrix;  //replica-dependent
  Vector3LD ClusterFMOP,ClusterCombinedOP;                 //replica-dependent
  Vector2LD ClusterStripyOP;                               //replica-dependent
  long double FMNorm, PerpNorm, ParNorm, CombinedNorm;     //replica-dependent
  long double FMNorm2, PerpNorm2, ParNorm2, CombinedNorm2; //replica-dependent
  long double FMNorm4, PerpNorm4, ParNorm4, CombinedNorm4; //replica-dependent


  TriangularLattice(const uint& l1, const uint& l2, const uint& num_defects,
                    const long double& jtau, const long double& lambda,
                    const long double& ising_y, const long double& defect,
                    const long double& h, Vector3LD& hdir);
  void CreateClusterPBC();
  Matrix3LD ReturnMPHamiltonian(const long double& angle);
  void FixMPHamiltonians();
  void AddDefectHamiltonia();
  void InitializeFMSpins(const long double& theta, const long double& phi);
  void InitializeRandomSpins();
  bool CheckIfPoisoned(uint lx, uint ly);
  void CalculateLocalEnergy(const uint& n1, const uint& n2, long double& energy);
  void MolecularField(const uint& n1, const uint& n2, Vector3LD& molec);
  void MetropolisSweep(const double& temperature, uint& accept);
  void DeterministicSweeps(const uint& max_sweeps);
  void PrintConfiguration(std::ostream &out);
  void SimulatedAnnealing(const uint& max_sweeps,
                          double& initial_T, double& final_T, double& rate);
  void ThermalizeConfiguration(double& temp, const uint& max_flips);
  void SampleConfiguration(double& temp, const uint& max_sweeps,
                           const uint& sampling_time);
  void CreateStripySignMatrices();
  void SelectStripyOP();
  void CalculateClusterEnergyandOP();
  void PrintThermalObservables(std::ostream &out);
  void CreateDefectPositions();
  void OverrelaxationSweep();
  void DoTheSweeps(double& temp, uint& accept);


  //private:
  //
};

#endif
