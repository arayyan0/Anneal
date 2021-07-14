///  @file     triangular.hpp
///  @author   Ahmed Rayyan
///  @date     June 7, 2021
///  @brief    defining the triangular Lattice class
#ifndef TRIANGULAR_HPP
#define TRIANGULAR_HPP

#include "common.hpp"
#include "spin.hpp"
// #include "hamiltonian.hpp"

struct Site
{
  //nn_1, nn_2, bond-dep Hamiltonian
  vector<std::tuple<int, int, Matrix3LD>> NearestNeighbours;
  Spin OnsiteSpin;
};

struct Sitee
{
  //nn_1, nn_2, bond-dep Hamiltonian
  vector<std::tuple<int, int, Matrix3LD>> NearestNeighbours;
  Vector3LD VectorXYZ;

};

class TriangularLattice
{
public:
  const uint L1, L2, NumSites, NumDefects;
  const long double JTau, Lambda, IsingY, Defect, HField;
  vector<vector<Site> > Cluster;
  vector<vector<uint> > Defects;
  Matrix3LD Hx, Hy, Hz, Hdefect1, Hdefect2;
  Vector3LD HDirection, ClusterFMOP,ClusterCombinedOP;
  Vector2LD ClusterStripyOP;
  Eigen::Matrix<long double, 3, 2> ClusterStripyOPMatrix;

  double FinalT;
  std::mt19937 RNG;
  std::uniform_real_distribution<double> unit_interval;
  std::uniform_int_distribution<uint> L1Dist, L2Dist;
  long double ClusterEnergy;
  long double EBar, E2Bar, E3Bar, E4Bar;
  long ActualDetSweeps;

  ArrayXXLD StripySignsX, StripySignsY, StripySignsZ;

  long double FMNorm, PerpNorm, ParNorm, CombinedNorm;
  long double FMNorm2, PerpNorm2, ParNorm2, CombinedNorm2;
  long double FMNorm4, PerpNorm4, ParNorm4, CombinedNorm4;

  uint overrelaxMCratio;

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
  void CalculateLocalEnergy(const Site& site, long double& energy);
  void MolecularField(const Site& site, Vector3LD& molec);
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
  void CalculateClusterEnergy();
  void PrintThermalObservables(std::ostream &out);
  void CreateDefectPositions();
  void OverrelaxationSweep();
  void DoTheSweeps(double& temp, uint& accept);


  //private:
  //
};

#endif
