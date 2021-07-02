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
  void MetropolisSweep(const double& temperature);
  void DeterministicSweeps(const uint& max_sweeps);
  void PrintConfiguration(std::ostream &out);

  void SimulatedAnnealing(const uint& max_sweeps,
                                            double& initial_T, double& final_T, double& rate);

  void ThermalizeConfiguration(double& temp, const uint& max_flips);
  void SampleConfiguration(double& temp, const uint& max_sweeps, const uint& sampling_time);

  void MetropolisFlip(
    uint& uc_x, uint& uc_y,
    Spin&  old_spin_at_chosen_site,
    long double &old_local_energy, long double &new_local_energy, long double &energy_diff, double &r,
    double &pd,
    const double& temperature
  );
  void CreateStripySignMatrices();

  void CalculateClusterEnergyandOP();
  void CalculateClusterEnergy();
  void CalculateClusterOP();
  void PrintThermalObservables(std::ostream &out);

  void CreateDefectPositions();
  void OverrelaxationFlip(uint& uc_x, uint& uc_y);
  void OverrelaxationSweep();

  //private:
  //
};

#endif
