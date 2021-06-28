///  @file     triangular.hpp
///  @author   Ahmed Rayyan
///  @date     June 7, 2021
///  @brief    defining the triangular Lattice class
#ifndef TRIANGULAR_HPP
#define TRIANGULAR_HPP

#include "common.hpp"
#include "spin.hpp"
#include "hamiltonian.hpp"

struct Site
{
  //nn_1, nn_2, bond-dep Hamiltonian
  vector<std::tuple<int, int, Eigen::Matrix3d>> NearestNeighbours;
  Spin OnsiteSpin;
};

class TriangularLattice
{
public:
  const uint L1, L2, NumSites, NumDefects;
  const double JTau, Lambda, IsingY, Defect, HField;
  vector<vector<Site> > Cluster;
  vector<vector<uint> > Defects;
  Eigen::Matrix3d Hx, Hy, Hz, Hdefect, Ham;
  Eigen::Vector3d HDirection, ClusterFMOP,ClusterCombinedOP;
  Eigen::Vector2d ClusterStripyOP;
  double FinalT;
  std::mt19937 RNG;
  std::uniform_real_distribution<double> unit_interval;
  std::uniform_int_distribution<uint> L1Dist, L2Dist;
  long double ClusterEnergy;
  long double EBar, E2Bar, E3Bar, E4Bar;
  long ActualDetFlips;

  bool PoisonedSite_Flag;

  Eigen::ArrayXXd StripySignsX, StripySignsY, StripySignsZ;

  long double FMNorm, PerpNorm, ParNorm, CombinedNorm;
  long double FMNorm2, PerpNorm2, ParNorm2, CombinedNorm2;
  long double FMNorm4, PerpNorm4, ParNorm4, CombinedNorm4;


  TriangularLattice(const uint& l1, const uint& l2, const uint& num_defects,
                    const double& jtau, const double& lambda,
                    const double& ising_y, const double& defect, const double& h,
                    Eigen::Vector3d& hdir);

  void CreateClusterPBC();
  Eigen::Matrix3d ReturnMPHamiltonian(const double& angle);
  void FixMPHamiltonians();
  void InitializeFMSpins(const double& theta, const double& phi);
  void InitializeRandomSpins();
  bool CheckIfPoisoned(uint lx, uint ly);
  void CalculateLocalEnergy(const Site& site, long double& energy);
  void MolecularField(const Site& site, Eigen::Vector3d& molec);
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
  void OverrelaxationFlip();

  //private:
  //
};

#endif
