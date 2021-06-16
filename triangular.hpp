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
  //nn_1, nn_2, bond_type
  vector<std::tuple<int, int, int>> NearestNeighbours;
  Spin OnsiteSpin;
};

class TriangularLattice
{
public:
  const uint L1, L2, NumSites, PoisonedX, PoisonedY;
  const double JTau, Lambda, IsingY, Defect, HField;
  vector<vector<Site> > Cluster;
  Eigen::Matrix3d Hx, Hy, Hz, Hdefect, Ham;
  Eigen::Vector3d HDirection;
  double FinalT;
  std::mt19937 RNG;
  std::uniform_real_distribution<double> unit_interval;
  std::uniform_int_distribution<uint> L1Dist, L2Dist;
  double ClusterEnergy, SpecificHeat;
  long ActualDetFlips;

  bool PoisonedSite_Flag;

  TriangularLattice(const uint& l1, const uint& l2, const double& jtau, const double& lambda,
                    const double& ising_y, const double& defect, const double& h,
                    Eigen::Vector3d& hdir);

  void CreateClusterPBC();
  Eigen::Matrix3d ReturnMPHamiltonian(const double& angle);
  void FixMPHamiltonians();
  void InitializeFMSpins(const double& theta, const double& phi);
  void InitializeRandomSpins();
  bool CheckIfPoisoned(uint lx, uint ly);
  void CalculateLocalEnergy(const Site& site, double& energy);
  void MolecularField(const Site& site, Eigen::Vector3d& molec);
  void CalculateClusterEnergy();
  void MetropolisSweep(const double& temperature);
  void SimulatedAnnealing1(const uint& max_sweeps, double& initial_T, double& final_T);
  void DeterministicSweeps(const uint& max_sweeps);
  void PrintConfiguration(std::ostream &out);

  void SimulatedAnnealing2(const uint& max_sweeps,
                                            double& initial_T, double& final_T, double& rate);

  void ThermalizeConfiguration(double& temp, const uint& max_flips);
  void SampleConfiguration(double& temp, const uint& max_sweeps);

  void MetropolisFlip(
    uint& uc_x, uint& uc_y,
    Spin&  old_spin_at_chosen_site,
    double &old_local_energy, double &new_local_energy, double &energy_diff, double &r,
    double &pd,
    const double& temperature
  );

  //private:
  //
};

#endif
