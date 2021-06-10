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
  const uint L1, L2, NumSites;
  const double JTau, Lambda;
  vector<vector<Site> > Cluster;
  Eigen::Matrix3d Hx, Hy, Hz;
  double FinalT;
  std::mt19937 RNG;
  std::uniform_real_distribution<double> unit_interval;
  std::uniform_int_distribution<uint> L1Dist, L2Dist;
  double ClusterEnergy;
  long ActualDetFlips;


  TriangularLattice(const uint& l1, const uint& l2, const double& jtau, const double& lambda);

  void CreateClusterPBC();
  Eigen::Matrix3d ReturnMPHamiltonian(const double& angle, const double& lambda);
  void FixMPHamiltonians(const double& jtau, const double& lambda);
  void InitializeFMSpins(const double& theta, const double& phi);
  void InitializeRandomSpins();
  void CalculateLocalEnergy(const Site& site, double& energy);
  void MolecularField(const Site& site, Eigen::Vector3d& molec);
  void CalculateClusterEnergy();
  void MetropolisSweep(const double& temperature);
  void SimulatedAnnealing(const uint& max_sweeps, double& initial_T, double& final_T);
  void DeterministicSweeps(const uint& max_sweeps);
  void PrintConfiguration(std::ostream &out);

  //private:
  //
};

#endif
