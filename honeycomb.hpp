///  @file     honeycomb.hpp
///  @author   Ahmed Rayyan
///  @date     December 2, 2019
///  @brief    defining the honeycomb Lattice class
#ifndef HONEYCOMB_HPP
#define HONEYCOMB_HPP

#include "common.hpp"
#include "spin.hpp"
#include "hamiltonian.hpp"

struct Site
{
  vector<std::tuple<int, int, int, int>> NearestNeighbours;
  Spin OnsiteSpin;
};

class HoneycombLattice
{
public:
  const uint NumSites;
  std::vector<std::vector<std::vector<Site> > > Cluster;
  //std::vector<Site> ClusterList;
  double ClusterEnergy;
  long ActualDetFlips;

  HoneycombLattice(const uint& hc_or_kek, const uint& type, const uint& number_of_sublattices,
          const uint& l1, const uint& l2);

  void InitializeRandomSpins();
  void InitializeFMSpins(const double& theta, const double& phi);
  void InitializeFromFile(const std::vector<std::string>& file_lines);
  void CalculateClusterEnergy(const Parameters& p);
  void PrintConfiguration(std::ostream &out);
  void SimulatedAnnealing(const Parameters& p, const uint& max_sweeps,
                          double& initial_T, double& final_T);
  void DeterministicSweeps(const Parameters& p, const uint& max_sweeps);

private:
  const uint L1, L2, NumSublattices, NumUnitCells, HcOrKekule, ClusterType;
  double FinalT;
  std::mt19937 RNG;
  std::uniform_real_distribution<double> unit_interval;
  std::uniform_int_distribution<uint> L1Dist, L2Dist, SubDist;


  void CreateRhombicCluster1();
  void CreateRhombicCluster2();
  void CreateRectangularCluster1();
  void CreateRectangularCluster2();
  void CreateC3Cluster();
  void CreateKekuleCluster();
  void CalculateLocalEnergy(const Site& site, const Parameters& p, double& energy);
  void MolecularField(const Site& site, const Parameters& p, Eigen::Vector3d& molec);
  void MetropolisSweep(const Parameters& p, const double& temperature);

};

#endif // LATTICE_HPP
