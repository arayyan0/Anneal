///  @file     common.cpp
///  @author   Ahmed Rayyan
///  @date     December 2, 2019
///  @brief    useful functions for the simulated annealing algorithm
#include "common.hpp"

namespace MyRandom
{
  auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  std::mt19937 RNG(seed);
  std::uniform_real_distribution<double> unit_interval(0,1);
}

void PrintSimulationData(std::ostream &out, const uint& type, const uint& l1, const uint& l2,
                         const uint& sublattice, const double& T_i, double& T_f,
                         const uint& max_mflips, const uint& max_daligns)
{
  out << "-------------------------Simulation Parameters-------------------------\n";
  out << "Type of lattice\n";
  out << type << "\n";
  out << "Number of Sublattices/Unit Cells (l1, l2, s)\n";
  out << l1 << " " << l2 << " " << sublattice << "\n";
  out << "Initial Temperature\n";
  out << T_i << "\n";
  out << "Final Temperature\n";
  out << T_f << "\n";
  out << "Number of Metropolis Flips\n";
  out << max_mflips << "\n";
  out << "Number of Deterministic Aligns\n";
  out << max_daligns << "\n";
}
