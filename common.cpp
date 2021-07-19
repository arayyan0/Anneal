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

void PrintSimulationData(std::ostream &out, const uint& hc_or_kek, const uint& type,
                         const uint& sublattice, const uint& l1, const uint& l2,
                         const double& T_i, double& T_f, const uint& max_mflips,
                         const uint& max_daligns)
{
  out << "-------------------------Simulation Parameters-------------------------\n";
  out << "Honeycomb/Kekule? Which cluster type?\n";
  out << hc_or_kek << " " << type << "\n";
  out << "Number of Sublattices/Unit Cells (s, l1, l2)\n";
  out << sublattice << " " << l1 << " " << l2 << " " << "\n";
  out << "Initial Temperature\n";
  out << T_i << "\n";
  out << "Final Temperature\n";
  out << T_f << "\n";
  out << "Number of Metropolis Flips\n";
  out << max_mflips << "\n";
  out << "Number of Deterministic Aligns\n";
  out << max_daligns << "\n";
}

void PrintTriangularSimulationData(std::ostream &out, const uint& type,
                         const uint& sublattice, const uint& l1, const uint& l2,
                         const double& T_i, double& T_f, const uint& max_sa_sweeps,
                         const uint& max_thermal_sweeps, const uint& max_measuring_sweeps,
                         const uint& sampling_time, const uint& max_dsweeps)
{
  out << "-------------------------Simulation Parameters-------------------------\n";
  out << "Which cluster type?\n";
  out << type << "\n";
  out << "Number of Sublattices/Unit Cells (s, l1, l2)\n";
  out << sublattice << " " << l1 << " " << l2 << " " << "\n";
  out << "Initial Temperature\n";
  out << T_i << "\n";
  out << "Final Temperature\n";
  out << T_f << "\n";
  out << "Number of Metropolis sweeps (SA, thermal, measuring/sampling time)\n";
  out << max_sa_sweeps << " " << max_thermal_sweeps << " " << max_measuring_sweeps << "/" << sampling_time << "\n";
  out << "Number of Deterministic sweeps\n";
  out << max_dsweeps << "\n";
}

void BravaisIndicesToFlat(const uint& x, const uint&y, const uint & period, uint& flat_index)
{
  flat_index = x + y*period;
}

long double Lorentzian(long double s, long double x, long double l){
  long double cutoff=10000;
  long double dist = (x-0.5)/l;
  long double result = dist<cutoff ? s /pow( 1 + pow(dist,2), 2.0/2.0 ) : 0;
  return result;
}

long double Gaussian(long double s, long double x, long double l){
  long double cutoff=10000;
  long double dist = (x-0.5)/l;
  long double result = dist<cutoff ? s*exp(-pow(dist,2)/2.0) : 0;
  return result;
}
