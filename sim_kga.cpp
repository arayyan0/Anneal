///  @file     sim.cpp
///  @author   Ahmed Rayyan
///  @date     December 2, 2019
///  @brief    main simulated annealing file
#include "common.hpp"
#include "lattice.hpp"
#include "hamiltonian.hpp"

int main(int argc, char *argv[])
{
  const double initial_T = 1;

  const uint type = 0;
  const uint sublattice = strtol(argv[1], NULL, 10);
  const uint l1 = strtol(argv[2], NULL, 10);
  const uint l2 = strtol(argv[3], NULL, 10);
  double final_T = pow(0.9,strtol(argv[4], NULL, 10));
  const uint max_metro_sweeps = pow(10,strtol(argv[5], NULL, 10));
  const uint max_det_sweeps = pow(10,strtol(argv[6], NULL, 10));

  const double phi = strtod(argv[7], NULL);
  const double g = strtod(argv[8], NULL);
  const double a = strtod(argv[9], NULL);

  const double G0 = sin(phi*pi);
  const double K0 = cos(phi*pi);
  
  const double Gp = strtod(argv[10], NULL);

  Parameters p(
    0, 0, 0,        //Kitaev
    0, 0, 0,        //Gamma
    Gp,              //Gamma'
    0,              //Heisenberg
    0,              //Field strength
    0,              //Field theta in ABC coordinate basis
    0               //Field phi in ABC coordinate basis
  );

  p.Anisotropy(K0, g, a, 0, -1);
  p.Anisotropy(G0, g, a, 1, +1);

  Lattice honeycomb(type, l1, l2, sublattice);
  honeycomb.InitializeRandomSpins();
  honeycomb.SimulatedAnnealing(p, max_metro_sweeps, initial_T, final_T);
  honeycomb.DeterministicSweeps(p, max_det_sweeps);
  honeycomb.CalculateClusterEnergy(p);

  cout << std::fixed << std::setprecision(14);
  PrintSimulationData(cout, type, l1, l2, sublattice, initial_T, final_T, honeycomb.NumSites*max_metro_sweeps, honeycomb.ActualDetFlips);
  p.PrintParameters(cout);
  honeycomb.PrintConfiguration(cout);
  return 0;
}
