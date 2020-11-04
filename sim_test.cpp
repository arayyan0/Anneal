///  @file     sim.cpp
///  @author   Ahmed Rayyan
///  @date     December 2, 2019
///  @brief    main simulated annealing file
#include "common.hpp"
#include "lattice.hpp"
#include "hamiltonian.hpp"

int main(int argc, char *argv[])
{

  uint type, sublattice, l1, l2;
  std::string filename = "huh1.out";

  std::ifstream file(filename);
  std::vector<std::string> lines;

  if (file.is_open()) {
      std::string line;
      while (std::getline(file, line)) {
          lines.push_back(line);
      }
    file.close();
  }

  type =  std::stol(lines[2], NULL, 10);

  std::stringstream lattice_ss(lines[4]);
  std::vector<int> lattice_int;
  uint number;
  while (lattice_ss >> number)
    lattice_int.push_back(number);
  l1 = lattice_int[0];
  l2 = lattice_int[1];
  sublattice = lattice_int[2];

  Lattice honeycomb(type, l1, l2, sublattice);

  //##################parameters

  const double phi = strtod(argv[1], NULL);
  const double g = strtod(argv[2], NULL);
  const double a = strtod(argv[3], NULL);

  const double G0 = sin(phi*pi);
  const double K0 = cos(phi*pi);

  Parameters p(
    0, 0, 0,        //Kitaev
    0, 0, 0,        //Gamma
    0,              //Gamma'
    0,              //Heisenberg
    0,              //Field strength
    0,              //Field theta in ABC coordinate basis
    0               //Field phi in ABC coordinate basis
  );

  p.Anisotropy(K0, g, a, 0, -1);
  p.Anisotropy(G0, g, a, 1, +1);

  honeycomb.InitializeFromFile(lines);

  double initial_T = pow(0.9,strtol(argv[4], NULL, 10));
  double final_T = pow(0.9,strtol(argv[5], NULL, 10));
  const uint max_metro_sweeps = pow(10,strtol(argv[6], NULL, 10));
  const uint max_det_sweeps = pow(10,strtol(argv[7], NULL, 10));

  honeycomb.SimulatedAnnealing(p, max_metro_sweeps, initial_T, final_T);
  honeycomb.DeterministicSweeps(p, max_det_sweeps);
  honeycomb.CalculateClusterEnergy(p);

  // cout << std::fixed << std::setprecision(14);
  // PrintSimulationData(cout, type, l1, l2, sublattice, initial_T, final_T, honeycomb.NumSites*max_metro_sweeps, honeycomb.ActualDetFlips);
  // p.PrintParameters(cout);
  // honeycomb.PrintConfiguration(cout);
  return 0;
}
