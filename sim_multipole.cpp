///  @file     sim.cpp
///  @author   Ahmed Rayyan
///  @date     December 2, 2019
///  @brief    main simulated annealing file for multipolar simulation
#include "common.hpp"
#include "triangular.hpp"
#include "hamiltonian.hpp"

int main(int argc, char *argv[])
{
  const uint l1 = strtol(argv[1], NULL, 10);
  const uint l2 = strtol(argv[2], NULL, 10);

  double jtau = 1;
  double lambda = 3.2/2.0;
  double initial_T = 1*abs(jtau);
  double final_T = pow(0.9,strtol(argv[3], NULL, 10));
  const uint max_metro_sweeps = pow(10,strtol(argv[4], NULL, 10));
  const uint max_det_sweeps = pow(10,strtol(argv[5], NULL, 10));

  TriangularLattice triangular(l1, l2, jtau, lambda);
  triangular.CreateClusterPBC();

  triangular.InitializeRandomSpins();
  triangular.SimulatedAnnealing(max_metro_sweeps, initial_T, final_T);
  triangular.DeterministicSweeps(max_det_sweeps);
  triangular.CalculateClusterEnergy();


  cout << std::fixed << std::setprecision(14);
  triangular.PrintConfiguration(cout);


  // prints nearest neighbours of each site
  // for (int y=0; y<l2; ++y)
  // {
  //   for (int x=0; x<l1; ++x)
  //   {
  //     cout << "(x,y)= (" << x << "," << y << ")" << endl;
  //     for (auto i : triangular.Cluster[x][y].NearestNeighbours){
  //     cout << "(nn_x,nn_y)= (" << std::get<0>(i) << "," << std::get<1>(i) << "); type: " << std::get<2>(i)<< endl;
  //     }
  //     cout << "........." << endl;
  //   }
  // }
  return 0;
}
