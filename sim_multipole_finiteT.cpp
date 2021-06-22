///  @file     sim.cpp
///  @author   Ahmed Rayyan
///  @date     December 2, 2019
///  @brief    main simulated annealing file for multipolar simulation
#include "common.hpp"
#include "triangular.hpp"
#include "hamiltonian.hpp"

int main(int argc, char *argv[])
{

  const uint type = 2; //v1 or v2
  const uint sublattice = 1; //1 is rhombic
  const uint l1 = strtol(argv[1], NULL, 10);
  const uint l2 = strtol(argv[2], NULL, 10);

  double jtau = 1;
  const double lambda = 0;
  const double ising_y = strtod(argv[3], NULL);
  const double defect = strtod(argv[4], NULL);
  const double h_field = 0.000;
  Eigen::Vector3d hdir = Eigen::Vector3d(0,0,1).normalized();
  // cout << l1<<l2<<ising_y<<defect<<endl;
  TriangularLattice triangular(l1, l2, jtau, lambda, ising_y, defect, h_field, hdir);

  const uint num_sweeps_SA = 5*1e4;
  const uint num_SA_steps = 30;
  double cooling_rate = 0.9;
  double final_T = strtod(argv[5], NULL)*abs(jtau);
  double initial_T = final_T/pow(cooling_rate,num_SA_steps-1);
  triangular.SimulatedAnnealing(num_sweeps_SA, initial_T, final_T,cooling_rate);

  const uint num_sweeps_thermal = 1e4;
  triangular.ThermalizeConfiguration(final_T, num_sweeps_thermal);

  triangular.CalculateClusterEnergyandOP();
  const uint num_sweeps_measurement = 1e6;
  const uint sampling_time = 1e1;

  triangular.SampleConfiguration(final_T, num_sweeps_measurement, sampling_time);

  triangular.CalculateClusterEnergyandOP();
  double actual_det_flips=0;
  cout << std::fixed << std::setprecision(14);
  PrintTriangularSimulationData(cout, type, sublattice, l1, l2,
                                           initial_T, final_T, num_sweeps_SA,
                                           num_sweeps_thermal,num_sweeps_measurement,
                                           actual_det_flips);
  cout << "------------------------Hamiltonian Parameters------------------------\n";
  cout << "J_Tau\n";
  cout << jtau << "\n";
  cout << "Lambda\n";
  cout << lambda << "\n";
  cout << "IsingY\n";
  cout << ising_y << "\n";
  cout << "defect\n";
  cout << defect << "\n";
  cout << "HField\n";
  cout << h_field << "\n";
  cout << "HDirection\n";
  cout << hdir.transpose() << "\n";
  triangular.PrintConfiguration(cout);

  triangular.PrintThermalObservables(cout);


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
