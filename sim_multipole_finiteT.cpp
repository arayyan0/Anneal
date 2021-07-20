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
  if (l1!=l2){
    cout << "CAUTION: L1 != L2" << endl;
    cout << "In principle this is OK if you are not including defects." << endl;
    cout << "If so, ensure that defect strength=0." << endl;
  }

  const uint num_defects = strtol(argv[3], NULL, 10); //should only be 1,3,9

  long double jtau = 1;
  const long double lambda = 0;
  const long double ising_y = strtod(argv[4], NULL);
  const long double defect = strtod(argv[5], NULL);
  const long double defect_lengthscale = strtod(argv[6], NULL);
  const long double h_field = 0.000;
  Vector3LD hdir = Vector3LD(0,0,1).normalized();
  // cout << l1<<l2<<ising_y<<defect<<endl;
  TriangularLattice triangular(l1, l2, num_defects,jtau, lambda, ising_y, defect, defect_lengthscale, h_field, hdir);

  double final_T = strtod(argv[7], NULL)*abs(jtau);

  uint num_sweeps_SA, num_sweeps_thermal;
  double initial_T, cooling_rate;
  // double abovecritical_T = 3*abs(jtau);
  // if (final_T > abovecritical_T){
  //   num_sweeps_SA = 0;
  //   initial_T = 0;
  //   //skip simulated annealing and do extra long thermalization instead
  //   num_sweeps_thermal = 1e5;
  // } else {
  //   num_sweeps_SA = 1*1e3;
  //   initial_T = abovecritical_T;
  //   cooling_rate = 0.9;
  //   triangular.SimulatedAnnealing(num_sweeps_SA, initial_T, final_T, cooling_rate);
  //
  //   num_sweeps_thermal = 1e4;
  // }
  num_sweeps_SA = 0;
  initial_T = 0;
  num_sweeps_thermal = 5*1e4;
  triangular.ThermalizeConfiguration(final_T, num_sweeps_thermal);
  //
  const uint sampling_time = 1e2;
  const uint num_sweeps_measurement = 1*(1e4)*sampling_time;
  triangular.SampleConfiguration(final_T, num_sweeps_measurement, sampling_time);
  //
  double actual_det_sweeps=0;
  std::ostream &which = std::cout;
  which << std::fixed << std::setprecision(14);
  PrintTriangularSimulationData(which, type, sublattice, l1, l2,
                                           initial_T, final_T, num_sweeps_SA,
                                           num_sweeps_thermal,num_sweeps_measurement,
                                           sampling_time, actual_det_sweeps);
  which << "------------------------Hamiltonian Parameters------------------------\n";
  which << "J_Tau\n";
  which << jtau << "\n";
  which << "Lambda\n";
  which << lambda << "\n";
  which << "IsingY\n";
  which << ising_y << "\n";
  which << "defect/length scale/number of defects\n";
  which << defect << "/" << defect_lengthscale << "/" <<num_defects<< "\n";
  which << "HField\n";
  which << h_field << "\n";
  which << "HDirection\n";
  which << hdir.transpose() << "\n";
  triangular.PrintConfiguration(which);

  triangular.PrintThermalObservables(which);


  // // prints nearest neighbours of each site
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
