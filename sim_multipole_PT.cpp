///  @file     sim.cpp
///  @author   Ahmed Rayyan
///  @date     December 2, 2019
///  @brief    main simulated annealing file for multipolar simulation
#include "common.hpp"
#include "triangular.hpp"
#include "hamiltonian.hpp"
#include <mpi.h>

int main(int argc, char *argv[])
{
  //lattice settings
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

  //hamiltonian settings
  long double jtau = 1;
  const long double lambda = 0.0;
  const long double ising_y = strtod(argv[4], NULL);
  const long double defect = strtod(argv[5], NULL);
  const long double h_field = 0.000;
  Vector3LD hdir = Vector3LD(0,1,0).normalized();

  //simulated annealing settings
  const uint num_SA_steps = strtol(argv[6], NULL, 10);
  const uint num_sweeps_SA = pow(10, strtol(argv[7], NULL, 10));
  double cooling_rate = 0.9;
  double initial_T = 1*abs(jtau);
  double final_T = pow(cooling_rate,num_SA_steps);

  //deterministic sweep settings
  const uint max_det_sweeps = pow(10,strtol(argv[8], NULL, 10));

  //no thermalization or measurements in pure simulated annealing
  const uint num_sweeps_thermal = 0;
  const uint num_sweeps_measurement = 0;
  const uint sampling_time = 0;

  //initiale MPI
  MPI_Init(&argc, &argv);
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  //pointer to array of size size which contains cluster energy
  //out of which we will find the lowest energy configuration
  long double* minarr = (long double *)malloc(sizeof(long double)*size);

  //initiate simulation on each process
  TriangularLattice triangular(l1, l2, num_defects, jtau, lambda, ising_y, defect, h_field, hdir);
  triangular.SimulatedAnnealing(num_sweeps_SA, initial_T, final_T,cooling_rate);
  triangular.DeterministicSweeps(max_det_sweeps);
  triangular.CalculateClusterEnergyandOP();

  //Allgather cluster energies into array, convert to STL vector, and find rank
  //of the lowest ground state
  MPI_Allgather(&triangular.ClusterEnergy,1,MPI_LONG_DOUBLE,
                minarr, 1, MPI_LONG_DOUBLE,MPI_COMM_WORLD);
  std::vector<long double> v(minarr,minarr+size);
  int firstlowestrank = std::min_element(v.begin(),v.end()) - v.begin();

  //output the information of the cluster with lowest energy
  if (rank == firstlowestrank){
    std::ostream &which = std::cout;
    which << std::fixed << std::setprecision(14);
    PrintTriangularSimulationData(which, type, sublattice, l1, l2,
                                             initial_T, final_T, num_sweeps_SA,
                                             num_sweeps_thermal, num_sweeps_measurement,
                                             sampling_time, max_det_sweeps);
    which << "------------------------Hamiltonian Parameters------------------------\n";
    which << "J_Tau\n";
    which << jtau << "\n";
    which << "Lambda\n";
    which << lambda << "\n";
    which << "IsingY\n";
    which << ising_y << "\n";
    which << "defect/number of defects\n";
    which << defect << "/" << num_defects<< "\n";
    which << "HField\n";
    which << h_field << "\n";
    which << "HDirection\n";
    which << hdir.transpose() << "\n";
    triangular.PrintConfiguration(which);
  }
  MPI_Finalize();

  return 0;
}
