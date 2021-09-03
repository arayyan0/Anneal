#include "MC.hpp"

int main(int argc, char *argv[]){
  //initiate MPI
  MPI_Init(&argc, &argv);
  int mpisize, mpirank;
  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);

  const uint l1 = strtol(argv[1], NULL, 10);
  const uint l2 = strtol(argv[2], NULL, 10);          //should be equal to l1
  const uint num_sublattices = 1;
  const uint num_defects = strtol(argv[3], NULL, 10); //should only be 1,3,9

  const long double jtau = 1; //J_tau = +1 sets energy unit
  const long double lambda = 0;
  const long double jquad = strtod(argv[4], NULL);
  const long double jocto = strtod(argv[5], NULL);
  const long double h_magnitude = 0.000;
  Vector3LD h_direction = Vector3LD(0,1,0).normalized();
  Hamiltonia hams(jtau, lambda, jquad, jocto, h_magnitude, h_direction);

  const long double defect_quad = strtod(argv[6], NULL);
  const long double defect_octo = strtod(argv[7], NULL);
  const long double defect_lengthscale = strtod(argv[8], NULL);

  Triangular tri(l1, l2, num_sublattices, num_defects, hams,
    defect_quad, defect_octo, defect_lengthscale);

  uint simulation = 1; //0 for simulated annealing, 1 for finite T

  if (simulation == 0){
    const double cooling_rate = 0.9;
    const uint num_SA_steps = strtol(argv[9], NULL, 10);
    const long double final_T = pow(cooling_rate,num_SA_steps);
    const uint num_overrelax_ratio = 0;

    const bool printstats = false;
    MonteCarlo mc(tri, final_T, num_overrelax_ratio, printstats, mpirank, mpisize);

    std::ostream &which = std::cout;
    which << std::fixed << std::setprecision(14);
    const long double initial_T = 1.0;
    const uint num_MC_sweeps = pow(10,strtol(argv[10], NULL, 10));
    const uint num_D_sweeps = pow(10,strtol(argv[11], NULL, 10));
    // const uint num_D_sweeps = 0;
    mc.PerformSimulatedAnnealing(which, cooling_rate, initial_T, num_SA_steps,
                                                                 num_MC_sweeps,
                                                                 num_D_sweeps);
  } else if (simulation == 1){
    const long double final_T = strtod(argv[9], NULL);
    const uint num_overrelax_ratio = 0;

    const bool printstats = false;
    MonteCarlo mc(tri, final_T, num_overrelax_ratio, printstats, mpirank, mpisize);

    std::ostream &which = std::cout;
    which << std::fixed << std::setprecision(14);

    const uint num_therm_sweeps = 5*1e4;
    const uint sampling_time = 1e2;
    const uint num_sweeps_measurement = 1*(1e4)*sampling_time;
    mc.PerformFiniteT(which, num_therm_sweeps, num_sweeps_measurement, sampling_time);
  }

  MPI_Finalize();

  return 0;
}
