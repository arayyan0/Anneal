#include "MC.hpp"

int main(int argc, char *argv[])
{
  uint simulation = 0; //0 for simulated annealing, 1 for finite T

  int mpisize, mpirank;
  if (simulation == 0){
    //initiate MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  } else if (simulation == 1){
    //do not initiate MPI
    mpisize = 1;
    mpirank = 0;
  }

  const uint hc_or_kek = 0; //Honeycomb or Kekule
  const uint cluster_type = strtol(argv[1], NULL, 10);
  const uint num_sublattices = strtol(argv[2], NULL, 10);
  const uint l1 = strtol(argv[3], NULL, 10);
  const uint l2 = strtol(argv[4], NULL, 10);

  const long double phi = strtod(argv[5], NULL);
  const long double g   = strtod(argv[6], NULL);
  const long double a   = strtod(argv[7], NULL);
  const long double h_magnitude = 0.000;
  Vector3LD h_direction = Vector3LD(1,1,1).normalized();
  Hamiltonia hams(phi, g, a, h_magnitude, h_direction);

  Honeycomb honey(hc_or_kek, cluster_type, num_sublattices, l1, l2, hams);

  if (simulation == 0){
    const double cooling_rate = 0.9;
    const uint num_SA_steps = strtol(argv[8], NULL, 10);
    const long double final_T = pow(cooling_rate,num_SA_steps);
    const uint num_overrelax_ratio = 0;

    const bool printstats = false;
    MonteCarlo mc(honey, final_T, num_overrelax_ratio, printstats, mpirank, mpisize);

    std::ostream &which = std::cout;
    which << std::fixed << std::setprecision(14);
    const long double initial_T = 1.0;
    const uint num_MC_sweeps = pow(10,strtol(argv[9], NULL, 10));
    const uint num_D_sweeps = pow(10,strtol(argv[10], NULL, 10));
    // const uint num_D_sweeps = 0;
    mc.PerformSimulatedAnnealing(which, cooling_rate, initial_T, num_SA_steps,
                                                                 num_MC_sweeps,
                                                                 num_D_sweeps);
    MPI_Finalize();
  } else if (simulation == 1){
    const long double final_T = strtod(argv[8], NULL);
    const uint num_overrelax_ratio = 5;

    const bool printstats = false;
    MonteCarlo mc(honey, final_T, num_overrelax_ratio, printstats, mpirank, mpisize);

    std::ostream &which = std::cout;
    which << std::fixed << std::setprecision(14);

    const uint num_therm_sweeps = 5*1e4;
    const uint sampling_time = 1e2;
    const uint num_sweeps_measurement = 1*(1e4)*sampling_time;
    mc.PerformFiniteT(which, num_therm_sweeps, num_sweeps_measurement, sampling_time);
  }

  return 0;
}
