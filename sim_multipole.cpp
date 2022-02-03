#include "MC.hpp"

int main(int argc, char *argv[]){

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
  const uint lattice_info = strtol(argv[1], NULL, 10); //0 for tri, 1 for hc, 2 for fcc

  const uint shape = strtol(argv[2], NULL, 10);
  //tri/hc: 0 for primitive, 1 for conventional (rect), 2 for conventional (sqrt)
  //   fcc: 0 for primitive, 1 for conventional

  //linear dimensions of cluster in T1, T2, T3 directions
  const uint    l1 = strtol(argv[3], NULL, 10);
  const uint    l2 = strtol(argv[4], NULL, 10);
  const uint    l3 = strtol(argv[5], NULL, 10);

  const long double a            = strtod(argv[6], NULL);
  const long double b            = strtod(argv[7], NULL);

  Eigen::Array<long double, 5, 1> params;
  long double jtau_unitless, jb_unitless, jquad_unitless, jocto_unitless, h_unitless;
  uint param_entry = 0;
  if (param_entry == 0){
    //scale = jtau2 + jq2 + jo2(=1), jb and h in units of scale
    long double theta = atan(1.0/2.0)/pi; //CAUTION HERE!
    long double phi   = b;

    jtau_unitless    = cos(theta*pi);
    jb_unitless      = strtod(argv[8], NULL);
    jquad_unitless   = sin(theta*pi)*cos(phi*pi);
    jocto_unitless   = sin(theta*pi)*sin(phi*pi);
    h_unitless       = strtod(argv[9], NULL);
  }
  else if (param_entry == 1){
    //scale = abs(jtau)(=1), jb jquad, jocto and h in units of scale
    jtau_unitless    = 1.0;
    jb_unitless      = strtod(argv[8], NULL)/abs(jtau_unitless);
    jquad_unitless   = a/abs(jtau_unitless);
    jocto_unitless   = b/abs(jtau_unitless);
    h_unitless       = strtod(argv[9], NULL)/abs(jtau_unitless);
  }
  params << jtau_unitless, jb_unitless, jquad_unitless, jocto_unitless, h_unitless;
  Vector3LD h_direction = Vector3LD(0.0,1.0,0.0).normalized();
  Hamiltonia hams(params, h_direction);
  Lattice lat(lattice_info, shape, l1, l2, l3, hams);

  // need to implement defect input at some point as well as in jtau_sweep
  // const long double defect_quad = strtod(argv[7], NULL);
  // const long double defect_octo = strtod(argv[8], NULL);
  // const long double defect_lengthscale = strtod(argv[9], NULL);
  // lat.CreateDefectPositions();
  // lat.AddDefectHamiltonia(defect_quad,
  //                         defect_octo,
  //                         defect_lengthscale);
  // lat.CheckHamiltonians();

  if (simulation == 0){
    const double cooling_rate = 0.9;
    const uint num_SA_steps = strtol(argv[10], NULL, 10);
    const long double final_T = pow(cooling_rate,num_SA_steps);
    const uint num_overrelax_ratio = 0;

    const bool printstats = false;
    MonteCarlo mc(lat, final_T, num_overrelax_ratio, printstats, mpirank, mpisize);

    std::ostream &which = std::cout;
    which << std::fixed << std::setprecision(14);
    const long double initial_T = 1.0;
    const uint num_MC_sweeps = pow(10,strtol(argv[11], NULL, 10));
    const uint num_D_sweeps = pow(10,strtol(argv[12], NULL, 10));
    // const uint num_D_sweeps = 0;
    mc.PerformSimulatedAnnealing(which, cooling_rate, initial_T, num_SA_steps,
                                                                 num_MC_sweeps,
                                                                 num_D_sweeps);
    MPI_Finalize();
  } else if (simulation == 1){
    // const long double final_T = strtod(argv[8], NULL);
    // const uint num_overrelax_ratio = 0;
    //
    // const bool printstats = false;
    // MonteCarlo mc(lat, final_T, num_overrelax_ratio, printstats, mpirank, mpisize);
    //
    // std::ostream &which = std::cout;
    // which << std::fixed << std::setprecision(14);
    //
    // const uint num_therm_sweeps = 5*1e2;
    // const uint sampling_time = 1e2;
    // const uint num_sweeps_measurement = 1*(1e2)*sampling_time;
    // mc.PerformFiniteT(which, num_therm_sweeps, num_sweeps_measurement, sampling_time);
  }

  return 0;
}
