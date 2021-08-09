#include "MC.hpp"

int main(int argc, char *argv[]){

  const uint l1 = strtol(argv[1], NULL, 10);
  const uint l2 = strtol(argv[2], NULL, 10);          //should be equal to l1
  const uint num_sublattices = 1;
  const uint num_defects = strtol(argv[3], NULL, 10); //should only be 1,3,9

  const long double jtau = 1; //J_tau = +1 sets energy unit
  const long double lambda = 0;
  const long double jquad = strtod(argv[4], NULL);
  const long double jocto = strtod(argv[5], NULL);
  const long double defect_strength = strtod(argv[6], NULL);
  const long double defect_lengthscale = strtod(argv[7], NULL);
  const long double h_magnitude = 0.000;
  Vector3LD h_direction = Vector3LD(0,1,0).normalized();

  TriLatt tri(l1, l2, num_sublattices, num_defects, jtau, lambda, jquad, jocto,
    defect_strength, defect_lengthscale, h_magnitude, h_direction);

  const long double final_T = strtod(argv[8], NULL);
  const uint num_overrelax_ratio = 5;

  const bool printstats = false;
  MonteCarlo mc(tri, final_T, num_overrelax_ratio, printstats);

  std::ostream &which = std::cout;
  which << std::fixed << std::setprecision(14);

  const uint num_therm_sweeps = 5*1e4;
  const uint sampling_time = 1e2;
  const uint num_sweeps_measurement = 1*(1e4)*sampling_time;
  mc.PerformFiniteT(which, num_therm_sweeps, num_sweeps_measurement, sampling_time);

  return 0;
}
