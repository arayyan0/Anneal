#include "MC.hpp"

int main(int argc, char *argv[])
{
  const uint hc_or_kek = 0; //Honeycomb or Kekule
  const uint cluster_type = strtol(argv[1], NULL, 10);
  const uint num_sublattices = strtol(argv[2], NULL, 10);
  const uint l1 = strtol(argv[3], NULL, 10);
  const uint l2 = strtol(argv[4], NULL, 10);

  const long double phi = strtod(argv[5], NULL);
  const long double g   = strtod(argv[6], NULL);
  const long double a   = strtod(argv[7], NULL);
  const long double h_magnitude = 0.000;
  Vector3LD h_direction = Vector3LD(0,1,0).normalized();

  HoneyLatt honey(hc_or_kek, cluster_type, num_sublattices, l1, l2, phi, g, a,
                  h_magnitude, h_direction);
  //
  const double cooling_rate = 0.9;
  const uint num_SA_steps = strtol(argv[8], NULL, 10);
  const long double final_T = pow(cooling_rate,num_SA_steps);
  const uint num_overrelax_ratio = 0;

  const bool printstats = false;
  MonteCarlo mc(honey, final_T, num_overrelax_ratio, printstats);

  std::ostream &which = std::cout;
  which << std::fixed << std::setprecision(14);
  const long double initial_T = 1.0;
  const uint num_MC_sweeps = pow(10,strtol(argv[9], NULL, 10));
  const uint num_D_sweeps = pow(10,strtol(argv[10], NULL, 10));
  // const uint num_D_sweeps = 0;
  mc.PerformSimulatedAnnealing(which, cooling_rate, initial_T, num_SA_steps,
                                                               num_MC_sweeps,
                                                               num_D_sweeps);
  return 0;
}
