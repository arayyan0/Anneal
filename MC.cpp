///  @file     triangular.hpp
///  @author   Ahmed Rayyan
///  @date     July 23, 2021
///  @brief    defining a separate Monte Carlo class
#include "MC.hpp"

MonteCarloStatistics::MonteCarloStatistics(const uint& num_temps, const uint& num_sweeps):
NumTemps(num_temps), NumSweeps(num_sweeps)
{
  EnergyDensity.resize(NumTemps*NumSweeps);
  AcceptanceRate.resize(NumTemps*NumSweeps);
}

void MonteCarloStatistics::WriteStatisticsFile()
{
  std::ofstream outfile;
  outfile.open("stats.out");
  outfile << std::fixed << std::setprecision(14);
  outfile << NumTemps << " " << NumSweeps <<  " " << 1 << endl;
  for (uint i=0; i<NumTemps*NumSweeps;++i){
    outfile << i << " " << EnergyDensity[i] << " " << AcceptanceRate[i] << endl;
  }
  outfile.close();
}

// MonteCarlo::MonteCarlo(Triangular& lattice, const double& final_T, const uint& num_overrelax_ratio,
MonteCarlo::MonteCarlo(Honeycomb& lattice, const double& final_T, const uint& num_overrelax_ratio,
                       const bool& recordstats, const int& mpirank, const int& mpisize):
Lattice(lattice), FinalT(final_T),OverrelaxMCRatio(num_overrelax_ratio), RecordStats(recordstats),
MPIRank(mpirank), MPISize(mpisize)
{
  auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> unit_interval(0,1);
  std::uniform_int_distribution<uint> sd(0, Lattice.NumSites-1);

  RNG = rng;
  UnitInterval = unit_interval;
  SiteDist = sd;
}

void MonteCarlo::InitializeFMSpins(const long double& theta, const long double& phi)
{
  Vector3LD v = SphericalAnglesToCubic(theta, phi);
  for (uint flat_index=0; flat_index<Lattice.NumSites; ++flat_index){
    Lattice.Cluster[flat_index]= v;
  }
}

void MonteCarlo::InitializeRandomSpins()
{
  for (uint flat_index=0; flat_index<Lattice.NumSites; ++flat_index){
    SpherePointPicker(Lattice.Cluster[flat_index]);
  }
}

void MonteCarlo::MolecularField(const uint& flat_index, Vector3LD& molec)
{
  Vector3LD v = Vector3LD::Zero();
  uint flat_index_nn;
  for (auto &j : Lattice.ClusterInfo[flat_index].NearestNeighbours){
    BravaisIndicesToFlat(get<2>(j), get<0>(j), Lattice.NumSublattices, get<1>(j), Lattice.L1*Lattice.NumSublattices, flat_index_nn);
    v += -Lattice.Cluster[flat_index_nn].transpose()*get<3>(j);
  }
  molec = v.transpose()+ Lattice.hField.transpose();
}

void MonteCarlo::CalculateLocalEnergy(const uint& flat_index, long double& energy)
{
  Vector3LD molec;
  MolecularField(flat_index, molec);
  energy = -Lattice.Cluster[flat_index].dot(molec+ Lattice.hField);
  // cout << energy<< endl;
}

void MonteCarlo::CalculateClusterEnergy()
{
  long double e=0;
  long double local_energy;
  uint flat_index;

  for (uint flat_index=0; flat_index<Lattice.NumSites; ++flat_index){
    // local_energy=0;
    //calculate energy
    CalculateLocalEnergy(flat_index, local_energy);
    e += local_energy;
  }
  Lattice.ClusterEnergy = e/2.0;
}

void MonteCarlo::PerformSimulatedAnnealing(std::ostream &out, const double& cooling_rate,
                                                              const double& initial_T,
                                                              const uint& num_SA_steps,
                                                              const uint& num_MC_sweeps,
                                                              const uint& num_D_sweeps)
{
  const uint max_thermal_sweeps = 0;
  const uint max_measuring_sweeps = 0;
  const uint sampling_time = 0;

  if (MPIRank == 0){
    Lattice.PrintLatticeParameters(out);
    PrintSimulationData(out, initial_T, FinalT, num_MC_sweeps,
                 max_thermal_sweeps, max_measuring_sweeps,
                 sampling_time, num_D_sweeps);
    Lattice.PrintHamiltonianParameters(out);
  }

  long double* minarr = (long double *)malloc(sizeof(long double)*MPISize);

  InitializeFMSpins(pi/2.0,pi/2.0);
  // InitializeRandomSpins();
  CalculateClusterEnergy();

  MonteCarloStatistics statistics(num_SA_steps+1, num_MC_sweeps);

  long double temp_T;
  uint single_sweep_accept, index;
  for (uint temp_counter=0;temp_counter<num_SA_steps+1;++temp_counter){
    temp_T = pow(cooling_rate,temp_counter)*initial_T;
    for (uint sweep = 0; sweep < num_MC_sweeps; sweep++){
      for (uint i=0;i<OverrelaxMCRatio;i++){
        OverrelaxationSweep(); //perform overrelaxation sweeps
      }
      single_sweep_accept=0; // reset total acceptance per sweep
      MetropolisSweep(temp_T,single_sweep_accept); //perform Metropolis sweep
      if (RecordStats == true){ //save statistics of every sweep
        index = sweep+temp_counter*num_MC_sweeps;
        statistics.EnergyDensity[index] = (double)Lattice.ClusterEnergy/(double)Lattice.NumSites;
        statistics.AcceptanceRate[index] = (double)single_sweep_accept/(double)Lattice.NumSites;
      }
    }
  }

  DeterministicSweeps(num_D_sweeps);
  CalculateClusterEnergy();

  MPI_Allgather(&Lattice.ClusterEnergy,1,MPI_LONG_DOUBLE,
                minarr, 1, MPI_LONG_DOUBLE,MPI_COMM_WORLD);
  std::vector<long double> v(minarr,minarr+MPISize);
  int firstlowestrank = std::min_element(v.begin(),v.end()) - v.begin();
  free(minarr);

  // if (MPIRank == 0){
  //   for (uint i=0; i<MPISize; i++){
  //     cout << i << " " << v[i]/(double)Lattice.NumSites <<  endl;
  //   }
  //   cout << firstlowestrank << endl;
  // }

  if (MPIRank == firstlowestrank){
    if (RecordStats == true){
      statistics.WriteStatisticsFile();
    }
    PrintConfiguration(out);
    Lattice.CalculateClusterOP();
    Lattice.PrintOP(out);
  }
}

void MonteCarlo::PerformFiniteT(std::ostream &out, const uint& num_thermal_sweeps,
                                                   const uint& num_measuring_sweeps,
                                                   const uint& sampling_time)
{
  const double initial_T = 0;
  const double num_MC_sweeps = 0;
  const uint num_D_sweeps = 0;
  Lattice.PrintLatticeParameters(out);
  PrintSimulationData(out, initial_T, FinalT, num_MC_sweeps,
               num_thermal_sweeps, num_measuring_sweeps,
               sampling_time, num_D_sweeps);
  Lattice.PrintHamiltonianParameters(out);

  InitializeFMSpins(pi/2.0,pi/2.0);
  // InitializeRandomSpins();
  CalculateClusterEnergy();

  MonteCarloStatistics statistics(1, num_thermal_sweeps);

  uint single_sweep_accept, index;          //sweep counter
  for (uint sweep = 0; sweep < num_thermal_sweeps; sweep++){
    for (uint i=0;i<OverrelaxMCRatio;i++){
      OverrelaxationSweep(); //perform overrelaxation sweeps
    }
    single_sweep_accept=0; // reset total acceptance per sweep
    MetropolisSweep(FinalT,single_sweep_accept); //perform Metropolis sweeps
    if (RecordStats == true){//save statistics
      index = sweep;
      statistics.EnergyDensity[index] = (double)Lattice.ClusterEnergy/(double)Lattice.NumSites;
      statistics.AcceptanceRate[index] = (double)single_sweep_accept/(double)Lattice.NumSites;
    }
  }

  if (RecordStats == true){
    statistics.WriteStatisticsFile();
  }

  PrintConfiguration(out);
  Lattice.CalculateClusterOP();
  Lattice.PrintOP(out);
}

void MonteCarlo::PrintConfiguration(std::ostream &out){
  uint flat_index;
  out << "--------------------------------Final Configuration--------------------------------\n";
  out << "Energy per site\n";
  out << std::setprecision(14) << Lattice.ClusterEnergy/Lattice.NumSites << "\n";
  out << "Spin configuration\n";
  for (uint y=0; y<Lattice.L2; ++y){
    for (uint x=0; x<Lattice.L1; ++x){
      for (uint sub=0; sub<Lattice.NumSublattices;++sub){
        BravaisIndicesToFlat(sub, x, Lattice.NumSublattices, y, Lattice.L1*Lattice.NumSublattices, flat_index);
        out << std::setprecision(14) << x << " " << y << " " << sub << " " << Lattice.Cluster[flat_index].transpose() << "\n";
      }
    }
  }
}

void MonteCarlo::PrintSimulationData(std::ostream &out,
                         const double& T_i, const double& T_f, const uint& max_sa_sweeps,
                         const uint& max_thermal_sweeps, const uint& max_measuring_sweeps,
                         const uint& sampling_time, const uint& max_dsweeps)
{
  out << "Initial Temperature\n";
  out << T_i << "\n";
  out << "Final Temperature\n";
  out << T_f << "\n";
  out << "Number of Metropolis sweeps (SA, thermal, measuring/sampling time)\n";
  out << max_sa_sweeps << " " << max_thermal_sweeps << " " << max_measuring_sweeps << "/" << sampling_time << "\n";
  out << "Number of Deterministic sweeps\n";
  out << max_dsweeps << "\n";
}

void MonteCarlo::SpherePointPicker(Vector3LD& some_spin)
{
  long double u = UnitInterval(RNG);
  long double v = UnitInterval(RNG);
  long double theta = acos(2*u-1); long double phi = 2*pi*v;
  some_spin =SphericalAnglesToCubic(theta,phi);
}

void MonteCarlo::OverrelaxationSweep(){
  Vector3LD *chosen_site_ptr;
  Vector3LD old_spin_vec,molec_field, spindiff;
  uint flat_index;

  uint flip=0;
  while (flip<Lattice.NumSites){
    flat_index = SiteDist(RNG);
    chosen_site_ptr = &(Lattice.Cluster[flat_index]); //pointer to spin vector
    old_spin_vec = *chosen_site_ptr;    //saving it to calculate the OP dynamically

    MolecularField(flat_index, molec_field);
    molec_field.normalize();
    //overrelaxation step
    *chosen_site_ptr = (-old_spin_vec + 2.0*molec_field.dot(old_spin_vec)*molec_field).normalized();
    flip++;
  }
}

void MonteCarlo::MetropolisSweep(const double& temperature, uint& single_sweep_accept)
{
  uint n1, n2, flat_index;
  long double old_local_energy, energydiff, new_local_energy;
  Vector3LD *chosen_site_ptr;
  Vector3LD spindiff,old_spin_at_chosen_site;

  uint flip = 0;
  while (flip < Lattice.NumSites){
    flat_index = SiteDist(RNG);

    CalculateLocalEnergy(flat_index, old_local_energy);

    chosen_site_ptr = &(Lattice.Cluster[flat_index]);
    old_spin_at_chosen_site = *chosen_site_ptr;
    SpherePointPicker(*chosen_site_ptr); //selection of angle, currently uniform update
    CalculateLocalEnergy(flat_index, new_local_energy);
    energydiff = new_local_energy-old_local_energy;

    if (UnitInterval(RNG) < std::min(exp(-energydiff/temperature),1.0)){
      ++single_sweep_accept;
      Lattice.ClusterEnergy+=energydiff;
    }else{
      *chosen_site_ptr = old_spin_at_chosen_site;
    }
    ++flip;
  }
}


void MonteCarlo::ThermalizeConfiguration(const double& temp, const uint& max_sweeps, uint& all_sweep_accept)
{
  uint single_sweep_accept;          //accepted moves per sweep

  uint sweep = 0;
  while (sweep < max_sweeps){
    for (uint i=0;i<OverrelaxMCRatio;i++){
      OverrelaxationSweep();
    }
    MetropolisSweep(temp,single_sweep_accept);
    all_sweep_accept += single_sweep_accept;
    ++sweep;
  }
}

void MonteCarlo::DeterministicSweeps(const uint& num_D_sweeps)
{
  uint align;
  uint uc_x, uc_y, flat_index;
  Vector3LD *chosen_site_ptr;

  Vector3LD old_spin_vec;
  Vector3LD molec_field;

  uint sweep = 0;
  while (sweep < num_D_sweeps){
    align = 0;
    while (align < Lattice.NumSites){
        flat_index = SiteDist(RNG);
        chosen_site_ptr = &(Lattice.Cluster[flat_index]);

        old_spin_vec = *chosen_site_ptr;
        MolecularField(flat_index, molec_field);
        molec_field.normalize();
        *chosen_site_ptr = molec_field;
        align++;
    }
    sweep++;
  }
}
