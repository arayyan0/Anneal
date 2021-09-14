///  @file     MC.hpp
///  @author   Ahmed Rayyan
///  @date     July 23, 2021
///  @brief    defining a separate Monte Carlo class
#include "MC.hpp"

MonteCarloStatistics::MonteCarloStatistics(const bool& recordstats, const uint& num_temps,
const uint& num_sweeps):
NumTemps(num_temps), NumSweeps(num_sweeps)
{
  if (recordstats==true){
    EnergyDensity.resize(NumTemps*NumSweeps+1);
    AcceptanceRate.resize(NumTemps*NumSweeps+1);
  }
}

void MonteCarloStatistics::WriteStatisticsFile()
{
  std::ofstream outfile;
  outfile.open("stats.out");
  outfile << std::fixed << std::setprecision(14);
  outfile << NumTemps << " " << NumSweeps <<  " " << 1 << endl;
  for (uint i=0; i<EnergyDensity.size();++i){
    outfile << i << " " << EnergyDensity[i] << " " << AcceptanceRate[i] << endl;
  }
  outfile.close();
}

MonteCarlo::MonteCarlo(Triangular& lattice, const double& final_T, const uint& num_overrelax_ratio,
// MonteCarlo::MonteCarlo(Honeycomb& lattice, const double& final_T, const uint& num_overrelax_ratio,
                       const bool& recordstats, const int& mpirank, const int& mpisize):
Lattice(lattice), FinalT(final_T),OverrelaxMCRatio(num_overrelax_ratio), RecordStats(recordstats),
MPIRank(mpirank), MPISize(mpisize)
{
  auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  std::mt19937 rng(seed);
  std::uniform_real_distribution<long double> unit_interval(0,1);
  std::uniform_int_distribution<uint> sd(0, Lattice.NumSites-1);

  RNG = rng;
  UnitInterval = unit_interval;
  SiteDist = sd;
}

void MonteCarlo::InitializeFMSpins(const long double& theta, const long double& phi)
{
  Vector3LD v;
  SphericalAnglesToCubic(theta, phi, v);
  for (uint flat_index=0; flat_index<Lattice.NumSites; ++flat_index){
    Lattice.Cluster.col(flat_index) = v;
  }
}

void MonteCarlo::InitializeRandomSpins()
{
  Vector3LD v;
  for (uint flat_index=0; flat_index<Lattice.NumSites; ++flat_index){
    SpherePointPicker(v);
    Lattice.Cluster.col(flat_index) = v;
  }
}

void MonteCarlo::MolecularField(const uint& flat_index, Vector3LD& molec)
{
  Vector3LD v = Vector3LD::Zero();
  for (auto &j : Lattice.ClusterInfo[flat_index].NearestNeighbours){
    v += -Lattice.Cluster.col(get<6>(j)).transpose()*get<3>(j);
  }
  molec = v.transpose()+ Lattice.hField.transpose();
}

void MonteCarlo::CalculateLocalEnergy(const uint& flat_index, long double& energy)
{
  Vector3LD molec;
  MolecularField(flat_index, molec);
  energy = -Lattice.Cluster.col(flat_index).dot(molec+ Lattice.hField);
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

  // InitializeFMSpins(pi/2.0,pi/2.0);
  InitializeRandomSpins();
  CalculateClusterEnergy();

  MonteCarloStatistics statistics(RecordStats, num_SA_steps+1, num_MC_sweeps);

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

  if (RecordStats == true){//save statistics
    index++;
    statistics.EnergyDensity[index] = (double)Lattice.ClusterEnergy/(double)Lattice.NumSites;
    statistics.AcceptanceRate[index] = 1;
  }

  MPI_Allgather(&Lattice.ClusterEnergy,1,MPI_LONG_DOUBLE,
                minarr, 1, MPI_LONG_DOUBLE,MPI_COMM_WORLD);
  std::vector<long double> v(minarr,minarr+MPISize);
  int firstlowestrank = std::min_element(v.begin(),v.end()) - v.begin();
  free(minarr);

  // if (MPIRank == 0){
  //   for (uint i=0; i<MPISize; i++){
  //     std::cerr << i << " " << v[i]/(double)Lattice.NumSites <<  endl;
  //   }
  //   std::cerr << firstlowestrank << endl;
  // }

  if (MPIRank == firstlowestrank){
    if (RecordStats == true){
      statistics.WriteStatisticsFile();
    }
    PrintConfiguration(out);
    // Lattice.CalculateClusterOP();
    // Lattice.PrintOP(out);

    // Vector3LD aaa = Vector3LD::Zero();
    // for (uint flat_index=0; flat_index<Lattice.NumSites; ++flat_index){
    //   MolecularField(flat_index, aaa);
    //   std::cerr << aaa.normalized().transpose() << endl;
    // }

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

  // InitializeFMSpins(pi/2.0,pi/2.0);
  InitializeRandomSpins();
  CalculateClusterEnergy();

  //thermalization
  uint single_sweep_accept, index;          //sweep counter
  for (uint sweep = 0; sweep < num_thermal_sweeps; sweep++){
    for (uint i=0;i<OverrelaxMCRatio;i++){
      OverrelaxationSweep(); //perform overrelaxation sweeps
    }
    single_sweep_accept=0; // reset total acceptance per sweep
    MetropolisSweep(FinalT,single_sweep_accept); //perform Metropolis sweeps
  }

  //measurement sweeps
  long double ebar = 0;
  long double e2bar = 0;
  // vector<Vector3LD> ssfbar(Lattice.SSFPoints.size(), Vector3LD::Zero());
  uint n_samples = num_measuring_sweeps/sampling_time;
  for (uint sweep = 0; sweep < num_measuring_sweeps; sweep++){
    for (uint i=0;i<OverrelaxMCRatio;i++){
      OverrelaxationSweep(); //perform overrelaxation sweeps
    }
    single_sweep_accept=0; // reset total acceptance per sweep
    MetropolisSweep(FinalT,single_sweep_accept); //perform Metropolis sweeps

    //measuring in this loop
    if (sweep%sampling_time==0){
      ebar  += Lattice.ClusterEnergy;
      e2bar += pow(Lattice.ClusterEnergy,2);

      Lattice.CalculateClusterOP();

      // for (uint kk=0; kk < Lattice.SSFPoints.size(); kk++){
      //   ssfbar[kk] += Lattice.ClusterSSf[kk];
      // }

    }
  }
  Lattice.EBar  = ebar/n_samples;
  Lattice.E2Bar = e2bar/n_samples;
  // for (uint jj=0; jj < Lattice.SSFPoints.size(); jj++){
  //   Lattice.SSfBar.push_back(ssfbar[jj]/n_samples);
  // }

  PrintConfiguration(out);
  Lattice.CalculateClusterOP();
  // Lattice.PrintOP(out);
  // Lattice.PrintThermalObservables(out);
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
        // out << std::setprecision(14) << x << " " << y << " " << sub << " " << Lattice.Cluster[flat_index].transpose() << "\n";
        out << std::setprecision(14) << x << " " << y << " " << sub << " " << Lattice.Cluster.col(flat_index).transpose() << "\n";
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
  long double theta = acos(2*UnitInterval(RNG)-1);
  long double phi = 2*pi*UnitInterval(RNG);
  SphericalAnglesToCubic(theta,phi, some_spin);
}

void MonteCarlo::OverrelaxationSweep(){
  Vector3LD v_before, molec_field;
  uint flat_index;

  uint flip=0;
  while (flip<Lattice.NumSites){
    flat_index = SiteDist(RNG);
    v_before = Lattice.Cluster.col(flat_index); //pointer to spin vector

    MolecularField(flat_index, molec_field);
    molec_field.normalize();
    //overrelaxation step
    Lattice.Cluster.col(flat_index) = (-v_before + 2.0*molec_field.dot(v_before)*molec_field).normalized();
    flip++;
  }
}

void MonteCarlo::MetropolisSweep(const double& temperature, uint& single_sweep_accept)
{
  long double old_local_energy, energydiff, new_local_energy;
  Vector3LD v_before, v_after;
  uint flat_index;

  uint flip = 0;
  while (flip < Lattice.NumSites){
    flat_index = SiteDist(RNG);
    v_before = Lattice.Cluster.col(flat_index);
    CalculateLocalEnergy(flat_index, old_local_energy);

    SpherePointPicker(v_after); //selection of angle, currently uniform update
    Lattice.Cluster.col(flat_index) = v_after;

    CalculateLocalEnergy(flat_index, new_local_energy);
    energydiff = new_local_energy-old_local_energy;

    if (UnitInterval(RNG) < std::min<long double>(exp(-energydiff/temperature),1.0)){
      ++single_sweep_accept;
      Lattice.ClusterEnergy+=energydiff;
    }else{
      Lattice.Cluster.col(flat_index) = v_before;
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
  uint flat_index;

  Vector3LD molec_field;

  uint sweep = 0;
  while (sweep < num_D_sweeps){
    align = 0;
    while (align < Lattice.NumSites){
        flat_index = SiteDist(RNG);
        MolecularField(flat_index, molec_field);
        molec_field.normalize();
        Lattice.Cluster.col(flat_index) = molec_field;
        align++;
    }
    sweep++;
  }
}
