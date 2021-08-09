///  @file     triangular.hpp
///  @author   Ahmed Rayyan
///  @date     July 23, 2021
///  @brief    defining a separate Monte Carlo class
#ifndef MC_HPP
#define MC_HPP

#include "triangular.hpp"

class MonteCarloStatistics{
public:
  MonteCarloStatistics(const uint& num_temps, const uint& num_sweeps);
  void WriteStatisticsFile();
  vector<long double> EnergyDensity, AcceptanceRate;
private:
  uint NumTemps, NumSweeps;
};

class MonteCarlo
{
  public:
    // TriLatt Lattice;
    // MonteCarlo(TriLatt& lattice, const double& final_T, const uint& num_overrelax_ratio,
    //                                                     const bool& recordstats);
    HoneyLatt Lattice;
    MonteCarlo(HoneyLatt& lattice, const double& final_T, const uint& num_overrelax_ratio,
                                                        const bool& recordstats);

    void PerformSimulation(std::ostream &out);
    void PerformSimulatedAnnealing(std::ostream &out, const double& cooling_rate,
                                 const double& initial_T, const uint& num_SA_steps,
                                                          const uint& num_MC_sweeps,
                                                          const uint& num_D_sweeps);
    void PerformFiniteT(std::ostream &out, const uint& num_thermal_sweeps,
                                           const uint& max_measuring_sweeps,
                                           const uint& sampling_time);

  private:
    std::mt19937 RNG;
    std::uniform_real_distribution<double> UnitInterval;
    std::uniform_int_distribution<uint> SiteDist;

    const double FinalT;
    const uint OverrelaxMCRatio;

    const bool RecordStats;

    void InitializeRandomSpins();
    void InitializeFMSpins(const long double& theta, const long double& phi);
    void MolecularField(const uint& flat_index, Vector3LD& molec);
    void CalculateLocalEnergy(const uint& flat_index, long double& energy);
    void CalculateClusterEnergy();
    void PrintConfiguration(std::ostream &out);
    void SpherePointPicker(Vector3LD& some_spin);

    void OverrelaxationSweep();
    void MetropolisSweep(const double& temperature, uint& single_sweep_accept);
    void DeterministicSweeps(const uint& num_D_sweeps);

    void ThermalizeConfiguration(const double& temp, const uint& max_sweeps, uint& all_sweep_accept);
    void PrintSimulationData(std::ostream &out,
         const double& T_i, const double& T_f, const uint& max_sa_sweeps,
         const uint& max_thermal_sweeps, const uint& max_measuring_sweeps,
         const uint& sampling_time, const uint& max_dsweeps);
};

#endif
