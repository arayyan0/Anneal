///  @file     common.hpp
///  @author   Ahmed Rayyan
///  @date     December 2, 2019
///  @brief    useful functions for the simulated annealing algorithm
#ifndef COMMON_HPP
#define COMMON_HPP

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include "Eigen/Dense"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <vector>

using std::cout; using std::endl; using std::vector;

typedef unsigned int uint;

#define pi M_PI

#define FIELD_DIR const Eigen::Vector3d A_Dir = Eigen::Vector3d(1/sqrt(6), 1/sqrt(6), -2/sqrt(6));\
                  const Eigen::Vector3d B_Dir = Eigen::Vector3d(-1/sqrt(2), 1/sqrt(2), 0);\
                  const Eigen::Vector3d C_Dir = Eigen::Vector3d(1/sqrt(3), 1/sqrt(3), 1/sqrt(3));
//
// #define AXIS_DIR const Eigen::Vector3d X(1, 0, 0);\
//                  const Eigen::Vector3d Y(0, 1, 0);\
//                  const Eigen::Vector3d Z(0, 0, 1);

namespace MyRandom
{
  extern std::mt19937 RNG;
  extern std::uniform_real_distribution<double> unit_interval;
}

void PrintSimulationData(std::ostream &out, const uint& hc_or_kek, const uint& type,
                         const uint& sublattice, const uint& l1, const uint& l2,
                         const double& T_i, double& T_f, const uint& max_mflips,
                         const uint& max_daligns);

void PrintTriangularSimulationData(std::ostream &out, const uint& type,
                        const uint& sublattice, const uint& l1, const uint& l2,
                        const double& T_i, double& T_f, const uint& max_sa_sweeps,
                        const uint& max_thermal_sweeps, const uint& max_measuring_sweeps,
                        const uint& sampling_time, const uint& max_daligns);

#endif // COMMON_HPP
