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
#include <fstream>
#include <mpi.h>
#include <random>
#include <sstream>
#include <vector>

using std::cout; using std::endl; using std::vector; using std::get;

typedef unsigned int uint;

typedef Eigen::Matrix<long double, 2, 1> Vector2LD;
typedef Eigen::Matrix<long double, 3, 1> Vector3LD;
typedef Eigen::Matrix<long double, 3, 3> Matrix3LD;

typedef Eigen::Array<long double, Eigen::Dynamic, 1> ArrayXLD;
typedef Eigen::Array<long double, Eigen::Dynamic, Eigen::Dynamic> ArrayXXLD;

#define pi M_PI

#define FIELD_DIR const Vector3LD A_Dir = Vector3LD(1/sqrt(6), 1/sqrt(6), -2/sqrt(6));\
                  const Vector3LD B_Dir = Vector3LD(-1/sqrt(2), 1/sqrt(2), 0);\
                  const Vector3LD C_Dir = Vector3LD(1/sqrt(3), 1/sqrt(3), 1/sqrt(3));

#define LATTICE_DIR const Vector2LD a1 = Vector2LD( 1.0/2.0, sqrt(3)/2.0);\
                    const Vector2LD a2 = Vector2LD(-1.0/2.0, sqrt(3)/2.0);


void BravaisIndicesToFlat(const uint& x, const uint&y, const uint & period1,
                                         const uint&z, const uint & period2,
                                                           uint& flat_index);

long double Lorentzian(const long double& s, const long double& x, const long double& l);

long double Gaussian(const long double& s, const long double& x, const long double& l);

Vector3LD SphericalAnglesToCubic(const long double& theta, const long double& phi);


#endif // COMMON_HPP
