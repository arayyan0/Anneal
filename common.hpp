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
#include <complex>
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

typedef Eigen::Matrix<long double,              3, 1> Vector3LD;

typedef Eigen::Matrix<long double,              3,              3> Matrix3LD;
typedef Eigen::Matrix<long double,              3, Eigen::Dynamic> Matrix3XLD;
typedef Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> MatrixXXLD;
typedef Eigen::Matrix<std::complex<long double>, Eigen::Dynamic, Eigen::Dynamic> MatrixXXLDc;

typedef Eigen::Array<             long double , Eigen::Dynamic, Eigen::Dynamic> ArrayXXLD;
typedef Eigen::Array<std::complex<long double>, Eigen::Dynamic, Eigen::Dynamic> ArrayXXLDc;

#define pi M_PI

#define ABC_DIR const Vector3LD A_Dir = Vector3LD( 1/sqrt(6), 1/sqrt(6), -2/sqrt(6));\
                  const Vector3LD B_Dir = Vector3LD(-1/sqrt(2), 1/sqrt(2),          0);\
                  const Vector3LD C_Dir = Vector3LD( 1/sqrt(3), 1/sqrt(3),  1/sqrt(3));

void ThreeDBravaisIndicesToFlat(const uint& w, const uint&x, const uint & periodx,
                                               const uint&y, const uint & periody,
                                               const uint&z, const uint & periodz,
                                               uint& flat_index);

long double Lorentzian(const long double& r, const long double& l);

void SphericalAnglesToCubic(const long double& theta, const long double& phi, Vector3LD& some_spin);

void PBCIndices(const uint& i, const uint& length, int& lower_i, int& higher_i);

bool EqualWithinEpsilon(const long double& a, const long double& b, const long double& eps);

bool IsThisFloatAnIntegerWithinEpsilon(const long double& a, const long double& eps);


#endif // COMMON_HPP
