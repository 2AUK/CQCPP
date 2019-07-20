#pragma once
#include <iostream>
#include <tuple>
#include <eigen3/Eigen/Core>
#include <array>

/**
 *   \file BasisFunc.hpp
 *   \brief BasisFunction class encapsulating information on basis function
 *
 *  The BasisFunction class contains all the information required to drive the integral evaluation maths.
 *
 */


class BasisFunction
{
public:
  std::array<int,3> shell; /**< Cartesian Angular Momentum Shell. */
  Eigen::ArrayXf norm; /**< Normalisation constant for orbital with values for each primitive gaussian function. */
  Eigen::ArrayXf origin; /**< Orbital origin (essentially coordinates of atom). */
  Eigen::ArrayXf exps; /**< Array of exponents of primitive gaussians of each basis function. */
  Eigen::ArrayXf coefs; /**< Array of coefficients of primitive gaussians of each basis function. */
  BasisFunction(std::array<int, 3>, Eigen::ArrayXf, Eigen::ArrayXf, Eigen::ArrayXf);
  void normalize();
};
