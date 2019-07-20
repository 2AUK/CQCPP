#pragma once
#include <iostream>
#include <tuple>
#include <eigen3/Eigen/Core>
#include <array>

class BasisFunction
{
public:
  std::array<int,3> shell; //Cartesian Angular Momentum Shell
  Eigen::ArrayXf norm; //Normalisation constant for orbital - calculated using normalize() function - with values for each primitive gaussian function
  Eigen::ArrayXf origin; //Orbital origin (essentially coordinates of atom)
  Eigen::ArrayXf exps; //Exponents read from basis-set file
  Eigen::ArrayXf coefs; //Coefficients read from basis-set file
  BasisFunction(std::array<int, 3>, Eigen::ArrayXf, Eigen::ArrayXf, Eigen::ArrayXf); //Constructor
  void normalize(); //Normalization function - this sets norm to a value
};
