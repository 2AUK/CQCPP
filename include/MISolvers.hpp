#pragma once
#include "Molecule.hpp"
#include <string>

class IntegralSolver{
public:
  virtual float one_electron_integral(std::string) = 0;
  virtual float two_electron_integral(std::string) = 0;
  virtual Eigen::ArrayXXf one_electron_integral_matrix() = 0;
  virtual Eigen::ArrayXXf two_electron_integral_matrix() = 0;
};

class THO : public IntegralSolver{
public:
  THO(Molecule);
  float one_electron_integral(std::string);
  EigenArrayXXf one_electron_integral_matrix();
private:
  float binomial_factor();
};
