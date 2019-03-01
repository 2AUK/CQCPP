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
  THO(Molecule&);
  float one_electron_integral(std::string);
  Eigen::ArrayXXf one_electron_integral_matrix();
  float two_electron_integral(std::string);
  Eigen::ArrayXXf two_electron_integral_matrix();
  Molecule system;
private:
  float f(float, float, float, float, float);
  float overlap(std::tuple<int, int, int>, Eigen::ArrayXf, float, std::tuple<int, int, int>, Eigen::ArrayXf, float); 
  float overlap_1d(int, int, float, float, float);
};
