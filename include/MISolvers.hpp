#pragma once
#include "Molecule.hpp"
#include "BasisFunc.hpp"
#include <string>

class IntegralSolver{
public:
  virtual Eigen::ArrayXXf integral(std::string) = 0;
};

class THO : public IntegralSolver{
public:
  THO(Molecule&);
  Eigen::ArrayXXf integral(std::string);
  Molecule system;
private:
  float f(float, float, float, float, float);
  float S(BasisFunction, BasisFunction);
  float overlap(std::tuple<int, int, int>, Eigen::ArrayXf, float, std::tuple<int, int, int>, Eigen::ArrayXf, float); 
  float overlap_1d(int, int, float, float, float);
};
