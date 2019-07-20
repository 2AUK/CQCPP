#pragma once
#include "Molecule.hpp"
#include "BasisFunc.hpp"
#include <string>

class Integrator{
public:
  Integrator(Molecule&);
  Molecule system;
  Eigen::ArrayXXf S(); //Overlap
  Eigen::ArrayXXf T(); //Kinetic
  Eigen::ArrayXXf V(); //Nuclear Attraction
  Eigen::ArrayXf ERI(); 
  Eigen::ArrayXXf J(); //Multipole moment
private:
  float E(int, int, int, float, float, float);
  float overlap(Eigen::ArrayXf, std::array<int, 3>, float, Eigen::ArrayXf, std::array<int, 3>, float);
};
