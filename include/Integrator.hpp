#pragma once
#include "Molecule.hpp"
#include "BasisFunc.hpp"
#include <string>

/**
 *   \file Integrator.hpp
 *   \brief Integral evaluation module
 *
 *   McMurchie-Davidson scheme for integral evaluation
 */

class Integrator{
public:
  Integrator(Molecule&);
  Molecule system; /**< System for which integrals are generated. */
  Eigen::ArrayXXf SMatrix(); 
  Eigen::ArrayXXf TMatrix(); 
  Eigen::ArrayXXf VMatrix(); 
  Eigen::ArrayXf ERIMatrix();  
  Eigen::ArrayXXf JMatrix(); 
private:
  float E(int, int, int, float, float, float); 
  float overlap(Eigen::ArrayXf, std::array<int, 3>, float, Eigen::ArrayXf, std::array<int, 3>, float); 
  float S(BasisFunction, BasisFunction);
};
