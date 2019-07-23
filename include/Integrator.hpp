#pragma once
#include "Molecule.hpp"
#include "BasisFunc.hpp"
#include <string>

/**
 *   \file Integrator.hpp
 *   \brief Integral evaluation module
 *
 *   McMurchie-Davidson scheme for integral evaluation. Helgaker, Trygve, and Peter R. Taylor. “Gaussian basis sets and molecular integrals.” Modern Electronic Structure (1995).
 */

class Integrator{
public:
  Integrator(Molecule&);
  Molecule system; /**< System for which integrals are generated. */
  Eigen::ArrayXXd SMatrix(); 
  Eigen::ArrayXXd TMatrix(); 
  Eigen::ArrayXXd VMatrix(); 
  Eigen::ArrayXd ERIMatrix();  
  Eigen::ArrayXXd JMatrix(); 
  double E(int, int, int, double, double, double); 
  double overlap(Eigen::ArrayXd, std::array<int, 3>, double, Eigen::ArrayXd, std::array<int, 3>, double); 
  double S(BasisFunction, BasisFunction);
  double kinetic(Eigen::ArrayXd, std::array<int, 3>, double, Eigen::ArrayXd, std::array<int, 3>, double);
  double T(BasisFunction, BasisFunction);
  double R(int, int, int, int, Eigen::ArrayXd, double, double);
  double nuclear(Eigen::ArrayXd, std::array<int, 3>, double, Eigen::ArrayXd, std::array<int, 3>, double, Eigen::ArrayXd);
  double V(BasisFunction, BasisFunction, Eigen::ArrayXd);
};
