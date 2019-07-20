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
  Integrator(Molecule&); /**< Constructor for Integrator object. */
  Molecule system; /*!< System for which integrals are generated. */
  Eigen::ArrayXXf SMatrix(); /*!< Routine for generating overlap matrix from CGFs. */
  Eigen::ArrayXXf TMatrix(); /**< Routine for generating kinetic matrix from CGFs. */
  Eigen::ArrayXXf VMatrix(); /**< Routine for generating nuclear attraction matrix from CGFs. */
  Eigen::ArrayXf ERIMatrix();  /**< Routine for generating electron repulsion matrix from CGFs. */
  Eigen::ArrayXXf JMatrix(); /**< Routine for generating multipole moment matrix from CGFs. */
private:
  float E(int, int, int, float, float, float); /*!< Hermite expansion coefficients for calculating overlap from hermite gaussian functions. Implementation details in Integrator.cpp */
  float overlap(Eigen::ArrayXf, std::array<int, 3>, float, Eigen::ArrayXf, std::array<int, 3>, float); /*!< Function to calculate overlap between primitive gaussian functions. Implementation details in Integrator.cpp */
};
