#include "Integrator.hpp"
#include <math.h>

/**
 *   \file Integrator.cpp
 *   \brief Implementation for Integrator
 *
 *  All the implementation details for all Integrator routines held in Integrator.hpp
 *
 */


/**
 *  \internal
 *
 *  \brief Constructor for Integrator object
 *
 *  \param Molecule & Object
 *  \return Integrator Object
 */
Integrator::Integrator(Molecule& input_molecule) : system(input_molecule){
  std::cout << "Molecule Object loaded in to Integrator system" << std::endl;
}

/**
 *  \brief Calculates Hermite expansion coefficients (E)
 *
 * This function calculates the Hermite-to-Cartesian expansion coefficients (E) used for expanding the overlap distribution in Hermite gaussians  
 *  \param i Angular momentum number taken from shell of first primitive gaussian
 *  \param j Angular momentum number taken from shell of second primitive gaussian
 *  \param t Number of nodes in Hermite gaussian
 *  \param Qx Difference between 1 dimension of coordinates of both primitive gaussians
 *  \param a Exponent of first primitive gaussian
 *  \param b Exponent of second primitive gaussian
 *  \return Full solution of the expansion coefficient generated recursively
 */

float E(int i, int j, int t, float Qx, float a, float b){
  float p = a + b;
  float q = (a * b) / p;
  if (i == 0 && j == 0 && t == 0){
    return(exp(-q * Qx * Qx));
  } else if (t < 0 || t > (i + j)){
    return 0.0;
  } else if (i == 0){
    return (1 / (2 * p)) * E(i, j-1, t-1, Qx, a, b) + ((q * Qx) / b) * E(i, j-1, t, Qx, a, b) + (t+1) * E(i, j-1, t+1, Qx, a, b);
  } else {
    return (1 / (2 * p)) * E(i-1, j, t-1, Qx, a, b) - ((q * Qx) / a) * E(i-1, j, t, Qx, a, b) + (t+1) * E(i-1, j, t+1, Qx, a, b);
  }
}
/**
 *  \brief Primitive gaussian overlap distribution
 *
 *  This function calculates the overlap distribution between two primitive gaussian functions.
 *  \param A Coordinates for first primitive gaussian
 *  \param lmn1 Angular momentum shell for first primitive gaussian
 *  \param a Exponent for first primitive gaussian
 *  \param B Coordinates for second primitive gaussian
 *  \param lmn2 Angular momentum shell for second primitive gaussian
 *  \param b Exponent for second primitive gaussian
 *  \return Value of the overlap between two primitive gaussians
 */
float overlap(Eigen::ArrayXf A, std::array<int, 3> lmn1, float a, Eigen::ArrayXf B, std::array<int, 3> lmn2, float b){
  
  float p = a + b;
  float l1 = lmn1[0]; float m1 = lmn1[1]; float n1 = lmn1[2];
  float l2 = lmn2[0]; float m2 = lmn2[1]; float n2 = lmn2[2];

  float Sx = E(l1, l2, 0, A[0] - B[0], a, b);
  float Sy = E(m1, m2, 0, A[1] - B[1], a, b);
  float Sz = E(n1, n2, 0, A[2] - B[2], a, b);

  return pow(M_PI/p, 1.5) * Sx * Sy * Sz;
}


