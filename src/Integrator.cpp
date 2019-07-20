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

float Integrator::E(int i, int j, int t, float Qx, float a, float b){
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
float Integrator::overlap(Eigen::ArrayXf A, std::array<int, 3> lmn1, float a, Eigen::ArrayXf B, std::array<int, 3> lmn2, float b){
  
  float p = a + b;
  float l1 = lmn1[0]; float m1 = lmn1[1]; float n1 = lmn1[2];
  float l2 = lmn2[0]; float m2 = lmn2[1]; float n2 = lmn2[2];

  float Sx = E(l1, l2, 0, A[0] - B[0], a, b);
  float Sy = E(m1, m2, 0, A[1] - B[1], a, b);
  float Sz = E(n1, n2, 0, A[2] - B[2], a, b);

  return pow(M_PI/p, 1.5) * Sx * Sy * Sz;
  std::cout << pow(M_PI/p, 1.5) * Sx * Sy * Sz << std::endl;
}

/**
 *  \brief Contracted gaussian overlap distribution
 *  
 *  This function contracts the primitive gaussians and calculates the contracted gaussian overlap distribution.
 *  \param bf1 first BasisFunction object containing all the relevant information for integral calculation
 *  \param bf2 second BasisFunction object containing all the relevant information for integral calculation
 *  \return Value of the overlap between two contracted gaussians.
 */

float Integrator::S(BasisFunction bf1, BasisFunction bf2){
  float total = 0;
  for (int i = 0; i < bf1.coefs.size(); i++){
    for (int j = 0; j < bf2.coefs.size(); j++){
      total += bf1.norm[i] * bf2.norm[j] * bf1.coefs[i] * bf2.coefs[j] * overlap(bf1.origin, bf1.shell, bf1.exps[i], bf2.origin, bf2.shell, bf2.exps[j]);
    }
  }
  return total;
}

/**
 *  \brief Generates overlap matrix from contracted gaussians
 * 
 *  Member function to compute overlap matrix from system attribute of Integrator class.
 *  \return Eigen Array who's elements are the overlap distribution of contracted gaussians.
 */
Eigen::ArrayXXf Integrator::SMatrix(){
  int aos = system.nCGFs;
  Eigen::ArrayXXf retmat = Eigen::ArrayXXf(aos, aos);
  for (int i = 0; i < aos; i++){
    for (int j = 0; j < aos; j++){
      retmat(i, j) = S(system.cgbfs[i], system.cgbfs[j]);
    }
  }
  return retmat;
}
