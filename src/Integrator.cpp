#include "Integrator.hpp"
#include <math.h>
#include <experimental/array>
#include "utils.hpp"

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

double Integrator::E(int i, int j, int t, double Qx, double a, double b){
  double p = a + b;
  double q = (a * b) / p;
  if (t < 0 || t > (i + j)){
    return 0.0;
  } else if (t == 0 && i == 0 && j == 0){
    return exp(-q * Qx * Qx);
  } else if (j == 0){
    return (1 / (2 * p)) * E(i-1, j, t-1, Qx, a, b) - ((q * Qx) / a) * E(i-1, j, t, Qx, a, b) + (t+1) * E(i-1, j, t+1, Qx, a, b);
  } else {
    return (1 / (2 * p)) * E(i, j-1, t-1, Qx, a, b) + ((q * Qx) / b) * E(i, j-1, t, Qx, a, b) + (t+1) * E(i, j-1, t+1, Qx, a, b);
    
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
double Integrator::overlap(Eigen::ArrayXd A, std::array<int, 3> lmn1, double a, Eigen::ArrayXd B, std::array<int, 3> lmn2, double b){
  
  double p = a + b;
  int l1 = lmn1[0]; int m1 = lmn1[1]; int n1 = lmn1[2];
  int l2 = lmn2[0]; int m2 = lmn2[1]; int n2 = lmn2[2];

  double Sx = E(l1, l2, 0, A[0] - B[0], a, b);
  double Sy = E(m1, m2, 0, A[1] - B[1], a, b);
  double Sz = E(n1, n2, 0, A[2] - B[2], a, b);

  return pow(M_PI/p, 1.5) * Sx * Sy * Sz;
}

/**
 *  \brief Contracted gaussian overlap distribution
 *  
 *  This function contracts the primitive gaussians and calculates the contracted gaussian overlap distribution.
 *  \param bf1 first BasisFunction object containing all the relevant information for integral calculation
 *  \param bf2 second BasisFunction object containing all the relevant information for integral calculation
 *  \return Value of the overlap between two contracted gaussians.
 */

double Integrator::S(BasisFunction bf1, BasisFunction bf2){
  double total = 0;
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
Eigen::ArrayXXd Integrator::SMatrix(){
  int aos = system.nCGFs;
  Eigen::ArrayXXd retmat = Eigen::ArrayXXd(aos, aos);
  for (int i = 0; i < aos; i++){
    for (int j = 0; j < aos; j++){
      retmat(i, j) = S(system.cgbfs[i], system.cgbfs[j]);
    }
  }
  return retmat;
}


double Integrator::kinetic(Eigen::ArrayXd A, std::array<int, 3> lmn1, double a, Eigen::ArrayXd B, std::array<int, 3> lmn2, double b){
  int l1 = lmn1[0]; int m1 = lmn1[1]; int n1 = lmn1[2];
  int l2 = lmn2[0]; int m2 = lmn2[1]; int n2 = lmn2[2];

  double p = a + b;

  double Sx = E(l1, l2, 0, A[0] - B[0], a, b);
  double Sy = E(m1, m2, 0, A[1] - B[1], a, b);
  double Sz = E(n1, n2, 0, A[2] - B[2], a, b);

  double Tx = (l2*(l2-1)*E(l1, l2-2, 0, A[0]-B[0], a, b)) + (-2*b*(2*l2 + 1)*Sx) + (4*b*b*E(l1, l2+2, 0, A[0]-B[0], a, b));
  Tx *= Sy;
  Tx *= Sz;

  double Ty = (m2*(m2-1)*E(m1, m2-2, 0, A[1]-B[1], a, b)) + (-2*b*(2*m2 + 1)*Sy) + (4*b*b*E(m1, m2+2, 0, A[1]-B[1], a, b));
  Ty *= Sx;
  Ty *= Sz;

  double Tz = (n2*(n2-1)*E(n1, n2-2, 0, A[2]-B[2], a, b)) + (-2*b*(2*n2 + 1)*Sz) + (4*b*b*E(n1, n2+2, 0, A[2]-B[2], a, b));
  Tz *= Sy;
  Tz *= Sx;
  
  return -0.5 * (Tx + Ty + Tz) *  pow(M_PI/p, 1.5);
}

double Integrator::T(BasisFunction bf1, BasisFunction bf2){
  double total = 0;
  for (int i = 0; i < bf1.coefs.size(); i++){
    for (int j = 0; j < bf2.coefs.size(); j++){
      total += bf1.norm[i] * bf2.norm[j] * bf1.coefs[i] * bf2.coefs[j] * kinetic(bf1.origin, bf1.shell, bf1.exps[i], bf2.origin, bf2.shell, bf2.exps[j]);
    }
  }
  return total;
}

Eigen::ArrayXXd Integrator::TMatrix(){
  int aos = system.nCGFs;
  Eigen::ArrayXXd retmat = Eigen::ArrayXXd(aos, aos);
  for (int i = 0; i < aos; i++){
    for (int j = 0; j < aos; j++){
      retmat(i, j) = T(system.cgbfs[i], system.cgbfs[j]);
    }
  }
  return retmat;
}

double Integrator::R(int t, int u, int v, int n, Eigen::ArrayXd PC, double p, double RPC){
  double value = 0.0
  if (t == 0 && u == 0 && v == 0){
    val += pow(-2 * p, n) * boys(n, p * RPC * RPC);
  }
}
