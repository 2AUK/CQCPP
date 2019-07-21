#include "BasisFunc.hpp"
#include "utils.hpp"
#include <cmath>
#include <eigen3/Eigen/Core>

/**
 *   \file BasisFunc.cpp
 *   \brief Implementation for BasisFunction member functions
 *
 *   All the implementation details for all Molecule routines held in BasisFunc.hpp 
 */

/**
 *  \brief Constructor for BasisFunction object
 *
 *  \param a Cartesian angular momentum shell
 *  \param b Origin of basis function
 *  \param c Array of exponents of basis functions
 *  \param d Array of coefficients of basis functions
 *  \return BasisFunction object
 */
BasisFunction::BasisFunction(std::array<int, 3> a, 
			     Eigen::ArrayXd b,
			     Eigen::ArrayXd c,
			     Eigen::ArrayXd d)
{
  shell = a;
  origin = b;
  exps = c;
  coefs = d;
}


/**
 *  \brief Normalises the basis functions
 *
 *  This function populates the norm attribute with the normalisation constants, calculated from the exponent and shell information. Taketa, H., Huzinaga, S. and O-ohata, K., 1966. Gaussian-expansion methods for molecular integrals. Journal of the physical society of Japan, 21(11), pp.2313-2324.
 *  \return void
 */
void BasisFunction::normalize(){
  int l = shell[0]; int m = shell[1]; int n = shell[2];
  //float L = l+m+n;
  Eigen::ArrayXd num = std::pow(2.0, 2.0*(l+m+n) + 3.0/2.0) * exps.pow((l+m+n)+3.0/2.0);
  float denom = double_factorial(2*l-1) * double_factorial(2*m-1) * double_factorial(2*n-1) * std::pow(M_PI, 3.0/2.0);
  Eigen::ArrayXd nd = num / denom;
  norm = nd.sqrt();
  //This code isn't fully working
  //Double check notes on CGBF normalisation
  // float N = 0.0;
  // float prefac = std::pow(M_PI, 1.5) * double_factorial(2*l-1) * double_factorial(2*m-1) * double_factorial(2*n-1) / std::pow(2.0, L);
  // int nexps = exps.size();
  // for (int i = 0; i < nexps; i++){
  //   for (int j = 0; j < nexps; j++){
  //     N += norm[i] * norm[j] * coefs[i] * coefs[j] * std::pow(exps[i] + exps[j], L+1.5);
  //   }
  // }
  // N *= prefac;
  // N = std::pow(N, -0.5);
  // for (int k = 0; k < nexps; k++){
  //   coefs[k] *= N;
  // }
}
