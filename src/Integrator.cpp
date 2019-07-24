#include "Integrator.hpp"
#include <math.h>
#include <experimental/array>
#include "utils.hpp"

#define INDEX(i, j) i>j ? i*(i+1)/2 + j : j*(j+1)/2 + i


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
    return ((1 / (2 * p)) * E(i-1, j, t-1, Qx, a, b)) - ((q * Qx / a) * E(i-1, j, t, Qx, a, b)) + ((t+1) * E(i-1, j, t+1, Qx, a, b));
  } else {
    return ((1 / (2 * p)) * E(i, j-1, t-1, Qx, a, b)) + ((q * Qx / b) * E(i, j-1, t, Qx, a, b)) + ((t+1) * E(i, j-1, t+1, Qx, a, b));
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
  const int l1 = lmn1[0]; const int m1 = lmn1[1]; const int n1 = lmn1[2];
  const int l2 = lmn2[0]; const int m2 = lmn2[1]; const int n2 = lmn2[2];
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
  Eigen::ArrayXXd retmat = Eigen::ArrayXXd::Zero(aos, aos);
  for (int i = 0; i < aos; i++){
    for (int j = 0; j < aos; j++){
      retmat(i, j) = S(system.cgbfs[i], system.cgbfs[j]);
    }
  }
  return retmat;
}


double Integrator::kinetic(Eigen::ArrayXd A, std::array<int, 3> lmn1, double a, Eigen::ArrayXd B, std::array<int, 3> lmn2, double b){
  const int l1 = lmn1[0]; const int m1 = lmn1[1]; const int n1 = lmn1[2];
  const int l2 = lmn2[0]; const int m2 = lmn2[1]; const int n2 = lmn2[2];

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

double Integrator::R(int t, int u, int v, int n, double PCx, double PCy, double PCz, double p, double RPC){
  double val = 0.0;
  // if (t == 0 && u == 0 && v == 0){
  //   val += pow(-2 * p, n) * boys(n, p * RPC * RPC);
  // } else if(t < 0 || u < 0 || v < 0){
  //   return 0.0;
  // } else if(t == 0 && u == 0){
  //   val += (v-1)*R(t, u, v-2, n+1, PC, p, RPC) + PC[2] * R(t, u, v, n+1, PC, p, RPC);
  // } else if(t == 0){
  //   val += (u-1)*R(t, u-2, v, n+1, PC, p, RPC) + PC[1] * R(t, u, v, n+1, PC, p, RPC);
  // } else {
  //   val += (t-1)*R(t-2, u, v, n+1, PC, p, RPC) + PC[0] * R(t, u, v, n+1, PC, p, RPC);
  // }
  // return val;

  if (t == 0 && u == 0 && v == 0){
    val += pow(-2 * p, n) * boys(n, p * RPC * RPC);
  } else if(t == 0 && u == 0){
    if (v > 1){
      val += (v-1)*R(t, u, v-2, n+1, PCx, PCy, PCz, p, RPC);
    }
    val += PCz * R(t, u, v-1, n+1, PCx, PCy, PCz, p, RPC);
  } else if(t == 0){
    if (u > 1){
      val +=(u-1)*R(t, u-2, v, n+1, PCx, PCy, PCz, p, RPC);
    }
    val += PCy * R(t, u-1, v, n+1, PCx, PCy, PCz, p, RPC);
  } else {
    if (t > 1){
      val += (t-1)*R(t-2, u, v, n+1, PCx, PCy, PCz, p, RPC);
    }
    val += PCx * R(t-1, u, v, n+1, PCx, PCy, PCz, p, RPC);
  }
  return val;
}

double Integrator::nuclear(Eigen::ArrayXd A, std::array<int, 3> lmn1, double a, Eigen::ArrayXd B, std::array<int, 3> lmn2, double b, Eigen::ArrayXd C){
  double p = a+b;
  Eigen::Array3d P = GPC(a, A, b, B);
  double RPC = (P-C).matrix().norm();
  const int l1 = lmn1[0]; const int m1 = lmn1[1]; const int n1 = lmn1[2];
  const int l2 = lmn2[0]; const int m2 = lmn2[1]; const int n2 = lmn2[2];
  double val = 0.0;
  for (int t = 0; t < l1 + l2 + 1; t++){
    for (int u = 0; u < m1 + m2 + 1; u++){
      for (int v = 0; v < n1 + n2 + 1; v++){
	val +=
	  E(l1, l2, t, A[0] - B[0], a, b) *
	  E(m1, m2, u, A[1] - B[1], a, b) *
	  E(n1, n2, v, A[2] - B[2], a, b) *
	  R(t, u, v, 0, P[0]-C[0], P[1]-C[1], P[2]-C[2], p, RPC);
      }
    }
  }
  val *= (2 * M_PI / p);
  return val;
}

double Integrator::V(BasisFunction bf1, BasisFunction bf2, Eigen::Array3d C){
  double total = 0.0;
  for (int i = 0; i < bf1.coefs.size(); i++){
    for (int j = 0; j < bf2.coefs.size(); j++){
      total += bf1.norm[i] * bf2.norm[j] * bf1.coefs[i] * bf2.coefs[j] * nuclear(bf1.origin, bf1.shell, bf1.exps[i], bf2.origin, bf2.shell, bf2.exps[j], C);
    }
  }
  return total;
}

Eigen::ArrayXXd Integrator::VMatrix(){
  int aos = system.nCGFs;
  Eigen::ArrayXXd retmat = Eigen::ArrayXXd::Zero(aos, aos);
  for (int i = 0; i < aos; i++){
    for (int j = 0; j < aos; j++){
      for(auto &atom: system.atoms){
	retmat(i, j) += -atom.z_val * V(system.cgbfs[i], system.cgbfs[j], atom.coord);
      }
    }
  }
  return retmat;
}

double Integrator::electron(Eigen::ArrayXd A, std::array<int, 3> lmn1, double a,
			    Eigen::ArrayXd B, std::array<int, 3> lmn2, double b,
			    Eigen::ArrayXd C, std::array<int, 3> lmn3, double c,
			    Eigen::ArrayXd D, std::array<int, 3> lmn4, double d){
  const int l1 = lmn1[0]; const int m1 = lmn1[1]; const int n1 = lmn1[2];
  const int l2 = lmn2[0]; const int m2 = lmn2[1]; const int n2 = lmn2[2];
  const int l3 = lmn3[0]; const int m3 = lmn3[1]; const int n3 = lmn3[2];
  const int l4 = lmn4[0]; const int m4 = lmn4[1]; const int n4 = lmn4[2];
  Eigen::Array3d P = GPC(a, A, b, B);
  Eigen::Array3d Q = GPC(c, C, d, D);
  double RPQ = (P-Q).matrix().norm();
  double p = a+b;
  double q = c+d;
  double alpha = (p*q) / (p+q);
  double val = 0.0;
  for (int t = 0; t < l1 + l2 + 1; t++){
    for (int u = 0; u < m1 + m2 + 1; u++){
      for (int v = 0; v < n1 + n2 + 1; v++){
	for (int td = 0; td < l3 + l4 + 1; td++){
	  for (int ud = 0; ud < m3 + m4 + 1; ud++){
	    for (int vd = 0; vd < n3 + n4 + 1; vd++){
	      val += E(l1, l2, t, A[0] - B[0], a, b) *
		E(m1, m2, u, A[1] - B[1], a, b) *
		E(n1, n2, v, A[2] - B[2], a, b) *
		E(l3, l4, td, C[0] - D[0], c, d) *
		E(m3, m4, ud, C[1] - D[1], c, d) *
		E(n3, n4, vd, C[2] - D[2], c, d) *
		std::pow(-1, td + ud + vd) *
		R(t+td, u+ud, v+vd, 0, P[0] - Q[0], P[1] - Q[1], P[2] - Q[2], alpha, RPQ);
	    }
	  }
	}
      }
    }
  }
  val *= (2 *std::pow(M_PI, 2.5)) / (p*q*sqrt(p+q));
  return val;
}

double Integrator::ERI(BasisFunction bf1, BasisFunction bf2, BasisFunction bf3, BasisFunction bf4){
  double total = 0.0;
  for (int i = 0; i < bf1.coefs.size(); i++){
    for (int j = 0; j < bf2.coefs.size(); j++){
      for (int k = 0; k < bf3.coefs.size(); k++){
	for (int l = 0; l < bf4.coefs.size(); l++){
	  total +=
	    bf1.norm[i] * bf2.norm[j] * bf3.norm[k] * bf4.norm[l] *
	    bf1.coefs[i] * bf2.coefs[j] * bf3.coefs[k] * bf4.coefs[l] *
	    electron(bf1.origin, bf1.shell, bf1.exps[i],
		     bf2.origin, bf2.shell, bf2.exps[j],
		     bf3.origin, bf3.shell, bf3.exps[k],
		     bf4.origin, bf4.shell, bf4.exps[l]);
	}
      }
    }
  }
  return total;
}

Eigen::ArrayXd Integrator::ERIMatrix(){
  int aos = system.nCGFs;
  int size = (aos * (aos+1) * (std::pow(aos, 2) + aos + 2)) / 8;
  Eigen::ArrayXd retmat = Eigen::ArrayXd::Zero(size);
  int lcount = 0;
  for (int i = 0; i < aos; i++){
    for (int j = 0; j <= i; j++){
      int ij = i*(i+1)/2 + j;
      for (int k = 0; k < aos; k++){
  	for (int l = 0; l <= k; l++){
  	  int kl = k * (k+1)/2+l;

  	  if (ij >= kl){
  	    int ijkl = te_index(i, j, k, l);
  	    // int jikl = te_index(j, i, k, l);
  	    // int ijlk = te_index(i, j, l, k);
  	    // int jilk = te_index(j, i, l, k);
  	    // int klij = te_index(k, l, i, j);
  	    // int lkij = te_index(l, k, i, j);
  	    // int klji = te_index(k, l, j, i);
  	    // int lkji = te_index(l, k, j, i);
	    // This set up pretty much gets you what the crawdad project gives you
	    double val = ERI(system.cgbfs[i],
			     system.cgbfs[j],
			     system.cgbfs[k],
			     system.cgbfs[l]);
	    std::cout << i+1 << " " << j+1 << " " << k+1 << " " << l+1 << " " << ijkl << " "  << val << std::endl;

	    lcount++;
  	    // retmat(jikl) = val;
  	    // retmat(ijlk) = val;
  	    // retmat(jilk) = val;
  	    // retmat(klij) = val;
  	    // retmat(lkij) = val;
  	    // retmat(klji) = val;
  	    // retmat(lkji) = val;
  	  }
  	}
      }
    }
  }
  std::cout << lcount << std::endl;
  return retmat;
}
