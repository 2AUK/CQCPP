#include "MISolvers.hpp"
#include "utils.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

THO::THO(Molecule& input_molecule) : system(input_molecule){
    std::cout << f(0, 0, 0, 0, 0) << std::endl;
}

float THO::f(float j, float l, float m, float PA, float PB){
    float total = 0;
    for (float k = std::max(static_cast<float>(0.0), j-m); k < std::min(j, l)+1; k++){
	total += binomial(l, k) * binomial(m, j-k) * std::pow(PA, l-k) * std::pow(PB, m+k-j);
    }
    return total;
}

float THO::one_electron_integral(std::string int_type){
    if (int_type == "overlap"){

    }
}

float THO::overlap(std::tuple<int, int, int> lmn1,
		   Eigen::ArrayXf A,
		   float a,
		   std::tuple<int, int, int> lmn2,
		   Eigen::ArrayXf B,
		   float b)
{
  float l1, m1, n1, l2, m2, n2, S_x, S_y, S_z;
  std::tie(l1, m1, n1) = lmn1;
  std::tie(l2, m2, n2) = lmn2;
  float gamma = a + b;
  Eigen::ArrayXf Q = a * b * A * B / gamma;
  float dist2 = (A - B).matrix().squaredNorm();
  S_x = overlap_1d(l1, l2, Q[0] - A[0], Q[0] - B[0], gamma);
  S_y = overlap_1d(m1, m2, Q[1] - A[1], Q[1] - B[1], gamma);
  S_z = overlap_1d(n1, n2, Q[2] - A[2], Q[2] - B[2], gamma);
  return std::exp(-1 * a * b * dist2 / gamma) * S_x * S_y * S_z;
}

float THO::overlap_1d(int l1, int l2, float PAx, float PBx, float gamma){
  float total = 0;
  for (float j = 0; j < std::floor((l1+l2) / 2) + 1; j++){
    total += f(2*j, l1, l2, PAx, PAx) * double_factorial(2*j - 1) / (std::pow(2*gamma, j));
  }
  return std::sqrt(M_PI/gamma) * total;
}

Eigen::ArrayXXf THO::one_electron_integral_matrix(){}

float THO::two_electron_integral(std::string int_type){}

Eigen::ArrayXXf THO::two_electron_integral_matrix(){}

