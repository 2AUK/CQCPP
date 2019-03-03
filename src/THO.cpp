#include "MISolvers.hpp"
#include "utils.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

THO::THO(Molecule& input_molecule) : system(input_molecule){
  std::cout << integral("overlap") << std::endl;
  std::cout << integral("kinetic") << std::endl;
}

float THO::f(float j, float l, float m, float PA, float PB){
  float total = 0;
  for (float k = std::max(static_cast<float>(0.0), j-m); k < std::min(j, l)+1; k++){
    total += binomial(l, k) * binomial(m, j-k) * std::pow(PA, l-k) * std::pow(PB, m+k-j);
  }
  return total;
}

Eigen::ArrayXXf THO::integral(std::string int_type){
  int aos = static_cast<int>(system.cgbfs.size());
  if (int_type == "overlap"){
    Eigen::ArrayXXf S_mat = Eigen::ArrayXXf(aos, aos);
    for(int i = 0; i < aos; i++){
      for(int j = 0; j < aos; j++){
	S_mat(i,j) = S(system.cgbfs[i], system.cgbfs[j]);
      }
    }
    return S_mat;
  }else if(int_type == "kinetic"){
    Eigen::ArrayXXf T_mat = Eigen::ArrayXXf(aos, aos);
    for(int i = 0; i < aos; i++){
      for(int j = 0; j < aos; j++){
	T_mat(i,j) = T(system.cgbfs[i], system.cgbfs[j]);
      }
    }
  }
}

float THO::S(BasisFunction a, BasisFunction b){
  float total = 0;
  for (int i = 0; i < a.coefs.size(); i++){
    for (int j = 0; j < b.coefs.size(); j++){
      total += a.norm[i] * b.norm[j] * a.coefs[i] * b.coefs[j] * overlap(a.shell, a.origin, a.exps[i], b.shell, b.origin, b.exps[j]);
    }
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

float THO::kinetic(std::tuple<int, int, int> lmn1,
		   Eigen::ArrayXf A,
		   float a,
		   std::tuple<int, int, int> lmn2,
		   Eigen::ArrayXf B,
		   float b)
{
  float l1, m1, n1, l2, m2, n2;
  std::tie(l1, m1, n1) = lmn1;
  std::tie(l2, m2, n2) = lmn2;
  return (b*(2*(l2 + m2 + n2)+3)*overlap(lmn1, A, a, lmn2, B, b))
    - 2*b*b*(overlap(lmn1, A, a, std::make_tuple(l2+2, m2, n2), B, b)
	     + overlap(lmn1, A, a, std::make_tuple(l2, m2+2, n2), B, b)
	     + overlap(lmn1, A, a, std::make_tuple(l2, m2, n2+2), B, b))
    - (1/2*(l2*(l2-1)*overlap(lmn1, A, a, std::make_tuple(l2-2, m2, n2), B, b))
       + (m2*(m2-1)*overlap(lmn1, A, a, std::make_tuple(l2, m2-2, n2), B, b))
       + (n2*(n2-2)*overlap(lmn1, A, a, std::make_tuple(l2, m2, n2-2), B, b)));
}

float THO::T(BasisFunction a, BasisFunction b){
  float total = 0;
  for (int i = 0; i < a.coefs.size(); i++){
    for (int j = 0; j < b.coefs.size(); j++){
      total += a.norm[i] * b.norm[j] * a.coefs[i] * b.coefs[j] * kinetic(a.shell, a.origin, a.exps[i], b.shell, b.origin, b.exps[j]);
    }
  }
}
