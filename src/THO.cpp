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
		   Eigen::ArrayXf a,
		   std::tuple<int, int, int> lmn2,
		   Eigen::ArrayXf B,
		   Eigen::ArrayXf b)
{
  float l1, m1, n1, l2, m2, n2, S_x, S_y, S_z;
  std::tie(l1, m1, n1) = lmn1;
  std::tie(l2, m2, n2) = lmn2;
  Eigen::ArrayXf gamma = a + b;
}

Eigen::ArrayXXf THO::one_electron_integral_matrix(){}

float THO::two_electron_integral(std::string int_type){}

Eigen::ArrayXXf THO::two_electron_integral_matrix(){}

