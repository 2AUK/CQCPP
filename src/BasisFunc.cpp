#include "BasisFunc.hpp"
#include "utils.hpp"
#include <cmath>
#include <eigen3/Eigen/Core>

BasisFunction::BasisFunction(std::array<int, 3> a, 
			     Eigen::ArrayXf b,
			     Eigen::ArrayXf c,
			     Eigen::ArrayXf d)
{
  shell = a;
  origin = b;
  exps = c;
  coefs = d;
}

void BasisFunction::normalize(){
  int l = shell[0]; int m = shell[1]; int n = shell[2];
  float L = l+m+n;
  Eigen::ArrayXf num = std::pow(2.0, 2.0*(l+m+n) + 3.0/2.0) * exps.pow((l+m+n)+3.0/2.0);
  float denom = double_factorial(2*l-1) * double_factorial(2*m-1) * double_factorial(2*n-1) * std::pow(M_PI, 3.0/2.0);
  Eigen::ArrayXf nd = num / denom;
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
