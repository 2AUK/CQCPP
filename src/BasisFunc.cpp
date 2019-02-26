#include "BasisFunc.hpp"
#include "utils.hpp"
#include <cmath>
#include <eigen3/Eigen/Core>

BasisFunction::BasisFunction(std::tuple<int, int, int> a, 
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
  float l, m, n;
  std::tie(l,m,n) = shell;
  Eigen::ArrayXf term1 = 2*exps/M_PI;
  Eigen::ArrayXf term2 = 4*exps;
  float term3 = std::pow(double_factorial(2*l-1) * double_factorial(2*m-1) * double_factorial(2*n-1), 1/2);
  Eigen::ArrayXf term12 = term1.pow(3.0/4.0) * term2.pow((l+m+n)/2.0);
  norm = term12 / term3;
}
