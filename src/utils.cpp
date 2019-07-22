#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <eigen3/Eigen/Core>
#include <gsl/gsl_sf_hyperg.h>

std::vector<std::string> split(const std::string& s, char delimiter){
  std::vector<std::string> tokens;
  std::string token;
  std::istringstream tokenStream(s);
  while (std::getline(tokenStream, token, delimiter)){
    if (token != "")
      tokens.push_back(token);
  }
  return tokens;
}

float factorial(float n){
  if (n <= 1){
    return 1;
  } else {
    return n*factorial(n-1);
  }
}

float double_factorial(float n){
  if (n <= 1){
    return 1;
  } else {
    return n*double_factorial(n-2);
  }
}

float binomial(float n, float k){
  if (n == k){
    return 1;
  } else {
    return factorial(n) / (factorial(k) * factorial(n - k));
  }
}

Eigen::ArrayXd GPC(double a, Eigen::ArrayXd A, double b, Eigen::ArrayXd B){
  return (a * A + b * B) / (a+b);
}

double boys(double n, double T){
  return gsl_sf_hyperg_1F1(n+0.5, n+1.5, -T) / (2.0*n + 1.0);
}
