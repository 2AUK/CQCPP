#include <iostream>
#include <string>
#include <vector>
#include <cmath>
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

Eigen::Array3d GPC(double a, Eigen::Array3d A, double b, Eigen::Array3d B){
  double p = a+b;
  double x = ((a * A[0]) + (B[0] * b)) / p;
  double y = ((a * A[1]) + (B[1] * b)) / p;
  double z = ((a * A[2]) + (B[2] * b)) / p;

  Eigen::Array3d ret;
  ret << x, y, z;
  return ret;
}

void swap(int &i, int &j){
  int m = i;
  i = j;
  j = m;
}

double boys(double n, double T){
  return gsl_sf_hyperg_1F1(n+0.5, n+1.5, -T) / (2.0*n + 1.0);
}

double norm(Eigen::Array3d A, Eigen::Array3d B){
  return std::sqrt(pow(A[0] - B[0], 2) + pow(A[1] - B[1], 2) + pow(A[2] - B[2], 2));
}

int te_index(int i, int j, int k, int l){
  if (i < j){
    swap(i, j);
  }
  if (k < l){
    swap(k, l);
  }

  int ij = i * (i + 1) / 2 + j;
  int kl = k * (k + 1) / 2 + l;

  if(ij < kl){
    swap(ij, kl);
  }

  return ij * (ij + 1) / 2 + kl;
}


