#pragma once
#include <string>
#include <vector>
#include <sstream>
#include <eigen3/Eigen/Core>

std::vector<std::string> split(const std::string& s,char delimiter); // Tokenize string based on a char delimiter

float factorial(float n); //Calculates n!

float double_factorial(float n); //Calculates n!!

float binomial(float n, float k); //Calculates binomials

int te_index(int i, int j, int k, int l);

Eigen::Array3d GPC(double, Eigen::Array3d, double, Eigen::Array3d);

double boys(double, double);

double norm(Eigen::Array3d, Eigen::Array3d);

int te_index(int, int, int, int);

void swap(int&, int&);
