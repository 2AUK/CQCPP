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
