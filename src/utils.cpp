#include <iostream>
#include <string>
#include <vector>
#include <sstream>

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
  return factorial(n) / (factorial(k) * factorial(n - k));
}
