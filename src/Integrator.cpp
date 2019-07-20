#include "Integrator.hpp"
#include <math.h>

float E(int i, int j, int t, float Qx, float a, float b){
  float p = a + b;
  float q = (a * b) / p;
  if (i == 0 && j == 0 && t == 0){
    return(exp(-q * Qx * Qx));
  } else if (t < 0 || t > (i + j)){
    return 0.0;
  } else if (i == 0){
    return (1 / (2 * p)) * E(i, j-1, t-1, Qx, a, b) + ((q * Qx) / b) * E(i, j-1, t, Qx, a, b) + (t + 1) * E(i-1, j, t + 1, Qx, a, b);
  } else {
    return (1 / (2 * p)) * E(i-1, j, t-1, Qx, a, b) - ((q * Qx) / a) * E(i-1, j, t, Qx, a, b) + (t + 1) * E(i-1, j, t + 1, Qx, a, b);
  }
}
