#include "MISolvers.hpp"
#include "utils.hpp"
#include <cmath>
#include <algorithm>

THO::THO(Molecule input_molecule){

}

float THO::binomial_factor(float j, float l, float m, float PA, float PB){
    for (float i = min(j, l); i < max(0, j-m); i++){
	binomial(l, k) * binomial(m, j-k) * std::pow(PA, l-k) * std::pow(PB, m+k-j);
    }
}
