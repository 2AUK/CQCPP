#pragma once
#include <iostream>
#include <tuple>
#include <eigen3/Eigen/Core>

class BasisFunction
{
    public:
        std::tuple<int, int, int> shell; //Cartesian Angular Momentum Shell
        Eigen::ArrayXf norm; //Normalisation constant for orbital - calculated using normalize() function - with values for each primitive gaussian function
        Eigen::ArrayXf origin; //Orbital origin (essentially coordinates of atom)
        Eigen::ArrayXf exps; //Exponents read from basis-set file
        Eigen::ArrayXf coefs; //Coefficients read from basis-set file
        BasisFunction(std::tuple<int, int, int>, Eigen::ArrayXf, Eigen::ArrayXf, Eigen::ArrayXf); //Constructor
        void normalize(); //Normalization function - this sets norm to a value
};
