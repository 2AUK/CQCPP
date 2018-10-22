#pragma once
#include "Molecule.hpp"
#include <string>

class IntegralSolver{
    public:
        virtual float one_electron_integral(std::string) = 0;
        virtual float two_electron_integral(std::string) = 0;
        virtual Eigen::ArrayXXf one_electron_integral_matrix() = 0;
        virtual Eigen::ArrayXXf two_electron_integral_matrix() = 0;
};

class MMD : public IntegralSolver{
    public:
        MMD(Molecule);
        virtual float one_electron_integral(std::string);
        virtual float two_electron_integral(std::string);
        virtual Eigen::ArrayXXf one_electron_integral_matrix();
        virtual Eigen::ArrayXXf two_electron_integral_matrix();
    private:
        float E();
        float R();
};