#pragma once
#include "Molecule.hpp"

enum Integral_Type{
    overlap,
    kinetic,
    nuclear,
    electron
};

class IntegralSolver{
    public:
        virtual float one_electron_integral(Integral_Type) = 0;
        virtual float two_electron_integral(Integral_Type) = 0;
        virtual Eigen::ArrayXXf one_electron_integral_matrix() = 0;
        virtual Eigen::ArrayXXf two_electron_integral_matrix() = 0;
};

class MMD : public IntegralSolver{
    public:
        MMD(Molecule);
        virtual float one_electron_integral(Integral_Type);
        virtual float two_electron_integral(Integral_Type);
        virtual Eigen::ArrayXXf one_electron_integral_matrix();
        virtual Eigen::ArrayXXf two_electron_integral_matrix();
    private:
        float E();
        float R();
};