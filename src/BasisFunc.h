#pragma once
#include <iostream>
#include <tuple>
#include <eigen3/Eigen/Core>

class BasisFunction
{
    public:
        std::tuple<int, int, int> shell;
        Eigen::Array3f norm;
        Eigen::Array3f origin;
        Eigen::Array3f exps;
        Eigen::Array3f coefs;
        BasisFunction(std::tuple<int, int, int>, Eigen::Array3f, Eigen::Array3f, Eigen::Array3f);
        void normalize();
};
