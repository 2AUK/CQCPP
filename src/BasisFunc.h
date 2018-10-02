#pragma once
#include <iostream>
#include <tuple>
#include <eigen3/Eigen/Core>

class BasisFunction
{
    public:
        std::tuple<int, int, int> shell;
        Eigen::ArrayXf norm;
        Eigen::ArrayXf origin;
        Eigen::ArrayXf exps;
        Eigen::ArrayXf coefs;
        BasisFunction(std::tuple<int, int, int>, Eigen::ArrayXf, Eigen::ArrayXf, Eigen::ArrayXf);
        void normalize();
};
