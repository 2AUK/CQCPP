#pragma once
#include <iostream>
#include <tuple>
#include <vector>
#include <string>
#include <eigen3/Eigen/Core>
#include "BasisFunc.h"

class Atom
{
    public:
        int z_val;
        Eigen::Array3f coord;
        std::vector<BasisFunction> basisfunctions;
        Atom(int, Eigen::Array3f);
        void populate_basis(std::string in_basis);

    private:
        std::vector<std::string> read_basis(std::string in_basis);
};
