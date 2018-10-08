#pragma once
#include <iostream>
#include <tuple>
#include <vector>
#include <string>
#include <eigen3/Eigen/Core>
#include "BasisFunc.hpp"
#include "utils.hpp"

class Atom
{
    public:
        int z_val; //Atomic number
        Eigen::Array3f coord; //Cartesian coordinates
        std::vector<BasisFunction> basisfunctions; //Vector containing associated basisfunctions of atom
        Atom(int, Eigen::Array3f); //Constructor
        void populate_basis(std::string in_basis); //Function that finds orbital exponent and coefficient values from a Gaussian94 basis-set file. This requires the file to be in order of atomic number with no gaps

    private:
        std::vector<std::string> read_basis(std::string in_basis); //Helper function to initialise the input basis file before processing it with populate_basis
};
