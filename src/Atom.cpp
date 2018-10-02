#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "Atom.h"
#include "BasisFunc.h"
#include "utils.h"

Atom::Atom(int a, Eigen::Array3f c)
{
    z_val = a;
    coord = c;
}

std::vector<std::string> Atom::read_basis(std::string in_basis)
{
    std::ifstream basis;
    std::stringstream buffer;
    basis.open(in_basis);
    if (basis.is_open()){
        buffer << basis.rdbuf();
        std::vector<std::string> atoms = split(buffer.str(), '*');
        atoms.pop_back();
        return atoms;
    } else {
        std::cout << "Could not read basis-set file!\n";
    }
}

void Atom::populate_basis(std::string in_basis)
{
    std::vector<std::string> atoms = read_basis(in_basis);
    std::vector<std::string> lines = split(atoms[z_val-1], '\n');
    bool firstLine = true;
    Eigen::Array3f expTemp;
    Eigen::Array3f coefTemp;
    for (auto it = lines.begin(); it != lines.end(); it++){
        int primNum = 0;
        expTemp = Eigen::Array3f::Zero();
        coefTemp = Eigen::Array3f::Zero();
        if (firstLine)
        { firstLine = false; continue; }
        std::vector<std::string> tokens  = split(*it, ' ');
        if (tokens[0] == "S"){
            primNum = std::stoi(tokens[1]);
            for (int i = 1; i <= primNum; i++){
                auto nx = std::next(it, i);
                std::vector<std::string> toks = split(*nx, ' ');
                expTemp[i-1] = std::stof(toks[0]);
                coefTemp[i-1] = std::stof(toks[1]);
            }
            basisfunctions.push_back(BasisFunction(std::make_tuple(0, 0, 0), coord, expTemp, coefTemp));
        }
    }
}
