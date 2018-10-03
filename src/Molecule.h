#pragma once
#include "Atom.h"
#include <string>

class Molecule{

    public:
        int natoms;
        std::string molecule_name;
        std::vector<Atom> atoms;
        Molecule(std::string inXYZ);
    private:
        void read_xyz(std::string inXYZ);
};
