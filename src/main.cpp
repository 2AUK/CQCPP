#include <iostream>
#include <string>
#include "Molecule.h"

int main(){
    std::string input_file = "/home/abdullah/Code/C++/SCF/test/h20.dat";
    std::string input_basis = "/home/abdullah/Code/C++/SCF/basis_sets/sto3g.dat";
    Molecule water(input_file, input_basis);
    std::cout << water.molecule_name << '\n';
    int l,m,n;
    for (const auto atom: water.atoms){
        std::cout << "Atomic Number:\t" << atom.z_val << "\nCoordinates:\t" << atom.coord << '\n';
        for (const auto bf: atom.basisfunctions){
            std::tie(l,m,n) = bf.shell;
            std::cout << "Basisfunction Shell:\t" << l << m << n << "\nExponents:\t" << bf.exps << "\nCoefficients:\t" << bf.coefs << "\nOrigin:\t" << bf.origin << "\nNorm:\t" << bf.norm << '\n';
        }
    }
}
