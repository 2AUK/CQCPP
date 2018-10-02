#include <iostream>
#include "Atom.h"
#include "BasisFunc.h"

int main(){
    Eigen::Array3f coordinates;
    coordinates << 0, 0, 0;
    Atom hydrogen(1, coordinates);
    std::cout << coordinates << '\n';
    hydrogen.populate_basis("/home/abdullah/Code/C++/SCF/basis_set/sto3g.dat");
    for (const auto i: hydrogen.basisfunctions){
        std::cout << i.exps << i.coefs << i.origin << '\n';
    }
}
