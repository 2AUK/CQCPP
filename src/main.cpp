#include <iostream>
#include "Atom.h"
#include "BasisFunc.h"

int main(){
    Eigen::Array3f coordinates;
    coordinates << 0, 0, 0;
    Atom hydrogen(3, coordinates);
    hydrogen.populate_basis("/home/abdullah/Code/C++/SCF/basis_sets/631g.dat");
    for (const auto i: hydrogen.basisfunctions){
        std::cout << "Exponents\n" << i.exps << "\nCoefficients\n" << i.coefs << "\nOrigin\n" <<  i.origin << '\n';
    }
}
