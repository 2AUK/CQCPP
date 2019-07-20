#include <iostream>
#include <iomanip>
#include <string>
#include "Molecule.hpp"
#include "Integrator.hpp"
#include <eigen3/Eigen/Core>

int main(){
  std::string input_file = "/home/abdullah/Code/C++/SCF/test/H2.xyz";
  std::string input_basis = "/home/abdullah/Code/C++/SCF/basis_sets/sto3g.dat";
  Molecule water(input_file, input_basis);
  std::cout << water.molecule_name << '\n';
  for (const auto atom: water.atoms){
    std::cout << "Atomic Number:\t" << atom.z_val << "\nCoordinates:\n" << atom.coord << '\n';
    for (const auto bf: atom.basisfunctions){
      int l = bf.shell[0]; int m = bf.shell[1]; int n = bf.shell[2];
      std::cout << "Basisfunction Shell:\t" << l << m << n << "\nExponents:\n" << bf.exps << "\nCoefficients:\n" << bf.coefs << "\nOrigin:\n" << bf.origin << "\nNorm:\n"<< std::setprecision(10)  << bf.norm << '\n';
    }
  }
  std::cout << "Number of GTOS:\t" << water.nGTOs << std::endl;
  std::cout << "Number of CGFs:\t" << water.nCGFs << std::endl;
}
