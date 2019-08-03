#include <iostream>
#include <iomanip>
#include <string>
#include "Molecule.hpp"
#include "Integrator.hpp"
#include "utils.hpp"
#include <eigen3/Eigen/Core>
#include "HF.hpp"

int main(){
  Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "[", "]");
  std::string input_file = "/home/abdullah/Code/C++/SCF/examples/h20.dat";
  std::string input_basis = "/home/abdullah/Code/C++/SCF/basis_sets/631g.dat";
  Molecule water(input_file, input_basis);
  // std::cout << water.molecule_name << '\n';
  // for (const auto atom: water.atoms){
  //   std::cout << atom.z_val << "\t"  << atom.coord.format(CommaInitFmt) << '\n';
  //   for (const auto bf: atom.basisfunctions){
  //     int l = bf.shell[0]; int m = bf.shell[1]; int n = bf.shell[2];
  //     std::cout << "[" << l << ", " << m << ", "  << n << "]\t" << std::setprecision(10) << bf.exps.format(CommaInitFmt) << '\t' << bf.coefs.format(CommaInitFmt) << '\t' << bf.origin.format(CommaInitFmt) << '\t' << bf.norm.format(CommaInitFmt) << '\n';
  //   }
  // }
  HF hf(water);
  hf.initialise();
  hf.run(1e-12, 1e-12);
  std::cout << "Number of GTOS:\t" << water.nGTOs << std::endl;
  std::cout << "Number of CGFs:\t" << water.nCGFs << std::endl;
}
