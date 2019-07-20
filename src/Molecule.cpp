#include "Molecule.hpp"
#include "utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>

Molecule::Molecule(std::string inXYZ, std::string basis){
  read_xyz(inXYZ);
  int gtos = 0;
  int cgfs = 0;
  for (auto &atom: atoms){
    atom.populate_basis(basis);
    for (auto &bf: atom.basisfunctions){
      gtos += bf.coefs.size();
      cgfs += 1;
      bf.normalize();
      cgbfs.push_back(bf);
    }
  }
  nGTOs = gtos;
  nCGFs = cgfs;
}

void Molecule::read_xyz(std::string inXYZ){
  /* Currently the format specified is very controlled:
   * No of atoms
   * Molecule Name
   * Z_number X Y Z
   * Z_number_2 X_2 Y_2 Z_2
   * etc
   * Meaning that the first two lines are always header lines and everything is else atom information
   */
  std::ifstream mol;
  std::stringstream buffer;
  mol.open(inXYZ);
  if (mol.is_open()){
    buffer << mol.rdbuf();
    std::vector<std::string> mollines = split(buffer.str(), '\n');
    natoms = std::stoi(mollines[0]);
    molecule_name = mollines[1];        
    Eigen::Array3f coordinates;
    for (auto it = mollines.begin()+2; it != mollines.end(); it++){
      coordinates = Eigen::Array3f::Zero();
      std::vector<std::string> tok = split(*it, ' ');
      coordinates << std::stof(tok[1]), std::stof(tok[2]), std::stod(tok[3]);
      atoms.push_back(Atom(std::stoi(tok[0]), coordinates));
    }
  }
}
