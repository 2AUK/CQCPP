#include "Molecule.hpp"
#include "utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>

/**
 *   \file Molecule.cpp
 *   \brief Implementation for Molecule member functions
 *
 *   All the implementation details for all Molecule routines held in Molecule.hpp
 */


/**
 *  \brief Constructor for Molecule object
 *
 *  \param inXYZ String containing file path to input .xyz file
 *  \param basis String containing choice of basis-set to use 
 *  \return Molecule object
 */
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


/**
 *  \brief Reads .xyz file
 *
 *  Reads .xyz files of standard format. First two lines are header information, the rest atom data.
 *
 *  \param inXYZ String containing file path to input .xyz file
 *  \return void
 */
void Molecule::read_xyz(std::string inXYZ){
  std::ifstream mol;
  std::stringstream buffer;
  mol.open(inXYZ);
  if (mol.is_open()){
    buffer << mol.rdbuf();
    std::vector<std::string> mollines = split(buffer.str(), '\n');
    natoms = std::stoi(mollines[0]);
    molecule_name = mollines[1];        
    Eigen::Array3d coordinates;
    for (auto it = mollines.begin()+2; it != mollines.end(); it++){
      coordinates = Eigen::Array3d::Zero();
      std::vector<std::string> tok = split(*it, ' ');
      coordinates << std::stod(tok[1]), std::stod(tok[2]), std::stod(tok[3]);
      //coordinates << std::stod(tok[1])/0.52917721092, std::stod(tok[2])/0.52917721092, std::stod(tok[3])/0.52917721092;
      atoms.push_back(Atom(std::stoi(tok[0]), coordinates));
    }
  }
}
