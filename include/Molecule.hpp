#pragma once
#include "Atom.hpp"
#include <string>

class Molecule{

public:
  int natoms; //Number of atoms - read from first line of .xyz file
  std::string molecule_name; //Name of molecule - read from second line of .xyz file
  std::vector<Atom> atoms; //vector of atom objects associated with molecule
  Molecule(std::string, std::string); //Constructor
  std::vector<BasisFunction> cgbfs;
private:
  void read_xyz(std::string inXYZ); //Function to read .xyz file - used in the constructor
};
