#pragma once
#include "Atom.hpp"
#include <string>

/**
 *   \file Molecule.hpp
 *   \brief Molecule class encapsulating system information
 *
 *  The Molecule class contains all the relevant information for evaluating the molecular integrals
 *  and calculating energy via Hartree-Fock method.
 *
 */

class Molecule{

public:
  int natoms; /**< Number of atoms in molecule. */
  int nGTOs; /**< Number of primitive gaussians (Gaussian Type Orbitals) in each atom of molecule. */
  int nCGFs; /**< Number of contracted gaussians in each atom of molecule. */
  std::string molecule_name; /**< String containing name of the molecule (.xyz file header). */
  std::vector<Atom> atoms; /**< Vector of Atom objects. */
  Molecule(std::string, std::string); 
  std::vector<BasisFunction> cgbfs; /**< Vector containing BasisFunction objects of all atoms in molecule. */
private:
  void read_xyz(std::string inXYZ); 
};
