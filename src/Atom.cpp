#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <array>
#include <experimental/array>
#include "Atom.hpp"

/**
 *   \file Atom.cpp
 *   \brief Implementation for Atom
 *
 *   All the implementation details for all Atom routines held in Atom.hpp
 */

/**
 *  \brief Constructor for Atom object
 *
 *  \param a Atom number
 *  \param c Eigen Array containing coordinates of Atom
 *  \return return type
 */
Atom::Atom(int a, Eigen::Array3f c){
  z_val = a;
  coord = c;
}


/**
 *  \brief Tokenises basis-set files
 *
 * Basis-set files are contained in basis_sets directory and are in Gaussian94 format taken from <a href="https://www.basissetexchange.org/">Basis Set Exchange</a>. Requires for file to contain information from the beginning of periodic table.
 *  \param in_basis Name of desired basis-set
 *  \return Vector of strings containing tokenized basis-set information.
 */
std::vector<std::string> Atom::read_basis(std::string in_basis){
  std::ifstream basis;
  std::stringstream buffer;
  std::vector<std::string> atoms;
  basis.open(in_basis);
  if (basis.is_open()){
    buffer << basis.rdbuf();
    atoms = split(buffer.str(), '*');
    atoms.pop_back();
    return atoms;
  } else {
    std::cout << "Could not read basis-set file!\n";
    return atoms;
  }
}

/**
 *  \brief Populates basisfunctions vector in Atom object.
 *
 *  Finds orbital and coefficient values from Gaussian94 basis-set file.
 *
 *  \param in_basis Name of desired basis-set.
 *  \return void
 */
void Atom::populate_basis(std::string in_basis){
  std::vector<std::string> atoms = read_basis(in_basis);
  std::vector<std::string> lines = split(atoms[z_val-1], '\n');
  bool firstLine = true;
  Eigen::ArrayXf expTemp;
  Eigen::ArrayXf expSPTemp;
  Eigen::ArrayXf coefTemp;
  for (auto it = lines.begin(); it != lines.end(); it++){
    int primNum = 0;
    if (firstLine)
      { firstLine = false; continue; }
    std::vector<std::string> tokens  = split(*it, ' ');
    if (tokens[0] == "S"){
      primNum = std::stoi(tokens[1]);
      expTemp = Eigen::ArrayXf::Zero(primNum);
      coefTemp = Eigen::ArrayXf::Zero(primNum);
      for (int i = 1; i <= primNum; i++){
	auto nx = std::next(it, i);
	std::vector<std::string> toks = split(*nx, ' ');
	expTemp[i-1] = std::stof(toks[0]);
	coefTemp[i-1] = std::stof(toks[1]);
      }
      basisfunctions.push_back(BasisFunction(std::experimental::make_array(0, 0, 0), coord, expTemp, coefTemp));
      std::advance(it, primNum);
    } else if (tokens[0] == "SP"){
      primNum = std::stoi(tokens[1]);
      expTemp = Eigen::ArrayXf::Zero(primNum);
      expSPTemp = Eigen::ArrayXf::Zero(primNum);
      coefTemp = Eigen::ArrayXf::Zero(primNum);
      for (int i = 1; i <= primNum; i++){
	auto nx = std::next(it, i);
	std::vector<std::string> toks = split(*nx, ' ');
	expTemp[i-1] = std::stof(toks[0]);
	coefTemp[i-1] = std::stof(toks[1]);
	expSPTemp[i-1] = std::stof(toks[2]);
      }
      basisfunctions.push_back(BasisFunction(std::experimental::make_array(0,0,0), coord, expTemp, coefTemp));
      basisfunctions.push_back(BasisFunction(std::experimental::make_array(1,0,0), coord, expSPTemp, coefTemp));
      basisfunctions.push_back(BasisFunction(std::experimental::make_array(0,1,0), coord, expSPTemp, coefTemp));
      basisfunctions.push_back(BasisFunction(std::experimental::make_array(0,0,1), coord, expSPTemp, coefTemp));
      std::advance(it, primNum);
    } else if (tokens[0] == "P"){
      primNum = std::stoi(tokens[1]);
      expTemp = Eigen::ArrayXf::Zero(primNum);
      coefTemp = Eigen::ArrayXf::Zero(primNum);
      for (int i = 1; i <= primNum; i++){
	auto nx = std::next(it, i);
	std::vector<std::string> toks = split(*nx, ' ');
	expTemp[i-1] = std::stof(toks[0]);
	coefTemp[i-1] = std::stof(toks[1]);
      }
      basisfunctions.push_back(BasisFunction(std::experimental::make_array(1,0,0), coord, expTemp, coefTemp));
      basisfunctions.push_back(BasisFunction(std::experimental::make_array(0,1,0), coord, expTemp, coefTemp));
      basisfunctions.push_back(BasisFunction(std::experimental::make_array(0,0,1), coord, expTemp, coefTemp));
      std::advance(it, primNum);
    } else if (tokens[0] == "D"){
      primNum = std::stoi(tokens[1]);
      expTemp = Eigen::ArrayXf::Zero(primNum);
      coefTemp = Eigen::ArrayXf::Zero(primNum);
      for (int i = 1; i <= primNum; i++){
	auto nx = std::next(it, i);
	std::vector<std::string> toks = split(*nx, ' ');
	expTemp[i-1] = std::stof(toks[0]);
	coefTemp[i-1] = std::stof(toks[1]);
      }
      basisfunctions.push_back(BasisFunction(std::experimental::make_array(2,0,0), coord, expTemp, coefTemp));
      basisfunctions.push_back(BasisFunction(std::experimental::make_array(1,1,0), coord, expTemp, coefTemp));
      basisfunctions.push_back(BasisFunction(std::experimental::make_array(1,0,1), coord, expTemp, coefTemp));
      basisfunctions.push_back(BasisFunction(std::experimental::make_array(0,2,0), coord, expTemp, coefTemp));
      basisfunctions.push_back(BasisFunction(std::experimental::make_array(0,1,1), coord, expTemp, coefTemp));
      basisfunctions.push_back(BasisFunction(std::experimental::make_array(0,0,2), coord, expTemp, coefTemp));
      std::advance(it, primNum);
    } else if (tokens[0] == "F"){
      primNum = std::stoi(tokens[1]);
      expTemp = Eigen::ArrayXf::Zero(primNum);
      coefTemp = Eigen::ArrayXf::Zero(primNum);
      for (int i = 1; i <= primNum; i++){
	auto nx = std::next(it, i);
	std::vector<std::string> toks = split(*nx, ' ');
	expTemp[i-1] = std::stof(toks[0]);
	coefTemp[i-1] = std::stof(toks[1]);
      }
      basisfunctions.push_back(BasisFunction(std::experimental::make_array(3,0,0), coord, expTemp, coefTemp));
      basisfunctions.push_back(BasisFunction(std::experimental::make_array(2,1,0), coord, expTemp, coefTemp));
      basisfunctions.push_back(BasisFunction(std::experimental::make_array(2,0,1), coord, expTemp, coefTemp));
      basisfunctions.push_back(BasisFunction(std::experimental::make_array(1,2,0), coord, expTemp, coefTemp));
      basisfunctions.push_back(BasisFunction(std::experimental::make_array(1,1,1), coord, expTemp, coefTemp));
      basisfunctions.push_back(BasisFunction(std::experimental::make_array(1,0,2), coord, expTemp, coefTemp));
      basisfunctions.push_back(BasisFunction(std::experimental::make_array(0,3,0), coord, expTemp, coefTemp));
      basisfunctions.push_back(BasisFunction(std::experimental::make_array(0,2,1), coord, expTemp, coefTemp));
      basisfunctions.push_back(BasisFunction(std::experimental::make_array(0,1,2), coord, expTemp, coefTemp));
      basisfunctions.push_back(BasisFunction(std::experimental::make_array(0,0,3), coord, expTemp, coefTemp));
      std::advance(it, primNum);
    } else {
      std::cout << "Fatal Error: Could not read basis-set for current atom OR Orbital not implemented!\n";
    }
  }
}
