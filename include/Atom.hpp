#pragma once
#include <iostream>
#include <tuple>
#include <vector>
#include <string>
#include <eigen3/Eigen/Core>
#include "BasisFunc.hpp"
#include "utils.hpp"

/**
 *   \file Atom.hpp
 *   \brief Atom class encapsulating information of components of system.
 *
 *  Existence of this class allows for fine-tuned manipulation of the system.
 *
 */


class Atom
{
public:
  int z_val; /**< Atom number. */
  Eigen::Array3d coord; /**< Cartesian coordinates. */
  std::vector<BasisFunction> basisfunctions; /**< Vector containing BasisFunction objects for specific Atom. */
  Atom(int, Eigen::Array3d);
  void populate_basis(std::string in_basis); 
private:
  std::vector<std::string> read_basis(std::string in_basis);
};
