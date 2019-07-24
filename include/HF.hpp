#pragma once
#include "Molecule.hpp"
#include "Integrator.hpp"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigenvalues>
#include <iostream>

class HF{
public:
  HF(Molecule&);
  Molecule system;
  double nuc_rep;
  double current_energy;
  Eigen::ArrayXXd S;
  Eigen::ArrayXXd T;
  Eigen::ArrayXXd V;
  Eigen::ArrayXXd O;
  Eigen::ArrayXd ERI;
  Eigen::ArrayXXd H;
  Eigen::ArrayXXd Pnew;
  Eigen::ArrayXXd F;
  Eigen::ArrayXXd Fnew;
  void initialise();
  void run();
private:
  void step();
  double nuclear_repulsion();
};
