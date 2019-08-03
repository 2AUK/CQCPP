#pragma once
#include "Molecule.hpp"
#include "Integrator.hpp"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigenvalues>
#include <iostream>
#include <iomanip>

class HF{
public:
  HF(Molecule&);
  Molecule system;
  double nuc_rep;
  double current_Ee;
  double current_Et;
  double new_Ee;
  double new_Et;
  double delta_E = 1E20;
  double RMS = 1E20;
  int Nelec;
  Eigen::ArrayXXd S;
  Eigen::ArrayXXd T;
  Eigen::ArrayXXd V;
  Eigen::ArrayXXd O;
  Eigen::ArrayXd ERI;
  Eigen::ArrayXXd H;
  Eigen::ArrayXXd P;
  Eigen::ArrayXXd Pnew;
  Eigen::ArrayXXd F;
  Eigen::ArrayXXd Fnew;
  void initialise();
  void run(double, double, int=100);
private:
  void step();
  double nuclear_repulsion();
};
