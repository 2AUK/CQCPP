#include "HF.hpp"

HF::HF(Molecule& input_molecule) : system(input_molecule) {
  std::cout << "System loaded in to HF object" << std::endl;
}

void HF::initialise(){
  Integrator isys(system);
  S = isys.SMatrix();
  T = isys.TMatrix();
  V = isys.VMatrix();
  H = T + V;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S.matrix());
  O = es.operatorInverseSqrt();
  F = O.matrix().transpose() * H.matrix() * O.matrix();
}
