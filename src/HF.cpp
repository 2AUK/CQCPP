#include "HF.hpp"

HF::HF(Molecule& input_molecule) : system(input_molecule) {
  std::cout << "System loaded in to HF object" << std::endl;
}

void HF::initialise(){
  Integrator isys(system);
  S = isys.SMatrix();
  T = isys.TMatrix();
  V = isys.VMatrix();
  ERI = isys.ERIMatrix();
  std::cout << S << '\n';
  H = T + V;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S.matrix());
  O = es.operatorInverseSqrt();
  F = O.matrix().transpose() * H.matrix() * O.matrix();
  int Nelec = 0;
  int Ncgfs = system.nCGFs;
  for (auto atom: system.atoms){
    Nelec += atom.z_val;
  }
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esF(F.matrix());
  Eigen::ArrayXXd C = esF.eigenvectors();
  Eigen::ArrayXXd C_0 = O.matrix() * C.matrix();
  P = Eigen::ArrayXXd::Zero(Ncgfs, Ncgfs);
  for (int i = 0; i < Ncgfs; i++){
    for (int j = 0; j < Ncgfs; j++){
      for (int k = 0; k < Nelec/2; k++){
	P(i,j) += C_0(i,k) * C_0(j,k);
      }
    }
  }
  double E_e = (P * (H + F)).sum();
  std::cout << E_e << '\n';
  nuc_rep = nuclear_repulsion();
  std::cout << nuc_rep << '\n';
  std::cout << E_e + nuc_rep << std::endl;
}

double HF::nuclear_repulsion(){
  double total = 0;
  for(int i = 0; i < system.natoms; i++){
    for(int j = i+1; j < system.natoms; j++){
      total += system.atoms[i].z_val * system.atoms[j].z_val / (system.atoms[i].coord - system.atoms[j].coord).matrix().norm();
    }
  }
  return total;
}

void HF::step(){

}
