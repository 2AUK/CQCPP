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
  H = T + V;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S.matrix());
  O = es.operatorInverseSqrt();
  F = O.matrix().transpose() * H.matrix() * O.matrix();
  Nelec = 0;
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
  current_Ee = (P * (H + F)).sum();
  nuc_rep = nuclear_repulsion();
  current_Et =current_Ee + nuc_rep;
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
  int Ncgfs = system.nCGFs;
  Fnew = Eigen::ArrayXXd::Zero(Ncgfs, Ncgfs);
  for(int i = 0; i < Ncgfs; i++){
    for(int j = 0; j < Ncgfs; j++){
      Fnew(i, j) = H(i, j);
      for(int k = 0; k < Ncgfs; k++){
	for(int l = 0; l < Ncgfs; l++){
	  int ijkl = te_index(i, j, k, l);
	  int ikjl = te_index(i, k, j, l);
	  Fnew(i,j) += P(k, l) * (2.0 * ERI(ijkl) - ERI(ikjl));
	}
      }
    }
  }
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esF(Fnew.matrix());
  Eigen::ArrayXXd C = esF.eigenvectors();
  Eigen::ArrayXXd C_0 = O.matrix() * C.matrix();
  Pnew = Eigen::ArrayXXd::Zero(Ncgfs, Ncgfs);
  for (int i = 0; i < Ncgfs; i++){
    for (int j = 0; j < Ncgfs; j++){
      for (int k = 0; k < Nelec/2; k++){
	Pnew(i,j) += C_0(i,k) * C_0(j,k);
      }
    }
  }
  new_Ee = (Pnew * (H + Fnew)).sum();
  new_Et = new_Ee + nuc_rep;

  delta_E = new_Ee - current_Ee;
  RMS = std::sqrt(std::pow((Pnew - P).sum(), 2));

  current_Ee = new_Ee;
  current_Et = new_Et;
  for(int i = 0; i < Ncgfs; i++){
    for(int j = 0; j < Ncgfs; j++){
      P(i,j) = Pnew(i,j);
      F(i,j) = Fnew(i,j);
    }
  }
}

void HF::run(double tol1, double tol2, int steps){
  int i = 0;
  std::cout << "Iter" << std::setw(20) << "E_e" << std::setw(20) << "E_t" << std::setw(20) << "delta_E" << std::setw(20) << "RMS"<< '\n';
  while (i < steps){
    std::cout << i << std::setw(20) << current_Ee << std::setw(20) << current_Et << std::setw(20) << delta_E << std::setw(20) << RMS << '\n';
    step();
    if(std::abs(delta_E) < tol1 || std::abs(RMS) < tol2){
      std::cout << "Converged!\n";
      break;
    }
    i++;
  }
}
