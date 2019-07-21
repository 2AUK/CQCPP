#include "Integrator.hpp"
#include "Molecule.hpp"
#include <assert.h>
#include <eigen3/Eigen/Core>

void test_overlap_hydrogen_molecule_sto3g(){
  Eigen::Array22d actual_values;
  actual_values <<  0.9999999909, 0.6593182001, 0.6593182001, 0.9999999909;
  Molecule water("../examples/H2.xyz", "../basis_sets/sto3g.dat");
  Integrator water_system(water);
  Eigen::ArrayXXd computed_values = water_system.SMatrix();
}
