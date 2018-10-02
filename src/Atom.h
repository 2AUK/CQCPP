#include <iostream>
#include <tuple>
#include <vector>
#include <Eigen/Core>
#include "BasisFunc.h"

class Atom
{
    public:
        std::tuple<int, int, int> shell;
        Eigen::Array3f coord;
        std::vector<BasisFunction> basisfunctions;
}
