#include <iostream>
#include <tuple>
#include <Eigen/Core>

class BasisFunction
{
    public:
        std::tuple<int, int, int> shell;
        float norm;
        Eigen::Array3f origin;
        Eigen::Array3f exps;
        Eigen::Array3f coefs;
        void normalize();
}

