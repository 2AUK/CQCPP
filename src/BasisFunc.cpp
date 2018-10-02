#include "BasisFunc.h"
#include <cmath>
#include <eigen3/Eigen/Core>

BasisFunction::BasisFunction(std::tuple<int, int, int> a, 
                    Eigen::ArrayXf b,
                    Eigen::ArrayXf c,
                    Eigen::ArrayXf d)
{
    shell = a;
    origin = b;
    exps = c;
    coefs = d;
}

void BasisFunction::normalize()
{
    //Implement Double Factorials!!
    int l, m, n;
    std::tie(l,m,n) = shell;
}
