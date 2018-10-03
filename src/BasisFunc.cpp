#include "BasisFunc.h"
#include "utils.h"
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
    norm = std::pow(2*exps/M_PI, 3/4)*(std::pow(4*exps, (l+m+n)/2))/(std::pow(double_factorial(2*l-1)*double_factorial(2*m-1)*double_factorial(2*n-1)), 1/2);
}
