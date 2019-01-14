#pragma once

#include <armadillo>
#include <vector>

using namespace arma; // For armadillo classes

class Airfoil
{
private:
    std::vector<double> angle;
    std::vector<double> CL;
    std::vector<double> CD;
    std::vector<double> CM;

public:
};

