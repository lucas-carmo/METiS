#pragma once

#include <armadillo>
#include <vector>

using namespace arma; // For armadillo classes

class Airfoil
{
private:
    std::vector<double> m_angle;
    std::vector<double> m_CL;
    std::vector<double> m_CD;
    std::vector<double> m_CM;

public:
	Airfoil();
};

