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

	/*****************************************************
		Setters
	*****************************************************/
	void addAirfoilLine(double angle, double CL, double CD, double CM);

	/*****************************************************
		Getters
	*****************************************************/
	unsigned int size() const;
	double angle(unsigned int index) const;
	double CL(unsigned int index) const;
	double CD(unsigned int index) const;
	double CM(unsigned int index) const;
};

