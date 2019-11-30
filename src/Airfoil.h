#pragma once

#include <armadillo>
#include <vector>
#include "auxFunctions.h"

using namespace arma; // For armadillo classes

class Airfoil
{
private:
    std::vector<double> m_angle;
    std::vector<double> m_CL;
    std::vector<double> m_CD;
    std::vector<double> m_CM;

	// Splines for interpolation
	tk::spline m_spl_CL;
	tk::spline m_spl_CD;
	tk::spline m_spl_CM;

public:
	Airfoil();

	/*****************************************************
		Setters
	*****************************************************/
	void addAirfoilData(double angle, double CL, double CD, double CM);

	/*****************************************************
		Getters
	*****************************************************/
	unsigned int size() const;
	double angle(unsigned int index) const;
	double getCL(unsigned int index) const;
	double getCD(unsigned int index) const;
	double getCM(unsigned int index) const;

	/*****************************************************
		Getters based on the angle of attack (in degrees), using
		cubic spline interpolation
	*****************************************************/
	double CL(double angle) const;
	double CD(double angle) const;
	double CM(double angle) const;
};
