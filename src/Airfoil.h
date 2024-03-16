#pragma once

#include <armadillo>
#include <vector>
#include "auxFunctions.h"

using namespace arma; // For armadillo classes

class Airfoil
{
private:
    std::vector<double> m_angle;
    std::vector<realType> m_CL;
    std::vector<realType> m_CD;
    std::vector<realType> m_CM;

	// Splines for interpolation
	tk::spline m_spl_CL;
	tk::spline m_spl_CD;
	tk::spline m_spl_CM;

public:
	Airfoil();

	/*****************************************************
		Setters
	*****************************************************/
	void addAirfoilData(double angle, realType CL, realType CD, realType CM);

	/*****************************************************
		Getters
	*****************************************************/
	unsigned int size() const;
	double angle(unsigned int index) const;
	realType getCL(unsigned int index) const;
	realType getCD(unsigned int index) const;
	realType getCM(unsigned int index) const;

	/*****************************************************
		Getters based on the angle of attack (in degrees), using
		cubic spline interpolation
	*****************************************************/
	realType CL(double angle) const;
	realType CD(double angle) const;
	realType CM(double angle) const;
};
