#pragma once

#include <string>
#include <armadillo>

using namespace arma; // For armadillo classes

class Wave{
private:
    double m_height;
    double m_period;
    double m_direction;

public:
	/*****************************************************
		Constructors
	*****************************************************/
	Wave(double height = 0, double period = 0, double direction = 0);

	Wave(const std::string &wholeWaveLine);


	/*****************************************************
		Getters
	*****************************************************/
	double height() const;
	double period() const;
	double direction() const;
	double freq() const;
	double angFreq() const;
	double waveNumber(ENVIR &envir) const;
	double waveLength(ENVIR &envir) const;

	/*****************************************************
		Other functions
	*****************************************************/
	vec::fixed<3> Wave::fluidVel(ENVIR &envir, vec &point) const;
	vec::fixed<3> fluidAcc(ENVIR &envir, vec &point) const;	
};
