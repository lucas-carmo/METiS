#pragma once

#include <string>
#include <armadillo>

using namespace arma; // For armadillo classes

class Wave{
private:
	double m_height;
	double m_period;
	double m_direction; // In degrees
	double m_phase; // In degrees

	// Since the wave number and length are directly related to the wave period by the
	// dispersion relation, having them as members is a violation of the DRY principle.
	// However, as calculating them every time step would be a large waste of time, I
	// decided to keep them and to avoid any setter, in such a way that it is necessary
	// to use the constructors to set the parameters.
	//
	// In any case, they receive NaN in the default constructor. If the user tries to call
	// ENVIR.waveNumber() or ENVIR.length() and they are NaN, an exception is thrown.
	// In this case it is still safe to call ENVIR.waveNumber(watDepth, gravity).
	double m_waveNumber;
	double m_length;

	// Absolute tolerance for the numerical solution of the dispersion equation
	static double constexpr m_epsWave = 1e-7;

	// Some quantities are computed several times in the evaluation of pressure, velocity, etc. 
	// I decided to calculate them only once and store them here.
	double m_omega_x_A;
	double m_cosBeta;
	double m_sinBeta;

public:
	/*****************************************************
		Constructors
	****************z************************************/
	Wave();
	Wave(const std::string &waveType, const double height, const double freqORperiod, const double direction, const double phase, const double watDepth, const double gravity);


	/*****************************************************
		Getters
	*****************************************************/
	double height() const;
	double amp() const;
	double period() const;
	double direction() const;
	double phase() const;
	double freq() const;
	double angFreq() const;

	double waveNumber() const;
	double length() const;
	double waveNumber(const double watDepth, const double gravity) const;
	double length(const double watDepth, const double gravity) const;

	double omega_x_A() const;
	double cosBeta() const;
	double sinBeta() const;
};
