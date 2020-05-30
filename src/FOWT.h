#pragma once

#include <string>
#include <armadillo>
#include "Floater.h"
#include "RNA.h"
#include "ENVIR.h"


using namespace arma; // For armadillo classes

class FOWT
{
private:
	/*
	Physical members of the FOWT (floater, rotor nacelle assembly (RNA), tower, mooring lines, ...)
	*/
	Floater m_floater;
	//Tower m_tower;
	RNA m_rna;
	mat::fixed<6,6> m_extLinStiff;
	vec::fixed<6> m_extConstForce;

	/*
	Specification of the analysis
	*/
	int m_hydroMode{ 0 };
	int m_aeroMode{ 0 };
	int m_moorMode{ 0 };

	// Flags to specify the active degrees of freedom
	std::array<bool, 6> m_dofs = { 1, 1, 1, 1, 1, 1 };

	/* 
	FOWT properties derived from its subsystems
	*/
	double m_mass;
	vec::fixed<3> m_CoG;

	/*
	FOWT condition
	*/
	vec::fixed<6> m_disp; // m_disp(0:2) = Position with respect to the initial CoG (i.e. CoG(t) - CoG(0)) --- m_disp(3:5) = Rotation with respect to initial configuration. For now, we are considering small rotations
	vec::fixed<6> m_vel;
	vec::fixed<6> m_acc;

	// Axis system that follows the mean and slow drift.
	// They are evaluated by filtering the instantaneous position with the following parameters.
	double m_filterSD_omega;
	double m_filterSD_zeta;
	vec::fixed<6> m_disp_sd;
	vec::fixed<6> m_vel_sd;
	vec::fixed<6> m_acc_sd;

public:
	FOWT();


	/*****************************************************
		Setters
	*****************************************************/
	void setHydroMode(const int hydroKinMode);
	void setAeroMode(const int aeroMode);
	void setMoorMode(const int moorMode);
	void setDoFs(std::array<bool, 6> &dofs);

	void setExtLinStiff(const mat::fixed<6,6> &extLinStiff);
	void setExtConstForce(const vec::fixed<6> &extConstForce);

	void setFilderSD(const double omega, const double zeta);

	void setFloater(Floater &floater);
	void setRNA(RNA &rna);

	/*****************************************************
		Getters
	*****************************************************/
	int hydroMode() const;
	int aeroMode() const;
	int moorMode() const;

	double filterSD_omega() const;
	double filterSD_zeta() const;

	vec::fixed<3> CoG();
	double mass();

	vec::fixed<6> disp() const;
	vec::fixed<6> vel() const;
	vec::fixed<6> acc() const;
 
	vec::fixed<6> constForce() const;
	std::string printLinStiff() const;
	std::string printFloater() const;
	std::string printRNA() const;
	std::string printHydroMode() const;
	std::string printAeroMode() const;
	std::string printMoorMode() const;
	std::string printDoF() const;

	/*****************************************************
		To add to the string line that is written to the output file at each time step
	*****************************************************/
	void print2outLine() const;

	/*****************************************************
		Forces, acceleration, displacement, etc
	*****************************************************/
	vec::fixed<6> calcAcceleration(const ENVIR &envir);
	void update(const vec::fixed<6> &disp, const vec::fixed<6> &vel, const vec::fixed<6> &acc, const double dt);

	vec::fixed<6> hydrodynamicForce(const ENVIR &envir);
	vec::fixed<6> hydrostaticForce(const ENVIR &envir);
	vec::fixed<6> aeroForce(const ENVIR &envir);
	vec::fixed<6> mooringForce();
	vec::fixed<6> weightForce(const double gravity);
	vec::fixed<6> totalForce(const ENVIR &envir);
};