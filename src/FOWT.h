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
	vec::fixed<3> m_extLinStiff;
	vec::fixed<6> m_extConstForce;

	/*
	Specification of the analysis
	*/
	int m_hydroMode{ 0 };
	int m_hydroPosMode{ 0 };
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

public:
	FOWT();


	/*****************************************************
		Setters
	*****************************************************/
	void setHydroMode(const int hydroKinMode);
	void setHydroPosMode(const int hydroPosMode);
	void setAeroMode(const int aeroMode);
	void setMoorMode(const int moorMode);
	void setDoFs(std::array<bool, 6> &dofs);

	void setExtLinStiff(const vec::fixed<3> &extLinStiff);
	void setExtConstForce(const vec::fixed<6> &extConstForce);

	void setFloater(Floater &floater);
	void setRNA(RNA &rna);

	/*****************************************************
		Getters
	*****************************************************/
	int hydroMode() const;
	int hydroPosMode() const;
	int aeroMode() const;
	int moorMode() const;

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
	std::string printHydroPosMode() const;
	std::string printAeroMode() const;
	std::string printMoorMode() const;
	std::string printDoF() const;

	/*****************************************************
		Forces, acceleration, displacement, etc
	*****************************************************/
	vec::fixed<6> calcAcceleration(const ENVIR &envir);
	void update(const vec::fixed<6> &disp, const vec::fixed<6> &vel, const vec::fixed<6> &acc);

	vec::fixed<6> hydrodynamicForce(const ENVIR &envir);
	vec::fixed<6> hydrostaticForce(const ENVIR &envir);
	vec::fixed<6> aeroForce(const ENVIR &envir);
	vec::fixed<6> mooringForce();
	vec::fixed<6> weightForce(const double gravity);
	vec::fixed<6> totalForce(const ENVIR &envir);
};