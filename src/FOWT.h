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
	//Physical members of the FOWT (floater, rotor nacelle assembly (RNA), tower, mooring lines, ...)
	Floater m_floater;
	//Tower m_tower;
	RNA m_rna;
	vec::fixed<3> m_linStiff;

	/*
	Forces included in the analysis
	*/
	int m_hydroMode{ 0 };
	int m_aeroMode{ 0 };
	int m_moorMode{ 0 };

	/*
	Flags to specify the active degrees of freedom
	*/
	std::array<bool, 6> m_dofs = { 1, 1, 1, 1, 1, 1 };

	// Properties derived from the other ones
	double m_mass;
	vec::fixed<3> m_CoG;

	// FOWT condition
	// m_disp(0:2) = Position with respect to the initial CoG (i.e. CoG(t) - CoG(0))
	// m_disp(3:5) = Rotation with respect to initial configuration. For now, we are considering small rotations.
	vec::fixed<6> m_disp;
	vec::fixed<6> m_vel;
	vec::fixed<6> m_acc;

public:
	FOWT();


	/*****************************************************
		Setters
	*****************************************************/
	void readHydroMode(const std::string &data);
	void readAeroMode(const std::string &data);
	void readMoorMode(const std::string &data);
	void readDOFs(const std::string &data);

	void readLinStiff(const std::string &data);
	void setFloater(Floater &floater);
	void setRNA(RNA &rna);

	/*****************************************************
		Getters
	*****************************************************/
	int hydroMode() const;
	int aeroMode() const;
	int moorMode() const;

	vec::fixed<3> CoG();
	double mass();

	vec::fixed<6> disp() const;
	vec::fixed<6> vel() const;
	vec::fixed<6> acc() const;
 
	std::string printLinStiff() const;
	std::string printFloater() const;
	std::string printRNA() const;
	std::string printHydroMode() const;
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