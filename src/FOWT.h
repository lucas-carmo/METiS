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

	// Properties derived from the other ones
	double m_mass;
	vec::fixed<3> m_CoG;

	// FOWT condition
	// m_pos(0:2) = Position with respect to the initial CoG (i.e. CoG(t) - CoG(0))
	// m_pos(3:5) = Rotation with respect to initial configuration. For now, we are considering small rotations.
	vec::fixed<6> m_pos;
	vec::fixed<6> m_vel;
	vec::fixed<6> m_acc;

public:
	FOWT();

	/*****************************************************
		Overloaded operators
	*****************************************************/
	FOWT& operator=(const FOWT &fowt);
	

	/*****************************************************
		Setters
	*****************************************************/
	void readLinStiff(const std::string &data);
	void setFloater(Floater &floater);
	void setRNA(RNA &rna);

	/*****************************************************
		Getters
	*****************************************************/
	vec::fixed<3> CoG();
	double mass();

	vec::fixed<6> pos() const;
	vec::fixed<6> vel() const;
	vec::fixed<6> acc() const;
 
	std::string printLinStiff() const;
	std::string printFloater() const;
	std::string printRNA() const;

	/*****************************************************
		Forces, acceleration, position, etc
	*****************************************************/
	vec::fixed<6> calcAcceleration(const ENVIR &envir);
	void update(const vec::fixed<6> &pos, const vec::fixed<6> &vel, const vec::fixed<6> &acc);

	vec::fixed<6> hydrodynamicForce(const ENVIR &envir);
	vec::fixed<6> hydrostaticForce(const double watDensity, const double gravity);
	////vec aeroForce(const ENVIR &envir);
	vec::fixed<6> mooringForce();
	vec::fixed<6> weightForce(const double gravity);
	vec::fixed<6> totalForce(const ENVIR &envir);
};

