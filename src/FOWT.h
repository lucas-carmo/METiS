#pragma once

#include <string>
#include <armadillo>
#include "Floater.h"
#include "ENVIR.h"


using namespace arma; // For armadillo classes

class FOWT
{
private:
	Floater m_floater;
	//Tower m_tower;
	//Rotor m_rotor;
	//Nacelle m_nacelle;
	vec::fixed<3> m_linStiff;

	// Properties derived from the ones others
	double m_mass;
	vec::fixed<3> m_CoG;

	// FOWT condition (position, velocity, acceleration, etc)
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

	/*****************************************************
		Forces, acceleration, position, etc
	*****************************************************/
	vec::fixed<6> acceleration(const ENVIR &envir);
	void update(const vec::fixed<6> &pos, const vec::fixed<6> &vel, const vec::fixed<6> &acc);

	vec::fixed<6> hydrodynamicForce(const ENVIR &envir);
	vec::fixed<6> hydrostaticForce(const ENVIR &envir);
	////vec aeroForce(const ENVIR &envir);
	vec::fixed<6> mooringForce();
	vec::fixed<6> weightForce(const ENVIR &envir);
	vec::fixed<6> totalForce(const ENVIR &envir);
};

