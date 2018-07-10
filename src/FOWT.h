#pragma once

#include <string>
#include <armadillo>
#include "Floater.h"


using namespace arma; // For armadillo classes

class FOWT
{
private:
	Floater m_floater;
	//Tower m_tower;
	//Rotor m_rotor;
	//Nacelle m_nacelle;
	vec::fixed<3> m_linStiff;
	
	vec::fixed<3> m_pos;
	vec::fixed<3> m_vel;
	vec::fixed<3> m_acc;

public:
	FOWT() : m_linStiff(fill::zeros)
	{}


	/*****************************************************
		Setters
	*****************************************************/
	void readLinStiff(const std::string &data);
	void setFloater(const Floater &floater);

	

	/*****************************************************
		Getters
	*****************************************************/
	std::string printLinStiff() const;
	std::string printFloater() const;


	/*****************************************************
		NEED TO NAME THIS SECTION
	*****************************************************/
	//vec hydroForce(const ENVIR &envir);
	////vec aeroForce(const ENVIR &envir);
	////vec mooringForce(const ENVIR &envir);
	//// vec weightForce(const ENVIR &envir);
	//vec totalForce(const ENVIR &envir);
	//vec calcAcc(); // Passar um flag como argumento pra dizer qual tipo de equação de movimento vai ser usada. Linear ou formulação completa

	//void setVel();
	//void setPos();

	//vec CoG(); // Calculate CoG of the aggregation (floater + tower + ....)
};

