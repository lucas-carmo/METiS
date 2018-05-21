#pragma once

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
	vec m_linStiff;
	
	vec m_pos;
	vec m_vel;
	vec m_acc;	

public:
	FOWT() : m_linStiff(3)
	{}

	/*
	Setters used to set the data read from the input file.
	*/
	void readLinStiff(const std::string &data);

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

