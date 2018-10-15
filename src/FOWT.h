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
	

	// Acho que essas coisas aqui nao devem ser membros, e sim funcoes.
	// Ademais, acho que eh melhor da=las com relacao ao centro de gravidade do flutuador, e nao do sistema todo.	
	vec::fixed<3> m_pos;
	vec::fixed<3> m_vel;
	vec::fixed<3> m_acc;

public:
	FOWT();


	/*****************************************************
		Setters
	*****************************************************/
	void readLinStiff(const std::string &data);

	// You should be aware that the floater that is passed to this function will have
	// its vector of MorisonElements (member m_MorisonElements) pointing to nullptr
	// after this function is terminated.
	void setFloater(Floater &floater);

	

	/*****************************************************
		Getters
	*****************************************************/
	std::string printLinStiff() const;
	std::string printFloater() const;


	/*****************************************************
		Forces, acceleration, position, etc
	*****************************************************/
	vec::fixed<6> hydrodynamicForce(const ENVIR &envir) const;
	//vec hydrostaticForce();
	////vec aeroForce(const ENVIR &envir);
	////vec mooringForce(const ENVIR &envir);
	//// vec weightForce(const ENVIR &envir);
	//vec totalForce(const ENVIR &envir);
	//vec calcAcc(); // Passar um flag como argumento pra dizer qual tipo de equação de movimento vai ser usada. Linear ou formulação completa

	//void setVel();
	//void setPos();

	//vec CoG(); // Calculate CoG of the aggregation (floater + tower + ....)
};

