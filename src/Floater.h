#pragma once

#include "MorisonElement.h"
#include "ENVIR.h"

#include <string>
#include <vector> // For std::vector
#include <memory> // For std::unique_ptr
#include <armadillo> // Linear algebra library with usage similar to MATLAB


using namespace arma; // For armadillo classes

class Floater
{
private:
	double m_mass;
	vec::fixed<3> m_CoG; // Coordinates of center of gravity.
	vec::fixed<6> m_inertia; // Moments and products of inertia. It is a 6x1 array. Actually, it is a symmetric 3x3 matrix, hence 3 elements are simply repeated.
	std::vector< std::unique_ptr<MorisonElement> > m_MorisonElements;	

public:
	Floater();


	/*****************************************************
		Setters
	*****************************************************/
	void readMass(const std::string &data);
	void readCoG(const std::string &data);
	void addMorisonCirc(const std::string &data, const ENVIR &envir); // Need envir class for nodes location
	void addMorisonRect(const std::string &data, const ENVIR &envir);


	/*****************************************************
		Getters
	*****************************************************/
	std::string printMass() const;
	std::string printCoG() const;
	std::string printMorisonElements() const;

	/*****************************************************
		Overloaded operators
	*****************************************************/
	Floater& operator= (Floater &floater);

	//mat rotatMat(const vec &FOWTpos);
	//vec hydrodynamicForce(const ENVIR &envir, const vec &FOWTpos, const vec &FOWTvel, const vec &FOWTacc);
};

