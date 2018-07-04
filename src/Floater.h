#pragma once

#include "MorisonElement.h"

#include <string>
#include <vector> // For std::vector
#include <functional> // For std::reference_wrapper
#include <armadillo> // Linear algebra library with usage similar to MATLAB

using namespace arma; // For armadillo classes

class Floater
{
private:
	double m_mass;
	vec m_CoG; // Coordinates of center of gravity. It is a 3x1 array, but since armadillo matrices do not have compile-time sizes, size will be defined latter.
	vec m_inertia; // Moments and products of inertia. It is a 6x1 array. Actually, it is a symmetric 3x3 matrix, hence 3 elements are simply repeated.
	std::vector< std::reference_wrapper<MorisonElement> > m_morisonElements;	
	std::vector< vec > m_nodes;

public:
	Floater() : m_CoG(3)
	{}


	/*****************************************************
		Setters
	*****************************************************/
	void readMass(const std::string &data);
	void readCoG(const std::string &data);


	/*****************************************************
		Getters
	*****************************************************/
	std::string printMass() const;
	std::string printCoG() const;


	//mat rotatMat(const vec &FOWTpos);
	//vec hydrodynamicForce(const ENVIR &envir, const vec &FOWTpos, const vec &FOWTvel, const vec &FOWTacc);
};

