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
	std::vector<std::unique_ptr<MorisonElement>> m_MorisonElements;	

	// Floater condition (position, velocity, acceleration, etc)
	vec::fixed<6> m_pos;
	vec::fixed<6> m_vel;
	vec::fixed<6> m_acc;


public:
	Floater();


	/*****************************************************
		Overloaded operators
	*****************************************************/
	Floater& operator=(const Floater &floater);


	/*****************************************************
		Setters
	*****************************************************/
	void readMass(const std::string &data);
	void readInertia(const std::string &data);
	void readCoG(const std::string &data);
	void addMorisonCirc(const std::string &data, const ENVIR &envir); // Need envir class for nodes location
	void addMorisonRect(const std::string &data, const ENVIR &envir);


	/*****************************************************
		Getters
	*****************************************************/
	vec::fixed<3> CoG() const;
	double mass() const;
	mat::fixed<6,6> inertiaMatrix() const;

	std::string printMass() const;
	std::string printInertia() const;
	std::string printCoG() const;
	std::string printMorisonElements() const;


	/*****************************************************
		Forces, acceleration, position, etc
	*****************************************************/
	void updatePosVelAcc(const vec::fixed<6> &FOWTpos, const vec::fixed<6> &FOWTvel, const vec::fixed<6> &FOWTacc);
	mat::fixed<6, 6> addedMass(const double density) const;
	vec::fixed<6> hydrodynamicForce(const ENVIR &envir) const;
	vec::fixed<6> hydrostaticForce(const double watDensity, const double gravity) const;
};

