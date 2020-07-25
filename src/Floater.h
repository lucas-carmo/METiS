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

	vec::fixed<6> m_disp; // Floater displacement (i.e. [instantaneous position] - [initial position])
	vec::fixed<6> m_disp_sd; // Keep track of the slow position of the floater as well

public:
	Floater();

	/*****************************************************
		Overloaded operators
	*****************************************************/
	Floater& operator=(const Floater &floater);

	/*****************************************************
		Setters
	*****************************************************/
	void setMass(const double mass);
	void setInertia(const vec::fixed<6> &inertia);
	void setCoG(const vec::fixed<3> &cog);

	void addMorisonCirc(vec::fixed<3> &node1_coord, vec::fixed<3> &node2_coord, const double diam, const double CD, const double CM, const unsigned int numIntPoints,
	const double botDiam, const double topDiam, const double axialCD, const double axialCa, const bool botPressFlag);
	void addMorisonRect(vec::fixed<3> &node1_coord, vec::fixed<3> &node2_coord, vec::fixed<3> &node3_coord, const double diam_X, const double diam_Y,
	const double CD_X, const double CD_Y, const double CM_X, const double CM_Y, const unsigned int numIntPoints,
	const double botArea, const double topArea, const double axialCD, const double axialCa, const bool botPressFlag);

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
		Forces, acceleration, displacement, etc
	*****************************************************/
	void update(const vec::fixed<6> &FOWTdisp, const vec::fixed<6> &FOWTvel, const vec::fixed<6> &FOWTacc,
				const vec::fixed<6> &FOWTdisp_SD, const vec::fixed<6> &FOWTvel_SD);
	mat::fixed<6, 6> addedMass(const double density) const;
	vec::fixed<6> hydrodynamicForce(const ENVIR &envir, const int hydroMode) const;
	vec::fixed<6> hydrostaticForce(const ENVIR &envir) const;
};

