#pragma once

#include <vector> // For std::vector
#include <functional> //fpr std::reference_wrapper
#include <armadillo> // Linear algebra library with usage similar to MATLAB
#include "MorisonElement.h"


using namespace arma; // For armadillo stuff

class Floater
{
private:
	double m_floaterMass;
	vec m_floaterCoG; // Coordinates of center of gravity. It is a 3x1 array, but since armadillo matrices do not have compile-time sizes, size will be defined latter.
	vec m_floaterInertia; // Moments and products of inertia. It is a 6x1 array. Actually, it is a symmetric 3x3 matrix, hence 3 elements are simply repeated.
	std::vector< std::reference_wrapper<MorisonElement> > m_morisonElements;	

public:
	//mat rotatMat(const vec &FOWTpos);
	//vec hydrodynamicForce(const ENVIR &envir, const vec &FOWTpos, const vec &FOWTvel, const vec &FOWTacc);
};

