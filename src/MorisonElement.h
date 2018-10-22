#pragma once

#include <memory> // For std::unique_ptr
#include <armadillo> // Linear algebra library with usage similar to MATLAB
#include "ENVIR.h"

using namespace arma;

class MorisonElement
{
protected:
	vec::fixed<3> m_cog2node1;
	vec::fixed<3> m_cog2node2;
	int m_numIntPoints;

	bool m_botPressFlag;
	double m_axialCD;
	double m_axialCa;		

public:
	MorisonElement(vec cog2node1, vec cog2node2, int numIntPoints, 
				   bool botPressFlag, double axialCD, double axialCa);

	MorisonElement(std::unique_ptr<MorisonElement> sourceElement);

	virtual vec::fixed<3> node1Pos(const vec::fixed<6> &floaterPos) const = 0;
	virtual vec::fixed<3> node2Pos(const vec::fixed<6> &floaterPos) const = 0;
	virtual vec::fixed<3> node1Vel(const vec::fixed<6> &floaterPos, const vec::fixed<6> &floaterVel) const = 0;
	virtual vec::fixed<3> node2Vel(const vec::fixed<6> &floaterPos, const vec::fixed<6> &floaterVel) const = 0;
	virtual vec::fixed<3> node1Acc(const vec::fixed<6> &floaterPos, const vec::fixed<6> &floaterVel, const vec::fixed<6> &floaterAcc) const = 0;
	virtual vec::fixed<3> node2Acc(const vec::fixed<6> &floaterPos, const vec::fixed<6> &floaterVel, const vec::fixed<6> &floaterAcc) const = 0;		
	//virtual void updateNodeVel(const mat &rotatMatrix, const vec &FOWTpos, const vec &FOWTvel) = 0;
	//virtual void updateNodeAcc(const mat &rotatMatrix, const vec &FOWTpos, const vec &FOWTvel, const vec &FOWTacc) = 0;

	virtual vec::fixed<6> hydrostaticForce(const ENVIR &envir, const vec::fixed<6> &floaterPos) const = 0;
	virtual vec::fixed<6> hydrodynamicForce(const ENVIR &envir, const vec::fixed<6> &floaterPos, const vec::fixed<6> &floaterVel, const vec::fixed<6> &floaterAcc) const = 0;

	virtual std::string print() const = 0;
};

