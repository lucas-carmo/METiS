#pragma once

#include <armadillo> // Linear algebra library with usage similar to MATLAB
#include "ENVIR.h"

using namespace arma;

class MorisonElement
{
protected:
	double m_numIntPoints;

	vec m_cog2node1;
	vec m_cog2node2;
	
	vec m_node1pos;
	vec m_node2pos;

	vec m_node1vel;
	vec m_node2vel;

	vec m_node1acc;
	vec m_node2acc;

	vec m_fluidVelPerp;
	vec m_fluidAccPerp;
	vec m_fluidVelTang;
	vec m_fluidAccTang;

public:
	virtual void updateNodePos(const mat &rotatMatrix, const vec &FOWTpos) = 0;
	virtual void updateNodeVel(const mat &rotatMatrix, const vec &FOWTpos, const vec &FOWTvel) = 0;
	virtual void updateNodeAcc(const mat &rotatMatrix, const vec &FOWTpos, const vec &FOWTvel, const vec &FOWTacc) = 0;

	void decomposeFluidVelAcc(const ENVIR &envir);

	virtual vec hydrostaticForce(const ENVIR &envir) = 0;
	virtual vec morisonForce(const ENVIR &envir) = 0;
	virtual vec heavePlateForce(const ENVIR &envir) = 0;
};

