#pragma once

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
	

	//Acho que essas coisas tem que ser calculadas dentro de uma funcao for�a, e nao como membros da classe
	//vec m_node1pos;
	//vec m_node2pos;

	//vec m_node1vel;
	//vec m_node2vel;

	//vec m_node1acc;
	//vec m_node2acc;

	//vec m_fluidVelPerp;
	//vec m_fluidAccPerp;
	//vec m_fluidVelTang;
	//vec m_fluidAccTang;

public:
	MorisonElement(vec cog2node1, vec cog2node2, int numIntPoints, 
				   bool botPressFlag, double axialCD, double axialCa);

	//virtual void updateNodePos(const mat &rotatMatrix, const vec &FOWTpos) = 0;
	//virtual void updateNodeVel(const mat &rotatMatrix, const vec &FOWTpos, const vec &FOWTvel) = 0;
	//virtual void updateNodeAcc(const mat &rotatMatrix, const vec &FOWTpos, const vec &FOWTvel, const vec &FOWTacc) = 0;

	//void decomposeFluidVelAcc(const ENVIR &envir);

	virtual vec::fixed<6> hydrostaticForce(const ENVIR &envir) = 0;
	virtual vec::fixed<6> hydrodynamicForce(const ENVIR &envir) = 0;

	virtual std::string print() const = 0;
};

