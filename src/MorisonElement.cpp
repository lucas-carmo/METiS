#include "MorisonElement.h"

/*****************************************************
	Constructors
*****************************************************/

MorisonElement::MorisonElement(vec cog2node1, vec cog2node2, int numIntPoints, 
							   bool botPressFlag, double axialCD, double axialCa)
	: m_cog2node1(cog2node1), m_cog2node2(cog2node2), 
	  m_botPressFlag(botPressFlag), m_axialCD(axialCD), m_axialCa(axialCa)
{
	// Since Simpson's rule is employed for the integration of the forces along the 
	// Morison's element, we need to make sure that the number of integration points is odd
	if (numIntPoints % 2)
	{
		m_numIntPoints = numIntPoints + 1;
	}
	else
	{
		m_numIntPoints = numIntPoints;
	}
}

/*
	Functions for node position / velocity / acceleration
*/

vec::fixed<3> MorisonElement::node1Pos(const vec::fixed<6> &floaterPos) const
{
	// Fazer uma funcao que calcula a matriz de rotacao
	// mat::fixed<3,3> rotatMatrix(const vec::fixed<6> &FOWTpos) const
	double x = floaterPos[0] + m_cog2node1[0];
	double y = floaterPos[1] + m_cog2node1[1];
	double z = floaterPos[2] + m_cog2node1[2];

	vec::fixed<3> pos = { x, y, z };
	return pos;
}

vec::fixed<3> MorisonElement::node2Pos(const vec::fixed<6> &floaterPos) const
{
	double x = floaterPos[0] + m_cog2node2[0];
	double y = floaterPos[1] + m_cog2node2[1];
	double z = floaterPos[2] + m_cog2node2[2];

	vec::fixed<3> pos = { x, y, z };
	return pos;
}

vec::fixed<3> MorisonElement::node1Vel(const vec::fixed<6> &floaterPos, const vec::fixed<6> &floaterVel) const
{
	double x = 0;
	double y = 0;
	double z = 0;

	vec::fixed<3> vel = { x, y, z };
	return vel;
}

vec::fixed<3> MorisonElement::node2Vel(const vec::fixed<6> &floaterPos, const vec::fixed<6> &floaterVel) const
{
	double x = 0;
	double y = 0;
	double z = 0;

	vec::fixed<3> vel = { x, y, z };
	return vel;
}

vec::fixed<3> MorisonElement::node1Acc(const vec::fixed<6> &floaterPos, const vec::fixed<6> &floaterVel, const vec::fixed<6> &floaterAcc) const
{
	double x = 0;
	double y = 0;
	double z = 0;

	vec::fixed<3> acc = { x, y, z };
	return acc;
}

vec::fixed<3> MorisonElement::node2Acc(const vec::fixed<6> &floaterPos, const vec::fixed<6> &floaterVel, const vec::fixed<6> &floaterAcc) const
{
	double x = 0;
	double y = 0;
	double z = 0;

	vec::fixed<3> acc = { x, y, z };
	return acc;
}
