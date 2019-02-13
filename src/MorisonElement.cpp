#include "MorisonElement.h"

using namespace arma;

/*****************************************************
	Constructors
*****************************************************/

MorisonElement::MorisonElement(const vec &node1Pos, const vec &node2Pos, const vec &cog, const int numIntPoints,
							   const bool botPressFlag, const double axialCD, const double axialCa)
	: m_node1Pos(node1Pos), m_node2Pos(node2Pos), 
	  m_botPressFlag(botPressFlag), m_axialCD(axialCD), m_axialCa(axialCa),
	  m_cog2node1(fill::zeros), m_cog2node2(fill::zeros), m_node1Vel(fill::zeros), m_node2Vel(fill::zeros), m_node1Acc(fill::zeros), m_node2Acc(fill::zeros)
{
	m_cog2node1 = m_node1Pos - cog;
	m_cog2node2 = m_node2Pos - cog;
		
	// Since Simpson's rule is employed for the integration of the forces along the 
	// Morison's element, we need to make sure that the number of integration points is odd
	if (numIntPoints % 2 == 0)
	{
		m_numIntPoints = numIntPoints + 1;
	}
	else
	{
		m_numIntPoints = numIntPoints;
	}
}

/*****************************************************
	Functions for node position / velocity / acceleration
*****************************************************/
void MorisonElement::updateNodesPosVelAcc(const vec::fixed<6> &floaterCoGpos, const vec::fixed<6> &floaterVel, const vec::fixed<6> &floaterAcc)
{
	mat::fixed<3, 3> RotatMatrix(MorisonElement::RotatMatrix(floaterCoGpos.rows(3,5))); // Calculate it here so we just need to calculate the matrix once
	vec::fixed<3> R1 = RotatMatrix * m_cog2node1; // R1 and R2 are the vectors that give the node1 and node2 positions with respect to the CoG of the structure
	vec::fixed<3> R2 = RotatMatrix * m_cog2node2;

	m_node1Pos = floaterCoGpos.rows(0, 2) + R1;
	m_node2Pos = floaterCoGpos.rows(0, 2) + R2;

	m_node1Vel = floaterVel.rows(0, 2) + arma::cross(floaterVel.rows(3, 5), R1);
	m_node2Vel = floaterVel.rows(0, 2) + arma::cross(floaterVel.rows(3, 5), R2);

	m_node1AccCentrip = arma::cross(floaterVel.rows(3, 5), arma::cross(floaterVel.rows(3, 5), R1));
	m_node2AccCentrip = arma::cross(floaterVel.rows(3, 5), arma::cross(floaterVel.rows(3, 5), R2));

	m_node1Acc = floaterAcc.rows(0, 2) + arma::cross(floaterAcc.rows(3, 5), R1) + m_node1AccCentrip;
	m_node2Acc = floaterAcc.rows(0, 2) + arma::cross(floaterAcc.rows(3, 5), R2) + m_node2AccCentrip;
}

mat::fixed<3, 3> MorisonElement::RotatMatrix(const vec::fixed<3> &rotation) const
{
	/* Rotation matrix is RotatMat = RotatX * RotatY * RotatZ, i.e. a rotation around the Z axis,
	  followed by a rotation around the new y axis, and a rotation around the new x axis. Each rotation matrix is given by:
	
	  RotatX = { { 1 ,        0        ,         0        },
			     { 0 , cos(rotation(3)) , -sin(rotation(3)) },
			     { 0 , sin(rotation(3)) ,  cos(rotation(3)) }
			   };

	  RotatY = { { cos(rotation(4))  , 0 ,  sin(rotation(4)) },
		         {        0         , 1 ,         0        },
			     { -sin(rotation(4)) , 0 , cos(rotation(4)) }
			   };

	  RotatZ = { { cos(rotation(5)) , -sin(rotation(5)) , 0 },			     
				 { sin(rotation(5)) ,  cos(rotation(5)) , 0 },
			     {        0        ,         0        , 1 },
			   };

	  And the resulting matrix is the one below
	*/
	mat::fixed<3, 3> RotatMatrix = { 
									{                          cos(rotation(1)) * cos(rotation(2))                               ,                          -cos(rotation(1)) * sin(rotation(2))                               ,            sin(rotation(1))          },
									{ cos(rotation(0)) * sin(rotation(2)) + sin(rotation(0)) * sin(rotation(1)) * cos(rotation(2))  ,  cos(rotation(0)) * cos(rotation(1)) - sin(rotation(1)) * sin(rotation(1)) * sin(rotation(2))  ,  -sin(rotation(0)) * cos(rotation(1)) },
									{ sin(rotation(0)) * sin(rotation(2)) - cos(rotation(0)) * sin(rotation(1)) * cos(rotation(2))  ,  sin(rotation(0)) * cos(rotation(2)) + cos(rotation(1)) * sin(rotation(1)) * sin(rotation(2))  ,   cos(rotation(0)) * cos(rotation(1)) }
								   };

	return RotatMatrix;
}


vec::fixed<3> MorisonElement::node1Pos() const
{
	return m_node1Pos;
}

vec::fixed<3> MorisonElement::node2Pos() const
{
	return m_node2Pos;
}

vec::fixed<3> MorisonElement::node1Vel() const
{
	return m_node1Vel;
}

vec::fixed<3> MorisonElement::node2Vel() const
{
	return m_node2Vel;
}

vec::fixed<3> MorisonElement::node1Acc() const
{
	return m_node1Acc;
}

vec::fixed<3> MorisonElement::node2Acc() const
{
	return m_node2Acc;;
}


vec::fixed<3> MorisonElement::node1AccCentrip() const
{
	return m_node1AccCentrip;
}

vec::fixed<3> MorisonElement::node2AccCentrip() const
{
	return m_node2AccCentrip;
}