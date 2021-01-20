#include "MorisonElement.h"
#include "auxFunctions.h"

using namespace arma;

/*****************************************************
	Constructors
*****************************************************/
MorisonElement::MorisonElement(const vec &node1Pos, const vec &node2Pos, const vec &cog, const int numIntPoints,
							   const bool botPressFlag, const double axialCD_1, const double axialCa_1, const double axialCD_2, const double axialCa_2)
	: m_node1Pos(node1Pos), m_node2Pos(node2Pos), m_numIntPoints(numIntPoints),
	  m_botPressFlag(botPressFlag), m_axialCD_1(axialCD_1), m_axialCa_1(axialCa_1), m_axialCD_2(axialCD_2), m_axialCa_2(axialCa_2),
	  m_cog2node1(fill::zeros), m_cog2node2(fill::zeros), m_node1Vel(fill::zeros), m_node2Vel(fill::zeros),
	  m_node1Pos_t0(node1Pos), m_node2Pos_t0(node2Pos), m_node1Pos_sd(node1Pos), m_node2Pos_sd(node2Pos)
{
	m_cog2node1 = m_node1Pos - cog;
	m_cog2node2 = m_node2Pos - cog;		
}

/*****************************************************
	Functions for node position / velocity / acceleration
*****************************************************/
void MorisonElement::calcPosVel(const vec::fixed<6> &pos, const vec::fixed<6> &vel,
								   vec::fixed<3> &node1Pos, vec::fixed<3> &node2Pos, vec::fixed<3> &node1Vel, vec::fixed<3> &node2Vel,
                                   vec::fixed<3> &xvec, vec::fixed<3> &yvec, vec::fixed<3> &zvec)
{
	mat::fixed<3, 3> RotatMatrix(rotatMatrix(pos.rows(3, 5)));

	vec::fixed<3> R1 = RotatMatrix * m_cog2node1; // R1 and R2 are the vectors that give the node1 and node2 positions with respect to the CoG of the structure
	vec::fixed<3> R2 = RotatMatrix * m_cog2node2;

	node1Pos = pos.rows(0,2) + R1;
	node2Pos = pos.rows(0, 2) + R2;

	node1Vel = vel.rows(0, 2) + arma::cross(vel.rows(3, 5), R1);
	node2Vel = vel.rows(0, 2) + arma::cross(vel.rows(3, 5), R2);

	xvec = RotatMatrix * m_xvec_t0;
	yvec = RotatMatrix * m_yvec_t0;
	zvec = RotatMatrix * m_zvec_t0;
}

void MorisonElement::updateMorisonElement(const ENVIR &envir, const vec::fixed<6> &floaterCoGpos, const vec::fixed<6> &floaterVel, const vec::fixed<6> &floaterCoGpos_SD, const vec::fixed<6> &floaterVel_SD, const vec::fixed<6> &floaterCoGpos_1stOrd, const vec::fixed<6> &floaterVel_1stOrd)
{
	// Considering total body motion
	calcPosVel(floaterCoGpos, floaterVel, m_node1Pos, m_node2Pos, m_node1Vel, m_node2Vel, m_xvec, m_yvec, m_zvec);
	mat::fixed<3, 3> RotatMatrix(rotatMatrix(floaterCoGpos.rows(3, 5)));
	m_node1AccCentrip = arma::cross(floaterVel.rows(3, 5), arma::cross(floaterVel.rows(3, 5), RotatMatrix * m_cog2node1));
	m_node2AccCentrip = arma::cross(floaterVel.rows(3, 5), arma::cross(floaterVel.rows(3, 5), RotatMatrix * m_cog2node2));

	// Considering only the mean and slow drift motions
	calcPosVel(floaterCoGpos_SD, floaterVel_SD, m_node1Pos_sd, m_node2Pos_sd, m_node1Vel_sd, m_node2Vel_sd, m_xvec_sd, m_yvec_sd, m_zvec_sd);

	// Considering only motions due to firt-order wave forces
	calcPosVel(floaterCoGpos_1stOrd, floaterVel_1stOrd, m_node1Pos_1stOrd, m_node2Pos_1stOrd, m_node1Vel_1stOrd, m_node2Vel_1stOrd, m_xvec_1stOrd, m_yvec_1stOrd, m_zvec_1stOrd);

	// Find the intersection with the instantaneous waterline
	// Used in Wheeler stretching and in the calculation of the added mass matrix.
	// If no intersection with the WL is found, this is an array of NaNs
	m_intersectWL = findIntersectWL(envir);
}

vec::fixed<3> MorisonElement::node1Pos_t0() const
{
	return m_node1Pos_t0;
}

vec::fixed<3> MorisonElement::node2Pos_t0() const
{
	return m_node2Pos_t0;
}

vec::fixed<3> MorisonElement::node1Pos_sd() const
{
	return m_node1Pos_sd;
}

vec::fixed<3> MorisonElement::node2Pos_sd() const
{
	return m_node2Pos_sd;
}

vec::fixed<3> MorisonElement::node1Pos_1stOrd() const
{
	return m_node1Pos_1stOrd;
}

vec::fixed<3> MorisonElement::node2Pos_1stOrd() const
{
	return m_node2Pos_1stOrd;
}

vec::fixed<3> MorisonElement::node1Pos() const
{
	return m_node1Pos;
}

vec::fixed<3> MorisonElement::node2Pos() const
{
	return m_node2Pos;
}

vec::fixed<3> MorisonElement::node1Vel_sd() const
{
	return m_node1Vel_sd;
}

vec::fixed<3> MorisonElement::node2Vel_sd() const
{
	return m_node2Vel_sd;
}

vec::fixed<3> MorisonElement::node1Vel_1stOrd() const
{
	return m_node1Vel_1stOrd;
}

vec::fixed<3> MorisonElement::node2Vel_1stOrd() const
{
	return m_node2Vel_1stOrd;
}

vec::fixed<3> MorisonElement::node1Vel() const
{
	return m_node1Vel;
}

vec::fixed<3> MorisonElement::node2Vel() const
{
	return m_node2Vel;
}

vec::fixed<3> MorisonElement::node1AccCentrip() const
{
	return m_node1AccCentrip;
}

vec::fixed<3> MorisonElement::node2AccCentrip() const
{
	return m_node2AccCentrip;
}

// Calculate the intersection between the cylinder and waterline, considering the wave elevation.
// Returns arma::datum::nan if no intersection is found.
vec::fixed<3> MorisonElement::findIntersectWL(const ENVIR &envir) const
{
	// Nodes position
	vec::fixed<3> n1 = node1Pos();
	vec::fixed<3> n2 = node2Pos();

	// Find the intersection between the cylinder and the waterline.
	//
	// The equation we want to solve is s[2] - eta(s[0], s[1]) = 0, i.e. we want to find
	// when the 'z' coordinate of a point 's' along the cylinder is equal to the wave elevation
	// at the 'x' and 'y' coordinates of the same point 's'.		
	//
	// The first thing we do is to check if there is a solution.
	// This is not perfectly acurate, since there can be more than one intersection with the waterline.
	// However, this would be a pathological case that we do not want to deal with (at least by now).
	if ((n1[2] - envir.waveElev(n1[0], n1[1])) * (n2[2] - envir.waveElev(n2[0], n2[1])) >= 0)
	{
		return vec::fixed<3> {datum::nan, datum::nan, datum::nan};
	}

	// Brackets for solving the equation using the coordinate along the cylinder length
	vec::fixed<3> a = n1;
	vec::fixed<3> b = n2;

	// x is the coordinate along the cylinder length, 
	// with x = 0 corresponding to s = n1 and x = 1 to s = n2
	// s_i is the guess in the previous step and s_j in the current
	vec::fixed<3> s_i = n1;
	vec::fixed<3> s_j = n1;

	double eps = 0.001;
	while (std::abs(s_j[2] - envir.waveElev(s_j[0], s_j[1])) > eps)
	{
		s_i = s_j;
		if ((a[2] - envir.waveElev(a[0], a[1])) * (s_i[2] - envir.waveElev(s_i[0], s_i[1])) < 0) // Test with limit a
		{
			s_j = (a + s_i) / 2;
			b = s_i;
		}

		else if ((b[2] - envir.waveElev(b[0], b[1])) * (s_i[2] - envir.waveElev(s_i[0], s_i[1])) < 0) // Test with limit b
		{
			s_j = (b + s_i) / 2;
			a = s_i;
		}

		else // It is very unlikely that s_i will be the zero, but this possibility will be covered anyway
		{
			return s_i;
		}
	}

	return s_j;
}