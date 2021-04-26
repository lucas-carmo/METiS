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
	// Make sure that node1 is below node2 (or at the same height, at least).
	// Otherwise, need to swap them.
	if (node1Pos[2] > node2Pos[2])
	{
		m_node1Pos.swap(m_node2Pos);
		m_node1Pos_t0.swap(m_node2Pos_t0);
		std::swap(m_axialCD_1, m_axialCD_2);
		std::swap(m_axialCa_1, m_axialCa_2);
	}
			
	m_cog2node1 = m_node1Pos - cog;
	m_cog2node2 = m_node2Pos - cog;

	m_nodeWL = m_node1Pos + (m_node2Pos - m_node1Pos) * (0 - m_node1Pos.at(2)) / (m_node2Pos.at(2) - m_node1Pos.at(2));
	m_cog2nodeWL = m_nodeWL - cog;
	m_Zwl = m_nodeWL.at(2);

	// Don't know why, but without declaring the following sizes the code crashes due to lack of memory when compiling
	// the debug version in Visual Studio 2017
	m_waveElevAtWL = zeros(0, 0);
	m_hydroForce_1st_Array = zeros(0, 0);
	m_hydroForce_2nd_Array = zeros(0, 0);
	m_nodesArray = zeros(0, 0);

	m_u1_Array_x = zeros(0, 0);
	m_u1_Array_y = zeros(0, 0);
	m_u1_Array_z = zeros(0, 0);

	m_du1dx_Array_x = zeros(0, 0);
	m_du1dy_Array_y = zeros(0, 0);
	m_du1dz_Array_z = zeros(0, 0);
	m_du1dx_Array_y = zeros(0, 0);
	m_du1dx_Array_z = zeros(0, 0);	
	m_du1dy_Array_z = zeros(0, 0);	
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

	make_local_base(xvec, yvec, zvec, node1Pos, node2Pos);
}

void MorisonElement::updateMorisonElement(const ENVIR &envir, const vec::fixed<6> &floaterCoGpos, const vec::fixed<6> &floaterVel, const vec::fixed<6> &floaterCoGpos_SD, const vec::fixed<6> &floaterVel_SD)
{
	updateMorisonElement(floaterCoGpos, floaterVel, floaterCoGpos_SD, floaterVel_SD);

	// Find the intersection with the instantaneous waterline
	// Used in Wheeler stretching and in the calculation of the added mass matrix.
	// If no intersection with the WL is found, this is an array of NaNs
	m_intersectWL = findIntersectWL(envir);
}

void MorisonElement::updateMorisonElement(const vec::fixed<6> &floaterCoGpos, const vec::fixed<6> &floaterVel, const vec::fixed<6> &floaterCoGpos_SD, const vec::fixed<6> &floaterVel_SD)
{
	// Considering total body motion
	calcPosVel(floaterCoGpos, floaterVel, m_node1Pos, m_node2Pos, m_node1Vel, m_node2Vel, m_xvec, m_yvec, m_zvec);
	m_RotatMatrix = rotatMatrix(floaterCoGpos.rows(3, 5));
	m_node1AccCentrip = arma::cross(floaterVel.rows(3, 5), arma::cross(floaterVel.rows(3, 5), m_RotatMatrix * m_cog2node1));
	m_node2AccCentrip = arma::cross(floaterVel.rows(3, 5), arma::cross(floaterVel.rows(3, 5), m_RotatMatrix * m_cog2node2));

	m_nodeWL = floaterCoGpos.rows(0, 2) + m_RotatMatrix * m_cog2nodeWL;
	m_Zwl = m_nodeWL.at(2);

	// Considering only the mean and slow drift motions
	calcPosVel(floaterCoGpos_SD, floaterVel_SD, m_node1Pos_sd, m_node2Pos_sd, m_node1Vel_sd, m_node2Vel_sd, m_xvec_sd, m_yvec_sd, m_zvec_sd);
	m_RotatMatrix_sd = rotatMatrix(floaterCoGpos_SD.rows(3, 5));	
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
// Returns arma::datum::nan if no intersection is found
// For computational speed, it assumes that the waves are all so long and the cylinder is almost vertical, so that 
// the wave elevation is uniform across the cylinder's length.
vec::fixed<3> MorisonElement::findIntersectWL(const ENVIR &envir) const
{
	// Nodes position
	vec::fixed<3> n1 = node1Pos();
	vec::fixed<3> n2 = node2Pos();

	double zwl = 0;
	if (!m_waveElevAtWL.is_empty())
	{
		uword ind1 = envir.getInd4interp1();		
		const vec &t = envir.getTimeArray();

		zwl = m_waveElevAtWL.at(envir.getInd4interp1());
		if (envir.shouldInterp())
		{
			uword ind2 = envir.getInd4interp2();
			zwl += (m_waveElevAtWL.at(ind2) - m_waveElevAtWL.at(ind1)) * (envir.time() - t(ind1)) / (t(ind2) - t(ind1));
		}
	}
	else
	{
		vec::fixed<3> nodeAtMeanWL = n1 + (std::abs(0 - n1[2]) / (n2[2] - n1[2])) * norm(n2 - n1) * m_zvec;
		zwl = envir.waveElev(nodeAtMeanWL.at(0), nodeAtMeanWL.at(1));
	}

	return n1 + (std::abs(zwl - n1[2]) / (n2[2] - n1[2])) * norm(n2 - n1) * m_zvec;
}

void MorisonElement::calculateImmersedLengthProperties(double &Lw, int &ncyl, double &dL) const
{
	double L{ norm(m_node2Pos - m_node1Pos) };
	Lw = L;
	ncyl = m_numIntPoints;
	if (m_node2Pos.at(2) > 0)
	{
		Lw = L * (0 - m_node1Pos.at(2)) / (m_node2Pos.at(2) - m_node1Pos.at(2));
		ncyl = static_cast<int>(std::ceil(Lw / L * (m_numIntPoints - 1))) + 1; // Number of points to discretize the part of the cylinder that is below the water		
	}
	if (ncyl % 2 == 0)
	{
		ncyl += 1;
	}
	dL = Lw / (ncyl - 1); // length of each interval between points
}

void MorisonElement::calculateImmersedLengthProperties_sd(double &Lw, int &ncyl, double &dL) const
{
	double L{ norm(m_node2Pos_sd - m_node1Pos_sd) };
	Lw = L;
	ncyl = m_numIntPoints;
	if (m_node2Pos_sd.at(2) > 0)
	{
		Lw = L * (0 - m_node1Pos_sd.at(2)) / (m_node2Pos_sd.at(2) - m_node1Pos_sd.at(2));
		ncyl = static_cast<int>(std::ceil(Lw / L * (m_numIntPoints - 1))) + 1; // Number of points to discretize the part of the cylinder that is below the water		
	}
	if (ncyl % 2 == 0)
	{
		ncyl += 1;
	}
	dL = Lw / (ncyl - 1); // length of each interval between points
}
