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
	m_accNodeAtWL_1stOrd.zeros();
	m_node1Acc_1stOrd.zeros();
	m_node2Acc_1stOrd.zeros();

	// Don't know why, but without declaring the following sizes the code crashes due to lack of memory when compiling
	// the debug version in Visual Studio 2017
	m_waveElevAtWL = zeros(0, 0);
	m_hydroForce_1st_Array = zeros(0, 0);
	m_hydroForce_2nd_Array = zeros(0, 0);
	m_nodesArray = zeros(0, 0);

	m_du1dt_Array_x = zeros(0, 0);
	m_du1dt_Array_y = zeros(0, 0);
	m_du1dt_Array_z = zeros(0, 0);

	m_u1_Array_x = zeros(0, 0);
	m_u1_Array_y = zeros(0, 0);
	m_u1_Array_z = zeros(0, 0);

	m_du1dx_Array_x = zeros(0, 0);
	m_du1dy_Array_y = zeros(0, 0);
	m_du1dz_Array_z = zeros(0, 0);
	m_du1dx_Array_y = zeros(0, 0);
	m_du1dx_Array_z = zeros(0, 0);	
	m_du1dy_Array_z = zeros(0, 0);

	m_da1dx_Array_x = zeros(0, 0);
	m_da1dy_Array_y = zeros(0, 0);
	m_da1dz_Array_z = zeros(0, 0);
	m_da1dx_Array_y = zeros(0, 0);
	m_da1dx_Array_z = zeros(0, 0);
	m_da1dy_Array_z = zeros(0, 0);

	m_gradP1_Array_x = zeros(0, 0);
	m_gradP1_Array_y = zeros(0, 0);
	m_gradP1_Array_z = zeros(0, 0);
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

void MorisonElement::updateMorisonElement(const ENVIR &envir, const vec::fixed<6> &floaterCoGpos, const vec::fixed<6> &floaterVel, 
										 const vec::fixed<6> &floaterCoGpos_1stOrd, const vec::fixed<6> &floaterVel_1stOrd,
										 const vec::fixed<6> &floaterCoGpos_SD, const vec::fixed<6> &floaterVel_SD)
{
	updateMorisonElement(floaterCoGpos, floaterVel, floaterCoGpos_1stOrd, floaterVel_1stOrd, floaterCoGpos_SD, floaterVel_SD);

	// Find the intersection with the instantaneous waterline
	// Used in Wheeler stretching and in the calculation of the added mass matrix.
	// If no intersection with the WL is found, this is an array of NaNs
	m_intersectWL = findIntersectWL(envir);
}

void MorisonElement::updateMorisonElement(const vec::fixed<6> &floaterCoGpos, const vec::fixed<6> &floaterVel, 
										  const vec::fixed<6> &floaterCoGpos_1stOrd, const vec::fixed<6> &floaterVel_1stOrd,
										  const vec::fixed<6> &floaterCoGpos_SD, const vec::fixed<6> &floaterVel_SD)
{
	// Considering total body motion
	calcPosVel(floaterCoGpos, floaterVel, m_node1Pos, m_node2Pos, m_node1Vel, m_node2Vel, m_xvec, m_yvec, m_zvec);	

	// Considering only the 1st order motions
	calcPosVel(floaterCoGpos_1stOrd, floaterVel_1stOrd, m_node1Pos_1stOrd, m_node2Pos_1stOrd, m_node1Vel_1stOrd, m_node2Vel_1stOrd, m_xvec_1stOrd, m_yvec_1stOrd, m_zvec_1stOrd);
	mat::fixed<3, 3> R = rotatMatrix(floaterCoGpos_1stOrd.rows(3, 5));
	m_node1AccCentrip = arma::cross(floaterVel_1stOrd.rows(3, 5), arma::cross(floaterVel_1stOrd.rows(3, 5), R * m_cog2node1));
	m_node2AccCentrip = arma::cross(floaterVel_1stOrd.rows(3, 5), arma::cross(floaterVel_1stOrd.rows(3, 5), R * m_cog2node2));

	m_nodeWL = floaterCoGpos_1stOrd.rows(0, 2) + R * m_cog2nodeWL;
	m_Zwl = m_nodeWL.at(2);

	// Considering only the mean and slow drift motions
	calcPosVel(floaterCoGpos_SD, floaterVel_SD, m_node1Pos_sd, m_node2Pos_sd, m_node1Vel_sd, m_node2Vel_sd, m_xvec_sd, m_yvec_sd, m_zvec_sd);

	calculateImmersedLengthProperties_sd(m_Lw, m_numNodesBelowWL, m_dL);
}

void MorisonElement::updateAcc1stOrdNodeAtWL(const vec::fixed<6> &floaterAcc_1stOrd)
{
	m_accNodeAtWL_1stOrd = floaterAcc_1stOrd.rows(0, 2) + arma::cross(floaterAcc_1stOrd.rows(3, 5), m_cog2nodeWL);
	m_node1Acc_1stOrd = floaterAcc_1stOrd.rows(0, 2) + arma::cross(floaterAcc_1stOrd.rows(3, 5), m_cog2node1);
	m_node2Acc_1stOrd = floaterAcc_1stOrd.rows(0, 2) + arma::cross(floaterAcc_1stOrd.rows(3, 5), m_cog2node2);
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

vec::fixed<3> MorisonElement::nodePos_sd(const int nodeIndex) const
{
	if (m_u1_Array_x.is_empty())
	{
		return m_node1Pos_sd + m_dL * nodeIndex * m_zvec_sd;
	}
	else
	{
		return m_nodesArray.col(nodeIndex);
	}
}

double MorisonElement::waveElevAtWL(const ENVIR &envir) const
{
	double eta;
	if (m_waveElevAtWL.is_empty())
	{
		vec::fixed<3> n_wl = nodePos_sd(m_numNodesBelowWL-1); // Coordinates of the intersection with the still water line
		eta = envir.waveElev(n_wl.at(0), n_wl.at(1));
	}
	else
	{
		const vec &t = envir.getTimeArray();
		uword ind1 = envir.getInd4interp1();
		uword ind2 = envir.getInd4interp2();

		eta = m_waveElevAtWL.at(ind1);
		if (envir.shouldInterp())
		{
			eta += (m_waveElevAtWL.at(ind2) - m_waveElevAtWL.at(ind1)) * (envir.time() - t(ind1)) / (t(ind2) - t(ind1));
		}
	}
	return eta;
}

vec::fixed<3> MorisonElement::u1(const ENVIR &envir, const int nodeIndex) const
{
	vec::fixed<3> u1;
	if (m_u1_Array_x.is_empty())
	{
		vec::fixed<3> node = nodePos_sd(nodeIndex);
		double eta{ 0 };
		if (envir.waveStret() == 2)
		{
			eta = envir.waveElev(node.at(0), node.at(1));
		}
		u1 = envir.u1(node, eta);
	}
	else
	{
		const vec &t = envir.getTimeArray();
		uword ind1 = envir.getInd4interp1();
		uword ind2 = envir.getInd4interp2();
		
		double u_x = m_u1_Array_x.at(ind1, nodeIndex);
		double u_y = m_u1_Array_y.at(ind1, nodeIndex);
		double u_z = m_u1_Array_z.at(ind1, nodeIndex);
		if (envir.shouldInterp())
		{
			u_x += (m_u1_Array_x.at(ind2, nodeIndex) - m_u1_Array_x.at(ind1, nodeIndex)) * (envir.time() - t(ind1)) / (t(ind2) - t(ind1));
			u_y += (m_u1_Array_y.at(ind2, nodeIndex) - m_u1_Array_y.at(ind1, nodeIndex)) * (envir.time() - t(ind1)) / (t(ind2) - t(ind1));
			u_z += (m_u1_Array_z.at(ind2, nodeIndex) - m_u1_Array_z.at(ind1, nodeIndex)) * (envir.time() - t(ind1)) / (t(ind2) - t(ind1));
		}
		u1 = { u_x, u_y, u_z };
	}
	return u1;
}

vec::fixed<3> MorisonElement::du1dt(const ENVIR & envir, const int nodeIndex) const
{
	vec::fixed<3> du1dt;
	if (m_du1dt_Array_x.is_empty())
	{
		vec::fixed<3> node = nodePos_sd(nodeIndex);
		double eta{ 0 };
		if (envir.waveStret() == 2)
		{
			eta = envir.waveElev(node.at(0), node.at(1));
		}
		du1dt = envir.du1dt(node, eta);
	}
	else
	{
		const vec &t = envir.getTimeArray();
		uword ind1 = envir.getInd4interp1();
		uword ind2 = envir.getInd4interp2();

		double du1dt_x = m_du1dt_Array_x.at(ind1, nodeIndex);
		double du1dt_y = m_du1dt_Array_y.at(ind1, nodeIndex);
		double du1dt_z = m_du1dt_Array_z.at(ind1, nodeIndex);
		if (envir.shouldInterp())
		{
			du1dt_x += (m_du1dt_Array_x.at(ind2, nodeIndex) - m_du1dt_Array_x.at(ind1, nodeIndex)) * (envir.time() - t(ind1)) / (t(ind2) - t(ind1));
			du1dt_y += (m_du1dt_Array_y.at(ind2, nodeIndex) - m_du1dt_Array_y.at(ind1, nodeIndex)) * (envir.time() - t(ind1)) / (t(ind2) - t(ind1));
			du1dt_z += (m_du1dt_Array_z.at(ind2, nodeIndex) - m_du1dt_Array_z.at(ind1, nodeIndex)) * (envir.time() - t(ind1)) / (t(ind2) - t(ind1));
		}
		du1dt = { du1dt_x, du1dt_y, du1dt_z };
	}
	return du1dt;
}

vec::fixed<3> MorisonElement::du1dx(const ENVIR &envir, const int nodeIndex) const
{
	vec::fixed<3> du1dx;
	if (m_du1dx_Array_x.is_empty())
	{
		vec::fixed<3> node = nodePos_sd(nodeIndex);
		double eta{ 0 };
		if (envir.waveStret() == 2)
		{
			eta = envir.waveElev(node.at(0), node.at(1));
		}
	
		du1dx = envir.du1dx(node, eta);
	}
	else
	{
		const vec &t = envir.getTimeArray();
		uword ind1 = envir.getInd4interp1();
		uword ind2 = envir.getInd4interp2();

		double du1dx_x = m_du1dx_Array_x.at(ind1, nodeIndex);
		double du1dx_y = m_du1dx_Array_y.at(ind1, nodeIndex);
		double du1dx_z = m_du1dx_Array_z.at(ind1, nodeIndex);
		if (envir.shouldInterp())
		{
			du1dx_x += (m_du1dx_Array_x.at(ind2, nodeIndex) - m_du1dx_Array_x.at(ind1, nodeIndex)) * (envir.time() - t(ind1)) / (t(ind2) - t(ind1));
			du1dx_y += (m_du1dx_Array_y.at(ind2, nodeIndex) - m_du1dx_Array_y.at(ind1, nodeIndex)) * (envir.time() - t(ind1)) / (t(ind2) - t(ind1));
			du1dx_z += (m_du1dx_Array_z.at(ind2, nodeIndex) - m_du1dx_Array_z.at(ind1, nodeIndex)) * (envir.time() - t(ind1)) / (t(ind2) - t(ind1));
		}
		du1dx = { du1dx_x, du1dx_y, du1dx_z };
	}
	return du1dx;
}

vec::fixed<3> MorisonElement::du1dy(const ENVIR &envir, const int nodeIndex) const
{
	vec::fixed<3> du1dy;
	if (m_du1dx_Array_y.is_empty())
	{
		vec::fixed<3> node = nodePos_sd(nodeIndex);
		double eta{ 0 };
		if (envir.waveStret() == 2)
		{
			eta = envir.waveElev(node.at(0), node.at(1));
		}	
		du1dy = envir.du1dy(node, eta);
	}
	else
	{
		const vec &t = envir.getTimeArray();
		uword ind1 = envir.getInd4interp1();
		uword ind2 = envir.getInd4interp2();

		// Remember of symmetries of the velocity gradient
		// e.g. du1dy_x = du1dx_y
		double du1dy_x = m_du1dx_Array_y.at(ind1, nodeIndex);
		double du1dy_y = m_du1dy_Array_y.at(ind1, nodeIndex);
		double du1dy_z = m_du1dy_Array_z.at(ind1, nodeIndex);
		if (envir.shouldInterp())
		{
			du1dy_x += (m_du1dx_Array_y.at(ind2, nodeIndex) - m_du1dx_Array_y.at(ind1, nodeIndex)) * (envir.time() - t(ind1)) / (t(ind2) - t(ind1));
			du1dy_y += (m_du1dy_Array_y.at(ind2, nodeIndex) - m_du1dy_Array_y.at(ind1, nodeIndex)) * (envir.time() - t(ind1)) / (t(ind2) - t(ind1));
			du1dy_z += (m_du1dy_Array_z.at(ind2, nodeIndex) - m_du1dy_Array_z.at(ind1, nodeIndex)) * (envir.time() - t(ind1)) / (t(ind2) - t(ind1));
		}
		du1dy = { du1dy_x, du1dy_y, du1dy_z };
	}
	return du1dy;
}

vec::fixed<3> MorisonElement::du1dz(const ENVIR &envir, const int nodeIndex) const
{
	vec::fixed<3> du1dz;

	if (m_du1dx_Array_z.is_empty())
	{
		vec::fixed<3> node = nodePos_sd(nodeIndex);
		double eta{ 0 };
		if (envir.waveStret() == 2)
		{
			eta = envir.waveElev(node.at(0), node.at(1));
		}		
		du1dz = envir.du1dz(node, eta);
	}
	else
	{
		const vec &t = envir.getTimeArray();
		uword ind1 = envir.getInd4interp1();
		uword ind2 = envir.getInd4interp2();

		// Remember of symmetries of the velocity gradient
		// e.g. du1dz_x = du1dx_z
		double du1dz_x = m_du1dx_Array_z.at(ind1, nodeIndex);
		double du1dz_y = m_du1dy_Array_z.at(ind1, nodeIndex);
		double du1dz_z = m_du1dz_Array_z.at(ind1, nodeIndex);
		if (envir.shouldInterp())
		{
			du1dz_x += (m_du1dx_Array_z.at(ind2, nodeIndex) - m_du1dx_Array_z.at(ind1, nodeIndex)) * (envir.time() - t(ind1)) / (t(ind2) - t(ind1));
			du1dz_y += (m_du1dy_Array_z.at(ind2, nodeIndex) - m_du1dy_Array_z.at(ind1, nodeIndex)) * (envir.time() - t(ind1)) / (t(ind2) - t(ind1));
			du1dz_z += (m_du1dz_Array_z.at(ind2, nodeIndex) - m_du1dz_Array_z.at(ind1, nodeIndex)) * (envir.time() - t(ind1)) / (t(ind2) - t(ind1));
		}
		du1dz = { du1dz_x, du1dz_y, du1dz_z };
	}
	return du1dz;
}

vec::fixed<3> MorisonElement::da1dx(const ENVIR &envir, const int nodeIndex) const
{
	vec::fixed<3> da1dx;
	if (m_da1dx_Array_x.is_empty())
	{
		vec::fixed<3> node = nodePos_sd(nodeIndex);
		double eta{ 0 };
		if (envir.waveStret() == 2)
		{
			eta = envir.waveElev(node.at(0), node.at(1));
		}

		da1dx = envir.da1dx(node, eta);
	}
	else
	{
		const vec &t = envir.getTimeArray();
		uword ind1 = envir.getInd4interp1();
		uword ind2 = envir.getInd4interp2();

		double da1dx_x = m_da1dx_Array_x.at(ind1, nodeIndex);
		double da1dx_y = m_da1dx_Array_y.at(ind1, nodeIndex);
		double da1dx_z = m_da1dx_Array_z.at(ind1, nodeIndex);
		if (envir.shouldInterp())
		{
			da1dx_x += (m_da1dx_Array_x.at(ind2, nodeIndex) - m_da1dx_Array_x.at(ind1, nodeIndex)) * (envir.time() - t(ind1)) / (t(ind2) - t(ind1));
			da1dx_y += (m_da1dx_Array_y.at(ind2, nodeIndex) - m_da1dx_Array_y.at(ind1, nodeIndex)) * (envir.time() - t(ind1)) / (t(ind2) - t(ind1));
			da1dx_z += (m_da1dx_Array_z.at(ind2, nodeIndex) - m_da1dx_Array_z.at(ind1, nodeIndex)) * (envir.time() - t(ind1)) / (t(ind2) - t(ind1));
		}
		da1dx = { da1dx_x, da1dx_y, da1dx_z };
	}
	return da1dx;
}

vec::fixed<3> MorisonElement::da1dy(const ENVIR &envir, const int nodeIndex) const
{
	vec::fixed<3> da1dy;
	if (m_da1dx_Array_y.is_empty())
	{
		vec::fixed<3> node = nodePos_sd(nodeIndex);
		double eta{ 0 };
		if (envir.waveStret() == 2)
		{
			eta = envir.waveElev(node.at(0), node.at(1));
		}
		da1dy = envir.da1dy(node, eta);
	}
	else
	{
		const vec &t = envir.getTimeArray();
		uword ind1 = envir.getInd4interp1();
		uword ind2 = envir.getInd4interp2();

		// Remember of symmetries of the velocity gradient
		// e.g. da1dy_x = da1dx_y
		double da1dy_x = m_da1dx_Array_y.at(ind1, nodeIndex);
		double da1dy_y = m_da1dy_Array_y.at(ind1, nodeIndex);
		double da1dy_z = m_da1dy_Array_z.at(ind1, nodeIndex);
		if (envir.shouldInterp())
		{
			da1dy_x += (m_da1dx_Array_y.at(ind2, nodeIndex) - m_da1dx_Array_y.at(ind1, nodeIndex)) * (envir.time() - t(ind1)) / (t(ind2) - t(ind1));
			da1dy_y += (m_da1dy_Array_y.at(ind2, nodeIndex) - m_da1dy_Array_y.at(ind1, nodeIndex)) * (envir.time() - t(ind1)) / (t(ind2) - t(ind1));
			da1dy_z += (m_da1dy_Array_z.at(ind2, nodeIndex) - m_da1dy_Array_z.at(ind1, nodeIndex)) * (envir.time() - t(ind1)) / (t(ind2) - t(ind1));
		}
		da1dy = { da1dy_x, da1dy_y, da1dy_z };
	}
	return da1dy;
}

vec::fixed<3> MorisonElement::da1dz(const ENVIR &envir, const int nodeIndex) const
{
	vec::fixed<3> da1dz;

	if (m_da1dx_Array_z.is_empty())
	{
		vec::fixed<3> node = nodePos_sd(nodeIndex);
		double eta{ 0 };
		if (envir.waveStret() == 2)
		{
			eta = envir.waveElev(node.at(0), node.at(1));
		}
		da1dz = envir.da1dz(node, eta);
	}
	else
	{
		const vec &t = envir.getTimeArray();
		uword ind1 = envir.getInd4interp1();
		uword ind2 = envir.getInd4interp2();

		// Remember of symmetries of the velocity gradient
		// e.g. da1dz_x = da1dx_z
		double da1dz_x = m_da1dx_Array_z.at(ind1, nodeIndex);
		double da1dz_y = m_da1dy_Array_z.at(ind1, nodeIndex);
		double da1dz_z = m_da1dz_Array_z.at(ind1, nodeIndex);
		if (envir.shouldInterp())
		{
			da1dz_x += (m_da1dx_Array_z.at(ind2, nodeIndex) - m_da1dx_Array_z.at(ind1, nodeIndex)) * (envir.time() - t(ind1)) / (t(ind2) - t(ind1));
			da1dz_y += (m_da1dy_Array_z.at(ind2, nodeIndex) - m_da1dy_Array_z.at(ind1, nodeIndex)) * (envir.time() - t(ind1)) / (t(ind2) - t(ind1));
			da1dz_z += (m_da1dz_Array_z.at(ind2, nodeIndex) - m_da1dz_Array_z.at(ind1, nodeIndex)) * (envir.time() - t(ind1)) / (t(ind2) - t(ind1));
		}
		da1dz = { da1dz_x, da1dz_y, da1dz_z };
	}
	return da1dz;
}

vec::fixed<3> MorisonElement::gradP1(const ENVIR &envir, const int nodeIndex) const
{
	vec::fixed<3> gradP1;
	int localNodeIndex = nodeIndex;

	// As gradP1 is used only at the end nodes of the cylinder, the array has only 2 columns
	if (nodeIndex == m_nodesArray.n_cols-1) localNodeIndex = 1;

	if (localNodeIndex < 0 || localNodeIndex > 1)
		throw std::runtime_error("MorisonElement::gradP1 can only be called for the bottom or end nodes of the cylinder.");

	if (m_gradP1_Array_x.is_empty())
	{
		vec::fixed<3> node = nodePos_sd(nodeIndex);
		double eta{ 0 };
		if (envir.waveStret() == 2)
		{
			eta = envir.waveElev(node.at(0), node.at(1));
		}
		gradP1 = envir.gradP1(node, eta);
	}
	else
	{
		const vec &t = envir.getTimeArray();
		uword ind1 = envir.getInd4interp1();
		uword ind2 = envir.getInd4interp2();

		double gradP1_x = m_gradP1_Array_x.at(ind1, localNodeIndex);
		double gradP1_y = m_gradP1_Array_y.at(ind1, localNodeIndex);
		double gradP1_z = m_gradP1_Array_z.at(ind1, localNodeIndex);
		if (envir.shouldInterp())
		{
			gradP1_x += (m_gradP1_Array_x.at(ind2, localNodeIndex) - m_gradP1_Array_x.at(ind1, localNodeIndex)) * (envir.time() - t(ind1)) / (t(ind2) - t(ind1));
			gradP1_y += (m_gradP1_Array_y.at(ind2, localNodeIndex) - m_gradP1_Array_y.at(ind1, localNodeIndex)) * (envir.time() - t(ind1)) / (t(ind2) - t(ind1));
			gradP1_z += (m_gradP1_Array_z.at(ind2, localNodeIndex) - m_gradP1_Array_z.at(ind1, localNodeIndex)) * (envir.time() - t(ind1)) / (t(ind2) - t(ind1));
		}
		gradP1 = { gradP1_x, gradP1_y, gradP1_z };
	}
	return gradP1;
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

bool MorisonElement::flagFixed() const
{
	return m_flagFixed;
}
