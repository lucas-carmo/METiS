#include "MorisonCirc.h"
#include <cmath>


using namespace arma;

/*****************************************************
	Constructors
*****************************************************/
MorisonCirc::MorisonCirc(const vec &node1Pos, const vec &node2Pos, const vec &cog, const int numIntPoints,
	const bool botPressFlag, const double axialCD, const double axialCa,
	const double diam, const double CD, double CM, const double botDiam, const double topDiam)
	: MorisonElement(node1Pos, node2Pos, cog, numIntPoints, botPressFlag, axialCD, axialCa),
	m_diam(diam), m_CD(CD), m_CM(CM), m_botDiam(botDiam), m_topDiam(topDiam)
{}



/*****************************************************
	Forces acting on the Morison Element
*****************************************************/
void MorisonCirc::make_local_base(arma::vec::fixed<3> &xvec, arma::vec::fixed<3> &yvec, arma::vec::fixed<3> &zvec) const
{
	xvec.zeros();
	yvec.zeros();
	zvec = (m_node2Pos - m_node1Pos) / arma::norm(m_node2Pos - m_node1Pos, 2);

	// if Zlocal == Zglobal, REFlocal = REFglobal
	arma::vec::fixed<3> z_global = { 0,0,1 };
	if (arma::approx_equal(zvec, z_global, "absdiff", arma::datum::eps))
	{
		xvec = { 1, 0, 0 };
		yvec = { 0, 1, 0 };
	}
	else
	{
		yvec = arma::cross(z_global, zvec);
		yvec = yvec / arma::norm(yvec, 2);

		xvec = arma::cross(yvec, zvec);
		xvec = xvec / arma::norm(xvec, 2);
	}
}

// Same as make_local_base, but considering the initial position of the cylinder
void MorisonCirc::make_local_base_t0(arma::vec::fixed<3> &xvec, arma::vec::fixed<3> &yvec, arma::vec::fixed<3> &zvec) const
{
	xvec.zeros();
	yvec.zeros();
	zvec = (m_node2Pos_t0 - m_node1Pos_t0) / arma::norm(m_node2Pos_t0 - m_node1Pos_t0, 2);

	// if Zlocal == Zglobal, REFlocal = REFglobal
	arma::vec::fixed<3> z_global = { 0,0,1 };
	if (arma::approx_equal(zvec, z_global, "absdiff", arma::datum::eps))
	{
		xvec = { 1, 0, 0 };
		yvec = { 0, 1, 0 };
	}
	else
	{
		yvec = arma::cross(z_global, zvec);
		yvec = yvec / arma::norm(yvec, 2);

		xvec = arma::cross(yvec, zvec);
		xvec = xvec / arma::norm(xvec, 2);
	}
}

// Same as make_local_base, but considering only the mean and slow drift position of the cylinder
void MorisonCirc::make_local_base_sd(arma::vec::fixed<3> &xvec, arma::vec::fixed<3> &yvec, arma::vec::fixed<3> &zvec) const
{
	xvec.zeros();
	yvec.zeros();
	zvec = (m_node2Pos_sd - m_node1Pos_sd) / arma::norm(m_node2Pos_sd - m_node1Pos_sd, 2);

	// if Zlocal == Zglobal, REFlocal = REFglobal
	arma::vec::fixed<3> z_global = { 0,0,1 };
	if (arma::approx_equal(zvec, z_global, "absdiff", arma::datum::eps))
	{
		xvec = { 1, 0, 0 };
		yvec = { 0, 1, 0 };
	}
	else
	{
		yvec = arma::cross(z_global, zvec);
		yvec = yvec / arma::norm(yvec, 2);

		xvec = arma::cross(yvec, zvec);
		xvec = xvec / arma::norm(xvec, 2);
	}
}


vec::fixed<6> MorisonCirc::hydrostaticForce(const double rho, const double g) const
{
	// Forces and moments acting at the Morison Element
	vec::fixed<6> force(fill::zeros);

	// Use a more friendly notation
	double D = m_diam;

	// Nodes position
	vec::fixed<3> n1 = node1Pos();
	vec::fixed<3> n2 = node2Pos();

	// Vectors of the local coordinate system vectors
	vec::fixed<3> xvec(fill::zeros);
	vec::fixed<3> yvec(fill::zeros);
	vec::fixed<3> zvec(fill::zeros);
	MorisonCirc::make_local_base(xvec, yvec, zvec);

	// Make sure that node1 is below node2 (or at the same height, at least).
	// Otherwise, need to swap them.
	if (n1[2] > n2[2])
	{
		n1.swap(n2);

		// If the nodes position is reversed, we need to reverse the local base as well
		xvec = -xvec;
		yvec = -yvec;
		zvec = -zvec;
	}

	// If the cylinder is above the waterline, then the hydrostatic force is zero
	if (n1[2] >= 0)
	{
		return force;
	}

	// Calculation of the inclination of the cylinder (with respect to the
	// vertical), which is used in the calculation of the center of buoyoancy
	double alpha = acos(arma::dot(zvec, arma::vec::fixed<3> {0, 0, 1})); // zvec and {0, 0, 1} are both unit vectors
	double tanAlpha{ 0 };

	// Check if the angle is 90 degrees
	if (std::abs(alpha - arma::datum::pi / 2) > datum::eps)
	{
		tanAlpha = tan(alpha);
	}
	else
	{
		tanAlpha = arma::datum::inf;
	}

	// Length of the cylinder
	double L{ 0 };

	double xb{ 0 };
	double yb{ 0 };
	double zb{ 0 };

	// If the cylinder is completely submerged, then the center of volume is at the center of the cylinder
	if (n2[2] < 0)
	{
		L = norm(n2 - n1);
		xb = 0;
		yb = 0;
		zb = L / 2;
	}

	else if (is_finite(tanAlpha)) // otherwise, if the cylinder is not horizontal, the formulas for an inclined cylinder are used
	{
		// If only one of the nodes is below the water line, the coordinates of the other node
		// are changed by those of the intersection between the cylinder axis and the static
		// water line (defined by z_global = 0)
		n2 = n1 + (std::abs(0 - n1[2]) / (n2[2] - n1[2])) * norm(n2 - n1) * zvec;
		L = norm(n2 - n1);

		xb = tanAlpha * pow(D / 2, 2) / (4 * L);
		yb = 0;
		zb = (pow(tanAlpha*D / 2, 2) + 4 * pow(L, 2)) / (8 * L);
	}

	else // if the cylinder is horizontal and not completely submerged, forces and moments are equal to zero
	{
		return force;
	}

	// Vector Xb - n1(i.e., vector between the center of buoyancy and the
	// bottom node of the cylinder) written in the global coordinate frame
	vec::fixed<3> Xb_global = xb * xvec + yb * yvec + zb * zvec;

	// Volume of fluid displaced by the cylinder
	double Vol = arma::datum::pi * pow(D / 2, 2) * L;

	// Calculation of hydrostatic force and moment
	force[2] = rho * g * Vol; // Fx = Fy = 0 and Fz = Buoyancy force
	force.rows(3, 5) = cross(Xb_global, force.rows(0, 2));

	// The moment was calculated with relation to n1, which may be different from node1.
	// We need to change the fulcrum to node1
	force.rows(3, 5) = force.rows(3, 5) + cross(n1 - node1Pos(), force.rows(0, 2));

	return force;
}




// The vectors passed as reference are used to return the different components of the hydrodynamic force acting on the cylinder,
// without the need of calling three different methods for each component of the hydrodynamic force
vec::fixed<6> MorisonCirc::hydrodynamicForce(const ENVIR &envir, const int hydroMode, const mat::fixed<3, 3> &rotat,
	const vec::fixed<3> &refPt, const vec::fixed<3> &refPt_sd,
	vec::fixed<6> &force_inertia, vec::fixed<6> &force_drag, vec::fixed<6> &force_froudeKrylov,
	vec::fixed<6> &force_inertia_2nd_part1, vec::fixed<6> &force_inertia_2nd_part2,
	vec::fixed<6> &force_inertia_2nd_part3, vec::fixed<6> &force_inertia_2nd_part4,
	vec::fixed<6> &force_inertia_2nd_part5) const
{
	// Forces and moments acting at the Morison Element
	vec::fixed<6> force(fill::zeros);

	// Make sure that the force components that are passed as reference are set to zero
	force_inertia.zeros();
	force_drag.zeros();
	force_froudeKrylov.zeros();
	force_inertia_2nd_part1.zeros();
	force_inertia_2nd_part2.zeros();
	force_inertia_2nd_part3.zeros();
	force_inertia_2nd_part4.zeros();
	force_inertia_2nd_part5.zeros();

	// Use a more friendly notation
	double D = m_diam;
	double Cd = m_CD;
	double Cm = m_CM;
	double Cd_V = m_axialCD;
	double Ca_V = m_axialCa;
	double rho = envir.watDensity();
	double botDiam = m_botDiam;
	double topDiam = m_topDiam;

	// Nodes position and vectors of the local coordinate system
	vec::fixed<3> n1 = node1Pos();
	vec::fixed<3> n2 = node2Pos();
	vec::fixed<3> xvec(fill::zeros);
	vec::fixed<3> yvec(fill::zeros);
	vec::fixed<3> zvec(fill::zeros);
	MorisonCirc::make_local_base(xvec, yvec, zvec);

	// Same thing, but considering only the mean and slow drift displacements of the FOWT.
	// These ones are used to evaluate forces due to second order quantities (second order potential,
	// quadratic drag, etc).
	vec::fixed<3> n1_sd = node1Pos_sd();
	vec::fixed<3> n2_sd = node2Pos_sd();
	vec::fixed<3> xvec_sd(fill::zeros);
	vec::fixed<3> yvec_sd(fill::zeros);
	vec::fixed<3> zvec_sd(fill::zeros);
	MorisonCirc::make_local_base_sd(xvec_sd, yvec_sd, zvec_sd);

	// Velocity and acceleration of the cylinder nodes
	vec::fixed<3> v1 = node1Vel();
	vec::fixed<3> v2 = node2Vel();
	vec::fixed<3> a1 = node1AccCentrip();
	vec::fixed<3> a2 = node2AccCentrip();

	// Bottom and top diameter
	botDiam = m_botDiam;
	topDiam = m_topDiam;

	// If only first-order forces are going to be calculated, consider the fixed (or slow) position.
	if (hydroMode == 1)
	{
		n1 = n1_sd;
		n2 = n2_sd;
		xvec = xvec_sd;
		yvec = yvec_sd;
		zvec = zvec_sd;
	}

	// Make sure that node1 is below node2 (or at the same height, at least).
	// Otherwise, need to swap them. This is useful for a check below that speeds up the calculation
	// by ignoring the nodes above the free surface.
	//
	// There are some (rare) cases where this procedure causes n2_sd to be BELOW n1_sd, while n1 is
	// guaranteed to be below n2, but this is OK because it is taken in account when the aforementioned check is done.
	if (n1[2] > n2[2])
	{
		n1.swap(n2);
		n1_sd.swap(n2_sd);
		v1.swap(v2);
		a1.swap(a2);

		// If the nodes position is reversed, we need to reverse the local base as well
		xvec = -xvec;
		yvec = -yvec;
		zvec = -zvec;

		xvec_sd = -xvec_sd;
		yvec_sd = -yvec_sd;
		zvec_sd = -zvec_sd;

		// Bottom diameter and top diameter must be swapped as well
		botDiam = m_topDiam;
		topDiam = m_botDiam;
	}

	double L = arma::norm(n2 - n1, 2); // Total cylinder length
	double eta = 0; // Wave elevation above each integration node. Useful for Wheeler stretching method.
	double zwl = 0;
	vec::fixed<3> intersectWL(fill::zeros);
	if (envir.waveStret() == 2)
	{
		intersectWL = findIntersectWL(envir);

		// If no intersection with the WL is found, continue considering zwl == 0
		if (intersectWL.is_finite())
		{
			zwl = intersectWL.at(2);
		}
	}
	
	// Initialization of some variables
	vec::fixed<3> n_ii(fill::zeros); // Coordinates of the integration point
	vec::fixed<3> n_ii_sd(fill::zeros); // Coordinates of the integration point - considering the body fixed at the initial position
	vec::fixed<3> vel_ii(fill::zeros); // Velocity of the integration point
	vec::fixed<3> acc_ii(fill::zeros); // Acceleration of the integration point - Only centripetal part

	vec::fixed<3> u1(fill::zeros); // Fluid velocity at the integration point
	vec::fixed<3> du1dt(fill::zeros); // Components of fluid acceleration at the integration point
	vec::fixed<3> du1dx(fill::zeros);
	vec::fixed<3> du1dy(fill::zeros);
	vec::fixed<3> du1dz(fill::zeros);
	vec::fixed<3> du2dt(fill::zeros);
	vec::fixed<3> a_c(fill::zeros); // Convective acceleration of the fluid
	vec::fixed<3> a_a(fill::zeros); // Axial-divergence acceleration of the fluid
	vec::fixed<3> a_r(fill::zeros); // Fluid acceleration associated with body rotation

	// Forces acting at the integration point and moment (with relation to n1) due to the force acting at the integration point
	vec::fixed<3> force_inertia_ii(fill::zeros); // Inertial component
	vec::fixed<3> moment_inertia_ii(fill::zeros);

	vec::fixed<3> force_drag_ii(fill::zeros); // Drag component
	vec::fixed<3> moment_drag_ii(fill::zeros);

	vec::fixed<3> force_inertia_2nd_part1_ii(fill::zeros); // Inertial component - Second order - Part that is due to the second-order difference-frequency potential
	vec::fixed<3> moment_inertia_2nd_part1_ii(fill::zeros);
	vec::fixed<3> force_inertia_2nd_part3_ii(fill::zeros); // Inertial component - Second order - Part that is due to the convective acceleration
	vec::fixed<3> moment_inertia_2nd_part3_ii(fill::zeros);
	vec::fixed<3> force_inertia_2nd_part4_ii(fill::zeros); // Inertial component - Second order - Part that is due to the axial-divergence acceleration
	vec::fixed<3> moment_inertia_2nd_part4_ii(fill::zeros);
	vec::fixed<3> force_inertia_2nd_part5_ii(fill::zeros); // Inertial component - Second order - Part that is due to body rotation
	vec::fixed<3> moment_inertia_2nd_part5_ii(fill::zeros);

	// Relative distance between the integration point and the bottom node
	double lambda{ 0 };

	// Component of the velocity and acceleration that is parallel to the axis of the cylinder
	vec::fixed<3> v_axial = arma::dot(v1, zvec_sd) * zvec_sd; // Velocities are included in second order terms only, hence are calculate at the sd position
	vec::fixed<3> a_axial = arma::dot(a1, zvec) * zvec; // The situation is the opposite for the acceleration

	/*=================================
		Forces along the length of the cylinder - Considering the instantaneous position
	==================================*/
	double Lw = L;
	int ncyl = m_numIntPoints;
	if (n2.at(2) > zwl)
	{
		Lw = L * (zwl - n1.at(2)) / (n2.at(2) - n1.at(2));
		ncyl = static_cast<int>(std::ceil(Lw / L * (m_numIntPoints - 1))) + 1; // Number of points to discretize the part of the cylinder that is below the water		
	}
	if (ncyl % 2 == 0)
	{
		ncyl += 1;
	}
	double dL = Lw / (ncyl - 1); // length of each interval between points

	for (int ii = 1; ii <= ncyl; ++ii)
	{
		n_ii = n1 + dL * (ii - 1) * zvec; // Coordinates of the integration point

		if (n_ii[2] >= zwl && ii == ncyl)
		{
			n_ii[2] = zwl;
		}

		if (envir.waveStret() == 2 && hydroMode == 2)
		{
			eta = envir.waveElev(n_ii.at(0), n_ii.at(1));
		}

		// Acceleration of the integration point
		lambda = norm(n_ii_sd - n1_sd, 2) / L;
		acc_ii = a1 + lambda * (a2 - a1);

		// Component of the acceleration of the integration point that is perpendicular to the axis of the cylinder
		acc_ii -= a_axial;

		// Component of the fluid acceleration at the integration point that is perpendicular to the axis of the cylinder.
		// Written in the GLOBAL reference frame.
		du1dt = envir.du1dt(n_ii, eta);  // Wheeler stretching method requires 'eta' as input
		du1dt = du1dt - arma::dot(du1dt, zvec) * zvec;

		lambda = norm(n_ii - n1, 2) / L;
		
		// Force due to first-order acceleration integrated considering the instantaneous position of the cylinder
		force_inertia_ii = datum::pi * D*D / 4. * rho * (Cm * du1dt - (Cm - 1) * acc_ii);
		moment_inertia_ii = cross(n_ii - refPt, force_inertia_ii);

		// Integrate the forces along the cylinder using Simpson's Rule
		if (ii == 1 || ii == ncyl)
		{
			force_inertia += (dL / 3.0) * join_cols(force_inertia_ii, moment_inertia_ii);
		}
		else if (ii % 2 == 0)
		{
			force_inertia += (4 * dL / 3.0) * join_cols(force_inertia_ii, moment_inertia_ii);
		}
		else
		{
			force_inertia += (2 * dL / 3.0) * join_cols(force_inertia_ii, moment_inertia_ii);
		}
	}

	/*=================================
		Forces along the length of the cylinder - Considering the slow (or fixed) position
	==================================*/	
	Lw = L;
	ncyl = m_numIntPoints;
	if (n2_sd.at(2) > 0)
	{
		Lw = L * (0 - n1_sd.at(2)) / (n2_sd.at(2) - n1_sd.at(2));
		ncyl = static_cast<int>(std::ceil(Lw / L * (m_numIntPoints - 1))) + 1; // Number of points to discretize the part of the cylinder that is below the water		
	}
	if (ncyl % 2 == 0)
	{
		ncyl += 1;
	}
	dL = Lw / (ncyl - 1); // length of each interval between points
	
	for (int ii = 1; ii <= ncyl; ++ii)
	{
		n_ii_sd = n1_sd + dL * (ii - 1) * zvec_sd;

		if (n_ii_sd[2] >= 0 && ii == ncyl)
		{
			n_ii_sd[2] = 0;
		}

		if (envir.waveStret() == 2 && hydroMode == 2)
		{
			eta = envir.waveElev(n_ii_sd.at(0), n_ii_sd.at(1));
		}
				
		// Velocity of the integration point
		lambda = norm(n_ii_sd - n1_sd, 2) / L;
		vel_ii = v1 + lambda * (v2 - v1);

		// Component of the velocity of the integration point that is perpendicular to the axis of the cylinder
		vel_ii -= v_axial;

		// Fluid velocity at the integration point.		
		u1 = envir.u1(n_ii_sd, eta);
		u1 -= arma::dot(u1, zvec_sd) * zvec_sd;
	
		// Quadratic drag force.
		// Integrated considering the fixed (or slow) position of the cylinder.
		if (n_ii_sd[2] <= 0)
		{
			force_drag_ii = 0.5 * rho * Cd * D * norm(u1 - vel_ii, 2) * (u1 - vel_ii);
		}
		else
		{
			force_drag_ii.zeros();
		}
		moment_drag_ii = cross(n_ii_sd - refPt_sd, force_drag_ii);

		// If required, calculate the other second-order inertial forces,
		// which are integrated considering the fixed (or slow) position of the cylinder.
		if (hydroMode == 2)
		{
			if (n_ii_sd[2] <= 0)
			{
				// 1st component: Force due to the second-order potential
				du2dt = envir.du2dt(n_ii_sd);
				du2dt -= arma::dot(du2dt, zvec_sd) * zvec_sd;
				force_inertia_2nd_part1_ii = (datum::pi * D*D / 4.) * rho * Cm * du2dt;

				// 2nd component: Due to the wave elevation.
				// Calculated after this loop of integration along the cylinder length if Taylor stretching is used.
				// Otherwise, this component is zero and its effects are included in the integration of the forces due to the first-order fluid acceleration.

				// 3rd component: Force due to convective acceleration
				u1 = envir.u1(n_ii_sd, eta);
				du1dx = envir.du1dx(n_ii_sd, eta);
				du1dy = envir.du1dy(n_ii_sd, eta);
				du1dz = envir.du1dz(n_ii_sd, eta);

				a_c.at(0) = u1.at(0) * du1dx.at(0) + u1.at(1) * du1dy.at(0) + u1.at(2) * du1dz.at(0);
				a_c.at(1) = u1.at(0) * du1dx.at(1) + u1.at(1) * du1dy.at(1) + u1.at(2) * du1dz.at(1);
				a_c.at(2) = u1.at(0) * du1dx.at(2) + u1.at(1) * du1dy.at(2) + u1.at(2) * du1dz.at(2);
				a_c -= arma::dot(a_c, zvec_sd) * zvec_sd;

				force_inertia_2nd_part3_ii = (datum::pi * D*D / 4.) * rho * Cm * a_c;

				// 4th component: Force due to axial-divergence acceleration
				double dwdz = arma::dot(du1dx, zvec_sd) * zvec_sd.at(0) + arma::dot(du1dy, zvec_sd) * zvec_sd.at(1) + arma::dot(du1dz, zvec_sd) * zvec_sd.at(2);
				a_a = dwdz * (u1 - arma::dot(u1, zvec_sd)*zvec_sd - vel_ii); // vel_ii was already projected in the direction perpendicular to the cylinder
				force_inertia_2nd_part4_ii = (datum::pi * D*D / 4.) * rho * (Cm - 1) * a_a;

				// 5th component: Force due to cylinder rotation				
				a_r = 2*arma::dot(u1-v_axial, zvec_sd) * (1 / L) * (arma::dot(v2 - v1, xvec_sd) * yvec_sd + arma::dot(v2 - v1, yvec_sd) * xvec_sd);
				force_inertia_2nd_part5_ii = -(datum::pi * D*D / 4.) * rho * (Cm - 1) * a_r;
			}
			else
			{
				force_inertia_2nd_part1_ii.zeros();
				force_inertia_2nd_part3_ii.zeros();
				force_inertia_2nd_part4_ii.zeros();
				force_inertia_2nd_part5_ii.zeros();
			}
			moment_inertia_2nd_part1_ii = cross(n_ii_sd - refPt_sd, force_inertia_2nd_part1_ii);
			moment_inertia_2nd_part3_ii = cross(n_ii_sd - refPt_sd, force_inertia_2nd_part3_ii);
			moment_inertia_2nd_part4_ii = cross(n_ii_sd - refPt_sd, force_inertia_2nd_part4_ii);
			moment_inertia_2nd_part5_ii = cross(n_ii_sd - refPt_sd, force_inertia_2nd_part5_ii);
		}

		// Integrate the forces along the cylinder using Simpson's Rule
		if (ii == 1 || ii == ncyl)
		{
			force_drag += (dL / 3.0) * join_cols(force_drag_ii, moment_drag_ii);
			force_inertia_2nd_part1 += (dL / 3.0) * join_cols(force_inertia_2nd_part1_ii, moment_inertia_2nd_part1_ii);
			force_inertia_2nd_part3 += (dL / 3.0) * join_cols(force_inertia_2nd_part3_ii, moment_inertia_2nd_part3_ii);
			force_inertia_2nd_part4 += (dL / 3.0) * join_cols(force_inertia_2nd_part4_ii, moment_inertia_2nd_part4_ii);
			force_inertia_2nd_part5 += (dL / 3.0) * join_cols(force_inertia_2nd_part5_ii, moment_inertia_2nd_part5_ii);
		}
		else if (ii % 2 == 0)
		{
			force_drag += (4 * dL / 3.0) * join_cols(force_drag_ii, moment_drag_ii);
			force_inertia_2nd_part1 += (4 * dL / 3.0) * join_cols(force_inertia_2nd_part1_ii, moment_inertia_2nd_part1_ii);
			force_inertia_2nd_part3 += (4 * dL / 3.0) * join_cols(force_inertia_2nd_part3_ii, moment_inertia_2nd_part3_ii);
			force_inertia_2nd_part4 += (4 * dL / 3.0) * join_cols(force_inertia_2nd_part4_ii, moment_inertia_2nd_part4_ii);
			force_inertia_2nd_part5 += (4 * dL / 3.0) * join_cols(force_inertia_2nd_part5_ii, moment_inertia_2nd_part5_ii);
		}
		else
		{
			force_drag += (2 * dL / 3.0) * join_cols(force_drag_ii, moment_drag_ii);
			force_inertia_2nd_part1 += (2 * dL / 3.0) * join_cols(force_inertia_2nd_part1_ii, moment_inertia_2nd_part1_ii);
			force_inertia_2nd_part3 += (2 * dL / 3.0) * join_cols(force_inertia_2nd_part3_ii, moment_inertia_2nd_part3_ii);
			force_inertia_2nd_part4 += (2 * dL / 3.0) * join_cols(force_inertia_2nd_part4_ii, moment_inertia_2nd_part4_ii);
			force_inertia_2nd_part5 += (2 * dL / 3.0) * join_cols(force_inertia_2nd_part5_ii, moment_inertia_2nd_part5_ii);
		}
	}

	// 2nd component of the second order force: Force due to the wave elevation.
	// It is computed here only if Taylor series for wave stretching was chosen and
	// if the cylinder intersects the water line.
	// If waveStret == 0, this effect is not included in the analysis, 
	// and if waveStret > 1, this force component is included in force_inertia.
	if (hydroMode == 2 && envir.waveStret() == 1 && (n2.at(2)*n1.at(2) < 0))
	{
		n_ii = (n2 - n1) * (0 - n1.at(2)) / (n2.at(2) - n1.at(2)) + n1; // Coordinates of the intersection with the still water line;				
		n_ii.at(2) = 0; // Since envir.du1dt returns 0 for z > 0, this line is necessary to make sure that the z coordinate of n_ii is exactly 0, and not slightly above due to roundoff errors.
		du1dt = envir.du1dt(n_ii, 0);
		du1dt -= arma::dot(du1dt, zvec_sd) * zvec_sd;
		eta = envir.waveElev(n_ii.at(0), n_ii.at(1));
		force_inertia_2nd_part2.rows(0, 2) = (datum::pi * D*D / 4.) * rho * Cm * du1dt * eta;
		double R_ii = norm(n_ii - n1, 2) + eta / 2;
		force_inertia_2nd_part2.rows(3, 5) = cross(n_ii + eta / 2 * zvec - refPt, force_inertia_2nd_part2.rows(0, 2));
	}


	/*=================================
		Forces on the bottom of the cylinder
	==================================*/
	// Component of the fluid velocity/acceleration at the bottom node that is parallel to the cylinder axis
	du1dt = arma::dot(envir.du1dt(n1, 0), zvec) * zvec;
	u1 = arma::dot(envir.u1(n1_sd, 0), zvec_sd) * zvec_sd;	

	// Calculate the force acting on the bottom of the cylinder
	vec::fixed<3> force_inertia_axial = rho * Ca_V * (4 / 3.) * datum::pi * (D*D*D / 8.) * (du1dt - a_axial);
	vec::fixed<3> force_drag_axial = 0.5 * rho * Cd_V * datum::pi * (D*D / 4.) * norm(u1 - v_axial, 2) * (u1 - v_axial);

	force_inertia.rows(0, 2) += force_inertia_axial;
	force_inertia.rows(3, 5) += cross(n1 - refPt, force_inertia_axial);

	force_drag.rows(0, 2) += force_drag_axial;
	force_drag.rows(3, 5) += cross(n1_sd - refPt_sd, force_drag_axial);

	if (m_botPressFlag)
	{
		force_froudeKrylov.rows(0, 2) = datum::pi * (botDiam*botDiam / 4. * envir.wavePressure(n1)
			- (botDiam*botDiam / 4. - topDiam * topDiam / 4.) * envir.wavePressure(n2)) * zvec_sd;

		force_froudeKrylov.rows(3, 5) = cross(n1 - refPt, force_froudeKrylov.rows(0, 2));
	}

	/*
		Total force
	*/
	force = force_inertia + force_drag + force_froudeKrylov +
		force_inertia_2nd_part1 + force_inertia_2nd_part2 +
		force_inertia_2nd_part3 + force_inertia_2nd_part4 + force_inertia_2nd_part5;

	return force;
}


mat::fixed<6, 6> MorisonCirc::addedMass_perp(const double rho, const vec::fixed<3> &refPt, const int hydroMode) const
{
	mat::fixed<6, 6> A(fill::zeros);

	// Use a more friendly notation
	double Lambda = datum::pi * pow(m_diam / 2., 2) * rho * (m_CM - 1);
	double ncyl = m_numIntPoints;

	// Nodes position and vectors of the local coordinate system vectors
	vec::fixed<3> n1 = node1Pos();
	vec::fixed<3> n2 = node2Pos();
	if (hydroMode == 1)
	{
		n1 = node1Pos_sd();
		n2 = node2Pos_sd();
	}
	vec::fixed<3> xvec(fill::zeros);
	vec::fixed<3> yvec(fill::zeros);
	vec::fixed<3> zvec(fill::zeros);
	MorisonCirc::make_local_base_sd(xvec, yvec, zvec);

	// Make sure that node1 is below node2 (or at the same height, at least).
	// Otherwise, need to swap them.
	if (n1[2] > n2[2])
	{
		n1.swap(n2);

		// If the nodes position is reversed, we need to reverse the local base as well
		xvec = -xvec;
		yvec = -yvec;
		zvec = -zvec;
	}


	// Since n2 is above n1, if n1[2] > 0, the cylinder is above the waterline
	if (n1[2] > 0)
	{
		return A;
	}

	// If only one of the nodes is above the water line, the coordinates of the other node
	// are changed by those of the intersection between the cylinder axis and the static
	// water line(defined by z_global = 0)
	if (n2[2] > 0)
	{
		n2 = n1 + (std::abs(0 - n1[2]) / (n2[2] - n1[2])) * norm(n2 - n1) * zvec;
	}

	// Length of the cylinder and of each interval between points
	double L = norm(n2 - n1, 2);
	double dL = norm(n2 - n1, 2) / (ncyl - 1);

	// The purely translational elements of the matrix (Aij for i,j = 1, 2, 3) are integrated analytically
	for (int pp = 0; pp < 3; ++pp)
	{
		for (int qq = pp; qq < 3; ++qq)
		{
			A(pp, qq) = Lambda * L * A_perp(pp, qq, { 0,0,0 }, refPt, xvec, yvec);
		}
	}

	vec::fixed<3> n_ii;
	double step{ 0 }; // Used for Simpson's rule. See below.
	for (int ii = 1; ii <= ncyl; ++ii)
	{
		n_ii = (n2 - n1) * (ii - 1) / (ncyl - 1) + n1; // Coordinates of the integration point

		if (n_ii[2] > 0)
		{
			break; // Since n2 is above n1, if we reach n_ii[2] > 0, all the next n_ii are also above the waterline
		}

		// Integrate along the cylinder using Simpon's rule
		if (ii == 1 || ii == ncyl)
		{
			step = dL / 3.;
		}
		else if (ii % 2 == 0)
		{
			step = 4 * dL / 3.;
		}
		else
		{
			step = 2 * dL / 3.;
		}

		for (int pp = 0; pp < 6; ++pp)
		{
			int q0 = pp;
			if (pp < 3) // The first square of the matrix was already filled above
			{
				q0 = 3;
			}
			for (int qq = q0; qq < 6; ++qq)
			{
				A(pp, qq) += Lambda * step * A_perp(pp, qq, n_ii, refPt, xvec, yvec);
			}
			// (const int ii, const int jj, const vec::fixed<3> &x, const vec::fixed<3> &xG, const vec::fixed<3> &xvec, const vec::fixed<3> &yvec)
		}
	}

	// The matrix is symmetrical. In the lines above, only the upper triangle was filled.
	// In the loop below, the lower triangle is filled with the values from the upper triangle.
	for (int pp = 0; pp < 6; ++pp)
	{
		for (int qq = 0; qq < pp; ++qq)
		{
			A(pp, qq) = A(qq, pp);
		}
	}

	return A;
}


double MorisonCirc::A_perp(const int ii, const int jj, const vec::fixed<3> &x, const vec::fixed<3> &xG, const vec::fixed<3> &xvec, const vec::fixed<3> &yvec) const
{

	if (ii < 3 && jj < 3)
	{
		return (xvec.at(ii) * xvec.at(jj) + yvec.at(ii) * yvec.at(jj));
	}

	if (ii == 3 && jj == 3)
	{
		return (pow(x.at(1) - xG.at(1), 2) * (pow(xvec.at(2), 2) + pow(yvec.at(2), 2))
			+ pow(x.at(2) - xG.at(2), 2) * (pow(xvec.at(1), 2) + pow(yvec.at(1), 2))
			- 2 * (x.at(1) - xG.at(1)) * (x.at(2) - xG.at(2))  * (xvec.at(1) * xvec.at(2) + yvec.at(1) * yvec.at(2))
			);
	}

	if (ii == 4 && jj == 4)
	{
		return (pow(x.at(0) - xG.at(0), 2) * (pow(xvec.at(2), 2) + pow(yvec.at(2), 2))
			+ pow(x.at(2) - xG.at(2), 2) * (pow(xvec.at(0), 2) + pow(yvec.at(0), 2))
			- 2 * (x.at(0) - xG.at(0)) * (x.at(2) - xG.at(2))  * (xvec.at(0) * xvec.at(2) + yvec.at(0) * yvec.at(2))
			);
	}

	if (ii == 5 && jj == 5)
	{
		return (pow(x.at(0) - xG.at(0), 2) * (pow(xvec.at(1), 2) + pow(yvec.at(1), 2))
			+ pow(x.at(1) - xG.at(1), 2) * (pow(xvec.at(0), 2) + pow(yvec.at(0), 2))
			- 2 * (x.at(0) - xG.at(0)) * (x.at(1) - xG.at(1))  * (xvec.at(0) * xvec.at(1) + yvec.at(0) * yvec.at(1))
			);
	}

	if (ii == 0 && jj == 3)
	{
		return ((x.at(1) - xG.at(1)) * (xvec.at(0) * xvec.at(2) + yvec.at(0) * yvec.at(2))
			- (x.at(2) - xG.at(2)) * (xvec.at(0) * xvec.at(1) + yvec.at(0) * yvec.at(1))
			);
	}

	if (ii == 0 && jj == 4)
	{
		return ((x.at(2) - xG.at(2)) * (pow(xvec.at(0), 2) + pow(yvec.at(0), 2))
			- (x.at(0) - xG.at(0)) * (xvec.at(0) * xvec.at(2) + yvec.at(0) * yvec.at(2))
			);
	}

	if (ii == 0 && jj == 5)
	{
		return ((x.at(0) - xG.at(0)) * (xvec.at(0) * xvec.at(1) + yvec.at(0) * yvec.at(1))
			- (x.at(1) - xG.at(1)) * (pow(xvec.at(0), 2) + pow(yvec.at(0), 2))
			);
	}

	if (ii == 1 && jj == 3)
	{
		return ((x.at(1) - xG.at(1)) * (xvec.at(1) * xvec.at(2) + yvec.at(1) * yvec.at(2))
			- (x.at(2) - xG.at(2)) * (pow(xvec.at(1), 2) + pow(yvec.at(1), 2))
			);
	}

	if (ii == 1 && jj == 4)
	{
		return ((x.at(2) - xG.at(2)) * (xvec.at(0) * xvec.at(1) + yvec.at(0) * yvec.at(1))
			- (x.at(0) - xG.at(0)) * (xvec.at(1) * xvec.at(2) + yvec.at(1) * yvec.at(2))
			);
	}

	if (ii == 1 && jj == 5)
	{
		return ((x.at(0) - xG.at(0)) * (pow(xvec.at(1), 2) + pow(yvec.at(1), 2))
			- (x.at(1) - xG.at(1)) * (xvec.at(0) * xvec.at(1) + yvec.at(0) * yvec.at(1))
			);
	}

	if (ii == 2 && jj == 3)
	{
		return ((x.at(1) - xG.at(1)) * (pow(xvec.at(2), 2) + pow(yvec.at(2), 2))
			- (x.at(2) - xG.at(2)) * (xvec.at(1) * xvec.at(2) + yvec.at(1) * yvec.at(2))
			);
	}

	if (ii == 2 && jj == 4)
	{
		return ((x.at(2) - xG.at(2)) * (xvec.at(0) * xvec.at(2) + yvec.at(0) * yvec.at(2))
			- (x.at(0) - xG.at(0)) * (pow(xvec.at(2), 2) + pow(yvec.at(2), 2))
			);
	}

	if (ii == 2 && jj == 5)
	{
		return ((x.at(0) - xG.at(0)) * (xvec.at(1) * xvec.at(2) + yvec.at(1) * yvec.at(2))
			- (x.at(1) - xG.at(1)) * (xvec.at(0) * xvec.at(2) + yvec.at(0) * yvec.at(2))
			);
	}

	if (ii == 3 && jj == 4)
	{
		return (-(x.at(0) - xG.at(0)) * (x.at(1) - xG.at(1)) * (pow(xvec.at(2), 2) + pow(yvec.at(2), 2))
			- (x.at(0) - xG.at(0)) * (x.at(2) - xG.at(2)) * (xvec.at(1) * xvec.at(2) + yvec.at(1) * yvec.at(2))
			+ (x.at(1) - xG.at(1)) * (x.at(2) - xG.at(2)) * (xvec.at(0) * xvec.at(2) + yvec.at(0) * yvec.at(2))
			- pow(x.at(2) - xG.at(2), 2) * (xvec.at(0) * xvec.at(1) + yvec.at(0) * yvec.at(1))
			);
	}

	if (ii == 3 && jj == 5)
	{
		return ((x.at(0) - xG.at(0)) * (x.at(1) - xG.at(1)) * (xvec.at(1) * xvec.at(2) + yvec.at(1) * yvec.at(2))
			- (x.at(0) - xG.at(0)) * (x.at(2) - xG.at(2)) * (pow(xvec.at(1), 2) + pow(yvec.at(1), 2))
			- pow(x.at(1) - xG.at(1), 2) * (xvec.at(0) * xvec.at(2) + yvec.at(0) * yvec.at(2))
			+ (x.at(1) - xG.at(1)) * (x.at(2) - xG.at(2)) * (xvec.at(0) * xvec.at(1) + yvec.at(0) * yvec.at(1))
			);
	}

	if (ii == 4 && jj == 5)
	{
		return (-pow(x.at(0) - xG.at(0), 2) * (xvec.at(1) * xvec.at(2) + yvec.at(1) * yvec.at(2))
			+ (x.at(0) - xG.at(0)) * (x.at(1) - xG.at(1)) * (xvec.at(0) * xvec.at(2) + yvec.at(0) * yvec.at(2))
			+ (x.at(0) - xG.at(0)) * (x.at(1) - xG.at(1)) * (xvec.at(0) * xvec.at(1) + yvec.at(0) * yvec.at(1))
			- (x.at(1) - xG.at(1)) * (x.at(2) - xG.at(2)) * (pow(xvec.at(0), 2) + pow(yvec.at(0), 2))
			);
	}

	return 0;
}


// TODO: depois de debugar direitinho, tirar os bound checks (usar [] ao inves de () pra acessar elementos das matrizes)
mat::fixed<6, 6> MorisonCirc::addedMass_paral(const double rho, const vec::fixed<3> &refPt, const int hydroMode) const
{
	mat::fixed<6, 6> A(fill::zeros);

	// Use a more friendly notation
	double Lambda = rho * m_axialCa * (4 / 3.) * datum::pi * pow(m_diam / 2, 3);
	double ncyl = m_numIntPoints;

	// Nodes position and vectors of the local coordinate system vectors
	vec::fixed<3> n1 = node1Pos();
	vec::fixed<3> n2 = node2Pos();
	if (hydroMode == 1)
	{
		n1 = node1Pos_sd();
		n2 = node2Pos_sd();
	}
	vec::fixed<3> xvec(fill::zeros);
	vec::fixed<3> yvec(fill::zeros);
	vec::fixed<3> zvec(fill::zeros);
	MorisonCirc::make_local_base_sd(xvec, yvec, zvec);

	// Center of Gravity
	double xG = refPt[0];
	double yG = refPt[1];
	double zG = refPt[2];

	// Make sure that node1 is below node2 (or at the same height, at least).
	// Otherwise, need to swap them.
	if (n1[2] > n2[2])
	{
		n1.swap(n2);

		// If the nodes position is reversed, we need to reverse the local base as well
		xvec = -xvec;
		yvec = -yvec;
		zvec = -zvec;
	}

	// Vectors of the global coordinate system - arranjed, they result in an eye matrix
	mat::fixed<3, 3> globalBase(fill::eye);

	// Coordinates of the analysed point, i.e. the bottom node
	double x_ii = n1[0];
	double y_ii = n1[1];
	double z_ii = n1[2];

	if (z_ii > 0)
	{
		return A;
	}

	for (int pp = 0; pp < 3; ++pp)
	{
		for (int qq = pp; qq < 3; ++qq)
		{
			A(pp, qq) = Lambda * arma::dot(zvec, globalBase.col(pp))*arma::dot(zvec, globalBase.col(qq));
		}
	}

	A(3, 3) += Lambda *
		(
			pow(y_ii - yG, 2) * pow(arma::dot(zvec, globalBase.col(2)), 2)
			+ pow(z_ii - zG, 2) * pow(arma::dot(zvec, globalBase.col(1)), 2)
			- 2 * (y_ii - yG) * (z_ii - zG) * arma::dot(zvec, globalBase.col(1)) * arma::dot(zvec, globalBase.col(2))
			);

	A(4, 4) += Lambda *
		(
			pow(x_ii - xG, 2) * pow(arma::dot(zvec, globalBase.col(2)), 2)
			+ pow(z_ii - zG, 2) * pow(arma::dot(zvec, globalBase.col(0)), 2)
			- 2 * (x_ii - xG) * (z_ii - zG)  * arma::dot(zvec, globalBase.col(0)) * arma::dot(zvec, globalBase.col(2))
			);

	A(5, 5) += Lambda *
		(
			pow(x_ii - xG, 2) * pow(arma::dot(zvec, globalBase.col(1)), 2)
			+ pow(y_ii - yG, 2) * pow(arma::dot(zvec, globalBase.col(0)), 2)
			- 2 * (x_ii - xG) * (y_ii - yG)  * arma::dot(zvec, globalBase.col(0)) * arma::dot(zvec, globalBase.col(1))
			);

	A(0, 3) += Lambda *
		(
		(y_ii - yG) * arma::dot(zvec, globalBase.col(0))*arma::dot(zvec, globalBase.col(2))
			- (z_ii - zG) * arma::dot(zvec, globalBase.col(0))*arma::dot(zvec, globalBase.col(1))
			);

	A(0, 4) += Lambda *
		(
		(z_ii - zG) * pow(arma::dot(zvec, globalBase.col(0)), 2)
			- (x_ii - xG) * arma::dot(zvec, globalBase.col(0))*arma::dot(zvec, globalBase.col(2))
			);

	A(0, 5) += Lambda *
		(
		(x_ii - xG) * arma::dot(zvec, globalBase.col(0))*arma::dot(zvec, globalBase.col(1))
			- (y_ii - yG) * pow(arma::dot(zvec, globalBase.col(0)), 2)
			);

	A(1, 3) += Lambda *
		(
		(y_ii - yG) * arma::dot(zvec, globalBase.col(1))*arma::dot(zvec, globalBase.col(2))
			- (z_ii - zG) * pow(arma::dot(zvec, globalBase.col(1)), 2)
			);

	A(1, 4) += Lambda *
		(
		(z_ii - zG) * arma::dot(zvec, globalBase.col(0))*arma::dot(zvec, globalBase.col(1))
			- (x_ii - xG) * arma::dot(zvec, globalBase.col(1))*arma::dot(zvec, globalBase.col(2))
			);

	A(1, 5) += Lambda *
		(
		(x_ii - xG) * pow(arma::dot(zvec, globalBase.col(1)), 2)
			- (y_ii - yG) * arma::dot(zvec, globalBase.col(0))*arma::dot(zvec, globalBase.col(1))
			);

	A(2, 3) += Lambda *
		(
		(y_ii - yG) * pow(arma::dot(zvec, globalBase.col(2)), 2)
			- (z_ii - zG) * arma::dot(zvec, globalBase.col(1))*arma::dot(zvec, globalBase.col(2))
			);

	A(2, 4) += Lambda *
		(
		(z_ii - zG) * arma::dot(zvec, globalBase.col(0))*arma::dot(zvec, globalBase.col(2))
			- (x_ii - xG) * pow(arma::dot(zvec, globalBase.col(2)), 2)
			);

	A(2, 5) += Lambda *
		(
		(x_ii - xG) * arma::dot(zvec, globalBase.col(1))*arma::dot(zvec, globalBase.col(2))
			- (y_ii - yG) * arma::dot(zvec, globalBase.col(0))*arma::dot(zvec, globalBase.col(2))
			);

	A(3, 4) += Lambda *
		(
			-(x_ii - xG) * (y_ii - yG) * pow(arma::dot(zvec, globalBase.col(2)), 2)
			- (x_ii - xG) * (z_ii - zG) * arma::dot(zvec, globalBase.col(1))*arma::dot(zvec, globalBase.col(2))
			+ (y_ii - yG) * (z_ii - zG) * arma::dot(zvec, globalBase.col(0))*arma::dot(zvec, globalBase.col(2))
			- pow(z_ii - zG, 2) * arma::dot(zvec, globalBase.col(0))*arma::dot(zvec, globalBase.col(1))
			);

	A(3, 5) += Lambda *
		(
		(x_ii - xG) * (y_ii - yG) * arma::dot(zvec, globalBase.col(1))*arma::dot(zvec, globalBase.col(2))
			- (x_ii - xG) * (z_ii - zG) * pow(arma::dot(zvec, globalBase.col(1)), 2)
			- pow(y_ii - yG, 2) * arma::dot(zvec, globalBase.col(0))*arma::dot(zvec, globalBase.col(2))
			+ (y_ii - yG) * (z_ii - zG) * arma::dot(zvec, globalBase.col(0))*arma::dot(zvec, globalBase.col(1))
			);

	A(4, 5) += Lambda *
		(
			-pow(x_ii - xG, 2) * arma::dot(zvec, globalBase.col(1))*arma::dot(zvec, globalBase.col(2))
			+ (x_ii - xG) * (y_ii - yG) * arma::dot(zvec, globalBase.col(0))*arma::dot(zvec, globalBase.col(2))
			+ (x_ii - xG) * (y_ii - yG) * arma::dot(zvec, globalBase.col(0))*arma::dot(zvec, globalBase.col(1))
			- (y_ii - yG) * (z_ii - zG) * pow(arma::dot(zvec, globalBase.col(0)), 2)
			);


	// The matrix is symmetrical. In the lines above, only the upper triangle was filled.
	// In the loop below, the lower triangle is filled with the values from the upper triangle.
	for (int pp = 0; pp < 6; ++pp)
	{
		for (int qq = 0; qq < pp; ++qq)
		{
			A(pp, qq) = A(qq, pp);
		}
	}

	return A;
}




/*****************************************************
	Printing
*****************************************************/
std::string MorisonCirc::print() const
{
	std::string output = "";

	output = output + "CoG_2_node1:\t(" + std::to_string(m_cog2node1(0)) + ", " + std::to_string(m_cog2node1(1)) + ", " + std::to_string(m_cog2node1(2)) + ")\n";
	output = output + "CoG_2_node1:\t(" + std::to_string(m_cog2node2(0)) + ", " + std::to_string(m_cog2node2(1)) + ", " + std::to_string(m_cog2node2(2)) + ")\n";
	output = output + "Diameter:\t" + std::to_string(m_diam) + '\n';
	output = output + "Drag Coeff.:\t" + std::to_string(m_CD) + '\n';
	output = output + "Inert. Coeff.:\t" + std::to_string(m_CM) + '\n';
	output = output + "Numb. of Int. Points:\t" + std::to_string(m_numIntPoints) + '\n';
	output = output + "Bottom diameter:\t" + std::to_string(m_botDiam) + '\n';
	output = output + "Top diameter:\t" + std::to_string(m_topDiam) + '\n';
	output = output + "Axial CD:\t" + std::to_string(m_axialCD) + '\n';
	output = output + "Axial Ca:\t" + std::to_string(m_axialCa) + '\n';
	output = output + "Bot. Press. Flag.:\t" + std::to_string(m_botPressFlag) + '\n';

	return output;
}


/*****************************************************
	Clone for creating copies of the Morison Element
*****************************************************/
MorisonCirc* MorisonCirc::clone() const
{
	return (new MorisonCirc(*this));
}
