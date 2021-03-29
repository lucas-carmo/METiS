#include "MorisonCirc.h"
#include <cmath>


using namespace arma;

/*****************************************************
	Constructors
*****************************************************/
MorisonCirc::MorisonCirc(const vec &node1Pos, const vec &node2Pos, const vec &cog, const int numIntPoints,
	const bool botPressFlag, const double axialCD_1, const double axialCa_1, const double axialCD_2, const double axialCa_2,
	const double diam, const double CD, double CM)
	: MorisonElement(node1Pos, node2Pos, cog, numIntPoints, botPressFlag, axialCD_1, axialCa_1, axialCD_2, axialCa_2),
	m_diam(diam), m_CD(CD), m_CM(CM)
{
	make_local_base(m_xvec_t0, m_yvec_t0, m_zvec_t0, node1Pos, node2Pos);

	m_xvec_sd = m_xvec_t0;
	m_yvec_sd = m_yvec_t0;
	m_zvec_sd = m_zvec_t0;

	m_xvec = m_xvec_t0;
	m_yvec = m_yvec_t0;
	m_zvec = m_zvec_t0;
}



/*****************************************************
	Forces acting on the Morison Element
*****************************************************/
// Make the local base of the cylinder at t = 0
void MorisonCirc::make_local_base(arma::vec::fixed<3> &xvec, arma::vec::fixed<3> &yvec, arma::vec::fixed<3> &zvec, const arma::vec::fixed<3> &n1, const arma::vec::fixed<3> &n2) const
{
	xvec.zeros();
	yvec.zeros();
	zvec = (n2 - n1) / arma::norm(n2 - n1, 2);

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

	vec::fixed<3> n1 = node1Pos();
	vec::fixed<3> n2 = node2Pos();
	vec::fixed<3> xvec(m_xvec), yvec(m_yvec), zvec(m_zvec);

	// Make sure that node1 is below node2 (or at the same height, at least).
	// Otherwise, need to swap them.
	if (n1[2] > n2[2])
	{
		n1.swap(n2);
		xvec *= -1; yvec *= -1; zvec *= -1;
	}	
	
	
	// If the cylinder is above the waterline, then the hydrostatic force is zero
	if (n1[2] >= 0)
	{
		return force;
	}

	// Calculation of the inclination of the cylinder (with respect to the
	// vertical), which is used in the calculation of the center of buoyoancy
	double alpha = std::acos(arma::dot(zvec, arma::vec::fixed<3> {0, 0, 1})); // zvec and {0, 0, 1} are both unit vectors
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
vec::fixed<6> MorisonCirc::hydrodynamicForce(const ENVIR &envir, const int hydroMode, const vec::fixed<3> &refPt, const vec::fixed<3> &refPt_sd,
	vec::fixed<6> &force_drag, vec::fixed<6> &force_1, vec::fixed<6> &force_2,
	vec::fixed<6> &force_3, vec::fixed<6> &force_4, vec::fixed<6> &force_eta, vec::fixed<6> &force_rem,
	vec::fixed<6> &force_drag_ext, vec::fixed<6> &force_1_ext, vec::fixed<6> &force_2_ext,
	vec::fixed<6> &force_3_ext, vec::fixed<6> &force_rem_ext) const
{
	// Forces and moments acting at the Morison Element
	vec::fixed<6> force(fill::zeros);

	// Make sure that the force components that are passed as reference are set to zero
	force_drag.zeros();
	force_1.zeros();
	force_2.zeros();
	force_3.zeros();
	force_4.zeros();
	force_eta.zeros();
	force_rem.zeros();

	force_drag_ext.zeros();
	force_1_ext.zeros();
	force_2_ext.zeros();
	force_3_ext.zeros();
	force_rem_ext.zeros();

	// Use a more friendly notation
	double D = m_diam;
	double Cd = m_CD;
	double Cm = m_CM;
	double CdV_1 = m_axialCD_1;
	double CaV_1 = m_axialCa_1;
	double CdV_2 = m_axialCD_2;
	double CaV_2 = m_axialCa_2;
	double rho = envir.watDensity();

	// Nodes position and vectors of the local coordinate system
	vec::fixed<3> n1 = node1Pos();
	vec::fixed<3> n2 = node2Pos();
	vec::fixed<3> xvec = m_xvec;
	vec::fixed<3> yvec = m_yvec;
	vec::fixed<3> zvec = m_zvec;

	// Same thing, but considering only the mean and slow drift displacements of the FOWT.
	// These ones are used to evaluate forces due to second order quantities (second order potential,
	// quadratic drag, etc).
	vec::fixed<3> n1_sd = node1Pos_sd();
	vec::fixed<3> n2_sd = node2Pos_sd();
	vec::fixed<3> xvec_sd = m_xvec_sd;
	vec::fixed<3> yvec_sd = m_yvec_sd;
	vec::fixed<3> zvec_sd = m_zvec_sd;

	// Velocity and acceleration of the cylinder nodes
	vec::fixed<3> v1 = node1Vel();
	vec::fixed<3> v2 = node2Vel();
	vec::fixed<3> a1 = node1AccCentrip();
	vec::fixed<3> a2 = node2AccCentrip();

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

		// Axial coefficients must be swapped as well
		CdV_1 = m_axialCD_2;
		CaV_1 = m_axialCa_2;
		CdV_2 = m_axialCD_1;
		CaV_2 = m_axialCa_1;
	}

	double L = arma::norm(n2 - n1, 2); // Total cylinder length
	double eta = 0; // Wave elevation above each integration node. Useful for Wheeler stretching method.
	double zwl = 0;
	if (hydroMode == 2 && envir.waveStret() == 2 && m_intersectWL.is_finite())
	{
		zwl = m_intersectWL.at(2);
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
	// TODO: quase certeza que todas essas variaveis moment sao desnecessarias
	vec::fixed<3> force_drag_ii(fill::zeros); // Drag component
	vec::fixed<3> moment_drag_ii(fill::zeros);
	vec::fixed<3> force_1_ii(fill::zeros); // First part - Forces due to 1st order potential
	vec::fixed<3> moment_1_ii(fill::zeros);
	vec::fixed<3> force_2_ii(fill::zeros); // Second part - Forces due to 2nd order potential
	vec::fixed<3> moment_2_ii(fill::zeros);
	vec::fixed<3> force_3_ii(fill::zeros); // Third part - Forces due to the convective acceleration
	vec::fixed<3> moment_3_ii(fill::zeros);
	vec::fixed<3> force_4_ii(fill::zeros); // Fourth part - Forces due to the axial-divergence acceleration
	vec::fixed<3> moment_4_ii(fill::zeros);
	vec::fixed<3> force_rem_ii(fill::zeros); // Remaining components
	vec::fixed<3> moment_rem_ii(fill::zeros);

	// Relative distance between the integration point and the bottom node
	double lambda{ 0 };	
	
	// Component of the velocity and acceleration that is parallel to the axis of the cylinder
	vec::fixed<3> v_axial = arma::dot(v1, zvec_sd) * zvec_sd; // Velocities are included in second order terms only, hence are calculate at the sd position
	
	// Useful auxilliary variables to avoid recalculating things
	vec::fixed<3> aux_force(fill::zeros);
	double aux{ 0 };

	/*=================================
		Forces along the length of the cylinder - Considering the instantaneous position
	==================================*/

	// Force due to first-order acceleration integrated considering the instantaneous position of the cylinder
	// Integrated analytically

	vec::fixed<6> auxForce1 = morisonForce_inertia(envir, hydroMode);
	force_1 += auxForce1;
	force_1.rows(3, 5) += cross(n1-refPt, auxForce1.rows(0, 2));
	

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

	aux = datum::pi * D*D / 4. * rho;
	if (hydroMode > 1)
	{
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

			// Component of the acceleration of the integration point
			// that is perpendicular to the axis of the cylinder
			lambda = norm(n_ii - n1, 2) / L;
			acc_ii = a1 + lambda * (a2 - a1);
			acc_ii = arma::dot(acc_ii, xvec)*xvec + arma::dot(acc_ii, yvec)*yvec;

			force_rem_ii = -aux * (Cm - 1) * acc_ii;
			moment_rem_ii = cross(n_ii - refPt, force_rem_ii);

			// Integrate the forces along the cylinder using Simpson's Rule
			if (ii == 1 || ii == ncyl)
			{
				force_rem += (dL / 3.0) * join_cols(force_rem_ii, moment_rem_ii);
			}
			else if (ii % 2 == 0)
			{
				force_rem += (4 * dL / 3.0) * join_cols(force_rem_ii, moment_rem_ii);
			}
			else
			{
				force_rem += (2 * dL / 3.0) * join_cols(force_rem_ii, moment_rem_ii);
			}
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

	force_rem_ii.zeros();
	moment_rem_ii.zeros();	
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

		// Fluid velocity at the integration point.		
		u1 = envir.u1(n_ii_sd, eta);
		u1 -= arma::dot(u1, zvec_sd) * zvec_sd;

		// Quadratic drag force.
		force_drag_ii = 0.5 * rho * Cd * D * norm(u1 - (vel_ii-v_axial), 2) * (u1 - (vel_ii - v_axial));

		// If required, calculate the other second-order forces,
		if (hydroMode == 2)
		{			
			// 2nd component: Force due to the second-order potential
			du2dt = envir.du2dt(n_ii_sd);
			du2dt -= arma::dot(du2dt, zvec_sd) * zvec_sd;
			force_2_ii = aux * Cm * du2dt; // aux == datum::pi * D*D / 4. * rho;

			// 3rd component: Force due to convective acceleration
			u1 = envir.u1(n_ii_sd, eta);
			du1dx = envir.du1dx(n_ii_sd, eta);
			du1dy = envir.du1dy(n_ii_sd, eta);
			du1dz = envir.du1dz(n_ii_sd, eta);

			a_c.at(0) = u1.at(0) * du1dx.at(0) + u1.at(1) * du1dy.at(0) + u1.at(2) * du1dz.at(0);
			a_c.at(1) = u1.at(0) * du1dx.at(1) + u1.at(1) * du1dy.at(1) + u1.at(2) * du1dz.at(1);
			a_c.at(2) = u1.at(0) * du1dx.at(2) + u1.at(1) * du1dy.at(2) + u1.at(2) * du1dz.at(2);
			a_c -= arma::dot(a_c, zvec_sd) * zvec_sd;

			force_3_ii = aux * Cm * a_c;

			// 4th component: Force due to axial-divergence acceleration
			double dwdz = arma::dot(du1dx, zvec_sd) * zvec_sd.at(0) + arma::dot(du1dy, zvec_sd) * zvec_sd.at(1) + arma::dot(du1dz, zvec_sd) * zvec_sd.at(2);
			a_a = dwdz * (arma::dot(u1 - vel_ii, xvec_sd)*xvec_sd + arma::dot(u1 - vel_ii, yvec_sd)*yvec_sd);
			force_4_ii = aux * (Cm - 1) * a_a;

			// Add to remaining forces: Force due to cylinder rotation				
			a_r = 2 * arma::dot(u1 - vel_ii, zvec_sd) * (1 / L) * (arma::dot(v2 - v1, xvec_sd) * xvec_sd + arma::dot(v2 - v1, yvec_sd) * yvec_sd);
			force_rem_ii = -aux * (Cm - 1) * a_r;
		}

		moment_drag_ii = cross(n_ii_sd - refPt_sd, force_drag_ii);
		moment_2_ii = cross(n_ii_sd - refPt_sd, force_2_ii);
		moment_3_ii = cross(n_ii_sd - refPt_sd, force_3_ii);
		moment_4_ii = cross(n_ii_sd - refPt_sd, force_4_ii);
		moment_rem_ii = cross(n_ii_sd - refPt_sd, force_rem_ii);

		// Integrate the forces along the cylinder using Simpson's Rule
		if (ii == 1 || ii == ncyl)
		{
			force_drag += (dL / 3.0) * join_cols(force_drag_ii, moment_drag_ii);
			force_2 += (dL / 3.0) * join_cols(force_2_ii, moment_2_ii);
			force_3 += (dL / 3.0) * join_cols(force_3_ii, moment_3_ii);
			force_4 += (dL / 3.0) * join_cols(force_4_ii, moment_4_ii);
			force_rem += (dL / 3.0) * join_cols(force_rem_ii, moment_rem_ii);
		}
		else if (ii % 2 == 0)
		{
			force_drag += (4 * dL / 3.0) * join_cols(force_drag_ii, moment_drag_ii);
			force_2 += (4 * dL / 3.0) * join_cols(force_2_ii, moment_2_ii);
			force_3 += (4 * dL / 3.0) * join_cols(force_3_ii, moment_3_ii);
			force_4 += (4 * dL / 3.0) * join_cols(force_4_ii, moment_4_ii);
			force_rem += (4 * dL / 3.0) * join_cols(force_rem_ii, moment_rem_ii);
		}
		else
		{
			force_drag += (2 * dL / 3.0) * join_cols(force_drag_ii, moment_drag_ii);
			force_2 += (2 * dL / 3.0) * join_cols(force_2_ii, moment_2_ii);
			force_3 += (2 * dL / 3.0) * join_cols(force_3_ii, moment_3_ii);
			force_4 += (2 * dL / 3.0) * join_cols(force_4_ii, moment_4_ii);
			force_rem += (2 * dL / 3.0) * join_cols(force_rem_ii, moment_rem_ii);
		}
	}

	// Eta: Force due to the wave elevation.
	// It is computed here only if Taylor series for wave stretching was chosen and
	// if the cylinder intersects the water line.
	// If waveStret == 0, this effect is not included in the analysis, 
	// and if waveStret > 1, this force component is included in force_1.
	if (hydroMode == 2 && envir.waveStret() == 1 && (n2.at(2)*n1.at(2) < 0))
	{
		n_ii = (n2 - n1) * (0 - n1.at(2)) / (n2.at(2) - n1.at(2)) + n1; // Coordinates of the intersection with the still water line;				
		n_ii.at(2) = 0; // Since envir.du1dt returns 0 for z > 0, this line is necessary to make sure that the z coordinate of n_ii is exactly 0, and not slightly above due to roundoff errors.
		du1dt = envir.du1dt(n_ii, 0);
		du1dt = arma::dot(du1dt, xvec) * xvec + arma::dot(du1dt, yvec) * yvec;
		eta = envir.waveElev(n_ii.at(0), n_ii.at(1));
		force_eta.rows(0, 2) = (datum::pi * D*D / 4.) * rho * Cm * du1dt * eta;
		force_eta.rows(3, 5) = cross(n_ii - refPt, force_eta.rows(0, 2));
	}


	/*=================================
		Forces at the extremities of the cylinder
	==================================*/	
	// Calculate the force acting on the bottom of the cylinder
	if (n1.at(2) <= zwl)
	{
		// Kinematics
		du1dt = arma::dot(envir.du1dt(n1, 0), zvec) * zvec;
		u1 = envir.u1(n1_sd, 0);
		a_c.zeros();
		du2dt.zeros();
		if (hydroMode == 2)
		{
			du1dx = envir.du1dx(n1_sd, 0);
			du1dy = envir.du1dy(n1_sd, 0);
			du1dz = envir.du1dz(n1_sd, 0);

			a_c.at(0) = u1.at(0) * du1dx.at(0) + u1.at(1) * du1dy.at(0) + u1.at(2) * du1dz.at(0);
			a_c.at(1) = u1.at(0) * du1dx.at(1) + u1.at(1) * du1dy.at(1) + u1.at(2) * du1dz.at(1);
			a_c.at(2) = u1.at(0) * du1dx.at(2) + u1.at(1) * du1dy.at(2) + u1.at(2) * du1dz.at(2);
			a_c = arma::dot(a_c, zvec_sd) * zvec_sd;

			du2dt = envir.du2dt(n1_sd);
			du2dt = arma::dot(du2dt, zvec_sd) * zvec_sd;
		}

		// Drag force
		u1 = arma::dot(u1, zvec_sd) * zvec_sd;		
		aux_force = 0.5 * rho * CdV_1 * datum::pi * (D*D / 4.) * norm(u1 - v_axial, 2) * (u1 - v_axial);
		force_drag_ext.rows(0, 2) += aux_force;
		force_drag_ext.rows(3, 5) += cross(n1_sd - refPt_sd, aux_force);

		// Inertial force - pt1
		aux = rho * CaV_1 * (4 / 3.) * datum::pi * (D*D*D / 8.);
		aux_force = aux * du1dt;
		force_1_ext.rows(0, 2) += aux_force;
		force_1_ext.rows(3, 5) += cross(n1 - refPt, aux_force);

		// Inertial force - pt2
		aux_force = aux * du2dt;
		force_2_ext.rows(0, 2) += aux_force;
		force_2_ext.rows(3, 5) += cross(n1_sd - refPt_sd, aux_force);

		// Inertial force - pt3
		aux_force = aux * a_c;
		force_3_ext.rows(0, 2) += aux_force;
		force_3_ext.rows(3, 5) += cross(n1_sd - refPt_sd, aux_force);

		// Remaining force components
		aux_force = -aux * arma::dot(a1, zvec_sd) * zvec_sd;
		force_rem_ext.rows(0, 2) += aux_force;
		force_rem_ext.rows(3, 5) += cross(n1 - refPt, aux_force);

		if (m_botPressFlag)
		{
			// Force due to 1st order Froude-Krylov pressure is part of the first component
			aux = datum::pi * 0.25 * D*D;
			aux_force = aux * envir.wavePressure(n1) * zvec;
			force_1_ext.rows(0, 2) += aux_force;
			force_1_ext.rows(3, 5) += cross(n1 - refPt, aux_force);

			if (hydroMode == 2)
			{
				// Inertial force - pt2 - 2nd order Froude-Krylov pressure
				aux_force = aux * envir.wavePressure_2ndOrd(n1_sd) * zvec_sd;
				force_2_ext.rows(0, 2) += aux_force;
				force_2_ext.rows(3, 5) += cross(n1_sd - refPt_sd, aux_force);

				// Inertial force - pt3 - Quadratic pressure drop from Rainey's formulation
				aux_force = -0.5 * aux * rho * (Cm - 1) * (pow(arma::dot(envir.u1(n1_sd, 0) - v1, xvec_sd), 2) + pow(arma::dot(envir.u1(n1_sd, 0) - v1, yvec_sd), 2)) * zvec_sd;
				force_3_ext.rows(0, 2) += aux_force;
				force_3_ext.rows(3, 5) += cross(n1_sd - refPt_sd, aux_force);

				// Inertial force - Remaining - Point load from Rainey's formulation that results in a Munk moment
				aux_force = aux * (Cm-1) * cdot(u1 - v_axial, zvec_sd) * ((envir.u1(n1_sd, 0) - u1) - (v1 - v_axial));
				force_rem_ext.rows(0, 2) += aux_force;
				force_rem_ext.rows(3, 5) += cross(n1_sd - refPt_sd, aux_force);
			}
		}
	}

	// Calculate the force acting on the top of the cylinder
	if (n2.at(2) <= zwl)
	{
		// Kinematics
		du1dt = arma::dot(envir.du1dt(n2, 0), zvec) * zvec;
		u1 = envir.u1(n2_sd, 0);
		a_c.zeros();
		du2dt.zeros();
		if (hydroMode == 2)
		{			
			du1dx = envir.du1dx(n2_sd, 0);
			du1dy = envir.du1dy(n2_sd, 0);
			du1dz = envir.du1dz(n2_sd, 0);

			a_c.at(0) = u1.at(0) * du1dx.at(0) + u1.at(1) * du1dy.at(0) + u1.at(2) * du1dz.at(0);
			a_c.at(1) = u1.at(0) * du1dx.at(1) + u1.at(1) * du1dy.at(1) + u1.at(2) * du1dz.at(1);
			a_c.at(2) = u1.at(0) * du1dx.at(2) + u1.at(1) * du1dy.at(2) + u1.at(2) * du1dz.at(2);
			a_c = arma::dot(a_c, zvec_sd) * zvec_sd;

			du2dt = envir.du2dt(n2_sd);
			du2dt = arma::dot(du2dt, zvec_sd) * zvec_sd;
		}

		// Drag force
		u1 = arma::dot(u1, zvec_sd) * zvec_sd;
		aux_force = 0.5 * rho * CdV_2 * datum::pi * (D*D / 4.) * norm(u1 - v_axial, 2) * (u1 - v_axial);
		force_drag_ext.rows(0, 2) += aux_force;
		force_drag_ext.rows(3, 5) += cross(n2_sd - refPt_sd, aux_force);

		// Inertial force - pt1
		aux = rho * CaV_2 * (4 / 3.) * datum::pi * (D*D*D / 8.);
		aux_force = aux * du1dt;
		force_1_ext.rows(0, 2) += aux_force;
		force_1_ext.rows(3, 5) += cross(n2 - refPt, aux_force);

		// Inertial force - pt2
		aux_force = aux * du2dt;
		force_2_ext.rows(0, 2) += aux_force;
		force_2_ext.rows(3, 5) += cross(n2_sd - refPt_sd, aux_force);

		// Inertial force - pt3
		aux_force = aux * a_c;
		force_3_ext.rows(0, 2) += aux_force;
		force_3_ext.rows(3, 5) += cross(n2_sd - refPt_sd, aux_force);

		// Remaining force components
		aux_force = -aux * arma::dot(a2, zvec_sd) * zvec;
		force_rem_ext.rows(0, 2) += aux_force;
		force_rem_ext.rows(3, 5) += cross(n2 - refPt, aux_force);

		if (m_botPressFlag)
		{
			// Force due to 1st order Froude-Krylov pressure is part of the first component
			aux = datum::pi * 0.25 * D*D;
			aux_force = -aux * envir.wavePressure(n2) * zvec;
			force_1_ext.rows(0, 2) += aux_force;
			force_1_ext.rows(3, 5) += cross(n2 - refPt, aux_force);

			if (hydroMode == 2)
			{
				// Inertial force - pt2 - 2nd order Froude-Krylov pressure
				aux_force = -aux * envir.wavePressure_2ndOrd(n2_sd) * zvec_sd;
				force_2_ext.rows(0, 2) += aux_force;
				force_2_ext.rows(3, 5) += cross(n2_sd - refPt_sd, aux_force);

				// Inertial force - pt3 - Quadratic pressure drop from Rainey's formulation
				aux_force = 0.5 * aux * rho * (Cm - 1) * pow(norm((envir.u1(n2_sd, 0) - u1) - (v1 - v_axial)), 2) * zvec_sd;
				force_3_ext.rows(0, 2) += aux_force;
				force_3_ext.rows(3, 5) += cross(n2_sd - refPt_sd, aux_force);

				// Inertial force - Remaining - Point load from Rainey's formulation that results in a Munk moment
				aux_force = -aux * (Cm - 1) * cdot(u1 - v_axial, zvec_sd) * ((envir.u1(n2_sd, 0) - u1) - (v1 - v_axial));
				force_rem_ext.rows(0, 2) += aux_force;
				force_rem_ext.rows(3, 5) += cross(n2_sd - refPt_sd, aux_force);
			}
		}
	}

	/*
		Total force
	*/
	force = force_drag + force_1 + force_2 + force_3 + force_4 + force_eta + force_rem +
		force_drag_ext + force_1_ext + force_2_ext + force_3_ext + force_rem_ext;

	return force;
}


vec::fixed<6> MorisonCirc::morisonForce_inertia(const ENVIR &envir, const int hydroMode) const
{
	vec::fixed<6> force(fill::zeros);

	double t = envir.time();
	double rho = envir.watDensity();
	double h = envir.watDepth();
	double g = envir.gravity();
	double R = m_diam/2.;
	double Cm = m_CM;
	double pi = arma::datum::pi;
	vec::fixed<3> Zvec{ 0, 0, 1 }; // Vectors of the global basis
	vec::fixed<3> Xvec{ 1, 0, 0 };


	vec::fixed<3> n1{ node1Pos() }, n2{ node2Pos() };
	vec::fixed<3> xvec{ m_xvec }, yvec{ m_yvec }, zvec{ m_zvec };

	// If only first-order forces are going to be calculated, consider the slow.
	if (hydroMode == 1)
	{
		n1 = node1Pos_sd(); n2 = node2Pos_sd();
		xvec = m_xvec_sd; yvec = m_yvec_sd; zvec = m_zvec_sd;
	}

	// Make sure that node1 is below node2 (or at the same height, at least).
	// Otherwise, need to swap them.
	if (n1[2] > n2[2])
	{
		n1.swap(n2);
		xvec *= -1; yvec *= -1; zvec *= -1;
	}

	if (n2[2] > 0)
	{
		n2 = n1 + (abs(0 - n1(2)) / (n2(2) - n1(2))) * arma::norm(n2 - n1) * zvec;
	}

	double x1(n1[0]), y1(n1[1]), z1(n1[2]);
	double x2(n2[0]), y2(n2[1]), z2(n2[2]);
	double L = norm(n2 - n1);

	// Calculation of the inclination of the cylinder with respect to the vertical
	double alpha = std::acos(arma::dot(zvec, arma::vec::fixed<3> {0, 0, 1})); // zvec and {0, 0, 1} are both unit vectors
	double tanAlpha{ 0 }, cosAlpha{ 0 }, sinAlpha{ 0 };

	// Check if the angle is 90 degrees
	if (std::abs(alpha - arma::datum::pi / 2) > datum::eps)
	{
		tanAlpha = tan(alpha);
		cosAlpha = cos(alpha);
		sinAlpha = sin(alpha);
	}
	else
	{
		tanAlpha = arma::datum::inf;
		sinAlpha = 1;
	}

	// Calculate the angle that the cylinder makes with the global x axis
	// Different method than the one used for tanAlpha because in here we need the sign of the angle
	double psi{ 0 };
	if (std::abs(alpha) > datum::eps)
	{
		vec::fixed<3> zproj = zvec - dot(zvec, Zvec)*Zvec;
		zproj = zproj / norm(zproj);

		psi = -atan2(dot(cross(zproj, Xvec), Zvec), dot(zproj, Xvec));
	}


	for (unsigned int ii = 0; ii < envir.numberOfWaveComponents(); ++ii)
	{
		const Wave &wave(envir.getWave(ii));

		double w{ wave.angFreq() }, k{ wave.waveNumber() }, A{wave.amp()};
		double cosBeta{wave.cosBeta()}, sinBeta{wave.sinBeta()};
		double beta = wave.direction() * pi / 180.;
		double phase = wave.phase() * pi / 180.;
		double cosT(cos(-w * t + phase)), sinT(sin(-w * t + phase));		

		// Avoid recalculting cossine and sine functions
		double cosBetaPsi = cos(beta - psi); 
		double sinBetaPsi = sin(beta - psi);

		// The integration is different for a horizontal cylinder, as the integration variable
		// along the cylinder is not related to the vertical position. This special case is treated separately.
		double factor1_x(0), factor2_x(0), factor1_z(0), factor2_z(0), mult_h(0), mult_v(0);

		/* 
			FORCES
		*/
		if (arma::is_finite(tanAlpha))
		{
			// Contribution of the horizontal acceleration
			factor1_x = k * sin(k*cosBeta*x1 + k * sinBeta*y1 + (z2 - z1)*k*tanAlpha*cosBetaPsi) * sinh(k*(z2 + h))
				- k * tanAlpha*cosBetaPsi*cos(k*cosBeta*x1 + k*sinBeta*y1 + (z2 - z1)*k*tanAlpha*cosBetaPsi) * cosh(k*(z2 + h));
			factor1_x = factor1_x - k * sin(k*cosBeta*x1 + k*sinBeta*y1) * sinh(k*(z1 + h))
				+ k * tanAlpha*cosBetaPsi*cos(k*cosBeta*x1 + k*sinBeta*y1) * cosh(k*(z1 + h));
			factor1_x = factor1_x / (k*tanAlpha*cosBetaPsi*k*tanAlpha*cosBetaPsi + k*k);

			factor2_x = k * cos(k*cosBeta*x1 + k*sinBeta*y1 + (z2 - z1)*k*tanAlpha*cosBetaPsi) * sinh(k*(z2 + h))
				+ k * tanAlpha*cosBetaPsi*sin(k*cosBeta*x1 + k * sinBeta*y1 + (z2 - z1)*k*tanAlpha*cosBetaPsi) * cosh(k*(z2 + h));
			factor2_x = factor2_x - k * cos(k*cosBeta*x1 + k*sinBeta*y1) * sinh(k*(z1 + h))
				- k * tanAlpha*cosBetaPsi*sin(k*cosBeta*x1 + k * sinBeta*y1) * cosh(k*(z1 + h));
			factor2_x = factor2_x / (k*tanAlpha*cosBetaPsi*k*tanAlpha*cosBetaPsi + k*k);

			mult_h = factor1_x * cosT + factor2_x * sinT;

			// Contribution of the vertical acceleration
			factor1_z = k * cos(k*cosBeta*x1 + k * sinBeta*y1 + (z2 - z1)*k*tanAlpha*cosBetaPsi) * cosh(k*(z2 + h))
				+ k * tanAlpha*cosBetaPsi*sin(k*cosBeta*x1 + k * sinBeta*y1 + (z2 - z1)*k*tanAlpha*cosBetaPsi) * sinh(k*(z2 + h));
			factor1_z = factor1_z - k * cos(k*cosBeta*x1 + k * sinBeta*y1) * cosh(k*(z1 + h))
				- k * tanAlpha*cosBetaPsi*sin(k*cosBeta*x1 + k * sinBeta*y1) * sinh(k*(z1 + h));
			factor1_z = factor1_z / (k*tanAlpha*cosBetaPsi*k*tanAlpha*cosBetaPsi + k*k);

			factor2_z = k * sin(k*cosBeta*x1 + k * sinBeta*y1 + (z2 - z1)*k*tanAlpha*cosBetaPsi) * cosh(k*(z2 + h))
				- k * tanAlpha*cosBetaPsi*cos(k*cosBeta*x1 + k*sinBeta*y1 + (z2 - z1)*k*tanAlpha*cosBetaPsi) * sinh(k*(z2 + h));
			factor2_z = factor2_z - k * sin(k*cosBeta*x1 + k*sinBeta*y1) * cosh(k*(z1 + h))
				+ k * tanAlpha*cosBetaPsi*cos(k*cosBeta*x1 + k*sinBeta*y1) * sinh(k*(z1 + h));
			factor2_z = factor2_z / (k*tanAlpha*cosBetaPsi*k*tanAlpha*cosBetaPsi + k*k);

			mult_v = factor1_z * cosT - factor2_z * sinT;

			// This is due to the change of integration variable.
			// We know cosAlpha !=0 because horizontal cylinders are treated separately
			mult_h = mult_h / cosAlpha;
			mult_v = mult_v / cosAlpha;
		}
		else
		{			
			// Contribution of the horizontal acceleration
			factor1_x = -cos(k*cosBeta*x2 + k*sinBeta*y2) + cos(k*cosBeta*x1 + k*sinBeta*y1);
			factor1_x = L * factor1_x / (k*cosBeta*(x2 - x1) + k*sinBeta*(y2 - y1));

			factor2_x = sin(k*cosBeta*x2 + k*sinBeta*y2) - sin(k*cosBeta*x1 + k*sinBeta*y1);
			factor2_x = L * factor2_x / (k*cosBeta*(x2 - x1) + k*sinBeta*(y2 - y1));

			mult_h = factor1_x * cos(-w * t + phase) + factor2_x * sin(-w * t + phase);

			// Contribution of the vertical acceleration
			mult_v = factor2_x * cos(-w * t + phase) - factor1_x * sin(-w * t + phase);

			// There is also the term related to the z position 
			mult_h = mult_h * cosh(k*(z1 + h));
			mult_v = mult_v * sinh(k*(z1 + h));
		}
			
		// Forces along the local x and y axes
		double fx = w * w * A * (cosBetaPsi * cosAlpha * mult_h + sinAlpha * mult_v) / sinh(k*h);
		double fy = w * w * A * sinBetaPsi * mult_h / sinh(k*h);


		/*
			MOMENTS
		*/
		if (arma::is_finite(tanAlpha))
		{
			double b = k * cosBeta*x1 + k * sinBeta*y1;
			double q = k * tanAlpha*cosBetaPsi;

			// Contribution of the horizontal acceleration
			factor1_x = k * sinh(k*(z2 + h)) * ((k*k + q*q)*(z2 - z1)*sin(q*(z2 - z1) + b) + 2*q*cos(q*(z2 - z1) + b))
				+ cosh(k*(z2 + h)) * (-q * (k*k + q*q)*(z2 - z1)*cos(q*(z2 - z1) + b) + (q*q - k*k)*sin(q*(z2 - z1) + b));
			factor1_x = factor1_x - k * sinh(k*(z1 + h)) * 2 * q*cos(b)
				- cosh(k*(z1 + h)) * (q*q - k*k)*sin(b);
			factor1_x = factor1_x / (q*q + k*k) / (q*q + k * k);

			factor2_x = k * sinh(k*(z2 + h)) * ((k*k + q*q)*(z2 - z1)*cos(q*(z2 - z1) + b) - 2*q*sin(q*(z2 - z1) + b))
				+ cosh(k*(z2 + h)) * (q*(k*k+ q*q)*(z2 - z1)*sin(q*(z2 - z1) + b) + (q*q - k*k)*cos(q*(z2 - z1) + b));
			factor2_x = factor2_x + k * sinh(k*(z1 + h)) * 2*q*sin(b)
				- cosh(k*(z1 + h)) * (q*q - k*k)*cos(b);
			factor2_x = factor2_x / (q*q + k * k) / (q*q + k * k);

			mult_h = factor1_x * cos(-w * t + phase) + factor2_x * sin(-w * t + phase);

			// Contribution of the vertical acceleration
			factor1_z = k * cosh(k*(z2 + h)) * ((k*k + q*q)*(z2 - z1)*cos(q*(z2 - z1) + b) - 2*q*sin(q*(z2 - z1) + b))
				+ sinh(k*(z2 + h)) * (q*(k*k + q*q)*(z2 - z1)*sin(q*(z2 - z1) + b) + (q*q - k*k)*cos(q*(z2 - z1) + b));
			factor1_z = factor1_z + k * cosh(k*(z1 + h)) * 2*q*sin(b)
				- sinh(k*(z1 + h)) * (q*q - k*k)*cos(b);
			factor1_z = factor1_z / (q*q + k * k) / (q*q + k * k);

			factor2_z = k * cosh(k*(z2 + h)) * ((k*k + q*q)*(z2 - z1)*sin(q*(z2 - z1) + b) + 2*q*cos(q*(z2 - z1) + b))
				+ sinh(k*(z2 + h)) * (-q * (k*k + q*q)*(z2 - z1)*cos(q*(z2 - z1) + b) + (q*q - k*k)*sin(q*(z2 - z1) + b));
			factor2_z = factor2_z - k * cosh(k*(z1 + h)) * 2 * q*cos(b)
				- sinh(k*(z1 + h)) * (q*q - k*k)*sin(b);
			factor2_z = factor2_z / (q*q + k * k) / (q*q + k * k);

			mult_v = factor1_z * cosT - factor2_z * sinT;

			// This is due to the change of integration variable
			mult_h = mult_h / cosAlpha / cosAlpha;
			mult_v = mult_v / cosAlpha/ cosAlpha;
		}
		else
		{
			double aux = (k*cosBeta*(x2 - x1) + k*sinBeta*(y2 - y1)) / L;

			// Contribution of the horizontal acceleration
			factor1_x = -aux * L*cos(k*cosBeta*x2 + k*sinBeta*y2) + sin(k*cosBeta*x2 + k*sinBeta*y2)
				- sin(k*cosBeta*x1 + k * sinBeta*y1);
			factor1_x = factor1_x / (aux*aux);

			factor2_x = aux * L*sin(k*cosBeta*x2 + k*sinBeta*y2) + cos(k*cosBeta*x2 + k*sinBeta*y2)
				- cos(k*cosBeta*x1 + k*sinBeta*y1);
			factor2_x = factor2_x / (aux*aux);;

			mult_h = factor1_x * cos(-w * t + phase) + factor2_x * sin(-w * t + phase);

			// Contribution of the vertical acceleration
			mult_v = factor2_x * cos(-w * t + phase) - factor1_x * sin(-w * t + phase);

			// There is also the term related to the z position 
			mult_h = mult_h * cosh(k*(z1 + h));
			mult_v = mult_v * sinh(k*(z1 + h));
		}

		// Forces along the local x and y axes
		double mx = -w * w * A * sinBetaPsi * mult_h / sinh(k*h);
		double my = w * w * A * (cosBetaPsi * cosAlpha * mult_h + sinAlpha * mult_v) / sinh(k*h);		


		force.at(0) += fx * xvec(0) + fy * yvec(0);
		force.at(1) += fx * xvec(1) + fy * yvec(1);
		force.at(2) += fx * xvec(2) + fy * yvec(2);

		force.at(3) += mx * xvec(0) + my * yvec(0);
		force.at(4) += mx * xvec(1) + my * yvec(1);
		force.at(5) += mx * xvec(2) + my * yvec(2);
	}




	return rho * pi* R * R * Cm * force * envir.ramp();
}


mat::fixed<6, 6> MorisonCirc::addedMass_perp(const double rho, const vec::fixed<3> &refPt, const int hydroMode) const
{
	mat::fixed<6, 6> A(fill::zeros);

	// Use a more friendly notation
	double Lambda = datum::pi * pow(m_diam / 2., 2) * rho * (m_CM - 1);
	int ncyl = m_numIntPoints;
	
	// Get vertical coordinate of the intersection with the waterline
	double zwl = 0;
	if (hydroMode == 2 && m_intersectWL.is_finite())
	{
		zwl = m_intersectWL.at(2);
	}

	// Nodes position and vectors of the local coordinate system vectors
	vec::fixed<3> n1 = node1Pos();
	vec::fixed<3> n2 = node2Pos();
	vec::fixed<3> xvec = m_xvec; // xvec and yvec are used to project the acceleration. They are analogous to the normal vector, which should consider the slow or fixed position
	vec::fixed<3> yvec = m_yvec;
	vec::fixed<3> zvec = m_zvec; // While zvec is used only to evaluate the nodes position, and hence should be first-order position
	if (hydroMode < 2)
	{
		n1 = node1Pos_sd();
		n2 = node2Pos_sd();
		xvec = m_xvec_sd;
		yvec = m_yvec_sd;
		zvec = m_zvec_sd;
	}

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

	// Since n2 is above n1, if n1[2] > zwl, the cylinder is above the waterline
	if (n1[2] > zwl)
	{
		return A;
	}

	// If only one of the nodes is above the water line, the coordinates of the other node
	// are changed by those of the intersection between the cylinder axis and the static
	// water line(defined by z_global = 0)
	if (n2[2] > zwl)
	{
		n2 = n1 + (std::abs(zwl - n1[2]) / (n2[2] - n1[2])) * norm(n2 - n1) * zvec;
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
	//omp_set_num_threads(4);
	//#pragma omp parallel for
	for (int ii = 1; ii <= ncyl; ++ii)
	{
		n_ii = (n2 - n1) * (ii - 1) / (ncyl - 1) + n1; // Coordinates of the integration point

		if (n_ii[2] > zwl)
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

	// Nodes position and vectors of the local coordinate system vectors
	vec::fixed<3> n1 = node1Pos();
	vec::fixed<3> n2 = node2Pos();
	vec::fixed<3> zvec = m_zvec;
	if (hydroMode < 2)
	{
		n1 = node1Pos_sd();
		n2 = node2Pos_sd();
		zvec = m_zvec_sd;
	}

	// Center of Gravity
	double xG = refPt[0];
	double yG = refPt[1];
	double zG = refPt[2];
	
	if (n1[2] < 0)
	{
		for (int pp = 0; pp < 6; ++pp)
		{
			int q0 = pp;
			for (int qq = q0; qq < 6; ++qq)
			{
				A(pp, qq) += rho * m_axialCa_1 * (4 / 3.) * datum::pi * pow(m_diam / 2, 3) * A_paral(pp, qq, n1, refPt, zvec);
			}
		}
	}

	if (n2[2] < 0)
	{
		for (int pp = 0; pp < 6; ++pp)
		{
			int q0 = pp;
			for (int qq = q0; qq < 6; ++qq)
			{
				A(pp, qq) += rho * m_axialCa_2 * (4 / 3.) * datum::pi * pow(m_diam / 2, 3) * A_paral(pp, qq, n2, refPt, zvec);
			}
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

double MorisonCirc::A_paral(const int ii, const int jj, const vec::fixed<3> &x, const vec::fixed<3> &xG, const vec::fixed<3> &zvec) const
{
	if (ii < 3 && jj < 3)
	{
		return zvec(ii) * zvec(jj);
	}	

	if (ii == 3 && jj == 3)
	{
		return (pow(x.at(1) - xG.at(1), 2) * pow(zvec.at(2), 2)
			+ pow(x.at(2) - xG.at(2), 2) * pow(zvec.at(1), 2)
			- 2 * (x.at(1) - xG.at(1)) * (x.at(2) - xG.at(2)) * zvec.at(1) * zvec.at(2)
			);
	}

	if (ii == 4 && jj == 4)
	{
		return (pow(x.at(0) - xG.at(0), 2) * pow(zvec.at(2), 2)
			+ pow(x.at(2) - xG.at(2), 2) * pow(zvec.at(0), 2)
			- 2 * (x.at(0) - xG.at(0)) * (x.at(2) - xG.at(2)) * zvec.at(0) * zvec.at(2)
			);
	}

	if (ii == 5 && jj == 5)
	{
		return (pow(x.at(0) - xG.at(0), 2) * pow(zvec.at(1), 2)
			+ pow(x.at(1) - xG.at(1), 2) * pow(zvec.at(0), 2)
			- 2 * (x.at(0) - xG.at(0)) * (x.at(1) - xG.at(1)) * zvec.at(0) * zvec.at(1)
			);
	}
	
	if (ii == 0 && jj == 3)
	{
		return ((x.at(1) - xG.at(1)) * zvec.at(0) * zvec.at(2) 
			- (x.at(2) - xG.at(2)) * zvec.at(0) * zvec.at(1));
	}

	if (ii == 0 && jj == 4)
	{
		return ((x.at(2) - xG.at(2)) * pow(zvec.at(0), 2) 
			- (x.at(0) - xG.at(0)) * zvec.at(0) * zvec.at(2));
	}

	if (ii == 0 && jj == 5)
	{
		return ((x.at(0) - xG.at(0)) * zvec.at(0) * zvec.at(1)
			- (x.at(1) - xG.at(1)) * pow(zvec.at(0), 2));
	}	

	if (ii == 1 && jj == 3)
	{
		return ((x.at(1) - xG.at(1)) * zvec.at(1) * zvec.at(2)
			- (x.at(2) - xG.at(2)) * pow(zvec.at(1), 2));
	}

	if (ii == 1 && jj == 4)
	{
		return ((x.at(2) - xG.at(2)) * zvec.at(0) * zvec.at(1)
			- (x.at(0) - xG.at(0)) * zvec.at(1) * zvec.at(2));
	}

	if (ii == 1 && jj == 5)
	{
		return ((x.at(0) - xG.at(0)) * pow(zvec.at(1), 2)
			- (x.at(1) - xG.at(1)) * zvec.at(0) * zvec.at(1));
	}	

	if (ii == 2 && jj == 3)
	{
		return ((x.at(1) - xG.at(1)) * pow(zvec.at(2), 2)
			- (x.at(2) - xG.at(2)) * zvec.at(1) * zvec.at(2));
	}

	if (ii == 2 && jj == 4)
	{
		return ((x.at(2) - xG.at(2)) * zvec.at(0) * zvec.at(2)
			- (x.at(0) - xG.at(0)) * pow(zvec.at(2), 2));
	}

	if (ii == 2 && jj == 5)
	{
		return ((x.at(0) - xG.at(0)) * zvec.at(1) * zvec.at(2)
			- (x.at(1) - xG.at(1)) * zvec.at(0) * zvec.at(2));
	}		

	if (ii == 3 && jj == 4)
	{
		return (-(x.at(0) - xG.at(0)) * (x.at(1) - xG.at(1)) * pow(zvec.at(2), 2)
			- (x.at(0) - xG.at(0)) * (x.at(2) - xG.at(2)) * zvec.at(1) * zvec.at(2)
			+ (x.at(1) - xG.at(1)) * (x.at(2) - xG.at(2)) * zvec.at(0) * zvec.at(2)
			- pow(x.at(2) - xG.at(2), 2) * zvec.at(0) * zvec.at(1));
	}

	if (ii == 3 && jj == 5)
	{
		return ((x.at(0) - xG.at(0)) * (x.at(1) - xG.at(1)) * zvec.at(1) * zvec.at(2)
			- (x.at(0) - xG.at(0)) * (x.at(2) - xG.at(2)) * pow(zvec.at(1), 2)
			- pow(x.at(1) - xG.at(1), 2) * zvec.at(0) * zvec.at(2)
			+ (x.at(1) - xG.at(1)) * (x.at(2) - xG.at(2)) * zvec.at(0) * zvec.at(1));
	}

	if (ii == 4 && jj == 5)
	{
		return (-pow(x.at(0) - xG.at(0), 2) * zvec.at(1) * zvec.at(2)
			+ (x.at(0) - xG.at(0)) * (x.at(1) - xG.at(1)) * zvec.at(0) * zvec.at(2)
			+ (x.at(0) - xG.at(0)) * (x.at(1) - xG.at(1)) * zvec.at(0) * zvec.at(1)
			- (x.at(1) - xG.at(1)) * (x.at(2) - xG.at(2)) * pow(zvec.at(0), 2));
	}

	return 0;
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
	output = output + "Axial CD - Node 1:\t" + std::to_string(m_axialCD_1) + '\n';
	output = output + "Axial Ca - Node 1:\t" + std::to_string(m_axialCa_1) + '\n';
	output = output + "Axial CD - Node 2:\t" + std::to_string(m_axialCD_2) + '\n';
	output = output + "Axial Ca - Node 2:\t" + std::to_string(m_axialCa_2) + '\n';
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
