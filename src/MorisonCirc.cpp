#include "MorisonCirc.h"
#include "auxFunctions.h"
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

// TODO: There is a lot of code repetition in this function that could be avoided.
void MorisonCirc::evaluateQuantitiesAtBegin(const ENVIR &envir, const int hydroMode)
{
	typedef std::vector<cx_double> cx_stdvec;

	const vec &t = envir.getTimeArray();
	int nWaves = envir.numberOfWaveComponents();

	if (m_node1Pos_sd.at(2) > 0)
		return;

	// Calculate immersed length and its discretization
	double Lw;
	int npts;
	calculateImmersedLengthProperties_sd(Lw, npts, m_dL);

	// Position of each node along the cylinder length
	m_nodesArray = zeros(3, npts);
	for (int iNode = 0; iNode < npts; ++iNode)
	{
		m_nodesArray.col(iNode) = m_node1Pos_sd + iNode * m_dL * m_zvec_sd;
	}

	// Wave forces or kinematics at each node at each time step
	m_hydroForce_1st_Array = zeros(t.size(), 6);
	m_hydroForce_2nd_Array = zeros(t.size(), 6);
	m_waveElevAtWL = zeros(t.size(), 1);
	m_du1dt_Array_x = zeros(t.size(), npts);
	m_du1dt_Array_y = zeros(t.size(), npts);
	m_du1dt_Array_z = zeros(t.size(), npts);
	m_u1_Array_x = zeros(t.size(), npts);
	m_u1_Array_y = zeros(t.size(), npts);
	m_u1_Array_z = zeros(t.size(), npts);

	// Need to store the frequency of each wave in case the IFFT won't be used
	vec w{ zeros(nWaves, 1) };

	// Complex amplitudes	
	cx_mat amp_hydroForce_1st(nWaves, 6, fill::zeros);	
	cx_vec amp_waveElevAtWL(nWaves, 1, fill::zeros);
	cx_mat amp_du1dt_x(nWaves, npts, fill::zeros);
	cx_mat amp_du1dt_y(nWaves, npts, fill::zeros);
	cx_mat amp_du1dt_z(nWaves, npts, fill::zeros);
	cx_mat amp_u1_x(nWaves, npts, fill::zeros);
	cx_mat amp_u1_y(nWaves, npts, fill::zeros);
	cx_mat amp_u1_z(nWaves, npts, fill::zeros);
	for (unsigned int iWave = 0; iWave < nWaves; ++iWave)
	{
		const Wave &wave(envir.getWave(iWave));
		w.at(iWave) = wave.angFreq();

		if (wave.amp() == 0) continue;

		amp_hydroForce_1st.row(iWave) = hydroForce_1st_coefs(wave, envir.watDensity(), envir.watDepth(), envir.gravity()).st();
		amp_waveElevAtWL.at(iWave) = envir.waveElev_coef(m_nodesArray.at(0, npts - 1), m_nodesArray.at(1, npts - 1), iWave); // Wave elevation is evaluated at the last node, which is the intersection with the water line

		for (int iNode = 0; iNode < npts; ++iNode)
		{
			cx_vec::fixed<3> aux = envir.u1_coef(m_nodesArray.at(0, iNode), m_nodesArray.at(1, iNode), m_nodesArray.at(2, iNode), iWave);
			amp_u1_x.at(iWave, iNode) = aux.at(0);
			amp_u1_y.at(iWave, iNode) = aux.at(1);
			amp_u1_z.at(iWave, iNode) = aux.at(2);

			aux = envir.du1dt_coef(m_nodesArray.at(0, iNode), m_nodesArray.at(1, iNode), m_nodesArray.at(2, iNode), iWave);
			amp_du1dt_x.at(iWave, iNode) = aux.at(0);
			amp_du1dt_y.at(iWave, iNode) = aux.at(1);
			amp_du1dt_z.at(iWave, iNode) = aux.at(2);
		}		
	}

	m_waveElevAtWL = envir.timeSeriesFromAmp(amp_waveElevAtWL, w);
	m_hydroForce_1st_Array = envir.timeSeriesFromAmp(amp_hydroForce_1st, w);
	m_du1dt_Array_x = envir.timeSeriesFromAmp(amp_du1dt_x, w);
	m_du1dt_Array_y = envir.timeSeriesFromAmp(amp_du1dt_y, w);
	m_du1dt_Array_z = envir.timeSeriesFromAmp(amp_du1dt_z, w);
	m_u1_Array_x = envir.timeSeriesFromAmp(amp_u1_x, w);
	m_u1_Array_y = envir.timeSeriesFromAmp(amp_u1_y, w);
	m_u1_Array_z = envir.timeSeriesFromAmp(amp_u1_z, w);

	// Quantities that are used only in second order analysis
	if (hydroMode == 2)
	{
		m_du1dx_Array_x = zeros(t.size(), npts);
		m_du1dy_Array_y = zeros(t.size(), npts);
		m_du1dz_Array_z = zeros(t.size(), npts);
		m_du1dx_Array_y = zeros(t.size(), npts);
		m_du1dx_Array_z = zeros(t.size(), npts);		
		m_du1dy_Array_z = zeros(t.size(), npts);		

		m_da1dx_Array_x = zeros(t.size(), npts);
		m_da1dy_Array_y = zeros(t.size(), npts);
		m_da1dz_Array_z = zeros(t.size(), npts);
		m_da1dx_Array_y = zeros(t.size(), npts);
		m_da1dx_Array_z = zeros(t.size(), npts);
		m_da1dy_Array_z = zeros(t.size(), npts);

		m_gradP1_Array_x = zeros(t.size(), 2);
		m_gradP1_Array_y = zeros(t.size(), 2);
		m_gradP1_Array_z = zeros(t.size(), 2);

		cx_mat amp_du1dx_x(nWaves, npts, fill::zeros);
		cx_mat amp_du1dy_y(nWaves, npts, fill::zeros);
		cx_mat amp_du1dz_z(nWaves, npts, fill::zeros);
		cx_mat amp_du1dx_y(nWaves, npts, fill::zeros);
		cx_mat amp_du1dx_z(nWaves, npts, fill::zeros);		
		cx_mat amp_du1dy_z(nWaves, npts, fill::zeros);

		cx_mat amp_da1dx_x(nWaves, npts, fill::zeros);
		cx_mat amp_da1dy_y(nWaves, npts, fill::zeros);
		cx_mat amp_da1dz_z(nWaves, npts, fill::zeros);
		cx_mat amp_da1dx_y(nWaves, npts, fill::zeros);
		cx_mat amp_da1dx_z(nWaves, npts, fill::zeros);
		cx_mat amp_da1dy_z(nWaves, npts, fill::zeros);

		cx_mat amp_gradP1_x(nWaves, 2, fill::zeros);
		cx_mat amp_gradP1_y(nWaves, 2, fill::zeros);
		cx_mat amp_gradP1_z(nWaves, 2, fill::zeros);
		
		for (unsigned int iWave = 0; iWave < nWaves; ++iWave)
		{
			const Wave &wave(envir.getWave(iWave));

			for (int iNode = 0; iNode < npts; ++iNode)
			{
				// Velocity gradient
				cx_vec::fixed<3> aux = envir.du1dx_coef(m_nodesArray.at(0, iNode), m_nodesArray.at(1, iNode), m_nodesArray.at(2, iNode), iWave);
				amp_du1dx_x.at(iWave, iNode) = aux.at(0);
				amp_du1dx_y.at(iWave, iNode) = aux.at(1);
				amp_du1dx_z.at(iWave, iNode) = aux.at(2);

				aux = envir.du1dy_coef(m_nodesArray.at(0, iNode), m_nodesArray.at(1, iNode), m_nodesArray.at(2, iNode), iWave);
				amp_du1dy_y.at(iWave, iNode) = aux.at(1);
				amp_du1dy_z.at(iWave, iNode) = aux.at(2);

				aux = envir.du1dz_coef(m_nodesArray.at(0, iNode), m_nodesArray.at(1, iNode), m_nodesArray.at(2, iNode), iWave);
				amp_du1dz_z.at(iWave, iNode) = aux.at(2);

				// Acceleration gradient
				aux = envir.da1dx_coef(m_nodesArray.at(0, iNode), m_nodesArray.at(1, iNode), m_nodesArray.at(2, iNode), iWave);
				amp_da1dx_x.at(iWave, iNode) = aux.at(0);
				amp_da1dx_y.at(iWave, iNode) = aux.at(1);
				amp_da1dx_z.at(iWave, iNode) = aux.at(2);

				aux = envir.da1dy_coef(m_nodesArray.at(0, iNode), m_nodesArray.at(1, iNode), m_nodesArray.at(2, iNode), iWave);
				amp_da1dy_y.at(iWave, iNode) = aux.at(1);
				amp_da1dy_z.at(iWave, iNode) = aux.at(2);

				aux = envir.da1dz_coef(m_nodesArray.at(0, iNode), m_nodesArray.at(1, iNode), m_nodesArray.at(2, iNode), iWave);
				amp_da1dz_z.at(iWave, iNode) = aux.at(2);

				if (iNode == 0 || iNode == npts-1)
				{
					int auxInd = 0;
					if (iNode == npts - 1) auxInd = 1;

					// Pressure gradient
					aux = envir.gradP1_coef(m_nodesArray.at(0, iNode), m_nodesArray.at(1, iNode), m_nodesArray.at(2, iNode), iWave);
					amp_gradP1_x.at(iWave, auxInd) = aux.at(0);
					amp_gradP1_y.at(iWave, auxInd) = aux.at(1);
					amp_gradP1_z.at(iWave, auxInd) = aux.at(2);
				}
			}
		}

		m_du1dx_Array_x = envir.timeSeriesFromAmp(amp_du1dx_x, w);
		m_du1dy_Array_y = envir.timeSeriesFromAmp(amp_du1dy_y, w);
		m_du1dx_Array_y = envir.timeSeriesFromAmp(amp_du1dx_y, w);
		m_du1dx_Array_z = envir.timeSeriesFromAmp(amp_du1dx_z, w);
		m_du1dz_Array_z = envir.timeSeriesFromAmp(amp_du1dz_z, w);		
		m_du1dy_Array_z = envir.timeSeriesFromAmp(amp_du1dy_z, w);

		m_da1dx_Array_x = envir.timeSeriesFromAmp(amp_da1dx_x, w);
		m_da1dy_Array_y = envir.timeSeriesFromAmp(amp_da1dy_y, w);
		m_da1dx_Array_y = envir.timeSeriesFromAmp(amp_da1dx_y, w);
		m_da1dx_Array_z = envir.timeSeriesFromAmp(amp_da1dx_z, w);
		m_da1dz_Array_z = envir.timeSeriesFromAmp(amp_da1dz_z, w);
		m_da1dy_Array_z = envir.timeSeriesFromAmp(amp_da1dy_z, w);

		m_gradP1_Array_x = envir.timeSeriesFromAmp(amp_gradP1_x, w);
		m_gradP1_Array_y = envir.timeSeriesFromAmp(amp_gradP1_y, w);
		m_gradP1_Array_z = envir.timeSeriesFromAmp(amp_gradP1_z, w);
	}

	// Special treatment for the force due to the second-order wave potential due to the double sum
	// See comments in envir.evaluateWaveKinematics() for details
	if (hydroMode == 2)
	{
		if (envir.getFlagIFFT())
		{
			cx_mat amp_hydroForce_2nd(nWaves, 6, fill::zeros);
			for (int iWave = 0; iWave < nWaves; ++iWave)
			{
				const Wave &wave_ii(envir.getWave(iWave));
				if (wave_ii.amp() == 0) continue;

				for (int jWave = 0; jWave <= iWave; ++jWave)
				{
					const Wave &wave_jj(envir.getWave(jWave));
					if (wave_jj.amp() == 0) continue;

					cx_rowvec::fixed<6> aux = hydroForce_2ndPot_coefs(wave_ii, wave_jj, envir.watDensity(), envir.watDepth(), envir.gravity()).st();
					if (iWave != jWave)
					{
						aux *= 2; // Because we are summing only the upper part of the matrix
					}
					amp_hydroForce_2nd.row(iWave - jWave) += aux;
				}
			}

			m_hydroForce_2nd_Array = nWaves * mkl_ifft_real(amp_hydroForce_2nd) % repmat(envir.getRampArray(), 1, 6);
		}
		else
		{
			for (unsigned int it = 0; it < t.size(); ++it)
			{
				for (int iWave = 0; iWave < nWaves; ++iWave)
				{
					const Wave &wave_ii(envir.getWave(iWave));
					if (wave_ii.amp() == 0) continue;

					for (int jWave = 0; jWave <= iWave; ++jWave)
					{
						const Wave &wave_jj(envir.getWave(jWave));
						if (wave_jj.amp() == 0) continue;

						double w_ii{ wave_ii.angFreq() }, w_jj{ wave_jj.angFreq() };
						cx_double sinCos{ cos((w_ii - w_jj) * t.at(it)), sin((w_ii - w_jj) * t.at(it)) };
						cx_rowvec::fixed<6> amp_dw = hydroForce_2ndPot_coefs(wave_ii, wave_jj, envir.watDensity(), envir.watDepth(), envir.gravity()).st();
						if (iWave != jWave)
						{
							amp_dw *= 2;
						}

						m_hydroForce_2nd_Array.row(it) += real(amp_dw * sinCos);
					}
				}
			}
			m_hydroForce_2nd_Array %= repmat(envir.getRampArray(), 1, 6);
		}
	}
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

vec::fixed<6> MorisonCirc::hydrostaticForce_helper(const double rho, const double g, const vec::fixed<3> &refPt, const vec::fixed<3> &n1, const vec::fixed<3> &n2_in, const vec::fixed<3> &xvec, const vec::fixed<3> &yvec, const vec::fixed<3> &zvec) const
{
	// Forces and moments acting at the Morison Element
	vec::fixed<6> force(fill::zeros);
	vec::fixed<3> n2(n2_in);

	// Use a more friendly notation
	double D = m_diam;

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

	// The moments acting on the cylinders were calculated with respect to the first node
	// We need to change the fulcrum to refPt
	force.rows(3, 5) = force.rows(3, 5) + cross(n1 - refPt, force.rows(0, 2));

	return force;
}

vec::fixed<6> MorisonCirc::hydrostaticForce(const double rho, const double g, const vec::fixed<3> &refPt) const
{
	// Forces and moments acting at the Morison Element
	vec::fixed<6> force(fill::zeros);

	// Use a more friendly notation
	double D = m_diam;
	
	return hydrostaticForce_helper(rho, g, refPt, node1Pos(), node2Pos(), m_xvec, m_yvec, m_zvec);
}

vec::fixed<6> MorisonCirc::hydrostaticForce_sd(const double rho, const double g, const vec::fixed<3> &refPt) const
{
	// Forces and moments acting at the Morison Element
	vec::fixed<6> force(fill::zeros);

	// Use a more friendly notation
	double D = m_diam;
	
	return hydrostaticForce_helper(rho, g, refPt, node1Pos_sd(), node2Pos_sd(), m_xvec_sd, m_yvec_sd, m_zvec_sd);
}

/*
Functions to evaluate force components 
*/

// TODO: Implement Wheeler stretching in force components
vec::fixed<6> MorisonCirc::hydroForce_1st(const ENVIR &envir, const vec::fixed<3> &refPt) const
{
	vec::fixed<6> force(fill::zeros);

	if (m_hydroForce_1st_Array.is_empty())
	{
		for (unsigned int ii = 0; ii < envir.numberOfWaveComponents(); ++ii)
		{
			const Wave &wave(envir.getWave(ii));
			double w{ wave.angFreq() };
			cx_double sinCos({ cos(w * envir.time()), sin(w * envir.time()) });

			force += real(hydroForce_1st_coefs(wave, envir.watDensity(), envir.watDepth(), envir.gravity()) * sinCos) * envir.ramp();
		}
	}
	else
	{
		uword ind1 = envir.getInd4interp1();
		force = m_hydroForce_1st_Array.row(ind1).t();

		if (envir.shouldInterp())
		{
			uword ind2 = envir.getInd4interp2();
			const vec &t = envir.getTimeArray();
			force += (m_hydroForce_1st_Array.row(ind2).t() - m_hydroForce_1st_Array.row(ind1).t()) * (envir.time() - t(ind1)) / (t(ind2) - t(ind1));
		}
	}

	force.rows(3, 5) += cross(m_node1Pos_sd - refPt, force.rows(0, 2));
	return force;
}

vec::fixed<6> MorisonCirc::hydroForce_drag(const ENVIR &envir, const vec::fixed<3> &refPt, bool flagUse1stOrd) const
{
	vec::fixed<6> force(fill::zeros);
	vec::fixed<3> zvec{ m_zvec_sd };
	vec::fixed<3> vel1 = m_node1Vel;
	vec::fixed<3> vel2 = m_node2Vel;
	if (flagUse1stOrd)
	{
		vel1 = m_node1Vel_1stOrd;
		vel2 = m_node2Vel_1stOrd;
	}
	vec::fixed<3> v_axial = arma::dot(vel1, zvec) * zvec; // Since the cylinder is a rigid body, this is the same for all the nodes

	bool evaluateTopNode{ true };
	if (m_node2Pos_sd.at(2) > 0)
	{
		evaluateTopNode = false;
	}

	int ncyl = m_numNodesBelowWL; // If the body is fixed, this is the same as m_nodesArray.n_cols. Otherwise, it is updated at each time step
	for (int ii = 0; ii < ncyl; ++ii)
	{
		vec::fixed<3> n_ii = nodePos_sd(ii);

		// Velocity of the integration point
		double lambda = norm(n_ii - m_node1Pos_sd, 2) / norm(m_node2Pos_sd - m_node1Pos_sd);		
		vec::fixed<3> vel_ii = vel1 + lambda * (vel2 - vel1);

		vec::fixed<3> u1(this->u1(envir, ii));
		vec::fixed<3> u1_axial = arma::dot(u1, zvec) * zvec;
		if (ii == 0)
		{
			force.rows(0, 2) += 0.5 * envir.watDensity() * m_axialCD_1 * datum::pi * (m_diam*m_diam / 4.) * norm(u1_axial - v_axial, 2) * (u1_axial - v_axial);
		}
		else if (ii == ncyl - 1 && evaluateTopNode)
		{
			force.rows(0, 2) += 0.5 * envir.watDensity() * m_axialCD_2 * datum::pi * (m_diam*m_diam / 4.) * norm(u1_axial - v_axial, 2) * (u1_axial - v_axial);
		}

		// Quadratic drag force along the length
		u1 -= u1_axial;
		vec::fixed<3> force_ii = 0.5 * envir.watDensity() * m_CD * m_diam * norm(u1 - (vel_ii - v_axial), 2) * (u1 - (vel_ii - v_axial));

		// Integrate the forces along the cylinder using Simpson's Rule
		if (ii == 0 || ii == ncyl - 1)
		{
			force += (m_dL / 3.0) * join_cols(force_ii, cross(n_ii - m_node1Pos_sd, force_ii));
		}
		else if (ii % 2 != 0)
		{
			force += (4 * m_dL / 3.0) * join_cols(force_ii, cross(n_ii - m_node1Pos_sd, force_ii));
		}
		else
		{
			force += (2 * m_dL / 3.0) * join_cols(force_ii, cross(n_ii - m_node1Pos_sd, force_ii));
		}
	}

	force.rows(3, 5) += cross(m_node1Pos_sd - refPt, force.rows(0, 2));
	return force;
}

vec::fixed<6> MorisonCirc::hydroForce_relWaveElev(const ENVIR &envir, const vec::fixed<3> &refPt) const
{
	vec::fixed<6> force(fill::zeros);

	// The cylinder must cross the mean waterline for this component to make sense
	if (m_node1Pos_sd.at(2)*m_node2Pos_sd.at(2) > 0)
		return force;
		
	vec::fixed<3> du1dt(this->du1dt(envir, m_numNodesBelowWL-1));	
	du1dt -= dot(du1dt, m_zvec_sd) * m_zvec_sd;
	vec::fixed<3> acc = m_accNodeAtWL_1stOrd;
	acc -= dot(acc, m_zvec_sd) * m_zvec_sd;
	double eta(this->waveElevAtWL(envir));
	vec::fixed<3> n_wl = this->nodePos_sd(m_numNodesBelowWL - 1);
		
	double g{ envir.gravity() };
	double zBody = m_Zwl;
	double L = norm(m_node2Pos_sd - m_node1Pos_sd);
	double ry(dot(m_node2Pos_1stOrd - m_node1Pos_1stOrd, m_xvec_sd) / L), rx(dot(m_node2Pos_1stOrd - m_node1Pos_1stOrd, m_yvec_sd) / L);

	force.rows(0, 2) = (datum::pi * m_diam*m_diam / 4.) * envir.watDensity() * (eta - zBody) * (m_CM * du1dt + g * rx*m_yvec_sd - g * ry*m_xvec_sd);
	force.rows(3, 5) = cross(n_wl - refPt, force.rows(0, 2));

	return force;	
}

vec::fixed<6> MorisonCirc::hydroForce_2ndPot(const ENVIR &envir, const vec::fixed<3> &refPt) const
{
	vec::fixed<6> force(fill::zeros);

	if (m_hydroForce_2nd_Array.is_empty())
	{
		// When i == j, du2dt = {0,0,0}, so it is safe to skip this part of the loop.
		// Besides, as only the real part of the second-order difference-frequency potential is used,
		// the acceleration due to a pair ij is equal to ji.
		for (unsigned int ii = 0; ii < envir.numberOfWaveComponents(); ++ii)
		{
			for (unsigned int jj = ii + 1; jj < envir.numberOfWaveComponents(); ++jj)
			{
				const Wave &wave_ii(envir.getWave(ii));
				const Wave &wave_jj(envir.getWave(jj));

				double w_ii{ wave_ii.angFreq() }, w_jj{ wave_jj.angFreq() };
				cx_double sinCos({ cos((w_ii - w_jj) * envir.time()), sin((w_ii - w_jj) * envir.time()) });

				force += 2 * real(hydroForce_2ndPot_coefs(wave_ii, wave_jj, envir.watDensity(), envir.watDepth(), envir.gravity()) * sinCos) * envir.ramp();
			}
		}
	}
	else
	{
		uword ind1 = envir.getInd4interp1();
		force = m_hydroForce_2nd_Array.row(ind1).t();

		if (envir.shouldInterp())
		{
			uword ind2 = envir.getInd4interp2();
			const vec &t = envir.getTimeArray();
			force += (m_hydroForce_2nd_Array.row(ind2).t() - m_hydroForce_2nd_Array.row(ind1).t()) * (envir.time() - t(ind1)) / (t(ind2) - t(ind1));
		}
	}

	force.rows(3, 5) += cross(m_node1Pos_sd - refPt, force.rows(0, 2));
	return force;
}

vec::fixed<6> MorisonCirc::hydroForce_convecAcc(const ENVIR &envir, const vec::fixed<3> &refPt) const
{
	vec::fixed<6> force(fill::zeros);
	vec::fixed<3> zvec{ m_zvec_sd };

	bool evaluateTopNode{ true };
	if (m_node2Pos_sd.at(2) > 0)
	{
		evaluateTopNode = false;
	}

	int ncyl = m_numNodesBelowWL; // If the body is fixed, this is the same as m_nodesArray.n_cols. Otherwise, it is updated at each time step
	vec::fixed<3> a_c(fill::zeros);
	for (int ii = 0; ii < ncyl; ++ii)
	{		
		vec::fixed<3> n_ii = nodePos_sd(ii);

		vec::fixed<3> u1(this->u1(envir, ii));
		vec::fixed<3> du1dx(this->du1dx(envir, ii));
		vec::fixed<3> du1dy(this->du1dy(envir, ii));
		vec::fixed<3> du1dz(this->du1dz(envir, ii));

		a_c.at(0) = u1.at(0) * du1dx.at(0) + u1.at(1) * du1dy.at(0) + u1.at(2) * du1dz.at(0);
		a_c.at(1) = u1.at(0) * du1dx.at(1) + u1.at(1) * du1dy.at(1) + u1.at(2) * du1dz.at(1);
		a_c.at(2) = u1.at(0) * du1dx.at(2) + u1.at(1) * du1dy.at(2) + u1.at(2) * du1dz.at(2);
		vec::fixed<3> a_c_axial = arma::dot(a_c, zvec) * zvec;

		if (ii == 0)
		{
			force.rows(0, 2) += (4 / 3.) * datum::pi * (m_diam*m_diam*m_diam / 8.)  * envir.watDensity() * m_axialCa_1 * a_c_axial;
		}
		else if (ii == ncyl - 1 && evaluateTopNode)
		{
			force.rows(0, 2) += (4 / 3.) * datum::pi * (m_diam*m_diam*m_diam / 8.)  * envir.watDensity() * m_axialCa_2 * a_c_axial;
		}

		// Componente that is perpendicular to the cylinder axis
		a_c -= a_c_axial;
		vec::fixed<3> force_ii = datum::pi * (m_diam*m_diam / 4.) * envir.watDensity() * m_CM  * a_c;

		// Integrate the forces along the cylinder using Simpson's Rule
		if (ii == 0 || ii == ncyl - 1)
		{
			force += (m_dL / 3.0) * join_cols(force_ii, cross(n_ii - m_node1Pos_sd, force_ii));
		}
		else if (ii % 2 != 0)
		{
			force += (4 * m_dL / 3.0) * join_cols(force_ii, cross(n_ii - m_node1Pos_sd, force_ii));
		}
		else
		{
			force += (2 * m_dL / 3.0) * join_cols(force_ii, cross(n_ii - m_node1Pos_sd, force_ii));
		}
	}

	force.rows(3, 5) += cross(m_node1Pos_sd - refPt, force.rows(0, 2));
	return force;
}

vec::fixed<6> MorisonCirc::hydroForce_axDiverg(const ENVIR &envir, const vec::fixed<3> &refPt) const
{
	vec::fixed<6> force(fill::zeros);
	vec::fixed<3> xvec{ m_xvec_sd }, yvec{ m_yvec_sd }, zvec{ m_zvec_sd };

	int ncyl = m_numNodesBelowWL;
	for (int ii = 0; ii < ncyl; ++ii)
	{
		vec::fixed<3> n_ii = nodePos_sd(ii);

		// Velocity of the integration point
		double lambda = norm(n_ii - m_node1Pos_sd, 2) / norm(m_node2Pos_sd - m_node1Pos_sd);
		vec::fixed<3> vel_ii = m_node1Vel_1stOrd + lambda * (m_node2Vel_1stOrd - m_node1Vel_1stOrd);
		
		vec::fixed<3> u(this->u1(envir, ii));
		vec::fixed<3> du1dx(this->du1dx(envir, ii));
		vec::fixed<3> du1dy(this->du1dy(envir, ii));
		vec::fixed<3> du1dz(this->du1dz(envir, ii));

		double dwdz = dot(du1dx, zvec) * zvec.at(0) + dot(du1dy, zvec) * zvec.at(1) + dot(du1dz, zvec) * zvec.at(2);
		vec::fixed<3> a_axdiv = dwdz * (dot(u - vel_ii, xvec)*xvec + dot(u - vel_ii, yvec)*yvec);
		vec::fixed<3> force_ii = datum::pi * (m_diam*m_diam / 4.) * envir.watDensity() * (m_CM - 1) * a_axdiv;

		// Integrate the forces along the cylinder using Simpson's Rule
		if (ii == 0 || ii == ncyl - 1)
		{
			force += (m_dL / 3.0) * join_cols(force_ii, cross(n_ii - refPt, force_ii));
		}
		else if (ii % 2 != 0)
		{
			force += (4 * m_dL / 3.0) * join_cols(force_ii, cross(n_ii - refPt, force_ii));
		}
		else
		{
			force += (2 * m_dL / 3.0) * join_cols(force_ii, cross(n_ii - refPt, force_ii));
		}
	}

	return force;
}

vec::fixed<6> MorisonCirc::hydroForce_accGradient(const ENVIR & envir, const vec::fixed<3>& refPt) const // TODO: include axial contribution
{
	vec::fixed<6> force(fill::zeros);

	bool evaluateTopNode{ true };
	if (m_node2Pos_sd.at(2) > 0)
	{
		evaluateTopNode = false;
	}

	int ncyl = m_numNodesBelowWL; // If the body is fixed, this is the same as m_nodesArray.n_cols. Otherwise, it is updated at each time step
	vec::fixed<3> a_g(fill::zeros);
	for (int ii = 0; ii < ncyl; ++ii)
	{
		vec::fixed<3> n_ii = nodePos_sd(ii);

		// Diference between the 1st order position and the fixed position
		double xii = m_node1Pos_1stOrd.at(0) + m_dL * ii * m_zvec_1stOrd.at(0) - n_ii.at(0);
		double yii = m_node1Pos_1stOrd.at(1) + m_dL * ii * m_zvec_1stOrd.at(1) - n_ii.at(1);
		double zii = m_node1Pos_1stOrd.at(2) + m_dL * ii * m_zvec_1stOrd.at(2) - n_ii.at(2);

		vec::fixed<3> a_g(this->da1dx(envir, ii)*xii + this->da1dy(envir, ii)*yii + this->da1dz(envir, ii)*zii);		
		vec::fixed<3> a_g_axial = arma::dot(a_g, m_zvec_sd) * m_zvec_sd;
		if (ii == 0)
		{			
			force.rows(0, 2) += datum::pi * (m_diam*m_diam*m_diam / 6.) * envir.watDensity() * m_axialCa_1 * a_g_axial;

			if (m_botPressFlag)
			{
				double gP1 = arma::dot(gradP1(envir, ii), arma::vec{ xii, yii, zii });
				force.rows(0, 2) += datum::pi * (m_diam*m_diam / 4.) * gP1 * m_zvec_sd;
			}
		}
		else if (ii == ncyl - 1 && evaluateTopNode)
		{
			double gP1 = arma::dot(gradP1(envir, ii), arma::vec{ xii, yii, zii });
			force.rows(0, 2) += datum::pi * (m_diam*m_diam*m_diam / 6.) * envir.watDensity() * m_axialCa_2 * a_g_axial;

			if (m_botPressFlag)
			{
				double gP1 = arma::dot(gradP1(envir, ii), arma::vec{ xii, yii, zii });
				force.rows(0, 2) -= datum::pi * (m_diam*m_diam / 4.) * gP1 * m_zvec_sd;
			}
		}

		a_g -= a_g_axial;
		vec::fixed<3> force_ii = datum::pi * (m_diam*m_diam / 4.) * envir.watDensity() * m_CM  * a_g;

		// Integrate the forces along the cylinder using Simpson's Rule
		if (ii == 0 || ii == ncyl - 1)
		{
			force += (m_dL / 3.0) * join_cols(force_ii, cross(n_ii - refPt, force_ii));
		}
		else if (ii % 2 != 0)
		{
			force += (4 * m_dL / 3.0) * join_cols(force_ii, cross(n_ii - refPt, force_ii));
		}
		else
		{
			force += (2 * m_dL / 3.0) * join_cols(force_ii, cross(n_ii - refPt, force_ii));
		}
	}

	return force;
}

vec::fixed<6> MorisonCirc::hydroForce_slendBodyRot(const ENVIR & envir, const vec::fixed<3> &refPt) const
{
	vec::fixed<6> force(fill::zeros);

	vec::fixed<3> xvec{ m_xvec_sd }, yvec{ m_yvec_sd }, zvec{ m_zvec_sd };

	int ncyl = m_numNodesBelowWL;
	double L = norm(m_node2Pos_sd - m_node1Pos_sd);
	double wy(dot(m_node2Vel_1stOrd - m_node1Vel_1stOrd, xvec) / L), wx(dot(m_node2Vel_1stOrd - m_node1Vel_1stOrd, yvec) / L);
	for (int ii = 0; ii < ncyl; ++ii)
	{
		vec::fixed<3> n_ii = nodePos_sd(ii);

		// Velocity of the integration point
		double lambda = norm(n_ii - m_node1Pos_sd, 2) / norm(m_node2Pos_sd - m_node1Pos_sd);
		vec::fixed<3> vel_ii = m_node1Vel + lambda * (m_node2Vel - m_node1Vel);

		vec::fixed<3> u1(this->u1(envir,ii));
		vec::fixed<3> a_r = 2 * dot(u1 - vel_ii, zvec) *  (wy * xvec + wx * yvec);
		vec::fixed<3> force_ii = -datum::pi * (m_diam*m_diam / 4.) * envir.watDensity() * (m_CM - 1) * a_r;

		// Integrate the forces along the cylinder using Simpson's Rule
		if (ii == 0 || ii == ncyl - 1)
		{
			force += (m_dL / 3.0) * join_cols(force_ii, cross(n_ii - refPt, force_ii));
		}
		else if (ii % 2 != 0)
		{
			force += (4 * m_dL / 3.0) * join_cols(force_ii, cross(n_ii - refPt, force_ii));
		}
		else
		{
			force += (2 * m_dL / 3.0) * join_cols(force_ii, cross(n_ii - refPt, force_ii));
		}
	}

	return force;
}

vec::fixed<6> MorisonCirc::hydroForce_rem(const ENVIR & envir, const vec::fixed<3>& refPt) const
{
	vec::fixed<6> force(fill::zeros);
	vec::fixed<3> xvec{ m_xvec_sd }, yvec{ m_yvec_sd }, zvec{ m_zvec_sd };

	bool evaluateTopNode{ true };
	if (m_node2Pos_sd.at(2) > 0)
	{
		evaluateTopNode = false;
	}

	int ncyl = m_numNodesBelowWL; // If the body is fixed, this is the same as m_nodesArray.n_cols. Otherwise, it is updated at each time step
	for (int ii = 0; ii < ncyl; ++ii)
	{
		vec::fixed<3> n_ii = nodePos_sd(ii);

		// Centripetal acceleration of the integration point
		double lambda = norm(n_ii - m_node1Pos_sd, 2) / norm(m_node2Pos_sd - m_node1Pos_sd);
		vec::fixed<3> acc_ii = m_node1AccCentrip + lambda * (m_node2AccCentrip - m_node1AccCentrip);
		vec::fixed<3> acc_ii_axial = dot(acc_ii, zvec) * zvec;

		if (ii == 0 && m_botPressFlag)
		{
			// For due to the convective acceleration - axial part
			force.rows(0, 2) += -(4 / 3.) * datum::pi * (m_diam*m_diam*m_diam / 8.)  * envir.watDensity() * m_axialCa_1 * acc_ii_axial;

			// Quadratic pressure drop from Rainey's formulation
			vec::fixed<3> u1(this->u1(envir, ii));
			force.rows(0, 2) += -0.5 * datum::pi * (m_diam*m_diam / 4.) * envir.watDensity() * (m_CM - 1) * (pow(dot(u1 - m_node1Vel_1stOrd, xvec), 2) + pow(dot(u1 - m_node1Vel_1stOrd, yvec), 2)) * zvec;

			// Point load from Rainey's formulation that results in a Munk moment
			force.rows(0,2) += datum::pi * (m_diam*m_diam / 4.)* (m_CM - 1) * dot(u1 - m_node1Vel_1stOrd, zvec) * (cdot(u1 - m_node1Vel_1stOrd, xvec)*xvec + cdot(u1 - m_node1Vel_1stOrd, yvec)*yvec);
		}
		else if (ii == ncyl - 1 && m_botPressFlag && evaluateTopNode)
		{
			force.rows(0, 2) += -(4 / 3.) * datum::pi * (m_diam*m_diam*m_diam / 8.)  * envir.watDensity() * m_axialCa_2 * acc_ii_axial;

			// At the top node, the direction of the pressure drop is the opposite of the bottom
			vec::fixed<3> u1(this->u1(envir, ii));
			force.rows(0, 2) += 0.5 * datum::pi * (m_diam*m_diam / 4.) * envir.watDensity() * (m_CM - 1) * (pow(dot(u1 - m_node2Vel_1stOrd, xvec), 2) + pow(dot(u1 - m_node2Vel_1stOrd, yvec), 2)) * zvec;

			// As this component does not act along the axis of the cylinder, need to include the moment with respect to node1
			// It is also the opposite of the bottom
			vec::fixed<3> auxForce = -datum::pi * (m_diam*m_diam / 4.)* (m_CM - 1) * dot(u1 - m_node2Vel_1stOrd, zvec) * (cdot(u1 - m_node2Vel_1stOrd, xvec)*xvec + cdot(u1 - m_node2Vel_1stOrd, yvec)*yvec);
			force.rows(0, 2) += auxForce;
			force.rows(3, 5) += cross(m_node2Pos_sd - m_node1Pos_sd, auxForce);
		}

		// For due to the centripetal acceleration - part that is perpendicular to the cylinder axis
		acc_ii -= acc_ii_axial;
		vec::fixed<3> force_ii = -datum::pi * (m_diam*m_diam / 4.) * envir.watDensity() * (m_CM - 1)  * acc_ii;

		// Integrate the forces along the cylinder using Simpson's Rule
		if (ii == 0 || ii == ncyl - 1)
		{
			force += (m_dL / 3.0) * join_cols(force_ii, cross(n_ii - m_node1Pos_sd, force_ii));
		}
		else if (ii % 2 != 0)
		{
			force += (4 * m_dL / 3.0) * join_cols(force_ii, cross(n_ii - m_node1Pos_sd, force_ii));
		}
		else
		{
			force += (2 * m_dL / 3.0) * join_cols(force_ii, cross(n_ii - m_node1Pos_sd, force_ii));
		}
	}

	force.rows(3, 5) += cross(m_node1Pos_sd - refPt, force.rows(0, 2));
	return force;
}

void MorisonCirc::quantities4hydrostaticMatrix(double &zb, double &V, double &Awl, double &xwl, double &ywl, double &Ixx, double &Iyy, double &Ixy) const
{
	// Use a more friendly notation
	double R = 0.5*m_diam;	
		
	double cosAlpha = m_zvec.at(2); // Inclination of the cylinder (with respect to the vertical)
	double tanAlpha{ tan(acos(cosAlpha)) }; // Alpha can not be pi/2 because this would require the vertical position of the nodes to be the same, but this part of the code is only reached if the cilinder crosses the waterline
	double cosBeta = m_yvec.at(2); // Inclination in the plane of the mean water level. The rotations must be around the local y axis, and this is taken into account when creating the local basis
	double sin2Beta = 2 * cosBeta*sqrt(1 - cosBeta * cosBeta);
	double inclFactor = 1 / cosAlpha;	
	double Iyy_local = 0.25 * arma::datum::pi * R * R * R * R;
	double Ixx_local = Iyy_local * inclFactor;
	double L = norm(m_node2Pos - m_node1Pos);

	// If the cylinder is above the waterline or if it is completely submerged, then 
	// this terms are all zero
	if (m_node1Pos.at(2) >= 0 || m_node2Pos.at(2) <= 0)
	{
		Awl = 0;
		xwl = 0;
		ywl = 0;
		Ixx = 0;
		Ixy = 0;
	}
	else
	{
		vec::fixed<3> nwl = m_node1Pos + (m_node2Pos - m_node1Pos) * std::abs(0 - m_node1Pos.at(2)) / (m_node2Pos.at(2) - m_node1Pos.at(2)); // Intersection with the waterline
		L = norm(nwl - m_node1Pos);
		Awl = arma::datum::pi * R * R * inclFactor;
		xwl = nwl.at(0);
		ywl = nwl.at(1);
		Ixx = Ixx_local * cosBeta * cosBeta + Iyy_local * (1 - cosBeta * cosBeta);
		Iyy = Iyy_local * cosBeta * cosBeta + Ixx_local * (1 - cosBeta * cosBeta);
		Ixy = 0.5 * (Ixx_local + Iyy_local) * sin2Beta;
	}

	V = arma::datum::pi * R * R * L;		
	zb = (pow(tanAlpha*R, 2) + 4 * pow(L, 2)) / (8 * L) * m_zvec.at(2) + m_node1Pos.at(2);
}

/*
Below are the functions with the force expressions
*/
cx_vec::fixed<6> MorisonCirc::hydroForce_1st_coefs(const Wave &wave, double watDensity, double watDepth, double gravity) const
{
	double rho = watDensity;
	double h = watDepth;
	double g = gravity;
	double R = m_diam / 2.;
	double pi = arma::datum::pi;
	double k{ wave.waveNumber() };

	if (k <= 0)
	{
		return fill::zeros;
	}

	vec::fixed<3> Zvec{ 0, 0, 1 }; // Vectors of the global basis
	vec::fixed<3> Xvec{ 1, 0, 0 };

	vec::fixed<3> n1{ node1Pos_sd() }, n2{ node2Pos_sd() };
	vec::fixed<3> xvec{ m_xvec_sd }, yvec{ m_yvec_sd }, zvec{ m_zvec_sd };

	if (n1.at(2) >= 0)
	{
		return fill::zeros;
	}

	bool evaluateTopNode{ true };
	if (n2.at(2) > 0)
	{
		n2 = n1 + (abs(0 - n1(2)) / (n2(2) - n1(2))) * arma::norm(n2 - n1) * zvec;
		evaluateTopNode = false;
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

	double w{ wave.angFreq() }, A{ wave.amp() };
	double cosBeta{ wave.cosBeta() }, sinBeta{ wave.sinBeta() };
	double beta = wave.direction() * pi / 180.;
	double phase = wave.phase() * pi / 180.;
	double cosP{ cos(phase) }, sinP{ sin(phase) };

	// Avoid recalculting cossine and sine functions
	double cosBetaPsi = cos(beta - psi);
	double sinBetaPsi = sin(beta - psi);

	// Auxiliar variables to make the expressions shorter
	double mult_h_cos(0), mult_h_sin(0), mult_v_cos(0), mult_v_sin(0);
	double p1 = k * cosBeta*x1 + k * sinBeta*y1 + phase;
	double p2 = k * cosBeta*x2 + k * sinBeta*y2 + phase;
	double I = k * cosBetaPsi * tanAlpha;
	double Ih = (k * cosBeta*(x2 - x1) + k * sinBeta*(y2 - y1)) / L;
	double Q{ m_CM * w * w * A }, Qz1{ m_axialCa_1 * w * w * A };
	double Qz2 = (evaluateTopNode ? m_axialCa_2 * w * w * A : 0);


	// Ratios of hyperbolic functions need a special treatment otherwise they yield inf in large water depth
	double cosh_sinh_1{ 0 }, cosh_sinh_2{ 0 }, cosh_cosh_1{ 0 }, cosh_cosh_2{ 0 }, sinh_sinh_1{ 0 }, sinh_sinh_2{ 0 };
	if (k*h >= 10)
	{
		cosh_sinh_1 = exp(k*z1);
		cosh_cosh_1 = cosh_sinh_1;
		sinh_sinh_1 = cosh_sinh_1;

		cosh_sinh_2 = exp(k*z2);
		cosh_cosh_2 = cosh_sinh_2;
		sinh_sinh_2 = cosh_sinh_2;
	}
	else
	{
		cosh_sinh_1 = cosh(k * (z1 + h)) / sinh(k*h);
		cosh_cosh_1 = cosh(k * (z1 + h)) / cosh(k*h);
		sinh_sinh_1 = sinh(k * (z1 + h)) / sinh(k*h);

		cosh_sinh_2 = cosh(k * (z2 + h)) / sinh(k*h);
		cosh_cosh_2 = cosh(k * (z2 + h)) / cosh(k*h);
		sinh_sinh_2 = sinh(k * (z2 + h)) / sinh(k*h);
	}

	/*
		FORCES
	*/
	// The integration is different for a horizontal cylinder, as the integration variable
	// along the cylinder is not related to the vertical position. This special case is treated separately.
	if (arma::is_finite(tanAlpha))
	{
		// Contribution of the horizontal acceleration
		mult_h_cos = k * sin(p2) * sinh_sinh_2 - I * cos(p2) * cosh_sinh_2;
		mult_h_cos += -k * sin(p1) * sinh_sinh_1 + I * cos(p1) * cosh_sinh_1;
		mult_h_cos *= 1 / (I*I + k * k);

		mult_h_sin = k * cos(p2) * sinh_sinh_2 + I * sin(p2) * cosh_sinh_2;
		mult_h_sin += -k * cos(p1) * sinh_sinh_1 - I * sin(p1) * cosh_sinh_1;
		mult_h_sin *= 1 / (I*I + k * k);

		// This cosAlpha is due to the change of integration variable.
		// We know cosAlpha !=0 because horizontal cylinders are treated separately
		mult_h_cos *= 1 / cosAlpha;
		mult_h_sin *= 1 / cosAlpha;

		// Contribution of the vertical acceleration
		mult_v_cos = k * cos(p2) * cosh_sinh_2 + I * sin(p2) * sinh_sinh_2;
		mult_v_cos += -k * cos(p1) * cosh_sinh_1 - I * sin(p1) * sinh_sinh_1;
		mult_v_cos *= 1 / (I*I + k * k);

		mult_v_sin = k * sin(p2) * cosh_sinh_2 - I * cos(p2) * sinh_sinh_2;
		mult_v_sin += -k * sin(p1) * cosh_sinh_1 + I * cos(p1) * sinh_sinh_1;
		mult_v_sin *= 1 / (I*I + k * k);

		mult_v_cos *= 1 / cosAlpha;
		mult_v_sin *= -1 / cosAlpha;
	}
	else
	{
		// Contribution of the horizontal acceleration
		if (Ih != 0)
		{
			mult_h_cos = -cos(p2) + cos(p1);
			mult_h_cos *= 1 / Ih;

			mult_h_sin = sin(p2) - sin(p1);
			mult_h_sin *= 1 / Ih;
		}
		else
		{
			mult_h_cos = sin(p1)*L;
			mult_h_sin = cos(p1)*L;
		}

		// Contribution of the vertical acceleration
		mult_v_cos = mult_h_sin * sinh_sinh_1;
		mult_v_sin = -mult_h_cos * sinh_sinh_1;

		// Need the factor of the water depth in the horizontal acceleration as well
		mult_h_cos *= cosh_sinh_1;
		mult_h_sin *= cosh_sinh_1;
	}

	// Forces along the local x and y axes
	double fx_cos = Q * (cosBetaPsi * cosAlpha * mult_h_cos + sinAlpha * mult_v_cos);
	double fx_sin = Q * (cosBetaPsi * cosAlpha * mult_h_sin + sinAlpha * mult_v_sin);
	double fy_cos = Q * sinBetaPsi * mult_h_cos;
	double fy_sin = Q * sinBetaPsi * mult_h_sin;

	// Forces along the local z axis. Factor 4*R/3. due to this force being acceleration * volume, while in the end 
	// of this function things are multiplied by the area of the circle
	double fz_cos = (4 * R / 3.) * cosBetaPsi * sinAlpha * (Qz1 * cosh_sinh_1 * sin(p1) - Qz2 * cosh_sinh_2 * sin(p2)); // Contribution of horizontal acceleration
	double fz_sin = (4 * R / 3.) * cosBetaPsi * sinAlpha * (-Qz1 * cosh_sinh_1 * cos(p1) + Qz2 * cosh_sinh_2 * cos(p2));
	fz_cos += -(4 * R / 3.) * cosAlpha * (Qz1 * sinh_sinh_1 * cos(p1) - Qz2 * sinh_sinh_2 * cos(p2)); // Contribution of vertical acceleration
	fz_sin += (4 * R / 3.) * cosAlpha * (Qz1 * sinh_sinh_1 * sin(p1) - Qz2 * sinh_sinh_2 * sin(p2));

	/*
		MOMENTS
	*/
	if (arma::is_finite(tanAlpha))
	{
		// Contribution of the horizontal acceleration
		mult_h_cos = k * sinh_sinh_2 * ((k*k + I * I)*(z2 - z1)*sin(p2) + 2 * I*cos(p2))
			+ cosh_sinh_2 * (-I * (k*k + I * I)*(z2 - z1)*cos(p2) + (I*I - k * k)*sin(p2));
		mult_h_cos += -k * sinh_sinh_1 * 2 * I*cos(p1)
			- cosh_sinh_1 * (I*I - k * k)*sin(p1);
		mult_h_cos *= 1 / (I*I + k * k) / (I*I + k * k);

		mult_h_sin = k * sinh_sinh_2 * ((k*k + I * I)*(z2 - z1)*cos(p2) - 2 * I*sin(p2))
			+ cosh_sinh_2 * (I*(k*k + I * I)*(z2 - z1)*sin(p2) + (I*I - k * k)*cos(p2));
		mult_h_sin += k * sinh_sinh_1 * 2 * I*sin(p1)
			- cosh_sinh_1 * (I*I - k * k)*cos(p1);
		mult_h_sin *= 1 / (I*I + k * k) / (I*I + k * k);

		mult_h_cos *= 1 / cosAlpha / cosAlpha;
		mult_h_sin *= 1 / cosAlpha / cosAlpha;

		// Contribution of the vertical acceleration
		mult_v_cos = k * cosh_sinh_2 * ((k*k + I * I)*(z2 - z1)*cos(p2) - 2 * I*sin(p2))
			+ sinh_sinh_2 * (I*(k*k + I * I)*(z2 - z1)*sin(p2) + (I*I - k * k)*cos(p2));
		mult_v_cos += k * cosh_sinh_1 * 2 * I*sin(p1)
			- sinh_sinh_1 * (I*I - k * k)*cos(p1);
		mult_v_cos *= 1 / (I*I + k * k) / (I*I + k * k);

		mult_v_sin = k * cosh_sinh_2 * ((k*k + I * I)*(z2 - z1)*sin(p2) + 2 * I*cos(p2))
			+ sinh_sinh_2 * (-I * (k*k + I * I)*(z2 - z1)*cos(p2) + (I*I - k * k)*sin(p2));
		mult_v_sin += -k * cosh_sinh_1 * 2 * I*cos(p1)
			- sinh_sinh_1 * (I*I - k * k)*sin(p1);
		mult_v_sin *= 1 / (I*I + k * k) / (I*I + k * k);

		mult_v_cos *= 1 / cosAlpha / cosAlpha;
		mult_v_sin *= -1 / cosAlpha / cosAlpha;
	}
	else
	{
		// Contribution of the horizontal acceleration
		if (Ih != 0)
		{			
			mult_h_cos = -Ih * L*cos(p2) + sin(p2) - sin(p1);
			mult_h_cos *= 1 / (Ih*Ih);

			mult_h_sin = Ih * L*sin(p2) + cos(p2) - cos(p1);
			mult_h_sin *= 1 / (Ih*Ih);
		}
		else
		{
			mult_h_cos = sin(p1)*L*L*0.5;
			mult_h_sin = cos(p1)*L*L*0.5;
		}

		// Contribution of the vertical acceleration
		mult_v_cos = mult_h_sin * sinh_sinh_1;
		mult_v_sin = -mult_h_cos * sinh_sinh_1;

		// Need the factor of the water depth in the horizontal acceleration as well
		mult_h_cos *= cosh_sinh_1;
		mult_h_sin *= cosh_sinh_1;
	}

	// Forces along the local x and y axes
	double mx_cos = -Q * sinBetaPsi * mult_h_cos;
	double mx_sin = -Q * sinBetaPsi * mult_h_sin;
	double my_cos = Q * (cosBetaPsi * cosAlpha * mult_h_cos + sinAlpha * mult_v_cos);
	double my_sin = Q * (cosBetaPsi * cosAlpha * mult_h_sin + sinAlpha * mult_v_sin);

	if (m_botPressFlag)
	{
		fz_cos += g * A * cos(p1) * cosh_cosh_1;
		fz_sin += -g * A * sin(p1) * cosh_cosh_1;
		if (evaluateTopNode)
		{
			fz_cos -= g * A * cos(p2) * cosh_cosh_2;
			fz_sin -= -g * A * sin(p2) * cosh_cosh_2;
		}
	}

	cx_vec::fixed<6> coef(fill::zeros);
	for (int ii = 0; ii < 3; ++ii)
	{
		coef.at(ii) = cx_double{ fx_cos * xvec(ii) + fy_cos * yvec(ii) + fz_cos * zvec(ii) , fx_sin * xvec(ii) + fy_sin * yvec(ii) + fz_sin * zvec(ii) };
	}
	for (int ii = 3; ii < 6; ++ii)
	{
		coef.at(ii) = cx_double{ mx_cos * xvec(ii - 3) + my_cos * yvec(ii - 3) , mx_sin * xvec(ii - 3) + my_sin * yvec(ii - 3) };
	}

	return rho * pi* R * R * coef;
}

// TODO: this function has many similarities with its 1st order counterpart. It would be better to group 
// common code in an auxiliary function.
cx_vec::fixed<6> MorisonCirc::hydroForce_2ndPot_coefs(const Wave &wave_ii, const Wave &wave_jj , double watDensity, double watDepth, double gravity) const
{
	vec::fixed<6> force(fill::zeros);

	double rho = watDensity;
	double h = watDepth;
	double g = gravity;
	double R = m_diam / 2.;
	double pi = arma::datum::pi;
	vec::fixed<3> Zvec{ 0, 0, 1 }; // Vectors of the global basis
	vec::fixed<3> Xvec{ 1, 0, 0 };

	vec::fixed<3> n1{ node1Pos_sd() }, n2{ node2Pos_sd() };
	vec::fixed<3> xvec{ m_xvec_sd }, yvec{ m_yvec_sd }, zvec{ m_zvec_sd };

	if (n1.at(2) >= 0)
	{
		return fill::zeros;
	}

	bool evaluateTopNode{ true };
	if (n2.at(2) > 0)
	{
		n2 = n1 + (abs(0 - n1(2)) / (n2(2) - n1(2))) * arma::norm(n2 - n1) * zvec;
		evaluateTopNode = false;
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
	
	double w_ii{ wave_ii.angFreq() }, k_ii{ wave_ii.waveNumber() }, A_ii{ wave_ii.amp() };
	double w_jj{ wave_jj.angFreq() }, k_jj{ wave_jj.waveNumber() }, A_jj{ wave_jj.amp() };
	double cosBeta_ii{ wave_ii.cosBeta() }, sinBeta_ii{ wave_ii.sinBeta() };
	double cosBeta_jj{ wave_jj.cosBeta() }, sinBeta_jj{ wave_jj.sinBeta() };
	double beta_ii{ wave_ii.direction() * pi / 180. }, beta_jj{ wave_jj.direction() * pi / 180. };
	double phase_ii{ wave_ii.phase() * pi / 180. }, phase_jj{ wave_jj.phase() * pi / 180. };

	if (k_ii == 0 || k_jj == 0)
	{
		return fill::zeros;
	}

	// Second-order potential does not contribute to the mean force
	if (w_ii == w_jj)
	{
		return fill::zeros;
	}

	// Avoid recalculting cossine and sine functions
	double cosBetaPsi_ii{ cos(beta_ii - psi) }, cosBetaPsi_jj{ cos(beta_jj - psi) };
	double sinBetaPsi_ii{ sin(beta_ii - psi) }, sinBetaPsi_jj{ sin(beta_jj - psi) };

	// The integration is different for a horizontal cylinder, as the integration variable
	// along the cylinder is not related to the vertical position. This special case is treated separately.
	double mult_h_cos(0), mult_h_sin(0), mult_v_cos(0), mult_v_sin(0);;

	vec::fixed<2> km = { k_ii * cosBeta_ii - k_jj * cosBeta_jj, k_ii * sinBeta_ii - k_jj * sinBeta_jj };
	double k = arma::norm(km);
	double w = w_ii - w_jj;
	double I = tanAlpha * (k_ii*cosBetaPsi_ii - k_jj * cosBetaPsi_jj);
	double p1 = (k_ii*cosBeta_ii - k_jj * cosBeta_jj)*x1 + (k_ii*sinBeta_ii - k_jj * sinBeta_jj)*y1 + phase_ii - phase_jj;
	double p2 = (k_ii*cosBeta_ii - k_jj * cosBeta_jj)*x2 + (k_ii*sinBeta_ii - k_jj * sinBeta_jj)*y2 + phase_ii - phase_jj;
	double Ih = ((k_ii*cosBeta_ii - k_jj * cosBeta_jj)*(x2 - x1) + (k_ii*sinBeta_ii - k_jj * sinBeta_jj)*(y2 - y1)) / L;

	double aux = ((w_jj - w_ii) / (w_ii * w_jj)) * k_ii * k_jj * (std::cos(beta_ii - beta_jj) + std::tanh(k_ii*h) * std::tanh(k_jj*h))
		- 0.5 * (k_ii*k_ii / (w_ii * pow(std::cosh(k_ii*h), 2)) - k_jj * k_jj / (w_jj * pow(std::cosh(k_jj*h), 2)));
	aux = aux / (g * k * std::tanh(k * h) - w * w);
	double Q{ m_CM * 0.5 * A_ii * A_jj * g*g * aux * w }, Qz1{ Q * m_axialCa_1 / m_CM };
	double Qz2 = (evaluateTopNode ? Q * m_axialCa_2 / m_CM : 0);

	// Ratios of hyperbolic functions need a special treatment otherwise they yield inf in large water depth
	double cosh_sinh_1{ 0 }, cosh_sinh_2{ 0 }, cosh_cosh_1{ 0 }, cosh_cosh_2{ 0 }, sinh_sinh_1{ 0 }, sinh_sinh_2{ 0 }, sinh_cosh_1{ 0 }, sinh_cosh_2{ 0 };
	if (k*h >= 10)
	{
		cosh_sinh_1 = exp(k*z1);
		cosh_cosh_1 = cosh_sinh_1;
		sinh_sinh_1 = cosh_sinh_1;
		sinh_cosh_1 = cosh_sinh_1;

		cosh_sinh_2 = exp(k*z2);
		cosh_cosh_2 = cosh_sinh_2;
		sinh_sinh_2 = cosh_sinh_2;
		sinh_cosh_2 = cosh_sinh_2;
	}
	else
	{
		cosh_sinh_1 = cosh(k * (z1 + h)) / sinh(k*h);
		cosh_cosh_1 = cosh(k * (z1 + h)) / cosh(k*h);
		sinh_sinh_1 = sinh(k * (z1 + h)) / sinh(k*h);
		sinh_cosh_1 = sinh(k * (z1 + h)) / cosh(k*h);


		cosh_sinh_2 = cosh(k * (z2 + h)) / sinh(k*h);
		cosh_cosh_2 = cosh(k * (z2 + h)) / cosh(k*h);
		sinh_sinh_2 = sinh(k * (z2 + h)) / sinh(k*h);
		sinh_cosh_2 = sinh(k * (z2 + h)) / cosh(k*h);
	}

	/*
		FORCES
	*/
	if (arma::is_finite(tanAlpha))
	{
		// Contribution of the horizontal acceleration
		mult_h_cos = k * sin(p2) * sinh_cosh_2 - I * cos(p2) * cosh_cosh_2;
		mult_h_cos += -k * sin(p1) * sinh_cosh_1 + I * cos(p1) * cosh_cosh_1;
		mult_h_cos *= 1 / (I*I + k * k);

		mult_h_sin = k * cos(p2) * sinh_cosh_2 + I * sin(p2) *cosh_cosh_2;
		mult_h_sin += -k * cos(p1) * sinh_cosh_1 - I * sin(p1) * cosh_cosh_1;
		mult_h_sin *= 1 / (I*I + k * k);

		// This is due to the change of integration variable.
		// We know cosAlpha !=0 because horizontal cylinders are treated separately
		mult_h_cos *= 1 / cosAlpha;
		mult_h_sin *= 1 / cosAlpha;

		// Contribution of the vertical acceleration
		mult_v_cos = k * cos(p2) * cosh_cosh_2 + I * sin(p2) * sinh_cosh_2;
		mult_v_cos += -k * cos(p1) * cosh_cosh_1 - I * sin(p1) * sinh_cosh_1;
		mult_v_cos *= 1 / (I*I + k * k);

		mult_v_sin = k * sin(p2) * cosh_cosh_2 - I * cos(p2) * sinh_cosh_2;
		mult_v_sin += -k * sin(p1) * cosh_cosh_1 + I * cos(p1) * sinh_cosh_1;
		mult_v_sin *= 1 / (I*I + k * k);

		mult_v_cos *= 1 / cosAlpha;
		mult_v_sin *= -1 / cosAlpha;
	}
	else
	{
		// Contribution of the horizontal acceleration
		if (Ih != 0)
		{			
			mult_h_cos = -cos(p2) + cos(p1);
			mult_h_cos *= 1 / Ih;

			mult_h_sin = sin(p2) - sin(p1);
			mult_h_sin *= 1 / Ih;
		}
		else
		{
			mult_h_cos = sin(p1)*L;
			mult_h_sin = cos(p1)*L;
		}

		// Contribution of the vertical acceleration
		mult_v_cos = mult_h_sin * sinh_cosh_1;
		mult_v_sin = -mult_h_cos * sinh_sinh_1;

		// Need the factor of the water depth in the horizontal acceleration as well
		mult_h_cos *= cosh_cosh_1;
		mult_h_sin *= cosh_cosh_1;
	}

	// Forces along the local x and y axes
	double fx_cos = Q * ((k_ii*cosBetaPsi_ii - k_jj * cosBetaPsi_jj) * cosAlpha * mult_h_cos + k * sinAlpha * mult_v_cos);
	double fx_sin = Q * ((k_ii*cosBetaPsi_ii - k_jj * cosBetaPsi_jj) * cosAlpha * mult_h_sin + k * sinAlpha * mult_v_sin);
	double fy_cos = Q * (k_ii*sinBetaPsi_ii - k_jj * sinBetaPsi_jj) * mult_h_cos;
	double fy_sin = Q * (k_ii*sinBetaPsi_ii - k_jj * sinBetaPsi_jj) * mult_h_sin;

	// Forces along the local z axis. Factor 4*R/3. due to this force being acceleration * volume, while in the end 
	// of this function things are multiplied by the area of the circle
	double fz_cos = (4 * R / 3.) * (k_ii*cosBetaPsi_ii - k_jj * cosBetaPsi_jj) * sinAlpha * (Qz1 * cosh_cosh_1 * sin(p1) - Qz2 * cosh_cosh_2 * sin(p2)); // Contribution of horizontal acceleration
	double fz_sin = (4 * R / 3.) * (k_ii*cosBetaPsi_ii - k_jj * cosBetaPsi_jj) * sinAlpha * (-Qz1 * cosh_cosh_1 * cos(p1) + Qz2 * cosh_cosh_2 * cos(p2));
	fz_cos += -(4 * R / 3.) * cosAlpha * k * (Qz1 * sinh_cosh_1 * cos(p1) - Qz2 * sinh_cosh_2 * cos(p2)); // Contribution of vertical acceleration
	fz_sin += (4 * R / 3.) * cosAlpha * k * (Qz1 * sinh_cosh_1 * sin(p1) - Qz2 * sinh_cosh_2 * sin(p2));

	/*
		MOMENTS
	*/
	if (arma::is_finite(tanAlpha))
	{
		// Contribution of the horizontal acceleration
		mult_h_cos = k * sinh_cosh_2 * ((k*k + I * I)*(z2 - z1)*sin(p2) + 2 * I*cos(p2))
			+ cosh_cosh_2 * (-I * (k*k + I * I)*(z2 - z1)*cos(p2) + (I*I - k * k)*sin(p2));
		mult_h_cos += -k * sinh_cosh_1 * 2 * I*cos(p1)
			- cosh_cosh_1 * (I*I - k * k)*sin(p1);
		mult_h_cos *= 1 / (I*I + k * k) / (I*I + k * k);

		mult_h_sin = k * sinh_cosh_2 * ((k*k + I * I)*(z2 - z1)*cos(p2) - 2 * I*sin(p2))
			+ cosh_cosh_2 * (I*(k*k + I * I)*(z2 - z1)*sin(p2) + (I*I - k * k)*cos(p2));
		mult_h_sin += k * sinh_cosh_1 * 2 * I*sin(p1)
			- cosh_cosh_1 * (I*I - k * k)*cos(p1);
		mult_h_sin *= 1 / (I*I + k * k) / (I*I + k * k);

		mult_h_cos *= 1 / cosAlpha / cosAlpha;
		mult_h_sin *= 1 / cosAlpha / cosAlpha;

		// Contribution of the vertical acceleration
		mult_v_cos = k * cosh_cosh_2 * ((k*k + I * I)*(z2 - z1)*cos(p2) - 2 * I*sin(p2))
			+ sinh_cosh_2 * (I*(k*k + I * I)*(z2 - z1)*sin(p2) + (I*I - k * k)*cos(p2));
		mult_v_cos += k * cosh_cosh_1 * 2 * I*sin(p1)
			- sinh_cosh_1 * (I*I - k * k)*cos(p1);
		mult_v_cos *= 1 / (I*I + k * k) / (I*I + k * k);

		mult_v_sin = k * cosh_cosh_2 * ((k*k + I * I)*(z2 - z1)*sin(p2) + 2 * I*cos(p2))
			+ sinh_cosh_2 * (-I * (k*k + I * I)*(z2 - z1)*cos(p2) + (I*I - k * k)*sin(p2));
		mult_v_sin += -k * cosh_cosh_1 * 2 * I*cos(p1)
			- sinh_cosh_1 * (I*I - k * k)*sin(p1);
		mult_v_sin *= 1 / (I*I + k * k) / (I*I + k * k);

		// This is due to the change of integration variable
		mult_v_cos *= 1 / cosAlpha / cosAlpha;
		mult_v_sin *= -1 / cosAlpha / cosAlpha;
	}
	else
	{
		// Contribution of the horizontal acceleration
		if (Ih != 0)
		{
			mult_h_cos = -Ih * L*cos(p2) + sin(p2) - sin(p1);
			mult_h_cos *= 1 / (Ih*Ih);

			mult_h_sin = Ih * L*sin(p2) + cos(p2) - cos(p1);
			mult_h_sin *= 1 / (Ih*Ih);
		}
		else
		{
			mult_h_cos = sin(p1)*L*L*0.5;
			mult_h_sin = cos(p1)*L*L*0.5;
		}

		// Contribution of the vertical acceleration
		mult_v_cos = mult_h_sin * sinh_cosh_1;
		mult_v_sin = -mult_h_cos * sinh_cosh_1;

		// Need the factor of the water depth in the horizontal acceleration as well
		mult_h_cos *= cosh_cosh_1;
		mult_h_sin *= cosh_cosh_1;
	}

	// Forces along the local x and y axes
	double mx_cos = -Q * (k_ii*sinBetaPsi_ii - k_jj * sinBetaPsi_jj) * mult_h_cos;
	double mx_sin = -Q * (k_ii*sinBetaPsi_ii - k_jj * sinBetaPsi_jj) * mult_h_sin;
	double my_cos = Q * ((k_ii*cosBetaPsi_ii - k_jj * cosBetaPsi_jj) * cosAlpha * mult_h_cos + k * sinAlpha * mult_v_cos);
	double my_sin = Q * ((k_ii*cosBetaPsi_ii - k_jj * cosBetaPsi_jj) * cosAlpha * mult_h_sin + k * sinAlpha * mult_v_sin);

	if (m_botPressFlag)
	{
		fz_cos += Q / m_CM * cosh_cosh_1 * cos(p1);
		fz_sin += Q / m_CM * cosh_cosh_1 * sin(p1);
		if (evaluateTopNode)
		{
			fz_cos -= Q / m_CM * cosh_cosh_2 * cos(p2);
			fz_sin -= Q / m_CM * cosh_cosh_2 * sin(p2);
		}
	}

	cx_vec::fixed<6> coef(fill::zeros);
	for (int ii = 0; ii < 3; ++ii)
	{
		coef.at(ii) = cx_double{ fx_cos * xvec(ii) + fy_cos * yvec(ii) + fz_cos * zvec(ii) , fx_sin * xvec(ii) + fy_sin * yvec(ii) + fz_sin * zvec(ii) };
	}
	for (int ii = 3; ii < 6; ++ii)
	{
		coef.at(ii) = cx_double{ mx_cos * xvec(ii - 3) + my_cos * yvec(ii - 3) , mx_sin * xvec(ii - 3) + my_sin * yvec(ii - 3) };
	}

	return rho * pi* R * R * coef;
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
	vec::fixed<3> n1 = m_node1Pos_1stOrd;
	vec::fixed<3> n2 = m_node2Pos_1stOrd;
	vec::fixed<3> xvec = m_xvec_1stOrd; // xvec and yvec are used to project the acceleration. They are analogous to the normal vector.
	vec::fixed<3> yvec = m_yvec_1stOrd;
	vec::fixed<3> zvec = m_zvec_1stOrd; // While zvec is used only to evaluate the nodes position
	if (hydroMode < 2)
	{
		n1 = node1Pos_sd();
		n2 = node2Pos_sd();
		xvec = m_xvec_sd;
		yvec = m_yvec_sd;
		zvec = m_zvec_sd;
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

mat::fixed<6, 6> MorisonCirc::addedMass_paral(const double rho, const vec::fixed<3> &refPt, const int hydroMode) const
{
	mat::fixed<6, 6> A(fill::zeros);

	// Nodes position and vectors of the local coordinate system vectors
	vec::fixed<3> n1 = m_node1Pos_1stOrd;
	vec::fixed<3> n2 = m_node1Pos_1stOrd;
	vec::fixed<3> zvec = m_zvec_1stOrd;
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
	output = output + "node1 at t=0:\t(" + std::to_string(m_node1Pos(0)) + ", " + std::to_string(m_node1Pos(1)) + ", " + std::to_string(m_node1Pos(2)) + ")\n";
	output = output + "node2 at t=0:\t(" + std::to_string(m_node2Pos(0)) + ", " + std::to_string(m_node2Pos(1)) + ", " + std::to_string(m_node2Pos(2)) + ")\n";
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
