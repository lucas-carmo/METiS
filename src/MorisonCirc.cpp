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
	vec::fixed<3> z_global = { 0,0,1 };
	if (approx_equal(zvec, z_global, "absdiff", datum::eps))
	{
		xvec = { 1, 0, 0 };
		yvec = { 0, 1, 0 };
	}
	else
	{
		yvec = arma::cross(z_global, zvec);
		yvec = yvec / norm(yvec, 2);

		xvec = cross(yvec, zvec);
		xvec = xvec / norm(xvec, 2);
	}

}


// TODO: depois de debugar direitinho, tirar os bound checks (usar [] ao inves de () pra acessar elementos das matrizes)
mat::fixed<6, 6> MorisonCirc::addedMass_perp(const double rho, const int hydroMode) const
{
	mat::fixed<6, 6> A(fill::zeros);

	// Use a more friendly notation
	double Lambda = datum::pi * pow(m_diam / 2., 2) * rho * (m_CM - 1);
	double ncyl = m_numIntPoints;

	// Nodes position
	vec::fixed<3> n1 = node1Pos();
	vec::fixed<3> n2 = node2Pos();
	if (hydroMode == 1) // Check if the hydrodynamic force should be calculated considering the initial position of the floater
	{
		n1 = m_node1Pos_t0;
		n2 = m_node2Pos_t0;
	}

	// Center of Gravity
	double xG = n1[0] - m_cog2node1[0];
	double yG = n1[1] - m_cog2node1[1];
	double zG = n1[2] - m_cog2node1[2];

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


	// Vectors of the global coordinate system - arranjed, they result in an eye matrix
	mat::fixed<3, 3> globalBase(fill::eye);

	// The purely translational elements of the matrix (Aij for i,j = 1, 2, 3) are integrated analytically
	for (int pp = 0; pp < 3; ++pp)
	{
		for (int qq = pp; qq < 3; ++qq)
		{
			A(pp, qq) = Lambda * L * ( dot(xvec, globalBase.col(pp))*dot(xvec, globalBase.col(qq))
							         + dot(yvec, globalBase.col(pp))*dot(yvec, globalBase.col(qq))
									 );
		}
	}

	vec::fixed<3> n_ii;
	double x_ii{ 0 };
	double y_ii{ 0 };
	double z_ii{ 0 };
	double step{ 0 }; // Used for Simpson's rule. See below.
	for (int ii = 1; ii <= ncyl; ++ii)
	{
		n_ii = (n2 - n1) * (ii - 1) / (ncyl - 1) + n1; // Coordinates of the integration point

		if (n_ii[2] > 0)
		{
			break; // Since n2 is above n1, if we reach n_ii[2] > 0, all the next n_ii are also above the waterline
		}

		x_ii = n_ii[0];
		y_ii = n_ii[1];
		z_ii = n_ii[2];

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

		A(3, 3) += Lambda * step *
				 (
				   pow(y_ii - yG, 2) * ( pow(dot(xvec, globalBase.col(2)), 2) + pow(dot(yvec, globalBase.col(2)), 2) )
			     + pow(z_ii - zG, 2) * ( pow(dot(xvec, globalBase.col(1)), 2) + pow(dot(yvec, globalBase.col(1)), 2) )
				 - 2 * (y_ii - yG) * (z_ii - zG)  * ( dot(xvec, globalBase.col(1)) * dot(xvec, globalBase.col(2)) + dot(yvec, globalBase.col(1)) * dot(yvec, globalBase.col(2)) )
				 );

		A(4, 4) += Lambda * step *
				(
				  pow(x_ii - xG, 2) * (pow(dot(xvec, globalBase.col(2)), 2) + pow(dot(yvec, globalBase.col(2)), 2))
				+ pow(z_ii - zG, 2) * (pow(dot(xvec, globalBase.col(0)), 2) + pow(dot(yvec, globalBase.col(0)), 2))
				- 2 * (x_ii - xG) * (z_ii - zG)  * (dot(xvec, globalBase.col(0)) * dot(xvec, globalBase.col(2)) + dot(yvec, globalBase.col(0)) * dot(yvec, globalBase.col(2)))
				);

		A(5, 5) += Lambda * step *
				(
				  pow(x_ii - xG, 2) * (pow(dot(xvec, globalBase.col(1)), 2) + pow(dot(yvec, globalBase.col(1)), 2))
				+ pow(y_ii - yG, 2) * (pow(dot(xvec, globalBase.col(0)), 2) + pow(dot(yvec, globalBase.col(0)), 2))
				- 2 * (x_ii - xG) * (y_ii - yG)  * (dot(xvec, globalBase.col(0)) * dot(xvec, globalBase.col(1)) + dot(yvec, globalBase.col(0)) * dot(yvec, globalBase.col(1)))
				);

		A(0, 3) += Lambda * step *
			    (
				  (y_ii - yG) * (dot(xvec, globalBase.col(0))*dot(xvec, globalBase.col(2)) + dot(yvec, globalBase.col(0))*dot(yvec, globalBase.col(2)))
				- (z_ii - zG) * (dot(xvec, globalBase.col(0))*dot(xvec, globalBase.col(1)) + dot(yvec, globalBase.col(0))*dot(yvec, globalBase.col(1)))
				);

		A(0, 4) += Lambda * step *
				(
				  (z_ii - zG) * (pow(dot(xvec, globalBase.col(0)), 2) + pow(dot(yvec, globalBase.col(0)), 2))
				- (x_ii - xG) * (dot(xvec, globalBase.col(0))*dot(xvec, globalBase.col(2)) + dot(yvec, globalBase.col(0))*dot(yvec, globalBase.col(2)))
				);

		A(0, 5) += Lambda * step *
				(
				  (x_ii - xG) * (dot(xvec, globalBase.col(0))*dot(xvec, globalBase.col(1)) + dot(yvec, globalBase.col(0))*dot(yvec, globalBase.col(1)))
				- (y_ii - yG) * (pow(dot(xvec, globalBase.col(0)), 2) + pow(dot(yvec, globalBase.col(0)), 2))
				);

		A(1, 3) += Lambda * step *
				(
				  (y_ii - yG) * (dot(xvec, globalBase.col(1))*dot(xvec, globalBase.col(2)) + dot(yvec, globalBase.col(1))*dot(yvec, globalBase.col(2)))
				- (z_ii - zG) * (pow(dot(xvec, globalBase.col(1)), 2) + pow(dot(yvec, globalBase.col(1)), 2))
				);

		A(1, 4) += Lambda * step *
				(
				  (z_ii - zG) * (dot(xvec, globalBase.col(0))*dot(xvec, globalBase.col(1)) + dot(yvec, globalBase.col(0))*dot(yvec, globalBase.col(1)))
				- (x_ii - xG) * (dot(xvec, globalBase.col(1))*dot(xvec, globalBase.col(2)) + dot(yvec, globalBase.col(1))*dot(yvec, globalBase.col(2)))
				);

		A(1, 5) += Lambda * step *
				(
				  (x_ii - xG) * (pow(dot(xvec, globalBase.col(1)), 2) + pow(dot(yvec, globalBase.col(1)), 2))
				- (y_ii - yG) * (dot(xvec, globalBase.col(0))*dot(xvec, globalBase.col(1)) + dot(yvec, globalBase.col(0))*dot(yvec, globalBase.col(1)))
				);

		A(2, 3) += Lambda * step *
				(
				  (y_ii - yG) * (pow(dot(xvec, globalBase.col(2)), 2) + pow(dot(yvec, globalBase.col(2)), 2))
				- (z_ii - zG) * (dot(xvec, globalBase.col(1))*dot(xvec, globalBase.col(2)) + dot(yvec, globalBase.col(1))*dot(yvec, globalBase.col(2)))
				);

		A(2, 4) += Lambda * step *
				(
				  (z_ii - zG) * (dot(xvec, globalBase.col(0))*dot(xvec, globalBase.col(2)) + dot(yvec, globalBase.col(0))*dot(yvec, globalBase.col(2)))
				- (x_ii - xG) * (pow(dot(xvec, globalBase.col(2)), 2) + pow(dot(yvec, globalBase.col(2)), 2))
				);

		A(2, 5) += Lambda * step *
				(
				  (x_ii - xG) * (dot(xvec, globalBase.col(1))*dot(xvec, globalBase.col(2)) + dot(yvec, globalBase.col(1))*dot(yvec, globalBase.col(2)))
				- (y_ii - yG) * (dot(xvec, globalBase.col(0))*dot(xvec, globalBase.col(2)) + dot(yvec, globalBase.col(0))*dot(yvec, globalBase.col(2)))
				);

		A(3, 4) += Lambda * step *
				(
				- (x_ii - xG) * (y_ii - yG) * (pow(dot(xvec, globalBase.col(2)), 2) + pow(dot(yvec, globalBase.col(2)), 2))
				- (x_ii - xG) * (z_ii - zG) * (dot(xvec, globalBase.col(1))*dot(xvec, globalBase.col(2)) + dot(yvec, globalBase.col(1))*dot(yvec, globalBase.col(2)))
				+ (y_ii - yG) * (z_ii - zG) * (dot(xvec, globalBase.col(0))*dot(xvec, globalBase.col(2)) + dot(yvec, globalBase.col(0))*dot(yvec, globalBase.col(2)))
				- pow(z_ii - zG, 2) * (dot(xvec, globalBase.col(0))*dot(xvec, globalBase.col(1)) + dot(yvec, globalBase.col(0))*dot(yvec, globalBase.col(1)))
				);

		A(3, 5) += Lambda * step *
				(
				  (x_ii - xG) * (y_ii - yG) * (dot(xvec, globalBase.col(1))*dot(xvec, globalBase.col(2)) + dot(yvec, globalBase.col(1))*dot(yvec, globalBase.col(2)))
				- (x_ii - xG) * (z_ii - zG) * (pow(dot(xvec, globalBase.col(1)), 2) + pow(dot(yvec, globalBase.col(1)), 2))
				- pow(y_ii - yG, 2) * (dot(xvec, globalBase.col(0))*dot(xvec, globalBase.col(2)) + dot(yvec, globalBase.col(0))*dot(yvec, globalBase.col(2)))
				+ (y_ii - yG) * (z_ii - zG) * (dot(xvec, globalBase.col(0))*dot(xvec, globalBase.col(1)) + dot(yvec, globalBase.col(0))*dot(yvec, globalBase.col(1)))
				);

		A(4, 5) += Lambda * step *
				(
				- pow(x_ii - xG, 2) * (dot(xvec, globalBase.col(1))*dot(xvec, globalBase.col(2)) + dot(yvec, globalBase.col(1))*dot(yvec, globalBase.col(2)))
				+ (x_ii - xG) * (y_ii - yG) * (dot(xvec, globalBase.col(0))*dot(xvec, globalBase.col(2)) + dot(yvec, globalBase.col(0))*dot(yvec, globalBase.col(2)))
				+ (x_ii - xG) * (y_ii - yG) * (dot(xvec, globalBase.col(0))*dot(xvec, globalBase.col(1)) + dot(yvec, globalBase.col(0))*dot(yvec, globalBase.col(1)))
				- (y_ii - yG) * (z_ii - zG) * (pow(dot(xvec, globalBase.col(0)), 2) + pow(dot(yvec, globalBase.col(0)), 2))
				);
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



// TODO: depois de debugar direitinho, tirar os bound checks (usar [] ao inves de () pra acessar elementos das matrizes)
mat::fixed<6, 6> MorisonCirc::addedMass_paral(const double rho, const int hydroMode) const
{
	mat::fixed<6, 6> A(fill::zeros);

	// Use a more friendly notation
	double Lambda = rho * m_axialCa * (4 / 3.) * datum::pi * pow(m_diam / 2, 3);
	double ncyl = m_numIntPoints;

	// Nodes position
	vec::fixed<3> n1 = node1Pos();
	vec::fixed<3> n2 = node2Pos();
	if (hydroMode == 1) // Check if the hydrodynamic force should be calculated considering the initial position of the floater
	{
		n1 = m_node1Pos_t0;
		n2 = m_node2Pos_t0;
	}

	// Center of Gravity
	double xG = n1[0] - m_cog2node1[0];
	double yG = n1[1] - m_cog2node1[1];
	double zG = n1[2] - m_cog2node1[2];

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
			A(pp, qq) = Lambda * dot(zvec, globalBase.col(pp))*dot(zvec, globalBase.col(qq));
		}
	}

	A(3, 3) += Lambda *
			(
			  pow(y_ii - yG, 2) * pow(dot(zvec, globalBase.col(2)), 2)
			+ pow(z_ii - zG, 2) * pow(dot(zvec, globalBase.col(1)), 2)
			- 2 * (y_ii - yG) * (z_ii - zG) * dot(zvec, globalBase.col(1)) * dot(zvec, globalBase.col(2))
			);

	A(4, 4) += Lambda *
			(
			  pow(x_ii - xG, 2) * pow(dot(zvec, globalBase.col(2)), 2)
			+ pow(z_ii - zG, 2) * pow(dot(zvec, globalBase.col(0)), 2)
			- 2 * (x_ii - xG) * (z_ii - zG)  * dot(zvec, globalBase.col(0)) * dot(zvec, globalBase.col(2))
			);

	A(5, 5) += Lambda *
			(
			  pow(x_ii - xG, 2) * pow(dot(zvec, globalBase.col(1)), 2)
			+ pow(y_ii - yG, 2) * pow(dot(zvec, globalBase.col(0)), 2)
			- 2 * (x_ii - xG) * (y_ii - yG)  * dot(zvec, globalBase.col(0)) * dot(zvec, globalBase.col(1))
			);

	A(0, 3) += Lambda *
			(
			  (y_ii - yG) * dot(zvec, globalBase.col(0))*dot(zvec, globalBase.col(2))
			- (z_ii - zG) * dot(zvec, globalBase.col(0))*dot(zvec, globalBase.col(1))
			);

	A(0, 4) += Lambda *
			(
			  (z_ii - zG) * pow(dot(zvec, globalBase.col(0)), 2)
			- (x_ii - xG) * dot(zvec, globalBase.col(0))*dot(zvec, globalBase.col(2))
			);

	A(0, 5) += Lambda *
			(
			  (x_ii - xG) * dot(zvec, globalBase.col(0))*dot(zvec, globalBase.col(1))
			- (y_ii - yG) * pow(dot(zvec, globalBase.col(0)), 2)
			);

	A(1, 3) += Lambda *
			(
			  (y_ii - yG) * dot(zvec, globalBase.col(1))*dot(zvec, globalBase.col(2))
			- (z_ii - zG) * pow(dot(zvec, globalBase.col(1)), 2)
			);

	A(1, 4) += Lambda *
			(
			  (z_ii - zG) * dot(zvec, globalBase.col(0))*dot(zvec, globalBase.col(1))
			- (x_ii - xG) * dot(zvec, globalBase.col(1))*dot(zvec, globalBase.col(2))
			);

	A(1, 5) += Lambda *
			(
			  (x_ii - xG) * pow(dot(zvec, globalBase.col(1)), 2)
			- (y_ii - yG) * dot(zvec, globalBase.col(0))*dot(zvec, globalBase.col(1))
			);

	A(2, 3) += Lambda *
			(
			  (y_ii - yG) * pow(dot(zvec, globalBase.col(2)), 2)
			- (z_ii - zG) * dot(zvec, globalBase.col(1))*dot(zvec, globalBase.col(2))
			);

	A(2, 4) += Lambda *
			(
			  (z_ii - zG) * dot(zvec, globalBase.col(0))*dot(zvec, globalBase.col(2))
			- (x_ii - xG) * pow(dot(zvec, globalBase.col(2)), 2)
			);

	A(2, 5) += Lambda *
			(
			  (x_ii - xG) * dot(zvec, globalBase.col(1))*dot(zvec, globalBase.col(2))
			- (y_ii - yG) * dot(zvec, globalBase.col(0))*dot(zvec, globalBase.col(2))
			);

	A(3, 4) += Lambda *
			(
			- (x_ii - xG) * (y_ii - yG) * pow(dot(zvec, globalBase.col(2)), 2)
			- (x_ii - xG) * (z_ii - zG) * dot(zvec, globalBase.col(1))*dot(zvec, globalBase.col(2))
			+ (y_ii - yG) * (z_ii - zG) * dot(zvec, globalBase.col(0))*dot(zvec, globalBase.col(2))
			- pow(z_ii - zG, 2) * dot(zvec, globalBase.col(0))*dot(zvec, globalBase.col(1))
			);

	A(3, 5) += Lambda *
			(
			  (x_ii - xG) * (y_ii - yG) * dot(zvec, globalBase.col(1))*dot(zvec, globalBase.col(2))
			- (x_ii - xG) * (z_ii - zG) * pow(dot(zvec, globalBase.col(1)), 2)
			- pow(y_ii - yG, 2) * dot(zvec, globalBase.col(0))*dot(zvec, globalBase.col(2))
			+ (y_ii - yG) * (z_ii - zG) * dot(zvec, globalBase.col(0))*dot(zvec, globalBase.col(1))
			);

	A(4, 5) += Lambda *
			(
			- pow(x_ii - xG, 2) * dot(zvec, globalBase.col(1))*dot(zvec, globalBase.col(2))
			+ (x_ii - xG) * (y_ii - yG) * dot(zvec, globalBase.col(0))*dot(zvec, globalBase.col(2))
			+ (x_ii - xG) * (y_ii - yG) * dot(zvec, globalBase.col(0))*dot(zvec, globalBase.col(1))
			- (y_ii - yG) * (z_ii - zG) * pow(dot(zvec, globalBase.col(0)), 2)
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
	double alpha = acos(dot(zvec, arma::vec::fixed<3> {0, 0, 1})); // zvec and {0, 0, 1} are both unit vectors
	double tanAlpha{ 0 };

	// Check if the angle is 90 degrees
	if (std::abs(alpha - arma::datum::pi/2) > datum::eps)
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

		xb = tanAlpha * pow(D/2, 2) / (4 * L);
		yb = 0;
		zb = (pow(tanAlpha*D/2, 2) + 4 * pow(L, 2)) / (8 * L);
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


// The three vectors passed as reference are used to return the different components of the hydrodynamic force acting on the cylinder,
// without the need of calling three different methods for each component of the hydrodynamic force
vec::fixed<6> MorisonCirc::hydrodynamicForce(const ENVIR &envir, const int hydroMode, vec::fixed<6> &force_inertia, vec::fixed<6> &force_drag, vec::fixed<6> &force_froudeKrylov) const
{
	// Forces and moments acting at the Morison Element
	vec::fixed<6> force(fill::zeros);

	// Use a more friendly notation
	double D = m_diam;
	double Cd = m_CD;
	double Cm = m_CM;
	double Cd_V = m_axialCD;
	double Ca_V = m_axialCa;
	double ncyl = m_numIntPoints;
	double rho = envir.watDensity();
	double botDiam = m_botDiam;
	double topDiam = m_topDiam;


	// Nodes position, velocity and acceleration
	vec::fixed<3> n1 = node1Pos();
	vec::fixed<3> n2 = node2Pos();
	if (hydroMode == 1) // Check if the hydrodynamic force should be calculated considering the initial position of the floater
	{
		n1 = m_node1Pos_t0;
		n2 = m_node2Pos_t0;
	}
	vec::fixed<3> v1 = node1Vel();
	vec::fixed<3> v2 = node2Vel();
	vec::fixed<3> a1 = node1AccCentrip();
	vec::fixed<3> a2 = node2AccCentrip();

	// Vectors of the local coordinate system vectors
	vec::fixed<3> xvec(fill::zeros);
	vec::fixed<3> yvec(fill::zeros);
	vec::fixed<3> zvec(fill::zeros);
	MorisonCirc::make_local_base(xvec, yvec, zvec);

	// Bottom and top diameter
	botDiam = m_botDiam;
	topDiam = m_topDiam;

	// Make sure that node1 is below node2 (or at the same height, at least).
	// Otherwise, need to swap them.
	if (n1[2] > n2[2])
	{
		n1.swap(n2);
		v1.swap(v2);
		a1.swap(a2);

		// If the nodes position is reversed, we need to reverse the local base as well
		xvec = -xvec;
		yvec = -yvec;
		zvec = -zvec;

		// Bottom diameter and top diameter must be swapped as well
		botDiam = m_topDiam;
		topDiam = m_botDiam;
	}

	// length of each interval between points
	double dL = norm(n2 - n1, 2) / (ncyl - 1);

	/*
		First part: forces on the length of the cylinder
	*/
	//  Loop to calculate the force/moment at each integration point
	vec::fixed<3> n_ii(fill::zeros); // Coordinates of the integration point
	vec::fixed<3> vel_ii(fill::zeros); // Velocity of the integration point
	vec::fixed<3> acc_ii(fill::zeros); // Acceleration of the integration point
	vec::fixed<3> u1(fill::zeros); // Fluid velocity at the integration point
	vec::fixed<3> du1dt(fill::zeros); // Fluid acceleration at the integration point

	// Forces acting at the integration point and moment (with relation to n1) due to the force acting at the integration point
	vec::fixed<3> force_inertia_ii(fill::zeros); // Inertial component
	vec::fixed<3> moment_inertia_ii(fill::zeros);

	vec::fixed<3> force_drag_ii(fill::zeros); // Drag component
	vec::fixed<3> moment_drag_ii(fill::zeros);

	vec::fixed<3> force_ii(fill::zeros); // Total force
	vec::fixed<3> moment_ii(fill::zeros);

	for (int ii = 1; ii <= ncyl; ++ii)
	{
		n_ii = (n2 - n1) * (ii-1)/(ncyl-1) + n1; // Coordinates of the integration point
		if (n_ii[2] > 0)
		{
			break; // Since n2 is above n1, if we reach n_ii[2] > 0, all the next n_ii are also above the waterline
		}

		/******
			Body velocity/acceleration
		******/
		// Absolute (R_ii) and relative (lambda) distance between the integration point and the bottom node
		double R_ii = norm(n_ii - n1, 2);
		double lambda = R_ii/norm(n2 - n1, 2);

		// Velocity and acceleration of the integration point
		vel_ii = v1 + lambda * ( v2 - v1 );
		acc_ii = a1 + lambda * ( a2 - a1 );

		// Component of the velocity and acceleration of the integration point that is perpendicular to the axis of the cylinder
		vel_ii = dot(vel_ii, xvec) * xvec + dot(vel_ii, yvec) * yvec;
		acc_ii = dot(acc_ii, xvec) * xvec + dot(acc_ii, yvec) * yvec;


		/*******
			Fluid velocity/acceleration
		******/
		// Fluid velocity and acceleration at the integration point
		u1 = envir.u1(n_ii);
		du1dt = envir.du1dt(n_ii);

		// Component of the fluid velocity and acceleration at the integration point that is perpendicular to the axis of the cylinder - written in the GLOBAL reference frame
		u1 = dot(u1, xvec) * xvec + dot(u1, yvec) * yvec;
		du1dt = dot(du1dt, xvec) * xvec + dot(du1dt, yvec) * yvec;

		// Calculation of the forces in the integration node using Morison's equation
		force_inertia_ii = (datum::pi * pow(D,2)/4) * rho * Cm * du1dt - (datum::pi * pow(D,2)/4) * rho * (Cm-1) * du1dt;
		force_drag_ii = 0.5 * rho * Cd * D * norm(u1 - vel_ii, 2) * (u1 - vel_ii);
		force_ii = force_inertia_ii + force_drag_ii;

		// Calculation of the moment (with respect to node 1) at the integration node using the force calculated above
		moment_inertia_ii = cross(R_ii * zvec, force_inertia_ii);
		moment_drag_ii = cross(R_ii * zvec, force_drag_ii);
		moment_ii = moment_inertia_ii + moment_drag_ii;


		// Integrate the forces along the cylinder using Simpson's Rule
		if (ii == 1 || ii == ncyl)
		{
			force += (dL/3.0) * join_cols(force_ii, moment_ii);
			force_inertia += (dL/3.0) * join_cols(force_inertia_ii, moment_inertia_ii);
			force_drag += (dL/3.0) * join_cols(force_drag_ii, moment_drag_ii);
		}
		else if (ii % 2 == 0)
		{
			force += (4*dL/3.0) * join_cols(force_ii, moment_ii);
			force_inertia += (4*dL/3.0) * join_cols(force_inertia_ii, moment_inertia_ii);
			force_drag += (4*dL/3.0) * join_cols(force_drag_ii, moment_drag_ii);
		}
		else
		{
			force += (2*dL/3.0) * join_cols(force_ii, moment_ii);
			force_inertia += (2*dL/3.0) * join_cols(force_inertia_ii, moment_inertia_ii);
			force_drag += (2*dL/3.0) * join_cols(force_drag_ii, moment_drag_ii);
		}
	}


	/*
		Second part: forces on the bottom of the cylinder
	*/
	// Component of the fluid velocity/acceleration at the bottom node that is parallel to the cylinder axis
	u1 = dot(envir.u1(n1), zvec) * zvec;
	du1dt = dot(envir.du1dt(n1), zvec) * zvec;

	// Component of the velocity and acceleration of the bottom that is perpendicular to the axis of the cylinder
	vec::fixed<3> v_axial = dot(v1, zvec) * zvec;
	vec::fixed<3> a_axial = dot(a1, zvec) * zvec;

	// Calculate the force acting on the bottom of the cylinder
	force.rows(0,2) += 0.5 * rho * Cd_V * datum::pi * pow(D/2, 2) * norm(u1 - v_axial, 2) * (u1 - v_axial)
					 + rho * Ca_V * (4/3.) * datum::pi * pow(D/2, 3) * (du1dt - a_axial);

	// These two are used to output the force components. They include the components from all the different
	// cylinders (note that they are passed by reference as an input).
	force_inertia.rows(0,2) += rho * Ca_V * (4 / 3.) * datum::pi * pow(D / 2, 3) * (du1dt - a_axial);
	force_drag.rows(0, 2) += 0.5 * rho * Cd_V * datum::pi * pow(D / 2, 2) * norm(u1 - v_axial, 2) * (u1 - v_axial);

	if (m_botPressFlag)
	{
		force.rows(0, 2) += datum::pi * ( pow(botDiam/2, 2) * envir.wavePressure(n1)
							- (pow(botDiam/2, 2) - pow(topDiam/2, 2)) * envir.wavePressure(n2) ) * zvec;

		force_froudeKrylov.rows(0, 2) += datum::pi * (pow(botDiam / 2, 2) * envir.wavePressure(n1)
										- (pow(botDiam / 2, 2) - pow(topDiam / 2, 2)) * envir.wavePressure(n2) ) * zvec;
	}

	// The moment was calculated with relation to n1, which may be different from node1.
	// We need to change the fulcrum to node1
	force.rows(3,5) = force.rows(3,5) + cross( n1 - node1Pos(), force.rows(0,2) );

	force_inertia.rows(3, 5) = force_inertia.rows(3, 5) + cross(n1 - node1Pos(), force_inertia.rows(0, 2));
	force_drag.rows(3, 5) = force_drag.rows(3, 5) + cross(n1 - node1Pos(), force_drag.rows(0, 2));

	return force;
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
