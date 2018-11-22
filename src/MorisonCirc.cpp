#include "MorisonCirc.h"


using namespace arma;

/*****************************************************
	Constructors
*****************************************************/
MorisonCirc::MorisonCirc(vec cog2node1, vec cog2node2, int numIntPoints, 
						 bool botPressFlag, double axialCD, double axialCa,
						 double diam, double CD, double CM, double botDiam, double topDiam)
	: MorisonElement(cog2node1, cog2node2, numIntPoints, botPressFlag, axialCD, axialCa), 
					 m_diam(diam), m_CD(CD), m_CM(CM), m_botDiam(botDiam), m_topDiam(topDiam)
{}



/*****************************************************
	Forces acting on the Morison Element and functions for node position/velocity/acceleration)
*****************************************************/
vec::fixed<3> MorisonCirc::node1Pos(const vec::fixed<6> &floaterPos) const
{	
	// Fazer uma funcao que calcula a matriz de rotacao
	// mat::fixed<3,3> rotatMatrix(const vec::fixed<6> &FOWTpos) const
	double x = floaterPos.at(0) + m_cog2node1.at(0);
	double y = floaterPos.at(1) + m_cog2node1.at(1);
	double z = floaterPos.at(2) + m_cog2node1.at(2);

	vec::fixed<3> pos = { x, y, z };
	return pos;
}

vec::fixed<3> MorisonCirc::node2Pos(const vec::fixed<6> &floaterPos) const
{	
	double x = floaterPos.at(0) + m_cog2node2.at(0);
	double y = floaterPos.at(1) + m_cog2node2.at(1);
	double z = floaterPos.at(2) + m_cog2node2.at(2);

	vec::fixed<3> pos = { x, y, z };
	return pos;
}

vec::fixed<3> MorisonCirc::node1Vel(const vec::fixed<6> &floaterPos, const vec::fixed<6> &floaterVel) const
{	
	double x = 0;
	double y = 0;
	double z = 0;

	vec::fixed<3> vel = { x, y, z };
	return vel;
}

vec::fixed<3> MorisonCirc::node2Vel(const vec::fixed<6> &floaterPos, const vec::fixed<6> &floaterVel) const
{	
	double x = 0;
	double y = 0;
	double z = 0;

	vec::fixed<3> vel = { x, y, z };
	return vel;
}

vec::fixed<3> MorisonCirc::node1Acc(const vec::fixed<6> &floaterPos, const vec::fixed<6> &floaterVel, const vec::fixed<6> &floaterAcc) const
{
	double x = 0;
	double y = 0;
	double z = 0;

	vec::fixed<3> acc = { x, y, z };
	return acc;
}

vec::fixed<3> MorisonCirc::node2Acc(const vec::fixed<6> &floaterPos, const vec::fixed<6> &floaterVel, const vec::fixed<6> &floaterAcc) const
{
	double x = 0;
	double y = 0;
	double z = 0;

	vec::fixed<3> acc = { x, y, z };
	return acc;
}

void MorisonCirc::make_local_base(arma::vec::fixed<3> &xvec, arma::vec::fixed<3> &yvec, arma::vec::fixed<3> &zvec) const
{
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

vec::fixed<6> MorisonCirc::hydrostaticForce(const ENVIR &envir, const vec::fixed<6> &floaterPos) const
{
	// Forces and moments acting at the Morison Element
	vec::fixed<6> force(fill::zeros);

	// Use a more friendly notation    
	double D = m_diam;
	double rho = envir.watDensity();
	double g = envir.gravity();

	// Nodes position, velocity and acceleration
	vec::fixed<3> n1 = node1Pos(floaterPos);
	vec::fixed<3> n2 = node2Pos(floaterPos);

	// If both nodes are above the water line, then the hydrostatic force is zero
	if (n1(3) >= 0 && n2(3) >= 0)
	{
		return force;
	}


	// Vectors of the local coordinate system vectors
	vec::fixed<3> xvec(fill::zeros);
	vec::fixed<3> yvec(fill::zeros);
	vec::fixed<3> zvec = (n2 - n1) / arma::norm(n2 - n1, 2);
	MorisonCirc::make_local_base(xvec, yvec, zvec);

	// Calculation of the inclination of the cylinder (with respect to the
	// vertical), which is used in the calculation of the center of buoyoancy
	double alpha = acos(dot(zvec, arma::vec::fixed<3> {0,0,1})); // zvec and {0, 0, 1} are both unit vectors
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


	// If only one of the nodes is above the water line, its coordinates are changed
	// by those of the intersection between the cylinder axis and the static
	// water line(defined by z_global = 0)
	if (n1(3) >= 0)
	{
		n1 = n2 + std::abs(0 - n2(3)) / (n1(3) - n2(3)) * norm(n2 - n1) * zvec;
	}

	if (n2(3) >= 0)
	{
		n2 = n1 + std::abs(0 - n1(3)) / (n2(3) - n1(3)) * norm(n2 - n1) * zvec;
	}

	// Length of the cylinder
	double L = norm(n2 - n1);

	// Volume of fluid displaced by the cylinder
	double Vol = arma::datum::pi * pow(D/2, 2) * L;

	// If the cylinder is below the water line, the center of volume is at the
	// center of the cylinder
	double xb{ 0 };
	double yb{ 0 };
	double zb{ 0 };
	if (n1(3) < 0 && n2(3) < 0)
	{
		xb = 0;
		yb = 0;
		zb = L / 2;
	}

	else if (is_finite(tanAlpha)) // otherwise, if the cylinder is not horizontal, use the formulas for an inclined cylinder are used
	{
		xb = tanAlpha * pow(D/2, 2) / (4 * L);
		yb = 0;
		zb = ( pow(tanAlpha*D/2, 2) + 4 * pow(L, 2) ) / (8 * L);
	}

	else // if the cylinder is horizontal and not completely submerged, forces and moments are equal to zero
	{
		return force;
	}


	// Matrix for transforming the coordinates from the local coordinate system
	// to the global coordinate system
	mat::fixed<3, 3> M_local2global = join_rows(join_rows(xvec, yvec), zvec);

	// Vector Xb - Xnode1(i.e., vector between the center of buoyancy and the
	// first node of the cylinder) written in the global coordinate frame
	arma::vec::fixed<3> Xb_global = M_local2global * arma::vec::fixed<3> { xb, yb, zb };

	// Calculation of hydrostatic force and moment
	force.at(3) = rho * g * Vol; // Fx = Fy = 0 and Fz = Buoyancy force
	force.rows(4, 6) = cross( Xb_global, force.rows(1,3) );

	return force;
}

vec::fixed<6> MorisonCirc::hydrodynamicForce(const ENVIR &envir, const vec::fixed<6> &floaterPos, const vec::fixed<6> &floaterVel, const vec::fixed<6> &floaterAcc) const
{
	// Forces and moments acting at the Morison Element
	vec::fixed<6> force(fill::zeros);

    // Use a more friendly notation    
	double D = m_diam;
	double Cd = m_CD;
	double Cm = m_CM;
	double ncyl = m_numIntPoints;	
	double rho = envir.watDensity();	

	// Nodes position, velocity and acceleration
	vec::fixed<3> n1 = node1Pos(floaterPos);
	vec::fixed<3> n2 = node2Pos(floaterPos);
	vec::fixed<3> v1 = node1Vel(floaterPos, floaterVel);
	vec::fixed<3> v2 = node2Vel(floaterPos, floaterVel);
	vec::fixed<3> a1 = node1Acc(floaterPos, floaterVel, floaterAcc);
	vec::fixed<3> a2 = node2Acc(floaterPos, floaterVel, floaterAcc);
	
	double dL = norm(n2 - n1, 2) / (ncyl - 1); // length of each interval between points

	// Vectors of the local coordinate system vectors
	vec::fixed<3> xvec(fill::zeros);
	vec::fixed<3> yvec(fill::zeros);	
    vec::fixed<3> zvec  = (n2 - n1) / arma::norm(n2 - n1, 2);
	MorisonCirc::make_local_base(xvec, yvec, zvec);

	// Matrix for transforming the coordinates from the local coordinate system
    // to the global coordinate system
    mat::fixed<3,3> M_local2global = join_rows(join_rows(xvec, yvec), zvec);
	

	//  Loop to calculate the force/moment at each integration point
	vec::fixed<3> n_ii(fill::zeros); // Coordinates of the integration point
	vec::fixed<3> vel_ii(fill::zeros); // Velocity of the integration point
	vec::fixed<3> acc_ii(fill::zeros); // Acceleration of the integration point
	vec::fixed<3> force_ii(fill::zeros); // Force acting at the integration point
	vec::fixed<3> moment_ii(fill::zeros); // Moment (with relation to node 1) due to the force acting at the integration point

	vec::fixed<3> velFluid(fill::zeros); // Fluid velocity at the integration point
	vec::fixed<3> accFluid(fill::zeros); // Fluid acceleration at the integration point
    for (int ii = 1; ii <= m_numIntPoints; ++ii) 
	{
        n_ii = (n2 - n1) * (ii-1)/(ncyl-1) + n1; // Coordinates of the integration point
        if (n_ii.at(3) > 0)        
		{
            continue;
		}

		/*******
			Fluid velocity/acceleration
		******/
		// Fluid velocity and acceleration at the integration point
		velFluid = envir.fluidVel(n_ii.at(0), n_ii.at(1), n_ii.at(2));
		accFluid = envir.fluidAcc(n_ii.at(0), n_ii.at(1), n_ii.at(2));

		// Fluid velocity/acceleration - written in the LOCAL reference frame
		// Component of the fluid velocity and acceleration at the integration point that is perpendicular to the axis of the cylinder - written in the LOCAL reference frame
		velFluid = { dot(velFluid, xvec), dot(velFluid, yvec), 0 };
		accFluid = { dot(accFluid, xvec), dot(accFluid, yvec), 0 };

		// Fluid velocity/acceleration - written in the GLOBALreference frame
		// Component of the velocity and acceleration of the integration point that is perpendicular to the axis of the cylinder - written in the GLOBAL reference frame
		velFluid = M_local2global * velFluid;
		accFluid = M_local2global * accFluid;


		/******
			Body velocity/acceleration
		******/
		// Absolute (R_ii) and relative (lambda) distance between the integration point and the first node                               
        double R_ii = norm(n_ii - n1, 2);
        double lambda = R_ii/norm(n2-n1, 2);

        // Velocity and acceleration of the integration point
        vel_ii = v1 + lambda * ( v2 - v1 );
        acc_ii = a1 + lambda * ( a2 - a1 );		

		// Component of the velocity and acceleration of the integration point that is perpendicular to the axis of the cylinder
		vel_ii = { dot(vel_ii, xvec), dot(vel_ii, yvec), 0 };
		acc_ii = { dot(acc_ii, xvec), dot(acc_ii, yvec), 0 };

		// Component of the velocity and acceleration of the integration point that is perpendicular to the axis of the cylinder
		vel_ii = M_local2global * vel_ii;
		acc_ii = M_local2global * acc_ii;

		// Calculation of the forces in the integration node using Morison's equation
		force_ii =  0.5 * rho * Cd * D * norm(velFluid - vel_ii, 2) * (velFluid - vel_ii)
				   + (datum::pi * pow(D,2)/4) * rho * Cm * accFluid
			       - (datum::pi * pow(D,2)/4) * rho * (Cm-1) * acc_ii;
	
		// Calculation of the moment (with respect to node 1) at the integration node using the force calculated above
		moment_ii = cross(R_ii * zvec, force_ii);


		if (ii == 1 || ii == ncyl)
		{
			force += (dL/3) * join_cols(force_ii, moment_ii);
		}		
		else if (ii % 2 == 0)
		{
			force += (4*dL/3) * join_cols(force_ii, moment_ii);
		}
		else
		{
			force += (2 * dL / 3) * join_cols(force_ii, moment_ii);
		}
	}
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