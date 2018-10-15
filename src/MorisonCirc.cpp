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

vec::fixed<6> MorisonCirc::hydrostaticForce(const ENVIR &envir, const vec::fixed<6> &floaterPos) const
{
	vec::fixed<6> force(fill::zeros); 	
	return force;
}

vec::fixed<6> MorisonCirc::hydrodynamicForce(const ENVIR &envir, const vec::fixed<6> &floaterPos, const vec::fixed<6> &floaterVel) const
{
	vec::fixed<6> force(fill::zeros);

    // Use a more friendly notation    
	vec::fixed<3> n1 = node1Pos(floaterPos);
	vec::fixed<3> n2 = node2Pos(floaterPos);
	double D = m_diam;
	double Cd = m_CD;
	double Cm = m_CM;
	double ncyl = m_numIntPoints;	
	double rho = envir.watDensity();

	// vec::fixed<3> n1_vel = 

	// Vectors of the local coordinate system vectors
	vec::fixed<3> xvec(fill::zeros);
	vec::fixed<3> yvec(fill::zeros);	
    vec::fixed<3> zvec  = (n2 - n1) / arma::norm(n2 - n1, 2);

	// if Zlocal == Zglobal, REFlocal = REFglobal
	vec::fixed<3> z_global = {0,0,1};
	if ( approx_equal(zvec, z_global, "absdiff", datum::eps) )
	{		
		xvec = {1, 0, 0};
		yvec = {0, 1, 0};
	}
	else
	{
		yvec = arma::cross(z_global, zvec);
		yvec = yvec / norm(yvec, 2);

		xvec = cross(yvec, zvec);
		xvec = xvec / norm(xvec, 2);
	}

	// Matrix for transforming the coordinates from the local coordinate system
    // to the global coordinate system
    mat::fixed<3,3> M_local2global = join_rows(join_rows(xvec, yvec), zvec);

	//  Loop to calculate the force/moment at each integration point
	vec::fixed<3> nii(fill::zeros);
    for (int ii = 1; ii <= m_numIntPoints; ++ii) 
	{
        nii = (n2 - n1) * (ii-1)/(ncyl-1) + n1; // Coordinates of the integration point
        if (nii.at(3) > 0)        
		{
            continue;
		}

		// Absolute (Rii) and relative (lambda) distance between the integration point and the first node                               
        double Rii = norm(nii - n1, 2);
        double lambda = Rii/norm(n2-n1, 2);

        // Velocity and acceleration of the integration point
        // velii = elem_t.vel1 + lambda * ( elem_t.vel2 - elem_t.vel1 );
        // accii = elem_t.acc1 + lambda * ( elem_t.acc2 - elem_t.acc1 );		


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