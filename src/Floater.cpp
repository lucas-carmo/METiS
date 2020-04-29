#include "Floater.h"
#include "MorisonCirc.h"
#include "MorisonRect.h"
#include "IO.h"

#include <iostream>
#include <vector>
#include <algorithm>    // std::binary_search
#include <utility> // For std::move
#include <stdexcept> // For std::exception

using namespace arma;

/*****************************************************
	Constructors
*****************************************************/
Floater::Floater()
{
	// Assign nan to these variable in order to check later whether they were read from the input file or not
	m_mass = datum::nan;
	m_CoG.fill(datum::nan);	
	m_inertia.fill(datum::nan);
}


/*****************************************************
	Overloaded operators
*****************************************************/
Floater& Floater::operator= (const Floater &floater)
{
	// Shallow copy for simple member variables
	m_CoG = floater.m_CoG;
	m_inertia = floater.m_inertia;
	m_mass = floater.m_mass;

	// The member variables that need deep copying are:
	// - m_MorisonElements;

	// Check for self-assignment
    if (this == &floater)
		return *this;

	// In case there is data in m_MorisonElements
	m_MorisonElements.clear();

	// Resize m_MorisonElement to match the size of the one in the input floater
	m_MorisonElements.resize( floater.m_MorisonElements.size() );

	// m_MorisonElements.at(ii) is a std::unique_ptr<MorisonElements>
	for (int ii = 0; ii < floater.m_MorisonElements.size(); ++ii)
	{
		MorisonElement* rawPtr = floater.m_MorisonElements.at(ii)->clone();
		std::unique_ptr<MorisonElement> smartPtr(rawPtr);
		m_MorisonElements.at(ii) = std::move(smartPtr);
	}
		
    return *this;	
}


/*****************************************************
	Setters
*****************************************************/
void Floater::setMass(const double mass)
{
	m_mass = mass;
}

void Floater::setInertia(const vec::fixed<6> &inertia)
{
	m_inertia = inertia;
}

void Floater::setCoG(const vec::fixed<3> &cog)
{
	m_CoG = cog;
}


/*
	Add circular cylinder Morison Element to m_MorisonElements
*/
void Floater::addMorisonCirc(vec::fixed<3> &node1_coord, vec::fixed<3> &node2_coord, const double diam, const double CD, const double CM, const unsigned int numIntPoints,
	const double botDiam, const double topDiam, const double axialCD, const double axialCa, const bool botPressFlag)
{
	// Since many times we need node1 to be below node2, it is better to swap them here in order to have less swaps in the future
	if (node1_coord[2] > node2_coord[2])
	{
		node1_coord.swap(node2_coord);
	}

	// Morison Elements need the position of the floater CoG
	// to define their relative position in the rigid body
	if (CoG().has_nan()) // If this is true, then the CoG of the floater wasn't read yet
	{
		throw std::runtime_error("Floater CoG must be specified before Morison Elements.");
	}

	// Create a circular cylinder Morison Element using the following constructor and add it to m_MorisonElements.
	m_MorisonElements.push_back(std::make_unique<MorisonCirc>(node1_coord, node2_coord, CoG(), numIntPoints,
		botPressFlag, axialCD, axialCa, diam, CD, CM, botDiam, topDiam));
}

/*
	Add rectangular cylinder Morison Element to m_MorisonElements
*/
void Floater::addMorisonRect(vec::fixed<3> &node1_coord, vec::fixed<3> &node2_coord, vec::fixed<3> &node3_coord, const double diam_X, const double diam_Y, 
	const double CD_X, const double CD_Y, const double CM_X, const double CM_Y, const unsigned int numIntPoints,
	const double botArea, const double topArea, const double axialCD, const double axialCa, const bool botPressFlag)
{
	// Morison Elements need the position of the floater CoG
	// to define their relative position in the rigid body
	if (CoG().has_nan()) // If this is true, then the CoG of the floater wasn't read yet
	{
		throw std::runtime_error("Floater CoG must be specified before Morison Elements.");
	}

	// Create a rectangular cylinder Morison Element using the following constructor and add it to m_MorisonElements.
	m_MorisonElements.push_back(std::make_unique<MorisonRect>(node1_coord, node2_coord, node3_coord, CoG(), numIntPoints,
		botPressFlag, axialCD, axialCa, diam_X, CD_X, CM_X, diam_Y, CD_Y, CM_Y, botArea, topArea));
}


/*****************************************************
	Getters
*****************************************************/
vec::fixed<3> Floater::CoG() const
{
	return m_CoG;
}

double Floater::mass() const
{
	return m_mass;
}

mat::fixed<6,6> Floater::inertiaMatrix() const
{
	mat::fixed<6, 6> A(fill::zeros);
	A(0, 0) = A(1, 1) = A(2, 2) = m_mass;
	A(3, 3) = m_inertia(0);
	A(4, 4) = m_inertia(1);
	A(5, 5) = m_inertia(2);

	A(3, 4) = A(4, 3) = m_inertia(3);
	A(3, 5) = A(5, 3) = m_inertia(4);
	A(4, 5) = A(5, 4) = m_inertia(5);

	return A;
}


std::string Floater::printMass() const
{
	return std::to_string(m_mass);
}

std::string Floater::printInertia() const
{
	return "Ixx = " + std::to_string(m_inertia(0)) + " Iyy = " + std::to_string(m_inertia(1)) + " Izz = " + std::to_string(m_inertia(2)) +
		  " Ixy = " + std::to_string(m_inertia(3)) + " Ixz = " + std::to_string(m_inertia(4)) + " Iyy = " + std::to_string(m_inertia(5));
}

std::string Floater::printCoG() const
{
	std::string output = "(" + std::to_string( m_CoG(0) );
	for (int ii = 1; ii < m_CoG.n_elem; ++ii)
	{
		output = output + "," + std::to_string( m_CoG(ii) );
	}
	
	return output + ")";
}


std::string Floater::printMorisonElements() const
{
	std::string output = "";	
	output = output + "Number of Morison Elements: " + std::to_string(m_MorisonElements.size()) + "\n";
	for (int ii = 0; ii < m_MorisonElements.size(); ++ii)
	{
		output = output + "Morison Element #" + std::to_string(ii) + '\n';
		output = output + m_MorisonElements.at(ii)->print() + '\n';
	}
	return output;
}



/*****************************************************
	Forces, acceleration, displacement, etc
*****************************************************/
void Floater::update(const vec::fixed<6> &FOWTdisp, const vec::fixed<6> &FOWTvel, const vec::fixed<6> &FOWTacc, const vec::fixed<6> &FOWTdisp_SD)
{
	m_disp = FOWTdisp;

	for (int ii = 0; ii < m_MorisonElements.size(); ++ii)
	{
		// The CoG position is added to the displacements because updateNodePosVelAcc requires the instantaneous position of the CoG of the floater
		m_MorisonElements.at(ii)->updateNodesPosVelAcc(FOWTdisp + join_cols(CoG(), zeros(3, 1)), FOWTvel, FOWTacc, FOWTdisp_SD + join_cols(CoG(), zeros(3, 1)));
	}
}

mat::fixed<6, 6> Floater::addedMass(const double density) const
{
	mat::fixed<6, 6> A(fill::zeros);

	for (int ii = 0; ii < m_MorisonElements.size(); ++ii)
	{
		A += m_MorisonElements.at(ii)->addedMass_perp(density) + m_MorisonElements.at(ii)->addedMass_paral(density);
	}

	return A;
}

vec::fixed<6> Floater::hydrodynamicForce(const ENVIR &envir, const int hydroMode) const
{	
	vec::fixed<6> force(fill::zeros); // Total hydrodynamic force acting on the floater
	vec::fixed<6> df(fill::zeros); // Total hydrodynamic force acting on each cylinder

	// These forces below are used to output the different components to the output files and internally in MorisonElement.hydrodynamicForce().
	vec::fixed<6> force_inertia(fill::zeros); // Total force acting on the floater
	vec::fixed<6> force_drag(fill::zeros);
	vec::fixed<6> force_froudeKrylov(fill::zeros);
	vec::fixed<6> force_inertia_2nd_part1(fill::zeros);
	vec::fixed<6> force_inertia_2nd_part2(fill::zeros);
	vec::fixed<6> force_inertia_2nd_part3(fill::zeros);

	vec::fixed<6> df_inertia(fill::zeros); // Forces acting on each Morison Element. They are set to zero inside MorisonElement.hydrodynamicForce().
	vec::fixed<6> df_drag(fill::zeros);
	vec::fixed<6> df_froudeKrylov(fill::zeros);
	vec::fixed<6> df_inertia_2nd_part1(fill::zeros);
	vec::fixed<6> df_inertia_2nd_part2(fill::zeros);
	vec::fixed<6> df_inertia_2nd_part3(fill::zeros);


	for (int ii = 0; ii < m_MorisonElements.size(); ++ii)
	{
		df = m_MorisonElements.at(ii)->hydrodynamicForce(envir, hydroMode, df_inertia, df_drag, df_froudeKrylov, df_inertia_2nd_part1, df_inertia_2nd_part2, df_inertia_2nd_part3);
		
		// The moments acting on the cylinders were calculated with respect to the first node
		// We need to change the fulcrum to the CoG
		df.rows(3,5) += cross( m_MorisonElements.at(ii)->node1Pos() - (m_disp.rows(0,2) + CoG()), df.rows(0,2) );		

		df_inertia.rows(3,5) += cross(m_MorisonElements.at(ii)->node1Pos() - (m_disp.rows(0, 2) + CoG()), df_inertia.rows(0, 2));
		df_drag.rows(3, 5) += cross(m_MorisonElements.at(ii)->node1Pos() - (m_disp.rows(0, 2) + CoG()), df_drag.rows(0, 2));
		df_froudeKrylov.rows(3, 5) += cross(m_MorisonElements.at(ii)->node1Pos() - (m_disp.rows(0, 2) + CoG()), df_froudeKrylov.rows(0, 2));
		df_inertia_2nd_part1.rows(3, 5) += cross(m_MorisonElements.at(ii)->node1Pos() - (m_disp.rows(0, 2) + CoG()), df_inertia_2nd_part1.rows(0, 2));
		df_inertia_2nd_part2.rows(3, 5) += cross(m_MorisonElements.at(ii)->node1Pos() - (m_disp.rows(0, 2) + CoG()), df_inertia_2nd_part2.rows(0, 2));
		df_inertia_2nd_part3.rows(3, 5) += cross(m_MorisonElements.at(ii)->node1Pos() - (m_disp.rows(0, 2) + CoG()), df_inertia_2nd_part3.rows(0, 2));

		// Add to the forces acting on the whole floater
		force += df;		
		force_inertia += df_inertia;
		force_drag += df_drag;
		force_froudeKrylov += df_froudeKrylov;
		force_inertia_2nd_part1 += df_inertia_2nd_part1;
		force_inertia_2nd_part2 += df_inertia_2nd_part2;
		force_inertia_2nd_part3 += df_inertia_2nd_part3;
	}
	
	IO::print2outLine(IO::OUTFLAG_HD_FORCE, force);
	IO::print2outLine(IO::OUTFLAG_HD_INERTIA_FORCE, force_inertia);
	IO::print2outLine(IO::OUTFLAG_HD_DRAG_FORCE, force_drag);
	IO::print2outLine(IO::OUTFLAG_HD_FK_FORCE, force_froudeKrylov);
	IO::print2outLine(IO::OUTFLAG_HD_2ND_FORCE_PART1, force_inertia_2nd_part1);
	IO::print2outLine(IO::OUTFLAG_HD_2ND_FORCE_PART2, force_inertia_2nd_part2);
	IO::print2outLine(IO::OUTFLAG_HD_2ND_FORCE_PART3, force_inertia_2nd_part3);
	
	return force;
}

vec::fixed<6> Floater::hydrostaticForce(const ENVIR &envir) const
{
	vec::fixed<6> force(fill::zeros);		
	vec::fixed<6> df(fill::zeros);
	
	for (int ii = 0; ii < m_MorisonElements.size(); ++ii)
	{
		df = m_MorisonElements.at(ii)->hydrostaticForce(envir.watDensity(), envir.gravity());

		// The moments acting on the cylinders were calculated with respect to the first node
		// We need to change the fulcrum to the CoG
		df.rows(3,5) = df.rows(3,5) + cross( m_MorisonElements.at(ii)->node1Pos() - (m_disp.rows(0, 2) + CoG()), df.rows(0,2) );

		force += df;
	}

	IO::print2outLine(IO::OUTFLAG_HS_FORCE, force);

	return force;
}


