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
	m_addedMass_t0.fill(datum::nan);
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
	const double axialCD_1, const double axialCa_1, const double axialCD_2, const double axialCa_2, const bool botPressFlag)
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
		botPressFlag, axialCD_1, axialCa_1, axialCD_2, axialCa_2, diam, CD, CM));
}

/*
	Add rectangular cylinder Morison Element to m_MorisonElements
*/
void Floater::addMorisonRect(vec::fixed<3> &node1_coord, vec::fixed<3> &node2_coord, vec::fixed<3> &node3_coord, const double diam_X, const double diam_Y, 
	const double CD_X, const double CD_Y, const double CM_X, const double CM_Y, const unsigned int numIntPoints,
	const double axialCD_1, const double axialCa_1, const double axialCD_2, const double axialCa_2, const bool botPressFlag)
{
	// Morison Elements need the position of the floater CoG
	// to define their relative position in the rigid body
	if (CoG().has_nan()) // If this is true, then the CoG of the floater wasn't read yet
	{
		throw std::runtime_error("Floater CoG must be specified before Morison Elements.");
	}

	// Create a rectangular cylinder Morison Element using the following constructor and add it to m_MorisonElements.
	m_MorisonElements.push_back(std::make_unique<MorisonRect>(node1_coord, node2_coord, node3_coord, CoG(), numIntPoints,
		botPressFlag, axialCD_1, axialCa_1, axialCD_2, axialCa_2, diam_X, CD_X, CM_X, diam_Y, CD_Y, CM_Y));
}

mat::fixed<6, 6> Floater::addedMass_t0() const
{
	if (m_addedMass_t0.is_finite())
	{
		return m_addedMass_t0;
	}
	else
	{
		throw std::runtime_error("Tried to call Floater::addedMass_t0(), but it was not calculated yet. Try calling Floater::addedMass(const double density, const int hydroMode) first.");
	}
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
void Floater::update(const ENVIR &envir, const vec::fixed<6> &FOWTdisp, const vec::fixed<6> &FOWTvel, const vec::fixed<6> &FOWTdisp_SD, const vec::fixed<6> &FOWTvel_SD)
{
	m_disp = FOWTdisp;
	m_disp_sd = FOWTdisp_SD;

	for (int ii = 0; ii < m_MorisonElements.size(); ++ii)
	{
		// The CoG position is added to the displacements because updateMorisonElement requires the instantaneous position of the CoG of the floater
		m_MorisonElements.at(ii)->updateMorisonElement(envir, FOWTdisp + join_cols(CoG(), zeros(3, 1)), FOWTvel, FOWTdisp_SD + join_cols(CoG(), zeros(3, 1)), FOWTvel_SD);
	}
}

mat::fixed<6, 6> Floater::addedMass(const double density, const int hydroMode) 
{
	mat::fixed<6, 6> A(fill::zeros);

	A = density*addedMass(hydroMode);

	if (!m_addedMass_t0.is_finite())
	{
		m_addedMass_t0 = A;
	}

	IO::print2outLine(IO::OUTFLAG_ADDED_MASS_DIAG, A.diag());
	return A;
}

mat::fixed<6, 6> Floater::addedMass(const int hydroMode) const
{
	mat::fixed<6, 6> A(fill::zeros);

	for (int ii = 0; ii < m_MorisonElements.size(); ++ii)
	{
		if (hydroMode == 1)
		{
			A += m_MorisonElements.at(ii)->addedMass_perp(1, m_disp_sd.rows(0, 2) + CoG(), hydroMode) + m_MorisonElements.at(ii)->addedMass_paral(1, m_disp_sd.rows(0, 2) + CoG(), hydroMode);
		}
		else
		{
			A += m_MorisonElements.at(ii)->addedMass_perp(1, m_disp.rows(0, 2) + CoG(), hydroMode) + m_MorisonElements.at(ii)->addedMass_paral(1, m_disp.rows(0, 2) + CoG(), hydroMode);
		}
	}
	IO::print2outLine(IO::OUTFLAG_ADDED_MASS_DIAG, A.diag());
	return A;
}

// Moments are output with respect to the CoG of the floater
vec::fixed<6> Floater::hydrodynamicForce(const ENVIR &envir, const int hydroMode) const
{	
	vec::fixed<6> force(fill::zeros); // Total hydrodynamic force acting on the floater
	vec::fixed<6> df(fill::zeros); // Total hydrodynamic force acting on each cylinder

	// These forces below are used to output the different components to the output files and internally in MorisonElement.hydrodynamicForce().	
	vec::fixed<6> force_drag(fill::zeros); // Drag force
	vec::fixed<6> force_1(fill::zeros);	// Component due to 1st order incoming flow
	vec::fixed<6> force_2(fill::zeros); // Component due to 2nd order incoming flow
	vec::fixed<6> force_3(fill::zeros); // Component due to axial acceleration
	vec::fixed<6> force_4(fill::zeros); // Component due to the axial-divergence acceleration
	vec::fixed<6> force_eta(fill::zeros); // Component due to the wave elevation
	vec::fixed<6> force_rem(fill::zeros); // Remaining force components

	// Same thing, but for the forces acting at the extremities of the cylinder
	vec::fixed<6> force_drag_ext(fill::zeros);
	vec::fixed<6> force_1_ext(fill::zeros);
	vec::fixed<6> force_2_ext(fill::zeros);
	vec::fixed<6> force_3_ext(fill::zeros);
	vec::fixed<6> force_eta_ext(fill::zeros);
	vec::fixed<6> force_rem_ext(fill::zeros);

	// Forces acting on each Morison Element. They are set to zero inside MorisonElement.hydrodynamicForce().
	vec::fixed<6> df_drag(fill::zeros);
	vec::fixed<6> df_1(fill::zeros);
	vec::fixed<6> df_2(fill::zeros);
	vec::fixed<6> df_3(fill::zeros);
	vec::fixed<6> df_4(fill::zeros);
	vec::fixed<6> df_eta(fill::zeros);
	vec::fixed<6> df_rem(fill::zeros);

	vec::fixed<6> df_drag_ext(fill::zeros);
	vec::fixed<6> df_1_ext(fill::zeros);
	vec::fixed<6> df_2_ext(fill::zeros);
	vec::fixed<6> df_3_ext(fill::zeros);
	vec::fixed<6> df_rem_ext(fill::zeros);

	for (int ii = 0; ii < m_MorisonElements.size(); ++ii)
	{
		// Force acting on each element
		if (hydroMode == 1)
		{
			df = m_MorisonElements.at(ii)->hydrodynamicForce(envir, hydroMode, (m_disp_sd.rows(0, 2) + CoG()), (m_disp_sd.rows(0, 2) + CoG()), df_drag, df_1, df_2, df_3, df_4, df_eta, df_rem, df_drag_ext, df_1_ext, df_2_ext, df_3_ext, df_rem_ext);
		}
		else
		{
			df = m_MorisonElements.at(ii)->hydrodynamicForce(envir, hydroMode, (m_disp.rows(0, 2) + CoG()), (m_disp_sd.rows(0, 2) + CoG()), df_drag, df_1, df_2, df_3, df_4, df_eta, df_rem, df_drag_ext, df_1_ext, df_2_ext, df_3_ext, df_rem_ext);
		}

		// Add to the forces acting on the whole floater
		force += df;

		force_drag += df_drag;
		force_1 += df_1;
		force_2 += df_2;
		force_3 += df_3;
		force_4 += df_4;
		force_eta += df_eta;
		force_rem += df_rem;

		force_drag_ext += df_drag_ext;
		force_1_ext += df_1_ext;
		force_2_ext += df_2_ext;
		force_3_ext += df_3_ext;
		force_rem_ext += df_rem_ext;
	}

	IO::print2outLine(IO::OUTFLAG_HD_FORCE, force);

	IO::print2outLine(IO::OUTFLAG_HD_FORCE_DRAG, force_drag);
	IO::print2outLine(IO::OUTFLAG_HD_FORCE_1, force_1);	
	IO::print2outLine(IO::OUTFLAG_HD_FORCE_2, force_2);
	IO::print2outLine(IO::OUTFLAG_HD_FORCE_3, force_3);
	IO::print2outLine(IO::OUTFLAG_HD_FORCE_4, force_4);
	IO::print2outLine(IO::OUTFLAG_HD_FORCE_ETA, force_eta);
	IO::print2outLine(IO::OUTFLAG_HD_FORCE_REM, force_rem);

	IO::print2outLine(IO::OUTFLAG_HD_FORCE_DRAG_EXT, force_drag);
	IO::print2outLine(IO::OUTFLAG_HD_FORCE_1_EXT, force_1);
	IO::print2outLine(IO::OUTFLAG_HD_FORCE_2_EXT, force_2);
	IO::print2outLine(IO::OUTFLAG_HD_FORCE_3_EXT, force_3);
	IO::print2outLine(IO::OUTFLAG_HD_FORCE_REM_EXT, force_rem);
	
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


