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
	m_mass = floater.m_mass;
	m_CoG = floater.m_CoG;
	m_inertia = floater.m_inertia;	
	m_addedMass_t0 = floater.m_addedMass_t0;
	m_disp = floater.m_disp;
	m_disp_sd = floater.m_disp_sd;

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

void Floater::setInstantSD(bool instantSD)
{
	m_instantSD = instantSD;
}

void Floater::evaluateQuantitiesAtBegin(const ENVIR &envir, const int hydroMode)
{
	for (int ii = 0; ii < m_MorisonElements.size(); ++ii)
	{
		m_MorisonElements.at(ii)->evaluateQuantitiesAtBegin(envir, hydroMode);
	}
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

vec::fixed<6> Floater::CoGPos() const
{
	return join_cols(m_CoG, zeros(3,1)) + m_disp;
}

vec::fixed<6> Floater::CoGPos_sd() const
{
	return join_cols(m_CoG, zeros(3, 1)) + m_disp_sd;
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
		m_MorisonElements.at(ii)->updateMorisonElement(envir, CoGPos(), FOWTvel, CoGPos_sd(), FOWTvel_SD);
	}
}

// Needed this function to evaluate the added mass at the beginning of the simulation
void Floater::setAddedMass_t0(const double density)
{
	m_addedMass_t0 = density*addedMass(1);
}

mat::fixed<6, 6> Floater::addedMass(const double density, const int hydroMode) 
{
	mat::fixed<6, 6> A(fill::zeros);

	A = density*addedMass(hydroMode);
	return A;
}

mat::fixed<6, 6> Floater::addedMass(const int hydroMode) const
{
	mat::fixed<6, 6> A(fill::zeros);
	vec::fixed<3> refPt(fill::zeros);

	for (int ii = 0; ii < m_MorisonElements.size(); ++ii)
	{
		// Evaluated at the fixed initial position
		if (hydroMode < 2)
		{
			refPt = CoGPos_sd().rows(0, 2);
		}
		else
		{
			refPt = CoGPos().rows(0,2);
		}

		A += m_MorisonElements.at(ii)->addedMass_perp(1, refPt, hydroMode) + m_MorisonElements.at(ii)->addedMass_paral(1, refPt, hydroMode);
	}
	return A;
}

// Moments are output with respect to the CoG of the floater
vec::fixed<6> Floater::hydrodynamicForce(const ENVIR &envir, const int hydroMode) const
{	
	vec::fixed<6> force(fill::zeros); // Total hydrodynamic force acting on the floater
	vec::fixed<6> df(fill::zeros); // Total hydrodynamic force acting on each cylinder

	// These forces below are used to output the different components to the output files and internally in MorisonElement.hydrodynamicForce().	
	vec::fixed<6> force_drag(fill::zeros); // Drag force
	vec::fixed<6> force_1stP(fill::zeros); // Due to first-order potential
	vec::fixed<6> force_eta(fill::zeros);  // Due to relative wave elevation
	vec::fixed<6> force_conv(fill::zeros); // Due to convective acceleration
	vec::fixed<6> force_axdv(fill::zeros); // Due to axial-divergence acceleration
	vec::fixed<6> force_acgr(fill::zeros); // Due to gradient of fluid acceleration
	vec::fixed<6> force_rotn(fill::zeros); // Due to the rotation of the normal vector
	vec::fixed<6> force_2ndP(fill::zeros); // Component due to 2nd order incoming flow
	vec::fixed<6> force_rslb(fill::zeros); // Due to the rotation term from slender-body approximation
	vec::fixed<6> force_rem(fill::zeros);  // Remaining force components

	vec::fixed<3> refPt = CoGPos_sd().rows(0, 2);
	mat::fixed<6, 6> R(fill::eye);
	R *= -1;
	R.rows(0, 2).cols(0, 2) += rotatMatrix(m_disp.rows(3, 5));
	R.rows(3, 5).cols(3, 5) = R.rows(0, 2).cols(0, 2);	

	// Force acting on each element
	for (int ii = 0; ii < m_MorisonElements.size(); ++ii)	
	{
		// If the cylinder is above the waterline, do not need to compute the hydrodynamic forces
		if (m_MorisonElements.at(ii)->node1Pos_sd().at(2) > 0)
			continue;		

		// Force due to first-order acceleration - part along the length of the cylinder is integrated analytically
		vec::fixed<6> auxForce = m_MorisonElements.at(ii)->hydroForce_1st(envir, hydroMode, refPt);
		force_1stP += auxForce;

		// Second order forces
		if (hydroMode == 2)
		{
			if (envir.waveStret() == 1)
			{
				force_eta += m_MorisonElements.at(ii)->hydroForce_relWaveElev(envir, refPt);
			}

			// If the body is not treated as fixed, these second order components are included in hydroForce_1st
			if (m_MorisonElements.at(ii)->flagFixed())
			{				
				force_rotn += R * auxForce;
				force_acgr += m_MorisonElements.at(ii)->hydroForce_accGradient(envir, refPt);
			}
			
			force_conv += m_MorisonElements.at(ii)->hydroForce_convecAcc(envir, refPt);
			force_axdv += m_MorisonElements.at(ii)->hydroForce_axDiverg(envir, refPt);
			force_2ndP+= m_MorisonElements.at(ii)->hydroForce_2ndPot(envir, refPt);
			force_rslb += m_MorisonElements.at(ii)->hydroForce_slendBodyRot(envir, refPt);						
		}
		force_drag += m_MorisonElements.at(ii)->hydroForce_drag(envir, refPt);
	}	

	force = force_drag + force_1stP + force_eta + force_conv + force_axdv + force_acgr + force_rotn + force_2ndP + force_rslb + force_rem;

	IO::print2outLine(IO::OUTFLAG_HD_FORCE, force);
	IO::print2outLine(IO::OUTFLAG_HD_FORCE_DRAG, force_drag);
	IO::print2outLine(IO::OUTFLAG_HD_FORCE_1STP, force_1stP);
	IO::print2outLine(IO::OUTFLAG_HD_FORCE_ETA, force_eta);
	IO::print2outLine(IO::OUTFLAG_HD_FORCE_CONV, force_conv);
	IO::print2outLine(IO::OUTFLAG_HD_FORCE_AXDV, force_axdv);
	IO::print2outLine(IO::OUTFLAG_HD_FORCE_ACGR, force_acgr);
	IO::print2outLine(IO::OUTFLAG_HD_FORCE_ROTN, force_rotn);
	IO::print2outLine(IO::OUTFLAG_HD_FORCE_2NDP, force_2ndP);
	IO::print2outLine(IO::OUTFLAG_HD_FORCE_RSLB, force_rslb);
	IO::print2outLine(IO::OUTFLAG_HD_FORCE_REM, force_rem);
	
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
		df.rows(3, 5) = df.rows(3, 5) + cross(m_MorisonElements.at(ii)->node1Pos() - CoGPos().rows(0, 2), df.rows(0, 2));
		force += df;
	}
	IO::print2outLine(IO::OUTFLAG_HS_FORCE, force);		

	return force;
}


