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
void Floater::readMass(const std::string &data)
{
	readDataFromString(data, m_mass);
}

void Floater::readInertia(const std::string &data)
{
	// The different components of the inertia matrix are separated by commas in the input string
	std::vector<std::string> input = stringTokenize(data, ",");

	if (input.size() != 6)
	{
		throw std::runtime_error("Unable to read the floater inertia matrix in input line " + std::to_string(IO::getInLineNumber()) + ". Wrong number of parameters.");
	}

	for (int ii = 0; ii < input.size(); ++ii)
	{
		readDataFromString(input.at(ii), m_inertia(ii));
	}
}

void Floater::readCoG(const std::string &data)
{
	// The coordinates of the center of gravity are separated by commas in the input string
	std::vector<std::string> input = stringTokenize(data, ",");

	if (input.size() != 3)
	{
		throw std::runtime_error("Unable to read the CoG in input line " + std::to_string(IO::getInLineNumber()) + ". Wrong number of parameters.");
	}

	for (int ii = 0; ii < input.size(); ++ii)
	{
		readDataFromString(input.at(ii), m_CoG(ii));
	}
}



/*
	Add circular cylinder Morison Element to m_MorisonElements
*/
void Floater::addMorisonCirc(const std::string &data, const ENVIR &envir)
{
	// Variables to handle the data read from the input file in a more user friendly way
	unsigned int node1_ID = 0;
	unsigned int node2_ID = 0;
	double diam = 0;
	double CD = 0;
	double CM = 0;
	int numIntPoints = 0;
	double botDiam = 0;
	double topDiam = 0;
	double axialCD = 0;
	double axialCa = 0;
	bool botPressFlag = false;

	// The eleven properties of a circular cylinder Morison's Element are separated by white spaces in the input string.
	std::vector<std::string> input = stringTokenize(data, " \t");

	// Check number of inputs
	if (input.size() != 11)
	{
		throw std::runtime_error("Unable to read the circular cylinder in input line " + std::to_string(IO::getInLineNumber()) + ". Wrong number of parameters.");
	}

	// Check whether nodes were specified
	if (envir.isNodeEmpty())
	{
		throw std::runtime_error( "Nodes should be specified before Morison Elements. Error in input line " + std::to_string(IO::getInLineNumber()) );
	}

	// Read data
	readDataFromString(input.at(0), node1_ID);
	readDataFromString(input.at(1), node2_ID);
	readDataFromString(input.at(2), diam);
	readDataFromString(input.at(3), CD);
	readDataFromString(input.at(4), CM);
	readDataFromString(input.at(5), numIntPoints);
	readDataFromString(input.at(6), botDiam);
	readDataFromString(input.at(7), topDiam);
	readDataFromString(input.at(8), axialCD);
	readDataFromString(input.at(9), axialCa);
	readDataFromString(input.at(10), botPressFlag);
	
	// Get coordinates of nodes based on their ID
	vec::fixed<3> node1_coord = envir.getNode(node1_ID);
	vec::fixed<3> node2_coord = envir.getNode(node2_ID);

	// Since many times we need node1 to be below node2, it is better to swap them here in order to have less swaps in the future
	if (node1_coord[2] > node2_coord[2])
	{
		node1_coord.swap(node2_coord);
	}

	// Morison Elements need the position of the floater CoG
	// to define their relative position in the rigid body
	if (CoG().has_nan()) // If this is true, then the CoG of the floater wasn't read yet
	{
		throw std::runtime_error("Floater CoG must be specified before Morison Elements. Error in input line " + std::to_string(IO::getInLineNumber()));
	}

	 // Create a circular cylinder Morison Element using the following constructor and add it to m_MorisonElements.
	 m_MorisonElements.push_back( std::make_unique<MorisonCirc>(node1_coord, node2_coord, CoG(), numIntPoints, 
	 							  botPressFlag, axialCD, axialCa, diam, CD, CM, botDiam, topDiam) );
}



/*
	Add rectangular cylinder Morison Element to m_MorisonElements
*/
void Floater::addMorisonRect(const std::string &data, const ENVIR &envir)
{
	// Variables to handle the data read from the input file in a more user friendly way
	unsigned int node1_ID = 0;
	unsigned int node2_ID = 0;
	unsigned int node3_ID = 0;
	double diam_X = 0;
	double diam_Y = 0;
	double CD_X = 0;
	double CD_Y = 0;
	double CM_X = 0;
	double CM_Y = 0;
	int numIntPoints = 0;
	double botArea = 0;
	double topArea = 0;
	double axialCD = 0;
	double axialCa = 0;
	bool botPressFlag = false;

	// The eleven properties of a circular cylinder Morison's Element are separated by white spaces in the input string.
	std::vector<std::string> input = stringTokenize(data, " \t");

	// Check number of inputs
	if (input.size() != 15)
	{
		throw std::runtime_error("Unable to read the rectangular cylinder in input line " + std::to_string(IO::getInLineNumber()) + ". Wrong number of parameters.");
	}

	// Check whether nodes were specified
	if (envir.isNodeEmpty())
	{
		throw std::runtime_error( "Nodes should be specified before Morison Elements. Error in input line " + std::to_string(IO::getInLineNumber()) );
	}

	// Read data
	readDataFromString(input.at(0), node1_ID);
	readDataFromString(input.at(1), node2_ID);
	readDataFromString(input.at(2), node3_ID);
	readDataFromString(input.at(3), diam_X);
	readDataFromString(input.at(4), CD_X);
	readDataFromString(input.at(5), CM_X);
	readDataFromString(input.at(6), diam_Y);
	readDataFromString(input.at(7), CD_Y);
	readDataFromString(input.at(8), CM_Y);
	readDataFromString(input.at(9), numIntPoints);
	readDataFromString(input.at(10), botArea);
	readDataFromString(input.at(11), topArea);
	readDataFromString(input.at(12), axialCD);
	readDataFromString(input.at(13), axialCa);
	readDataFromString(input.at(14), botPressFlag);
	
	// Get coordinates of nodes based on their ID
	vec::fixed<3> node1_coord = envir.getNode(node1_ID);
	vec::fixed<3> node2_coord = envir.getNode(node2_ID);
	vec::fixed<3> node3_coord = envir.getNode(node3_ID);


	// Morison Elements need the position of the floater CoG
	// to define their relative position in the rigid body
	if ( CoG().has_nan() ) // If this is true, then the CoG of the floater wasn't read yet
	{
		throw std::runtime_error( "Floater CoG must be specified before Morison Elements. Error in input line " + std::to_string(IO::getInLineNumber()) );
	}
			
	// Create a rectangular cylinder Morison Element using the following constructor and add it to m_MorisonElements.
	m_MorisonElements.push_back( std::make_unique<MorisonRect>(node1_coord, node2_coord, node3_coord, CoG(), numIntPoints,
	  							  botPressFlag, axialCD, axialCa, diam_X, CD_X, CM_X, diam_Y, CD_Y, CM_Y, botArea, topArea) );
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
void Floater::update(const vec::fixed<6> &FOWTdisp, const vec::fixed<6> &FOWTvel, const vec::fixed<6> &FOWTacc)
{
	m_disp = FOWTdisp;

	for (int ii = 0; ii < m_MorisonElements.size(); ++ii)
	{
		m_MorisonElements.at(ii)->updateNodesPosVelAcc(FOWTdisp + join_cols(CoG(), zeros(3, 1)), FOWTvel, FOWTacc);
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

	// These forces below are only used to output the different components to the output files.	
	vec::fixed<6> force_inertia(fill::zeros); // Total force acting on the floater
	vec::fixed<6> force_drag(fill::zeros);
	vec::fixed<6> force_froudeKrylov(fill::zeros);

	vec::fixed<6> df_inertia(fill::zeros); // Forces acting on each Morison Element
	vec::fixed<6> df_drag(fill::zeros);
	vec::fixed<6> df_froudeKrylov(fill::zeros);


	for (int ii = 0; ii < m_MorisonElements.size(); ++ii)
	{
		// Make sure that the force components acting on each Morison Element were set to zero
		df_inertia.zeros();
		df_drag.zeros();
		df_froudeKrylov.zeros();

		df = m_MorisonElements.at(ii)->hydrodynamicForce(envir, hydroMode, df_inertia, df_drag, df_froudeKrylov);
		
		// The moments acting on the cylinders were calculated with respect to the first node
		// We need to change the fulcrum to the CoG
		df.rows(3,5) += cross( m_MorisonElements.at(ii)->node1Pos() - (m_disp.rows(0,2) + CoG()), df.rows(0,2) );		

		df_inertia.rows(3,5) += cross(m_MorisonElements.at(ii)->node1Pos() - (m_disp.rows(0, 2) + CoG()), df_inertia.rows(0, 2));
		df_drag.rows(3, 5) += cross(m_MorisonElements.at(ii)->node1Pos() - (m_disp.rows(0, 2) + CoG()), df_drag.rows(0, 2));
		df_froudeKrylov.rows(3, 5) += cross(m_MorisonElements.at(ii)->node1Pos() - (m_disp.rows(0, 2) + CoG()), df_froudeKrylov.rows(0, 2));

		// Add to the forces acting on the whole floater
		force += df;		

		force_inertia += df_inertia;
		force_drag += df_drag;
		force_froudeKrylov += df_froudeKrylov;
	}
	
	IO::print2outLine(IO::OUTFLAG_HD_FORCE, force);
	IO::print2outLine(IO::OUTFLAG_HD_INERTIA_FORCE, force_inertia);
	IO::print2outLine(IO::OUTFLAG_HD_DRAG_FORCE, force_drag);
	IO::print2outLine(IO::OUTFLAG_HD_FK_FORCE, force_froudeKrylov);
	
	return force;
}

vec::fixed<6> Floater::hydrostaticForce(const double watDensity, const double gravity, const int hydroMode) const
{
	vec::fixed<6> force(fill::zeros);		
	vec::fixed<6> df(fill::zeros);

	// Z coordinate of cylinder intersection with water line. 
	// If hydroMode = 1, z_wl = 0, i.e. wave elevation is not considered. This works just like traditional linear hydrostatics if displacements are indeed small)
	// If hydroMode = 2, z_wl is set to the z coordinate of the intersection of the cylinder with the wave elevation.
	double z_wl{ 0 }; 

	for (int ii = 0; ii < m_MorisonElements.size(); ++ii)
	{
		//if (hydroMode == 2)
		//{
		//	// Ainda não faz nada
		//}

		df = m_MorisonElements.at(ii)->hydrostaticForce(watDensity, gravity, z_wl);

		// The moments acting on the cylinders were calculated with respect to the first node
		// We need to change the fulcrum to the CoG
		df.rows(3,5) = df.rows(3,5) + cross( m_MorisonElements.at(ii)->node1Pos() - (m_disp.rows(0, 2) + CoG()), df.rows(0,2) );

		force += df;
	}

	IO::print2outLine(IO::OUTFLAG_HS_FORCE, force);

	return force;
}


