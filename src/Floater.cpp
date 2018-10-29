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


	// However, Morison Elements are not actually specified by the nodes coordinates, 
	// but rather by the distance between the nodes and the floater CoG		
	if ( m_CoG.has_nan() ) // If this is true, then the CoG of the floater wasn't read yet
	{
		throw std::runtime_error( "Floater CoG must be specified before Morison Elements. Error in input line " + std::to_string(IO::getInLineNumber()) );
	}
	
	node1_coord = node1_coord - m_CoG;
	node2_coord = node2_coord - m_CoG;


	 // Create a circular cylinder Morison Element using the following constructor and add it to m_MorisonElements.
	 m_MorisonElements.push_back( std::make_shared<MorisonCirc>(node1_coord, node2_coord, numIntPoints, 
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


	// However, Morison Elements are not actually specified by the nodes coordinates, 
	// but rather by the distance between the nodes and the floater CoG		
	if ( m_CoG.has_nan() ) // If this is true, then the CoG of the floater wasn't read yet
	{
		throw std::runtime_error( "Floater CoG must be specified before Morison Elements. Error in input line " + std::to_string(IO::getInLineNumber()) );
	}
		
	node1_coord = node1_coord - m_CoG;
	node2_coord = node2_coord - m_CoG;
	node3_coord = node3_coord - m_CoG;

	
	 // Create a rectangular cylinder Morison Element using the following constructor and add it to m_MorisonElements.
	m_MorisonElements.push_back( std::make_shared<MorisonRect>(node1_coord, node2_coord, node3_coord, numIntPoints,
	  							  botPressFlag, axialCD, axialCa, diam_X, CD_X, CM_X, diam_Y, CD_Y, CM_Y, botArea, topArea) );
}


/*****************************************************
	Getters
*****************************************************/

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
	Overloaded operators
*****************************************************/
Floater& Floater::operator= (Floater &floater)
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

	// Resize m_MorisonElemen to match the size of the one in the input floater
	m_MorisonElements.resize( floater.m_MorisonElements.size() );

	for (int ii = 0; ii < floater.m_MorisonElements.size(); ++ii)
	{		
		m_MorisonElements.at(ii) = floater.m_MorisonElements.at(ii);
	}		
		
    return *this;	
}



/*****************************************************
	Forces, acceleration, position, etc
*****************************************************/
vec::fixed<6> Floater::hydrodynamicForce(const ENVIR &envir, const vec &FOWTpos, const vec &FOWTvel, const vec &FOWTacc) const
{
	vec::fixed<6> force(fill::zeros);
	for (int ii = 0; ii < m_MorisonElements.size(); ++ii)
	{
		force += m_MorisonElements.at(ii)->hydrodynamicForce(envir, FOWTpos, FOWTvel, FOWTacc);
	}
	return force;
}

//mat Floater::rotatMat(const vec &FOWTpos)
//{}
//
//


