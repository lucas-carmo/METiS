#include "Floater.h"
#include "MorisonCirc.h"
#include "MorisonRect.h"

#include "IO.h"
#include <iostream>
#include <vector>
#include <algorithm>    // std::binary_search

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



void Floater::readCoG(const std::string &data)
{
	// The coordinates of the center of gravity are separated by commas in the input string
	std::vector<std::string> input = stringTokenize(data, ",");

	if (input.size() != 3)
	{
		std::cout << "Deu ruim \n";
	}

	for (int ii = 0; ii < input.size(); ++ii)
	{
		readDataFromString(input.at(ii), m_CoG(ii));
	}
}



void Floater::addNode(const std::string &data)
{
	// Nodes are specified by a vec with four components: ID, X coord, Y coord, and Z coord. 
	// They are separated by commas in the input string.
	std::vector<std::string> input = stringTokenize(data, ",");

	if (input.size() != 4)
	{
		std::cout << "Deu ruim na leitura dos nos. Linha " << IO::getInLineNumber() << ".\n";
	}

	// Read node ID
	unsigned int nodeID{0};
	readDataFromString( input.at(0), nodeID );

	if (m_nodesID.size() != 0) // If this is not the first node that will be added to m_nodesID
	{
		if (nodeID <= m_nodesID.back()) // Then verify if its ID is larger than the previous one, thus garanteeing that m_nodesID is in ascending order (this is needed to use binary search to find nodes IDs)
		{
			std::cout << "Nodes tem que estar organizados em ordem crescente. Erro na linha " << IO::getInLineNumber() << ".\n";
		}
	}

	m_nodesID.push_back( nodeID );


	// Read node coord
	vec::fixed<3> nodeCoord(fill::zeros);
	for (int ii = 0; ii < nodeCoord.n_elem; ++ii)
	{
		readDataFromString( input.at(ii+1), nodeCoord(ii) );
	}

	m_nodesCoord.push_back( nodeCoord );
}



void Floater::addMorisonCirc(const std::string &data)
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
		std::cout << "Deu ruim na leitura do cilindro circular da linha " << IO::getInLineNumber() << ".\n";
		return;
	}

	// Check whether nodes were specified
	if (m_nodesID.empty())
	{
		std::cout << "Nodes should be specified before Morison's Elements. Error in line " << IO::getInLineNumber() << '\n';
		return;
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
	
	// Get coordinates of the FIRST node based on its ID
	std::vector<unsigned int>::iterator iter1 = std::find(m_nodesID.begin(), m_nodesID.end(), node1_ID); // Find first node by its ID.
	vec::fixed<3> node1_coord(fill::zeros);
	if (iter1 != m_nodesID.end())
	{
		auto index = std::distance(m_nodesID.begin(), iter1); // Get index by the distance between the iterator found above and m_nodes.begin()
		node1_coord = m_nodesCoord.at(index);		
	}
	else
	{
		std::cout << "Unable to find the first node of the Morison Element specified in line " << IO::getInLineNumber() << '\n';
	}

	// Get coordinates of the SECOND node based on its ID
	std::vector<unsigned int>::iterator iter2 = std::find(m_nodesID.begin(), m_nodesID.end(), node2_ID); // Find first node by its ID.
	vec::fixed<3> node2_coord(fill::zeros);
	if (iter2 != m_nodesID.end())
	{
		auto index = std::distance(m_nodesID.begin(), iter2); // Get index by the distance between the iterator found above and m_nodes.begin()
		node2_coord = m_nodesCoord.at(index);
	}
	else
	{
		std::cout << "Unable to find the second node of the Morison Element specified in line " << IO::getInLineNumber() << '\n';
	}


	// However, Morison Elements are not actually specified by the nodes coordinates, but rather by the distance between the nodes and the floater CoG

	node1_coord = node1_coord - m_CoG;
	node2_coord = node2_coord - m_CoG;

	
	// Create a circular cylinder Morison Element using the following constructor and add it to m_MorisonElements.
	MorisonCirc morison_element(node1_coord, node2_coord, numIntPoints, botPressFlag, diam, CD, CM, botDiam, topDiam, axialCD, axialCa); // As we are using std::reference_wrapper for m_MorisonElements, we can't use anonymous objects
	m_MorisonElements.push_back( morison_element );
}



/*****************************************************
	Getters
*****************************************************/

std::string Floater::printMass() const
{
	return std::to_string(m_mass);
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



std::string Floater::printNodes() const
{
	std::string output = "";
	for (int ii = 0; ii < m_nodesID.size(); ++ii)
	{
		output = output + "( " + std::to_string( m_nodesID.at(ii) ) + 
						  ", " + std::to_string( m_nodesCoord.at(ii)(0) ) +
			              ", " + std::to_string( m_nodesCoord.at(ii)(1) ) +
			              ", " + std::to_string( m_nodesCoord.at(ii)(2) ) + " )\n";
	}
	return output;
}



std::string Floater::printMorisonElements() const
{
	std::string output = "";
	int a = m_MorisonElements.size();
	int b = 1;
	for (int ii = 0; ii < m_MorisonElements.size(); ++ii)
	{
		std::cout << "Morison Element #" << ii << '\n';
		m_MorisonElements.at(ii).get().print();
	}
	return output;
}

//mat Floater::rotatMat(const vec &FOWTpos)
//{}
//
//
//vec Floater::hydrodynamicForce(const ENVIR &envir, const vec &FOWTpos, const vec &FOWTvel, const vec &FOWTacc)
//{}


