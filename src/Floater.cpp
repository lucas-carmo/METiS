#include "Floater.h"

#include "IO.h"
#include <iostream>
#include <vector>


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


//mat Floater::rotatMat(const vec &FOWTpos)
//{}
//
//
//vec Floater::hydrodynamicForce(const ENVIR &envir, const vec &FOWTpos, const vec &FOWTvel, const vec &FOWTacc)
//{}


