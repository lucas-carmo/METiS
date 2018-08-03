#include "ENVIR.h"
#include "IO.h"

#include <iostream>
#include <vector>
#include <algorithm>    // std::binary_search
#include <utility> // For std::move


/*****************************************************
	Setters
*****************************************************/
void ENVIR::readTimeStep(const std::string &data)
{
    readDataFromString(data, m_timeStep);
}



void ENVIR::readTimeTotal(const std::string &data)
{
    readDataFromString(data, m_timeTotal);
}



void ENVIR::readTimeRamp(const std::string &data)
{
    readDataFromString(data, m_timeRamp);
}



void ENVIR::readGrav(const std::string &data)
{
	readDataFromString(data, m_watDens);
}



void ENVIR::readWatDens(const std::string &data)
{
	readDataFromString(data, m_gravity);
}



void ENVIR::readWatDepth(const std::string &data)
{
	readDataFromString(data, m_watDepth);
}


void ENVIR::addWave(const Wave &wave)
{
	m_wave.push_back( wave );
}

void ENVIR::addNode(const std::string &data)
{
	// Nodes are specified by a vec with four components: ID, X coord, Y coord, and Z coord. 
	// They are separated by commas in the input string.
	std::vector<std::string> input = stringTokenize(data, ",");

	if (input.size() != 4)
	{
		throw std::runtime_error("Unable to read the node in line " + std::to_string(IO::getInLineNumber()) + ". Wrong number of parameters.");
		return;
	}

	// Read node ID
	unsigned int nodeID{0};
	readDataFromString( input.at(0), nodeID );

	if (m_nodesID.size() != 0) // If this is not the first node that will be added to m_nodesID
	{
		if (nodeID <= m_nodesID.back()) // Then verify if its ID is larger than the previous one, thus garanteeing that m_nodesID is in ascending order (this is needed to use binary search to find nodes IDs)
		{
			throw std::runtime_error( "Nodes must be organized in ascending order. Error in line " + std::to_string(IO::getInLineNumber()) + ".");
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
std::string ENVIR::printTimeStep() const
{
	return std::to_string(m_timeStep);
}

std::string ENVIR::printTimeTotal() const
{
	return std::to_string(m_timeTotal);
}

std::string ENVIR::printTimeRamp() const
{
	return std::to_string(m_timeRamp);
}

std::string ENVIR::printNodes() const
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

std::string ENVIR::printGrav() const
{
	return std::to_string(m_gravity);
}

std::string ENVIR::printWatDens() const
{
	return std::to_string(m_watDens);
}

std::string ENVIR::printWatDepth() const
{
	return std::to_string(m_watDepth);
}

std::string ENVIR::printWave() const
{
	std::string output = "";
	for (int ii = 0; ii < m_wave.size(); ++ii)
	{
		output = output + "Wave #" + std::to_string(ii) + "\n";
		output = output + "Height: " + std::to_string( m_wave.at(ii).height() ) + "\n";
		output = output + "Period: " + std::to_string( m_wave.at(ii).period() ) + "\n";
		output = output + "Direction: " + std::to_string( m_wave.at(ii).direction() ) + "\n\n";
	}
	return output;
}


/*****************************************************
	Other functions
*****************************************************/
bool ENVIR::isNodeEmpty() const
{
	return m_nodesID.empty();
}


// Return coordinates of a node based on its ID
// Throws a std::runtime_error if the node could not be found
arma::vec::fixed<3> ENVIR::getNode(unsigned int ID) const
{
	std::vector<unsigned int>::const_iterator iter = std::find(m_nodesID.begin(), m_nodesID.end(), ID); // Find node by its ID.
	vec::fixed<3> node_coord(fill::zeros);
	if (iter != m_nodesID.end())
	{
		auto index = std::distance(m_nodesID.begin(), iter); // Get index by the distance between the iterator found above and m_nodes.begin()
		node_coord = m_nodesCoord.at(index);		
	}
	else
	{
		throw std::runtime_error( "Unable to find node with ID " + std::to_string(ID) + " in line " + std::to_string(IO::getInLineNumber()) + "." );
	}

	return node_coord;
}