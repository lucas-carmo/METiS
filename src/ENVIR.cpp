#include "ENVIR.h"
#include "IO.h"

#include <iostream>
#include <vector>
#include <algorithm>    // std::binary_search
#include <utility> // For std::move


/*****************************************************
	Constructors
*****************************************************/
ENVIR::ENVIR()
{
	// Initialize with NaN so we can check whether they were defined later
	m_gravity = arma::datum::nan;
	m_watDepth = arma::datum::nan;
}

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
	readDataFromString(data, m_gravity);
}



void ENVIR::readWatDens(const std::string &data)
{
	readDataFromString(data, m_watDens);
}



void ENVIR::readWatDepth(const std::string &data)
{
	readDataFromString(data, m_watDepth);
}


void ENVIR::addWave(const Wave &wave)
{
	// Check whether the water depth was defined
	if ( !is_finite(m_watDepth) )
	{
		throw std::runtime_error("You should specify the water depth before the waves. Error in input line " + std::to_string(IO::getInLineNumber()) + ".");
	}

	// Check whether the acceleration of gravity was defined
	if (!is_finite(m_gravity))
	{
		throw std::runtime_error("You should specify the gravity before the waves. Error in input line " + std::to_string(IO::getInLineNumber()) + ".");
	}

	m_wave.push_back(Wave(wave.height(), wave.period(), wave.direction(), m_watDepth, m_gravity));
}


void ENVIR::addWaveLocation(const std::string &data)
{
	// The wave locations are specified by node IDs separated by tabs or white-spaces	
	std::vector<std::string> input = stringTokenize(data, " \t");	
	
	// Check whether input is not empty
	if (input.empty())
	{
		throw std::runtime_error("You should specify at least one node ID for defining a wave location. Error in input line " + std::to_string(IO::getInLineNumber()) + ".");
	}

	// Check whether nodes were specified
	if (this->isNodeEmpty())
	{		
		throw std::runtime_error("Nodes should be specified before adding wave locations. Error in input line " + std::to_string(IO::getInLineNumber()) + ".");
	}

	// For each of the node IDs:
	for (int ii = 0; ii < input.size(); ++ii)
	{
		unsigned int nodeID(0); // Initialize a variable to read the node ID
		readDataFromString(input.at(ii), nodeID); // Read the node ID specified as a string to the nodeID variable
		m_waveLocation.push_back( this->getNode(nodeID) ); // Get the node coordinate and add it to m_waveLocation
		m_waveLocation.back().at(2) = 0; // Set z=0
	}
}

void ENVIR::addNode(const std::string &data)
{
	// Nodes are specified by a vec with four components: ID, X coord, Y coord, and Z coord. 
	// They are separated by commas in the input string.
	std::vector<std::string> input = stringTokenize(data, ",");

	if (input.size() != 4)
	{
		throw std::runtime_error("Unable to read the node in input line " + std::to_string(IO::getInLineNumber()) + ". Wrong number of parameters.");
		return;
	}

	// Read node ID
	unsigned int nodeID{0};
	readDataFromString( input.at(0), nodeID );

	if (m_nodesID.size() != 0) // If this is not the first node that will be added to m_nodesID
	{
		if (nodeID <= m_nodesID.back()) // Then verify if its ID is larger than the previous one, thus garanteeing that m_nodesID is in ascending order (this is needed to use binary search to find nodes IDs)
		{
			throw std::runtime_error( "Nodes must be organized in ascending order. Error in input line " + std::to_string(IO::getInLineNumber()) + ".");
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
		throw std::runtime_error("Unable to find node with ID " + std::to_string(ID) + ". Error in input line " + std::to_string(IO::getInLineNumber()) + ".");
	}

	return node_coord;
}

double ENVIR::timeStep() const
{
	return m_timeStep;
}

double ENVIR::timeTotal() const
{
	return m_timeTotal;
}

double ENVIR::time() const
{
	return m_time;
}

double ENVIR::gravity() const
{
	return m_gravity;
}

double ENVIR::watDepth() const
{
	return m_watDepth;
}

double ENVIR::watDensity() const
{
	return m_watDens;
}

/*****************************************************
	Printing
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
		output = output + "Wave number: " + std::to_string( m_wave.at(ii).waveNumber() ) + "\n";
		output = output + "Length: " + std::to_string(m_wave.at(ii).length()) + "\n";
		output = output + "Direction: " + std::to_string(m_wave.at(ii).direction()) + "\n\n";
	}
	return output;
}

std::string ENVIR::printWaveLocation() const
{
	std::string output = "";
	for (int ii = 0; ii < m_waveLocation.size(); ++ii)
	{
		output = output + "Location #" + std::to_string(ii) + ": (" + std::to_string(m_waveLocation.at(ii).at(0)) 
						+ "," + std::to_string(m_waveLocation.at(ii).at(1)) + "," + std::to_string(m_waveLocation.at(ii).at(2)) + ")\n";
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

bool ENVIR::isWaveLocationEmpty() const
{
	return m_waveLocation.empty();
}

void ENVIR::stepTime()
{
	m_time += m_timeStep;
}

arma::vec::fixed<3> ENVIR::fluidVel(double x, double y, double z) const
{
	arma::vec::fixed<3> vel = {0,0,0};
	for (int ii = 0; ii < m_wave.size(); ++ii)
	{
		vel += m_wave.at(ii).fluidVel(x, y, z, m_time, m_watDepth);
	}

	return vel;
}

arma::vec::fixed<3> ENVIR::fluidAcc(double x, double y, double z) const
{
	arma::vec::fixed<3> acc = { 0,0,0 };
	for (int ii = 0; ii < m_wave.size(); ++ii)
	{
		acc += m_wave.at(ii).fluidAcc(x, y, z, m_time, m_watDepth);
	}

	return acc;
}
