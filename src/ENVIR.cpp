#include "ENVIR.h"
#include "IO.h"

#include <iostream>
#include <vector>
#include <algorithm>    // std::binary_search
#include <utility> // For std::move

using namespace arma;

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
void ENVIR::readTypeAnalysis(const std::string &data)
{
	readDataFromString(data, m_typeAnalysis);
}

void ENVIR::readDOFs(const std::string &data)
{
	// The flags for each of the six degrees of freedom are separated by white spaces in the input string (whitespace or tab)
	std::vector<std::string> input = stringTokenize(data, " \t");

	// Check number of inputs
	if (input.size() != 6)
	{
		throw std::runtime_error("Unable to read the DoFs in input line " + std::to_string(IO::getInLineNumber()) + ". Wrong number of parameters.");
	}	

	// Read data
	readDataFromString(input.at(0), m_dofs[0]);
	readDataFromString(input.at(1), m_dofs[1]);
	readDataFromString(input.at(2), m_dofs[2]);
	readDataFromString(input.at(3), m_dofs[3]);
	readDataFromString(input.at(4), m_dofs[4]);
	readDataFromString(input.at(5), m_dofs[5]);
}


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

void ENVIR::readUseBEMT(const std::string &data)
{
	readDataFromString(data, m_useBEMT);
}

void ENVIR::readUseTipLoss(const std::string &data)
{
	readDataFromString(data, m_useTipLoss);
}

void ENVIR::readUseHubLoss(const std::string &data)
{
	readDataFromString(data, m_useHubLoss);
}

void ENVIR::readUseSkewCorr(const std::string &data)
{
	readDataFromString(data, m_useSkewCorr);
}

void ENVIR::readGrav(const std::string &data)
{
	readDataFromString(data, m_gravity);
}

void ENVIR::readWatDens(const std::string &data)
{
	readDataFromString(data, m_watDens);
}

void ENVIR::readAirDens(const std::string &data)
{
	readDataFromString(data, m_airDens);
}

void ENVIR::readWindRefVel(const std::string &data)
{
	readDataFromString(data, m_windRefVel);
}

void ENVIR::readWindRefHeight(const std::string &data)
{
	readDataFromString(data, m_windRefHeight);
}

void ENVIR::readWindExp(const std::string &data)
{
	readDataFromString(data, m_windExp);
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

bool ENVIR::useBEMT() const
{
	return m_useBEMT;
}

bool ENVIR::useTipLoss() const
{
	return m_useTipLoss;
}

bool ENVIR::useHubLoss() const
{
	return m_useHubLoss;
}

bool ENVIR::useSkewCorr() const
{
	return m_useSkewCorr;
}

double ENVIR::gravity() const
{
	return m_gravity;
}

double ENVIR::watDensity() const
{
	return m_watDens;
}

double ENVIR::watDepth() const
{
	return m_watDepth;
}

double ENVIR::airDensity() const
{
	return m_airDens;
}

double ENVIR::windRefVel() const
{
	return m_windRefVel;
}

double ENVIR::windRefHeight() const
{
	return m_windRefHeight;
}

double ENVIR::windExp() const
{
	return m_windExp;
}

/*****************************************************
	Printing
*****************************************************/
std::string ENVIR::printTypeAnalysis() const
{
	return m_typeAnalysis;
}

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
bool ENVIR::isSurgeActive() const
{
	return m_dofs[0];
}

bool ENVIR::isSwayActive() const
{
	return m_dofs[1];
}

bool ENVIR::isHeaveActive() const
{
	return m_dofs[2];
}

bool ENVIR::isRollActive() const
{
	return m_dofs[3];
}

bool ENVIR::isPitchActive() const
{
	return m_dofs[4];
}

bool ENVIR::isYawActive() const
{
	return m_dofs[5];
}

bool ENVIR::isTypeFOWT() const
{
	return (caseInsCompare(m_typeAnalysis, "FOWT"));
}

bool ENVIR::isTypeFixedOffshore() const
{
	return (caseInsCompare(m_typeAnalysis, "FixedOffshore"));
}

bool ENVIR::isTypeOnshore() const
{
	return (caseInsCompare(m_typeAnalysis, "Onshore"));
}

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

void ENVIR::stepTime(double const step)
{
	m_time += step;
}


double ENVIR::ramp() const
{
	double ramp{1};

	if (m_time < m_timeRamp)
	{
		ramp = 0.5 * ( 1 - cos(datum::pi * m_time / m_timeRamp) );
	}

	return ramp;
}


vec::fixed<3> ENVIR::fluidVel(const vec::fixed<3> &coord) const
{
	arma::vec::fixed<3> vel = {0,0,0};
	for (int ii = 0; ii < m_wave.size(); ++ii)
	{
		vel += m_wave.at(ii).fluidVel(coord, m_time, m_watDepth);
	}

	vel = vel*ramp();
	return vel;
}

vec::fixed<3> ENVIR::fluidAcc(const vec::fixed<3> &coord) const
{
	arma::vec::fixed<3> acc = { 0,0,0 };
	for (int ii = 0; ii < m_wave.size(); ++ii)
	{
		acc += m_wave.at(ii).fluidAcc(coord, m_time, m_watDepth);
	}

	acc = acc*ramp();
	return acc;
}

double ENVIR::wavePressure(const vec::fixed<3> &coord) const
{
	double p{ 0 };
	for (int ii = 0; ii < m_wave.size(); ++ii)
	{
		p += m_wave.at(ii).pressure(coord, m_time, m_watDens, m_gravity, m_watDepth);
	}

	p = p*ramp();
	return p;
}


double ENVIR::windVel_X(const vec::fixed<3> &coord) const
{
	return ( windRefVel() * std::pow(coord[2] / windRefHeight(), windExp()) );
}