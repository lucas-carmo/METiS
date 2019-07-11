#include "ENVIR.h"
#include "IO.h"
#include "auxFunctions.h"

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

	m_wave.push_back(Wave(wave.height(), wave.period(), wave.direction(), wave.phase(), m_watDepth, m_gravity));
}

void ENVIR::addWave(const std::string &wholeWaveLine)
{
	// Test if the keyword is related to regular waves
	if (caseInsCompare(getKeyword(wholeWaveLine), "TRWave")
		|| caseInsCompare(getKeyword(wholeWaveLine), "FRWave")
		|| caseInsCompare(getKeyword(wholeWaveLine), "WRWave"))
	{
		addWave(Wave(wholeWaveLine));
	}

	// Check if it is a JONSWAP spectrum
	else if (caseInsCompare(getKeyword(wholeWaveLine), "JONSW"))
	{
		std::cout << "\n\n\nIsso ainda nao faz nada\n\n\n";
	}

	else
	{
		throw std::runtime_error("Unknown keyword '" + getKeyword(wholeWaveLine) + "' in input line " + std::to_string(IO::getInLineNumber()) + ".");
	}
	// Check whether the water depth was defined
	if (!is_finite(m_watDepth))
	{
		throw std::runtime_error("You should specify the water depth before the waves. Error in input line " + std::to_string(IO::getInLineNumber()) + ".");
	}

	// Check whether the acceleration of gravity was defined
	if (!is_finite(m_gravity))
	{
		throw std::runtime_error("You should specify the gravity before the waves. Error in input line " + std::to_string(IO::getInLineNumber()) + ".");
	}	
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
		m_waveLocationID.push_back(nodeID);
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
		output = output + "Direction: " + std::to_string(m_wave.at(ii).direction()) + "\n";
		output = output + "Phase: " + std::to_string(m_wave.at(ii).phase()) + "\n\n";
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


void ENVIR::printWaveCharact() const
{
	for (int ii = 0; ii < m_waveLocation.size(); ++ii)
	{
		IO::print2outLine(IO::OUTFLAG_WAVE_ELEV, m_waveLocationID[ii], waveElev(m_waveLocation[ii][0], m_waveLocation[ii][1]));
		IO::print2outLine(IO::OUTFLAG_WAVE_VEL, m_waveLocationID[ii], u1(m_waveLocation[ii]));
		IO::print2outLine(IO::OUTFLAG_WAVE_ACC, m_waveLocationID[ii], du1dt(m_waveLocation[ii]));
		IO::print2outLine(IO::OUTFLAG_WAVE_PRES, m_waveLocationID[ii], wavePressure(m_waveLocation[ii]));
	}
}


/*****************************************************
	Other functions
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


double ENVIR::waveElev(const double x, const double y) const
{
	double elev{ 0 };

	// Variables that are used for each wave
	double w(0), A(0), k(0), beta(0), phase(0);

	// We consider linear Airy waves, with velocity potential:
	// phi = g*A/w * cosh(k(z+h))/cosh(k*h) * sin(k*x - w*t)
	for (int ii = 0; ii < m_wave.size(); ++ii)
	{
		w = m_wave.at(ii).angFreq();
		A = m_wave.at(ii).amp();
		k = m_wave.at(ii).waveNumber();
		beta = m_wave.at(ii).direction() * arma::datum::pi / 180.;
		phase = m_wave.at(ii).phase();

		elev += A * cos(k*cos(beta)*x + k * sin(beta)*y - w * m_time + phase);
	}

	return elev * ramp();
}


vec::fixed<3> ENVIR::u1(const vec::fixed<3> &coord) const
{
	arma::vec::fixed<3> vel = {0,0,0};

	double x = coord[0];
	double y = coord[1];
	double z = coord[2];
	double h = m_watDepth;
	double t = m_time;

	// Variables that are used for each wave
	double w(0), A(0), k(0), beta(0), phase(0);
	double khz_xy(0), khz_z(0);

	// This formulation is valid only below the mean water level, i.e. z <= 0
	if (z <= 0)
	{
		for (int ii = 0; ii < m_wave.size(); ++ii)
		{
			w = m_wave.at(ii).angFreq();
			A = m_wave.at(ii).amp();
			k = m_wave.at(ii).waveNumber();
			beta = m_wave.at(ii).direction() * arma::datum::pi / 180.;
			phase = m_wave.at(ii).phase();

			// When k*h is too high, which happens for deep water/short waves, sinh(k*h) and cosh(k*h) become too large and are considered "inf".
			// Hence, we chose a threshold of 10, above which the deep water approximation is employed.
			if (k*h >= 10)
			{
				khz_xy = exp(k*z);
				khz_z = khz_xy;
			}
			else
			{
				khz_xy = cosh(k * (z + h)) / sinh(k*h);
				khz_z = sinh(k * (z + h)) / sinh(k*h);
			}

			vel[0] += w * A * khz_xy * cos(beta) * cos(k*cos(beta)*x + k * sin(beta)*y - w * t + phase);
			vel[1] += w * A * khz_xy * sin(beta) * cos(k*cos(beta)*x + k * sin(beta)*y - w * t + phase);
			vel[2] += w * A * khz_z * sin(k*cos(beta)*x + k * sin(beta)*y - w * t + phase);
		}
	}

	return vel * ramp();
}

vec::fixed<3> ENVIR::du1dt(const vec::fixed<3> &coord) const
{
	arma::vec::fixed<3> acc = {0,0,0};

	double x = coord[0];
	double y = coord[1];
	double z = coord[2];
	double h = m_watDepth;
	double t = m_time;

	// Variables that are used for each wave
	double w(0), A(0), k(0), beta(0), phase(0);
	double khz_xy(0), khz_z(0);

	// We consider linear Airy waves, with velocity potential:
	// phi = g*A/w * cosh(k(z+h))/cosh(k*h) * sin(k*x - w*t)
	// This formulation is valid only below the mean water level, i.e. z <= 0
	if (z <= 0)
	{
		for (int ii = 0; ii < m_wave.size(); ++ii)
		{
			w = m_wave.at(ii).angFreq();
			A = m_wave.at(ii).amp();
			k = m_wave.at(ii).waveNumber();
			beta = m_wave.at(ii).direction() * arma::datum::pi / 180.;
			phase = m_wave.at(ii).phase();

			// When k*h is too high, which happens for deep water/short waves, sinh(k*h) and cosh(k*h) become too large and are considered "inf".
			// Hence, we chose a threshold of 10, above which the deep water approximation is employed.			
			if (k*h >= 10)
			{
				khz_xy = exp(k*z);
				khz_z = khz_xy;
			}
			else
			{
				khz_xy = cosh(k * (z + h)) / sinh(k*h);
				khz_z = sinh(k * (z + h)) / sinh(k*h);
			}

			acc[0] += pow(w, 2) * A * khz_xy * cos(beta) * sin(k*cos(beta)*x + k * sin(beta)*y - w * t);
			acc[1] += pow(w, 2) * A * khz_xy * sin(beta) * sin(k*cos(beta)*x + k * sin(beta)*y - w * t);
			acc[2] += -pow(w, 2) * A * khz_z * cos(k*cos(beta)*x + k * sin(beta)*y - w * t + phase);
		}
	}

	return acc * ramp();
}


double ENVIR::wavePressure(const vec::fixed<3> &coord) const
{
	double p(0);

	double x = coord[0];
	double y = coord[1];
	double z = coord[2];
	double h = m_watDepth;
	double t = m_time;
	double rho = m_watDens;
	double g = m_gravity;

	// Variables that are used for each wave
	double w(0), A(0), k(0), beta(0), phase(0);

	// We consider linear Airy waves, with velocity potential:
	// phi = g*A/w * cosh(k(z+h))/cosh(k*h) * sin(k*x - w*t)
	// This formulation is valid only below the mean water level, i.e. z <= 0
	if (z <= 0)
	{
		for (int ii = 0; ii < m_wave.size(); ++ii)
		{
			w = m_wave.at(ii).angFreq();
			A = m_wave.at(ii).amp();
			k = m_wave.at(ii).waveNumber();
			beta = m_wave.at(ii).direction() * arma::datum::pi / 180.;
			phase = m_wave.at(ii).phase();

			// When k*h is too high, which happens for deep water/short waves, sinh(k*h) and cosh(k*h) become too large and are considered "inf".
			// Hence, we chose a threshold of 10, above which the deep water approximation is employed.
			double khz(0);
			if (k*h >= 10)
			{
				khz = exp(k*z);
			}
			else
			{
				khz = cosh(k * (z + h)) / cosh(k*h);
			}

			p = rho * g * A * khz * cos(k*cos(beta)*x + k * sin(beta)*y - w * t + phase);
		}
	}

	return p * ramp();
}



double ENVIR::windVel_X(const vec::fixed<3> &coord) const
{
	return ( windRefVel() * pow(coord[2] / windRefHeight(), windExp()) );
}