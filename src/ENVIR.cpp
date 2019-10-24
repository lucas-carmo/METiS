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
void ENVIR::setTimeStep(const double timeStep)
{
	m_timeStep = timeStep;
}

void ENVIR::setTimeTotal(const double timeTotal)
{
	m_timeTotal = timeTotal;
}

void ENVIR::setTimeRamp(const double timeRamp)
{
    m_timeRamp = timeRamp;
}

void ENVIR::setGravity(const double gravity)
{
	m_gravity = gravity;
}

void ENVIR::setWatDens(const double watDens)
{
	m_watDens = watDens;
}

void ENVIR::setWatDepth(const double watDepth)
{
	m_watDepth = watDepth;
}

void ENVIR::setAirDens(const double airDens)
{
	m_airDens = airDens;
}

void ENVIR::setWindRefVel(const double windRefVel)
{
	m_windRefVel = windRefVel;
}

void ENVIR::setWindRefHeight(const double windRefHeight)
{
	m_windRefHeight = windRefHeight;
}

void ENVIR::setWindExp(const double windExp)
{
	m_windExp = windExp;
}

void ENVIR::addNode(const unsigned int nodeID, const double nodeCoordX, const double nodeCoordY, const double nodeCoordZ)
{
	if (m_nodesID.size() != 0) // If this is not the first node that will be added to m_nodesID
	{
		if (nodeID <= m_nodesID.back()) // Then verify if its ID is larger than the previous one, thus garanteeing that m_nodesID is in ascending order (this is needed to use binary search to find nodes IDs)
		{
			throw std::runtime_error( "Nodes must be organized in ascending order. Error in input line " + std::to_string(IO::getInLineNumber()) + ".");
		}
	}

	m_nodesID.push_back( nodeID );	
	m_nodesCoord.push_back(vec::fixed<3> {nodeCoordX, nodeCoordY, nodeCoordZ});
}

void ENVIR::addWave(const std::string &wholeWaveLine)
{
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

	// Test if the keyword is related to regular waves
	if (caseInsCompare(getKeyword(wholeWaveLine), "TRWave")
		|| caseInsCompare(getKeyword(wholeWaveLine), "FRWave")
		|| caseInsCompare(getKeyword(wholeWaveLine), "WRWave"))
	{
		Wave tempWave(wholeWaveLine);
		m_wave.push_back(Wave(tempWave.height(), tempWave.period(), tempWave.direction(), tempWave.phase(), m_watDepth, m_gravity));
	}

	// Check if it is a JONSWAP spectrum
	else if (caseInsCompare(getKeyword(wholeWaveLine), "JONSW"))
	{
		jonswap(wholeWaveLine);
	}

	else
	{
		throw std::runtime_error("Unknown keyword '" + getKeyword(wholeWaveLine) + "' in input line " + std::to_string(IO::getInLineNumber()) + ".");
	}
}

void ENVIR::jonswap(const std::string &wholeWaveLine)
{
	// JONSWAP characteristics are divided by a space or a tab.
	// The characteristics are 
	// 1) Wave significant height
	// 2) Wave peark period
	// 3) Gamma
	// 4) Direction (in degrees)
	// 5) Lowest frequency for the spectrum (rad/s)
	// 6) Highest frequency for the spectrum (rad/s)
	std::vector<std::string> input = stringTokenize(getData(wholeWaveLine), " \t");

	// Check if there are exactly six outputs
	if (input.size() != 6)
	{
		throw std::runtime_error("Unable to read JONSWAP spectrum in input line " + std::to_string(IO::getInLineNumber()) + ". Wrong number of parameters.");
	}
	
	double Hs;
	double Tp;
	double gamma;
	double dir;
	double wlow;
	double whigh;

	// Finally, read data
	readDataFromString(input.at(0), Hs);
	readDataFromString(input.at(1), Tp);
	readDataFromString(input.at(2), gamma);
	readDataFromString(input.at(3), dir);
	readDataFromString(input.at(4), wlow);
	readDataFromString(input.at(5), whigh);

	// Frequency resolution and maximum possible frequency are determined by the total time simulation and time step
	// See Merigaud, A. and Ringwood, John V. - Free-Surface Time-Series Generation for Wave Energy Applications - IEEE J. of Oceanic Eng. - 2018
	double dw = 2 * arma::datum::pi / m_timeTotal;
	double wmax = arma::datum::pi / m_timeStep;

	// Wave parameters that will be calculated using the JONSWAP spectrum
	double height(0);
	double phase(0);

	// JONSWAP parameters
	double wp = 2* arma::datum::pi /Tp;
	double sigma(0);
	double A(0);
	double Sw(0);

	// Loop to create the waves
	int ii = 0;
	for (double w = wlow; w <= whigh; w += dw)
	{
		if (w <= wp)
		{
			sigma = 0.07;
		}
		else
		{
			sigma = 0.09;
		}

		A = std::exp(-0.5 * std::pow((w / wp - 1) / sigma, 2));
		Sw = 0.3125 * std::pow(Hs, 2) * Tp * std::pow(w / wp, -5) * std::exp(-1.25 * std::pow(wp / w, 4)) * (1 - 0.287*std::log(gamma)) * std::pow(gamma, A);

		height = 2 * std::sqrt(2 * Sw * dw);
		phase = -360 + (360 + 360) * randu(1, 1).at(0, 0);
		m_wave.push_back(Wave(height, 2*arma::datum::pi/w, dir, phase, m_watDepth, m_gravity));
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

			acc[0] += pow(w, 2) * A * khz_xy * cos(beta) * sin(k*cos(beta)*x + k * sin(beta)*y - w * t + phase);
			acc[1] += pow(w, 2) * A * khz_xy * sin(beta) * sin(k*cos(beta)*x + k * sin(beta)*y - w * t + phase);
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

			p += rho * g * A * khz * cos(k*cos(beta)*x + k * sin(beta)*y - w * t + phase);
		}
	}

	return p * ramp();
}



double ENVIR::windVel_X(const vec::fixed<3> &coord) const
{
	return ( windRefVel() * pow(coord[2] / windRefHeight(), windExp()) );
}