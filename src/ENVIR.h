#pragma once

#include <vector>
#include <array>
#include <string>
#include <armadillo>

#include "Wave.h"

using namespace arma;

class ENVIR
{
private:
    /*
    Data to specify the environment
    */
	std::vector< unsigned int > m_nodesID; // Nodes provide a spatial description of the environment
	std::vector< vec::fixed<3> > m_nodesCoord;   

    double m_gravity;
    double m_watDens;
    double m_watDepth;
    std::vector<Wave> m_wave;
	std::vector<unsigned int> m_waveLocationID;
    std::vector<vec::fixed<3>> m_waveLocation; // Coordinates of the points where the wave characteristics (elevation, velocity, etc) are calculated for output
    double m_airDens;
    double m_windRefVel;
	double m_windRefHeight;
    double m_windExp;

    /*
    Data to specify the numerical analysis
    */
    double m_timeStep;
    double m_timeTotal;
    double m_timeRamp;
    double m_time = 0;

public:
	ENVIR();

	/*****************************************************
		Setters
	*****************************************************/
	void setTimeStep(const double timeStep);
	void setTimeTotal(const double timeTotal);
	void setTimeRamp(const double timeRamp);
	void setGravity(const double gravity);
	void setWatDens(const double watDens);
	void setWatDepth(const double watDepth);
	void setAirDens(const double airDens);	
	void setWindRefVel(const double windRefVel);
	void setWindRefHeight(const double windRefHeight);
	void setWindExp(const double windExp);

	void addNode(const unsigned int nodeID, const double nodeCoordX, const double nodeCoordY, const double nodeCoordZ);
	
	void addWave(const std::string &wholeWaveLine);
	void jonswap(const std::string &wholeWaveLine);

	void addWaveLocation(const std::string &data);

	/*****************************************************
		Getters
	*****************************************************/    
    double timeStep() const;
    double timeTotal() const;
    double time() const;

    double gravity() const;
	double watDensity() const;
    double watDepth() const;
	double airDensity() const;
	double windRefVel() const;
	double windRefHeight() const;
	double windExp() const;
    

	/*****************************************************
		Printing
	*****************************************************/
	std::string printTimeStep() const;
	std::string printTimeTotal() const;
	std::string printTimeRamp() const;
    std::string printNodes() const;    
	std::string printGrav() const;
	std::string printWatDens() const;
	std::string printWatDepth() const;
	std::string printWave() const;
	std::string printWaveLocation() const;

	void printWaveCharact() const; // Print the wave characteristics (elevation, velocity, etc) specified for output in the locations given by m_waveLocation

	/*****************************************************
		Other functions
	*****************************************************/
    bool isNodeEmpty() const;
    bool isWaveLocationEmpty() const;
    arma::vec::fixed<3> getNode(unsigned int ID) const;

    void stepTime();
    void stepTime(double const step);

    double ramp() const;
	double waveElev(const double x, const double y) const;
	vec::fixed<3> u1(const vec::fixed<3> &coord) const;
	vec::fixed<3> du1dt(const vec::fixed<3> &coord) const;
	double wavePressure(const vec::fixed<3> &coord) const;

	double windVel_X(const vec::fixed<3> &coord) const;
};

