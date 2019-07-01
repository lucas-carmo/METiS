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

    bool m_useTipLoss;
    bool m_useHubLoss;
    // bool m_IncTIFac;
    // bool m_IncDragAIFac;
    // bool m_IncDragTIFac;
    bool m_useSkewCorr;
    bool m_TwrLoads;

public:
	ENVIR();

	/*****************************************************
		Setters
	*****************************************************/   
    void readTimeStep(const std::string &data);
    void readTimeTotal(const std::string &data);
    void readTimeRamp(const std::string &data);
	void readUseTipLoss(const std::string &data);
	void readUseHubLoss(const std::string &data);
	void readUseSkewCorr(const std::string &data);
	void addNode(const std::string &data);

    void readGrav(const std::string &data);
    void readWatDens(const std::string &data);
    void readWatDepth(const std::string &data);
	void readAirDens(const std::string &data);
	void readWindRefVel(const std::string &data);
	void readWindRefHeight(const std::string &data);
	void readWindExp(const std::string &data);

	void addWave(const Wave &wave);
	void addWaveLocation(const std::string &data);

	/*****************************************************
		Getters
	*****************************************************/    
    double timeStep() const;
    double timeTotal() const;
    double time() const;
	bool useTipLoss() const;
	bool useHubLoss() const;
	bool useSkewCorr() const;

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
	vec::fixed<3> fluidVel(const vec::fixed<3> &coord) const;
	vec::fixed<3> fluidAcc(const vec::fixed<3> &coord) const;
	double wavePressure(const vec::fixed<3> &coord) const;

	double windVel_X(const vec::fixed<3> &coord) const;
};

