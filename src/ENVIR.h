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
  std::vector<Wave> m_wave; // The sea is specified by a vector with regular wave components. You can add a wave component using the method addRegularWave(), or many components using addJonswap().
	unsigned int m_waveStret;

	double m_airDens;
  double m_windRefVel;
	double m_windRefHeight;
  double m_windExp;

	std::vector<unsigned int> m_waveProbeID;
	std::vector<vec::fixed<3>> m_waveProbe; // Coordinates of the points where the wave characteristics (elevation, velocity, etc) are calculated for output

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
	void setWaveStret(const unsigned int waveStret);
	void setAirDens(const double airDens);
	void setWindRefVel(const double windRefVel);
	void setWindRefHeight(const double windRefHeight);
	void setWindExp(const double windExp);

	void addNode(const unsigned int nodeID, const double nodeCoordX, const double nodeCoordY, const double nodeCoordZ);
	void addRegularWave(const std::string &waveType, const double height, const double freqORperiod, const double direction, const double phase);
	void addJonswap(const double Hs, const double Tp, const double gamma, const double direction, const double wlow, const double whigh);

	void addWaveProbe(const unsigned int ID);

	/*****************************************************
		Getters
	*****************************************************/
	double timeStep() const;
	double timeTotal() const;
	double time() const;

	double gravity() const;
	double watDensity() const;
	double watDepth() const;
	unsigned int waveStret() const;
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
	std::string printWaveProbe() const;

	void printWaveCharact() const; // Print the wave characteristics (elevation, velocity, etc) specified for output in the locations given by m_waveLocation

	/*****************************************************
		Main functions for calculation
	*****************************************************/
	bool isNodeEmpty() const;
	bool isWaveProbeEmpty() const;
	arma::vec::fixed<3> getNode(unsigned int ID) const;

	void stepTime();
	void stepTime(double const step);

	double ramp() const;
	double waveElev(const double x, const double y, const unsigned int waveIndex) const;
	double waveElev(const double x, const double y) const;
	double wavePressure(const vec::fixed<3> &coord, const unsigned int waveIndex) const;
	double wavePressure(const vec::fixed<3> &coord) const;
	vec::fixed<3> u1(const vec::fixed<3> &coord, const double zwl, const unsigned int waveIndex) const;
	vec::fixed<3> u1(const vec::fixed<3> &coord, const double zwl) const;
	vec::fixed<3> du1dt(const vec::fixed<3> &coord, const double zwl, const unsigned int waveIndex) const;
	vec::fixed<3> du1dt(const vec::fixed<3> &coord, const double zwl) const;
	vec::fixed<3> du2dt(const vec::fixed<3> &coord, const unsigned int waveIndex1, const unsigned int waveIndex2) const;
	vec::fixed<3> du2dt(const vec::fixed<3> &coord) const;

	double windVel_X(const vec::fixed<3> &coord) const;
};
