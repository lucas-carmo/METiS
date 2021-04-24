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
	double m_windDir;
	double m_windRefHeight;
	double m_windExp;

	std::vector<unsigned int> m_waveProbeID;
	std::vector<vec::fixed<3>> m_waveProbe; // Coordinates of the points where the wave characteristics (elevation, velocity, etc) are calculated for output

	/*
	Data to specify the numerical analysis
	*/
	double m_timeStep;
	double m_printStep;
	double m_timeTotal;
	double m_timeRamp;
	double m_time = 0;

	// Members that store variables evaluated at the beginning of the simulation, using IFFT or simple summation.
	bool m_flagIFFT{ false };
	bool m_shouldInterp{ false };
	uword m_ind4interp1;
	uword m_ind4interp2;
	vec m_timeArray;
	vec m_timeRampArray;
	mat m_waveElevArray;
	mat m_waveVel1stArray_x; // Component x of the first-order wave velocity
	mat m_waveVel1stArray_y;
	mat m_waveVel1stArray_z;	

public:
	ENVIR();

	/*****************************************************
		Setters
	*****************************************************/
	void setCurrentTime(const double time);
	void setTimeStep(const double timeStep);
	void setPrintStep(const double printStep);
	void setTimeTotal(const double timeTotal);
	void setTimeRamp(const double timeRamp);
	void setGravity(const double gravity);
	void setWatDens(const double watDens);
	void setWatDepth(const double watDepth);
	void setWaveStret(const unsigned int waveStret);
	void setAirDens(const double airDens);
	void setWindRefVel(const double windRefVel);
	void setWindDir(const double windDir);
	void setWindRefHeight(const double windRefHeight);
	void setWindExp(const double windExp);

	void addNode(const unsigned int nodeID, const double nodeCoordX, const double nodeCoordY, const double nodeCoordZ);
	void addRegularWave(const std::string &waveType, const double height, const double freqORperiod, const double direction, const double phase);
	void addJonswap(const double Hs, const double Tp, const double gamma, const double direction, const double wlow, const double whigh, const int numberOfRegularWaves, const double dwMax);
	void addWaveElevSeries(const std::string &elevFlPath, const double direction);

	void addWaveProbe(const unsigned int ID);
	void evaluateWaveKinematics();

	/*****************************************************
		Getters
	*****************************************************/
	double timeStep() const;
	double printStep() const;
	double timeTotal() const;
	double time() const;

	double gravity() const;
	double watDensity() const;
	double watDepth() const;
	unsigned int waveStret() const;
	double airDensity() const;
	double windRefVel() const;
	double windRefHeight() const;
	double windDir() const;
	double windExp() const;

	unsigned int numberOfWaveComponents() const;
	const Wave& getWave(unsigned int waveIndex) const;

	double waveElevAtProbe(const unsigned int ID) const;
	vec::fixed<3> waveVelAtProbe(const unsigned int ID) const;
	
	bool getFlagIFFT() const;	
	bool shouldInterp() const;
	uword getInd4interp1() const;
	uword getInd4interp2() const;
	const vec& getTimeArray() const;
	const vec& getRampArray() const;

	/*****************************************************
		Printing
	*****************************************************/
	std::string printTimeStep() const;
	std::string printPrintStep() const;
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
	double ramp(double time) const;	
	double waveElev(const double x, const double y) const;
	double waveElev(const double x, const double y, const double time) const;
	double wavePressure(const vec::fixed<3> &coord) const;
	double wavePressure(const vec::fixed<3> &coord, const double time) const;
	double wavePressure_2ndOrd(const vec::fixed<3> &coord) const;
	double wavePressure_2ndOrd(const vec::fixed<3> &coord, const double time) const;
	vec::fixed<3> u1(const vec::fixed<3> &coord, const double zwl) const;
	vec::fixed<3> u1(const vec::fixed<3> &coord, const double zwl, const double time) const;
	vec::fixed<3> du1dt(const vec::fixed<3> &coord, const double zwl) const;
	vec::fixed<3> du1dt(const vec::fixed<3> &coord, const double zwl, const double time) const;
	vec::fixed<3> du2dt(const vec::fixed<3> &coord) const;
	vec::fixed<3> du2dt(const vec::fixed<3> &coord, const double time) const;
	vec::fixed<3> du1dx(const vec::fixed<3> &coord, const double zwl) const;
	vec::fixed<3> du1dx(const vec::fixed<3> &coord, const double zwl, const double time) const;
	vec::fixed<3> du1dy(const vec::fixed<3> &coord, const double zwl) const;
	vec::fixed<3> du1dy(const vec::fixed<3> &coord, const double zwl, const double time) const;
	vec::fixed<3> du1dz(const vec::fixed<3> &coord, const double zwl) const;
	vec::fixed<3> du1dz(const vec::fixed<3> &coord, const double zwl, const double time) const;
	vec::fixed<3> da1dx(const vec::fixed<3> &coord, const double zwl) const;
	vec::fixed<3> da1dx(const vec::fixed<3> &coord, const double zwl, const double time) const;
	vec::fixed<3> da1dy(const vec::fixed<3> &coord, const double zwl) const;
	vec::fixed<3> da1dy(const vec::fixed<3> &coord, const double zwl, const double time) const;
	vec::fixed<3> da1dz(const vec::fixed<3> &coord, const double zwl) const;
	vec::fixed<3> da1dz(const vec::fixed<3> &coord, const double zwl, const double time) const;

	// Functions to evaluate properties of each wave
	cx_double waveElev_coef(const double x, const double y, const unsigned int waveIndex) const;
	cx_double wavePressure_coef(const double x, const double y, const double z, const unsigned int waveIndex) const;
	cx_double wavePressure_2ndOrd_coef(const vec::fixed<3> &coord, const unsigned int waveIndex1, const unsigned int waveIndex2) const;
	cx_vec::fixed<3> u1_coef(const double x, const double y, const double z, const unsigned int waveIndex) const;
	cx_vec::fixed<3> du1dt_coef(const double x, const double y, const double z, const unsigned int waveIndex) const;
	cx_vec::fixed<3> du2dt_coef(const vec::fixed<3> &coord, const unsigned int waveIndex1, const unsigned int waveIndex2) const;	
	cx_vec::fixed<3> du1dx_coef(const double x, const double y, const double z, const unsigned int waveIndex) const;
	cx_vec::fixed<3> du1dy_coef(const double x, const double y, const double z, const unsigned int waveIndex) const;
	cx_vec::fixed<3> du1dz_coef(const double x, const double y, const double z, const unsigned int waveIndex) const;
	cx_vec::fixed<3> da1dx_coef(const double x, const double y, const double z, const unsigned int waveIndex) const;
	cx_vec::fixed<3> da1dy_coef(const double x, const double y, const double z, const unsigned int waveIndex) const;
	cx_vec::fixed<3> da1dz_coef(const double x, const double y, const double z, const unsigned int waveIndex) const;

	double windVel_X(const vec::fixed<3> &coord) const;
	double windVel_Y(const vec::fixed<3> &coord) const;
};

// JONSWAP wave spectrum considering frequency in rad/s
double JONSWAP(const double w, const double Tp, const double Hs, const double gamma);
