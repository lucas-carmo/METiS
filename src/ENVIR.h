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

	// Wind Characteristics
	float m_airDens;
	float m_windDir;
	float m_windDirCos;
	float m_windDirSin;
	float m_windRefVel;	
	float m_windRefHeight;
	float m_windExp;

	// For uniform wind from .wnd file
	// Decided to include '_wnd' in the names of the variables
	// to avoid confusion with other members with similar names
	bool m_flagUnifWind = false;
	std::string m_wnd_fileName = "";
	float m_wnd_refLength{ 1 };
	arma::fvec m_wnd_time;
	arma::fvec m_wnd_windSpeed;
	arma::fvec m_wnd_windDir;
	arma::fvec m_wnd_vertSpeed;
	arma::fvec m_wnd_horizSheer;
	arma::fvec m_wnd_windExp;
	arma::fvec m_wnd_linVertSheer;
	arma::fvec m_wnd_gustSpeed;

	// For turbulent wind from .bts file
	bool m_flagTurbWind = false;
	bool m_flagPeriodicTurb = false;
	std::string m_turbFileName = "";
	arma::field<arma::fcube> m_windVelocity;
	arma::fvec m_windGrid_x; // Coordinates of the turbulent grid points
	arma::fvec m_windGrid_y;
	arma::fvec m_windGrid_z;
	float m_windDt; // Time discretization of the turbulent input file
	float m_windTimeTotal; // Max simulation time in TurbSim file
	std::vector<unsigned int> m_windProbeID;
	std::vector<vec::fixed<3>> m_windProbe; // Coordinates of the points where the wind velocity is calculated for output

	// End of wind characteristics

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
	mat m_wavePress2ndArray;

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

	void setWindFromTurbFile(const std::string &fileName);
	void setWindFromUnifFile(const std::string &fileName);
	void setWindRefLength(const double windRefLength);
	void setWindRefVel(const double windRefVel);
	void setWindDir(const double windDir);
	void setWindRefHeight(const double windRefHeight);
	void setWindExp(const double windExp);

	void addNode(const unsigned int nodeID, const double nodeCoordX, const double nodeCoordY, const double nodeCoordZ);
	void addRegularWave(const std::string &waveType, const double height, const double freqORperiod, const double direction, const double phase);
	void addJonswap(const double Hs, const double Tp, const double gamma, const double direction, const double wlow, const double whigh, const int numberOfRegularWaves, const double dwMax);
	void addWaveElevSeries(const std::string &elevFlPath, const double direction, const double wlow, const double whigh);

	void addWindProbe(const unsigned int ID);
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
	bool getFlagWindTurb() const;
	bool getFlagWindUnif() const;

	unsigned int numberOfWaveComponents() const;
	const Wave& getWave(unsigned int waveIndex) const;

	double waveElevAtProbe(const unsigned int ID) const;
	vec::fixed<3> waveVelAtProbe(const unsigned int ID) const;
	double wavePres2ndAtProbe(const unsigned int ID) const;
	
	bool getFlagIFFT() const;	
	bool shouldInterp() const;
	uword getInd4interp1() const;
	uword getInd4interp2() const;
	const vec& getTimeArray() const;
	const vec& getRampArray() const;
	std::string getTurbFileName() const;

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


	void printWaveCharact() const; // Print the wave characteristics (elevation, velocity, etc) specified for output in the locations given by m_waveProbe
	void printWindVelocity() const; // Same thing for wind velocity at m_windProbe

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
	vec::fixed<3> gradP1(const vec::fixed<3> &coord, const double zwl) const;
	vec::fixed<3> gradP1(const vec::fixed<3> &coord, const double zwl, const double time) const;

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
	cx_vec::fixed<3> gradP1_coef(const double x, const double y, const double z, const unsigned int waveIndex) const;

	void windVel(vec::fixed<3> &windVel, const vec::fixed<3> &coord) const;

	// Function that generates time series from complex amplitudes
	mat timeSeriesFromAmp(cx_mat &in, const vec &w) const;
};

// JONSWAP wave spectrum considering frequency in rad/s
double JONSWAP(const double w, const double Tp, const double Hs, const double gamma);
