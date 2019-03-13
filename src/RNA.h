#pragma once

#include "Blade.h"
#include "Airfoil.h"

#include <vector>

class RNA
{
private:
    double m_rotorSpeed;
    double m_rotorTilt;
    double m_rotorYaw;
    std::vector< Blade > m_blades;
    std::vector< Airfoil > m_airfoils;
    double m_hubRadius;
    double m_hubHeight;
	double m_hubHeight2CoG; // z coordinate of the hub in the FOWT coordinate system. It is equal to the relative height in t=0.
    double m_overhang;

public:
	RNA();

	/*****************************************************
		Setters
	*****************************************************/
	void readRotorSpeed(const std::string &data);
	void readRotorTilt(const std::string &data);
	void readRotorYaw(const std::string &data);

	void readNumBlades(const std::string &data);
	void readBladePrecone(const std::string &data);
	void readBladePitch(const std::string &data);
	void readBladeAeroLine(const std::string &data);
	void addAirfoil();
	void readAirfoilLine(const std::string &data);

	void readHubRadius(const std::string &data);
	void readHubHeight(const std::string &data);
	void readOverhang(const std::string &data);
	void setHubHeight2CoG(const double zCoG);

	/*****************************************************
		Getters
	*****************************************************/
	double rotorSpeed() const;
	double rotorTilt() const;
	double rotorYaw() const;
	double hubRadius() const;
	double hubHeight() const;
	double overhang() const;
	double hubHeight2CoG() const;
	unsigned int numBlades() const;
	double bladePrecone(const unsigned int ii) const;
	double bladePitch(const unsigned int ii) const;

	std::string printBladeAero() const;
	std::string printAirfoils() const;

	/*****************************************************
		Caculation functions
	*****************************************************/	
	double dAzimuth(const double time) const;
	vec::fixed<6> aeroForce(const ENVIR &envir, const vec::fixed<6> &FOWTpos, const vec::fixed<6> &FOWTvel) const;
};