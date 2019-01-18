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

    int m_numBlades;
	double m_bladePrecone;
	double m_bladePitch;
    Blade m_blade;
    std::vector< Airfoil > m_airfoils;

    double m_hubRadius;
    double m_hubHeight;
    double m_overhang;

public:
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

	/*****************************************************
		Getters
	*****************************************************/
	double rotorSpeed() const;
	double rotorTilt() const;
	double rotorYaw() const;

	int numBlades() const;
	double bladePrecone() const;
	double bladePitch() const;

	std::string printBladeAero() const;
	std::string printAirfoils() const;

	double hubRadius() const;
	double hubHeight() const;
	double overhang() const;
};

