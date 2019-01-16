#include "IO.h"
#include "RNA.h"

/*****************************************************
	Setters
*****************************************************/
void RNA::readRotorSpeed(const std::string &data)
{
	readDataFromString(data, m_rotorSpeed);
}

void RNA::readRotorTilt(const std::string &data)
{
	readDataFromString(data, m_rotorTilt);
}

void RNA::readRotorYaw(const std::string &data)
{
	readDataFromString(data, m_rotorYaw);
}

void RNA::readNumBlades(const std::string &data)
{
	readDataFromString(data, m_numBlades);
}

void RNA::readBladePrecone(const std::string &data)
{
	readDataFromString(data, m_bladePrecone);
}

void RNA::readBladePitch(const std::string &data)
{
	readDataFromString(data, m_bladePitch);
}

void RNA::readBladeLine(const std::string &data)
{
	double span = 0;
	double crvAC = 0;
	double swpAC = 0;
	double crvAng = 0;
	double twist = 0;
	double chord = 0;
	int airfoilID = 0;

	// The seven properties provided by each line are separated by white spaces in the input string (whitespace or tab)
	std::vector<std::string> input = stringTokenize(data, " \t");

	// Check number of inputs
	if (input.size() != 7)
	{
		throw std::runtime_error("Unable to read the blade aerodynamic properties in input line " + std::to_string(IO::getInLineNumber()) + ". Wrong number of parameters.");
	}

	// Read data
	readDataFromString(input.at(0), span);
	readDataFromString(input.at(1), crvAC);
	readDataFromString(input.at(2), swpAC);
	readDataFromString(input.at(3), crvAng);
	readDataFromString(input.at(4), twist);
	readDataFromString(input.at(5), chord);
	readDataFromString(input.at(6), airfoilID);

	// Add the data read in the current line to the blade element
	RNA.m_blade.addBladeLine(span, crvAC, swpAC, crvAng, twist, chord, airfoilID);
}

void RNA::addAirfoil()
{
	m_airfoils.push_back(Airfoil());
}

void RNA::readAirfoilLine(const std::string &data)
{

}

void RNA::readHubRadius(const std::string &data)
{
	readDataFromString(data, m_hubRadius);
}

void RNA::readHubHeight(const std::string &data)
{
	readDataFromString(data, m_hubHeight);
}

void RNA::readOverhang(const std::string &data)
{
	readDataFromString(data, m_overhang);
}


/*****************************************************
	Getters
*****************************************************/
double RNA::rotorSpeed() const
{
	return m_rotorSpeed;
}

double RNA::rotorTilt() const
{
	return m_rotorTilt;
}

double RNA::rotorYaw() const
{
	return m_rotorYaw;
}

int RNA::numBlades() const
{
	return m_numBlades;
}

double RNA::bladePrecone() const
{
	return m_bladePrecone;
}

double RNA::bladePitch() const
{
	return m_bladePitch;
}

/*std::string printBlade() const;
std::string printAirofoils() const;*/

double RNA::hubRadius() const
{
	return m_hubRadius;
}

double RNA::hubHeight() const
{
	return m_hubHeight;
}

double RNA::overhang() const
{
	return m_overhang;
}