#include "Blade.h"

using namespace arma;

/*****************************************************
	Setters
*****************************************************/
// Blade(const double precone, const double azimuth, const double pitch)
// {
// 	m_precone = precone;
// 	m_azimuth = azimuth;
// 	m_pitch = pitch;
// Como eles são lidos em linhas diferentes, acho melhor deixar o constructor com NaN e 
// fazer uma função set pra cada um
// }


void Blade::addBladeAeroLine(const double span, const double crvAC, const double swpAC, 
							 const double crvAng, const double twist, const double chord, 
							 const int airfoilID)
{
	m_span.push_back(span);
	m_crvAC.push_back(crvAC);
	m_swpAC.push_back(swpAC);
	m_crvAng.push_back(crvAng);
	m_twist.push_back(twist);
	m_chord.push_back(chord);
	m_airfoilID.push_back(airfoilID);
}

void Blade::setPrecone(const double precone)
{
	m_precone = precone;
}

void Blade::setPitch(const double pitch)
{
	m_pitch = pitch;
}

void Blade::setInitialAzimuth(const double initialAzimuth)
{
	m_initialAzimuth = initialAzimuth;
}


/*****************************************************
	Getters
*****************************************************/
unsigned int Blade::size() const
{
	return static_cast<unsigned int>(m_airfoilID.size());
}

double Blade::span(const unsigned int index) const
{
	return m_span.at(index);
}

double Blade::crvAC(const unsigned int index) const
{
	return m_crvAC.at(index);
}

double Blade::swpAC(const unsigned int index) const
{
	return m_swpAC.at(index);
}

double Blade::crvAng(const unsigned int index) const
{
	return m_crvAng.at(index);
}

double Blade::twist(const unsigned int index) const
{
	return m_twist.at(index);
}

double Blade::chord(const unsigned int index) const
{
	return m_chord.at(index);
}

int Blade::airoilID(const unsigned int index) const
{
	return m_airfoilID.at(index);
}

double Blade::precone() const
{
	return m_precone;
}

double Blade::pitch() const
{
	return m_pitch;
}

double Blade::initialAzimuth() const
{
	return m_initialAzimuth;
}


/*****************************************************
	Calculate node position in different coordinate systems
*****************************************************/	
vec::fixed<3> Blade::hubCoord(unsigned int index, double hubRadius) const
{
	vec::fixed<3> hubCoord;
	double r = hubRadius + span(index);		

	hubCoord[0] = r * tan(precone());
	hubCoord[1] = -r * sin(initialAzimuth()) * cos(precone());
	hubCoord[2] = r * cos(initialAzimuth()) * cos(precone());
}