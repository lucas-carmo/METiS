#include "Blade.h"

void Blade::addBladeAeroLine(double span, double crvAC, double swpAC, double crvAng, double twist, double chord, int airfoilID)
{
	m_span.push_back(span);
	m_crvAC.push_back(crvAC);
	m_swpAC.push_back(swpAC);
	m_crvAng.push_back(crvAng);
	m_twist.push_back(twist);
	m_chord.push_back(chord);
	m_airfoilID.push_back(airfoilID);
}


unsigned int Blade::size() const
{
	return static_cast<unsigned int>(m_airfoilID.size());
}

double Blade::span(unsigned int index) const
{
	return m_span.at(index);
}

double Blade::crvAC(unsigned int index) const
{
	return m_crvAC.at(index);
}

double Blade::swpAC(unsigned int index) const
{
	return m_swpAC.at(index);
}

double Blade::crvAng(unsigned int index) const
{
	return m_crvAng.at(index);
}

double Blade::twist(unsigned int index) const
{
	return m_twist.at(index);
}

double Blade::chord(unsigned int index) const
{
	return m_chord.at(index);
}

int Blade::airoilID(unsigned int index) const
{
	return m_airfoilID.at(index);
}