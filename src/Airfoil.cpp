#include "Airfoil.h"

Airfoil::Airfoil()
{}


/*****************************************************
	Setters
*****************************************************/
void Airfoil::addAirfoilLine(double angle, double CL, double CD, double CM)
{
	m_angle.push_back(angle);
	m_CL.push_back(CL);
	m_CD.push_back(CD);
	m_CM.push_back(CM);
}


/*****************************************************
	Getters
*****************************************************/
unsigned int Airfoil::size() const
{
	return static_cast<unsigned int>(m_angle.size());
}

double Airfoil::angle(unsigned int index) const
{
	return m_angle.at(index);
}

double Airfoil::CL(unsigned int index) const
{
	return m_CL.at(index);
}

double Airfoil::CD(unsigned int index) const
{
	return m_CD.at(index);
}

double Airfoil::CM(unsigned int index) const
{
	return m_CM.at(index);
}