#include "Airfoil.h"
#include "IO.h"
#include <string>

Airfoil::Airfoil()
{}


/*****************************************************
	Setters
*****************************************************/
void Airfoil::addAirfoilData(double angle, realType CL, realType CD, realType CM)
{
	if (m_angle.size() != 0) // If this is not the first angle that will be added to m_angle
	{
		if (angle <= m_angle.back()) // Then verify if the angle is larger than the previous one, thus garanteeing that m_angle is in ascending order. This is needed to use the spline interpolation implemented in auxFunctions.h
		{
			throw std::runtime_error("Airfoil angles must be organized in ascending order, but tried to input "
                               + std::to_string(angle) + "deg after " + std::to_string(m_angle.back()) + "deg.");
		}
	}

	m_angle.push_back(angle);
	m_CL.push_back(CL);
	m_CD.push_back(CD);
	m_CM.push_back(CM);

	if (m_angle.size() > 2)
	{
		m_spl_CL.set_points(m_angle, m_CL);
		m_spl_CD.set_points(m_angle, m_CD);
		m_spl_CM.set_points(m_angle, m_CM);
	}
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

realType Airfoil::getCL(unsigned int index) const
{
	return m_CL.at(index);
}

realType Airfoil::getCD(unsigned int index) const
{
	return m_CD.at(index);
}

realType Airfoil::getCM(unsigned int index) const
{
	return m_CM.at(index);
}

realType Airfoil::CL(double angle) const
{
	if (size() > 2)
	{
		return m_spl_CL(angle);
	}
	else
	{
		// - Talvez fosse bom dar um warning se for interpolado pra fora dos limites dos dados
		// - Dar um warning e fazer interpolacao linear se for 2.
		// - Retornar erro se for 1 so.
		return 0;
	}
}

realType Airfoil::CD(double angle) const
{
	if (size() > 2)
	{
		return m_spl_CD(angle);
	}
	else
	{
		return 0;
	}
}

realType Airfoil::CM(double angle) const
{
	if (size() > 2)
	{
		return m_spl_CM(angle);
	}
	else
	{
		return 0;
	}
}
