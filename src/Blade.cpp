#include "Blade.h"

void Blade::addBladeLine(double span, double crvAC, double swpAC, double crvAng, double twist, double chord, int airfoilID)
{
	m_span.push_back(span);
	m_crvAC.push_back(crvAC);
	m_swpAC.push_back(swpAC);
	m_crvAng.push_back(crvAng);
	m_twist.push_back(twist);
	m_chord.push_back(chord);
	m_airfoilID.push_back(airfoilID);
}