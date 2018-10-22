#include "MorisonElement.h"

/*****************************************************
	Constructors
*****************************************************/

MorisonElement::MorisonElement(vec cog2node1, vec cog2node2, int numIntPoints, 
							   bool botPressFlag, double axialCD, double axialCa)
	: m_cog2node1(cog2node1), m_cog2node2(cog2node2), 
	  m_botPressFlag(botPressFlag), m_axialCD(axialCD), m_axialCa(axialCa)
{
	// Since Simpson's rule is employed for the integration of the forces along the 
	// Morison's element, we need to make sure that the number of integration points is odd
	if (numIntPoints % 2)
	{
		m_numIntPoints = numIntPoints + 1;
	}
	else
	{
		m_numIntPoints = numIntPoints;
	}
}

