#include "MorisonElement.h"

/*****************************************************
	Constructors
*****************************************************/

MorisonElement::MorisonElement(vec cog2node1, vec cog2node2, int numIntPoints, bool botPressFlag)
	: m_numIntPoints(numIntPoints), m_cog2node1(cog2node1), m_cog2node2(cog2node2), m_botPressFlag(botPressFlag)
{}
