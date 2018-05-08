#pragma once

#include <armadillo>
#include "MorisonElement.h"

using namespace arma;

class MorisonCyl : public MorisonElement
{
private:
	double m_diam;
	double m_CD;
	double m_CM;

	double m_botDiam;
	double m_topDiam;
	double m_axialCD;
	double m_axialCM;
	bool m_botPressFlag;

public:
	virtual vec hydrostaticForce(const ENVIR &envir);
	virtual vec morisonForce(const ENVIR &envir);
	virtual vec heavePlateForce(const ENVIR &envir);	
};