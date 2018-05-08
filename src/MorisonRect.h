#pragma once

#include <armadillo>
#include "MorisonElement.h"

using namespace arma;

class MorisonRect : public MorisonElement
{
private:
	vec m_cog2node3;

	double m_diam_X;
	double m_CD_X;
	double m_CM_X;

	double m_diam_Y;
	double m_CD_Y;
	double m_CM_Y;

	double m_botArea;
	double m_topArea;
	double m_axialCD;
	double m_axialCM;
	bool m_botPressFlag;

public:
	virtual vec hydrostaticForce(const ENVIR &envir);
	virtual vec morisonForce(const ENVIR &envir);
	virtual vec heavePlateForce(const ENVIR &envir);
};