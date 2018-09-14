#pragma once

#include <armadillo>
#include <string>
#include "MorisonElement.h"


using namespace arma;

class MorisonCirc : public MorisonElement
{
private:
	double m_diam;
	double m_CD;
	double m_CM;

	double m_botDiam;
	double m_topDiam;

public:
	/*****************************************************
		Constructors
	*****************************************************/
	MorisonCirc(vec cog2node1, vec cog2node2, int numIntPoints, bool botPressFlag,
				double axialCD, double axialCa, double diam, double CD, double CM,
				double botDiam, double topDiam);

	/*****************************************************
		Forces acting on the Morison Element
	*****************************************************/
	virtual vec::fixed<6> hydrostaticForce(const ENVIR &envir) override;
	virtual vec::fixed<6> hydrodynamicForce(const ENVIR &envir) override;

	/*****************************************************
		Printing
	*****************************************************/
	virtual std::string print() const override;
};

