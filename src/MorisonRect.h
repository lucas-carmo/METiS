#pragma once

#include <armadillo>
#include "MorisonElement.h"

using namespace arma;

class MorisonRect : public MorisonElement
{
private:
	vec::fixed<3> m_cog2node3;

	double m_diam_X;
	double m_CD_X;
	double m_CM_X;

	double m_diam_Y;
	double m_CD_Y;
	double m_CM_Y;

	double m_botArea;
	double m_topArea;

public:
	/*****************************************************
		Constructors
	*****************************************************/
	MorisonRect(vec cog2node1, vec cog2node2, vec cog2node3, int numIntPoints, bool botPressFlag,
				double axialCD, double axialCa, double diam_X, double CD_X, double CM_X,
				double diam_Y, double CD_Y, double CM_Y,
				double botArea, double topArea);


	/*****************************************************
		Forces acting on the Morison Element and functions for node position/velocity/acceleration)
	*****************************************************/
	virtual vec::fixed<6> hydrostaticForce(const ENVIR &envir) const override;	
	virtual vec::fixed<6> hydrodynamicForce(const ENVIR &envir) const override;

	/*****************************************************
		Printing
	*****************************************************/
	virtual std::string print() const override;

	/*****************************************************
		Clone for creating copies of the Morison Element
	*****************************************************/
	virtual MorisonRect* clone() const override;
};