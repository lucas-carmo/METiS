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
		Forces acting on the Morison Element and functions for node position/velocity/acceleration)
	*****************************************************/
	virtual vec::fixed<3> node1Pos(const vec::fixed<6> &FOWTpos) const override;
	virtual vec::fixed<3> node2Pos(const vec::fixed<6> &FOWTpos) const override;
	virtual vec::fixed<3> node1Vel(const vec::fixed<6> &floaterPos, const vec::fixed<6> &floaterVel) const override;
	virtual vec::fixed<3> node2Vel(const vec::fixed<6> &floaterPos, const vec::fixed<6> &floaterVel) const override;
	virtual vec::fixed<3> node1Acc(const vec::fixed<6> &floaterPos, const vec::fixed<6> &floaterVel, const vec::fixed<6> &floaterAcc) const override;
	virtual vec::fixed<3> node2Acc(const vec::fixed<6> &floaterPos, const vec::fixed<6> &floaterVel, const vec::fixed<6> &floaterAcc) const override;		

	virtual vec::fixed<6> hydrostaticForce(const ENVIR &envir, const vec::fixed<6> &floaterPos) const override;	
	virtual vec::fixed<6> hydrodynamicForce(const ENVIR &envir, const vec::fixed<6> &floaterPos, const vec::fixed<6> &floaterVel) const override;

	/*****************************************************
		Printing
	*****************************************************/
	virtual std::string print() const override;
};

