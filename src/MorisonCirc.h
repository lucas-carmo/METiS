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
	void make_local_base(arma::vec::fixed<3> &xvec, arma::vec::fixed<3> &yvec, arma::vec::fixed<3> &zvec) const;
	virtual vec::fixed<6> hydrostaticForce(const ENVIR &envir) const override;	
	virtual vec::fixed<6> hydrodynamicForce(const ENVIR &envi) const override;

	/*****************************************************
		Printing
	*****************************************************/
	virtual std::string print() const override;

	/*****************************************************
		Clone for creating copies of the Morison Element
	*****************************************************/
	virtual MorisonCirc* clone() const override;
};