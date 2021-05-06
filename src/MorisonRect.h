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

public:
	/*****************************************************
		Constructors
	*****************************************************/
	MorisonRect(const vec &node1Pos, const vec &node2Pos, const vec &node3Pos, const vec &cog, const int numIntPoints,
		const bool botPressFlag, const double axialCD_1, const double axialCa_1, const double axialCD_2, const double axialCa_2,
		const double diam_X, const double CD_X, const double CM_X,
		const double diam_Y, const double CD_Y, const double CM_Y);

	virtual void evaluateQuantitiesAtBegin(const ENVIR &envir, const int hydroMode) override;

	/*****************************************************
		Forces acting on the Morison Element and functions for node position/velocity/acceleration)
	*****************************************************/
	virtual void make_local_base(arma::vec::fixed<3> &xvec, arma::vec::fixed<3> &yvec, arma::vec::fixed<3> &zvec, const arma::vec::fixed<3> &n1, const arma::vec::fixed<3> &n2) const override;
	virtual mat::fixed<6, 6> addedMass_perp(const double rho, const vec::fixed<3> &refPt, const int hydroMode) const override;
	virtual double A_perp(const int ii, const int jj, const vec::fixed<3> &x, const vec::fixed<3> &xG, const vec::fixed<3> &xvec, const vec::fixed<3> &yvec) const override;
	virtual mat::fixed<6, 6> addedMass_paral(const double rho, const vec::fixed<3> &refPt, const int hydroMode) const override;
	virtual double A_paral(const int ii, const int jj, const vec::fixed<3> &x, const vec::fixed<3> &xG, const vec::fixed<3> &zvec) const override;


	// Forces up to second order - to be used in the evaluation of the total acceleration
	virtual vec::fixed<6> hydrostaticForce(const double rho, const double g) const override;
	virtual vec::fixed<6> hydrostaticForce_sd(const double rho, const double g) const override;
	virtual vec::fixed<6> hydroForce_1st(const ENVIR &envir, const int hydroMode, const vec::fixed<3> &refPt) const override;
	virtual vec::fixed<6> hydroForce_drag(const ENVIR &envir, const vec::fixed<3> &refPt) const override;
	virtual vec::fixed<6> hydroForce_relWaveElev(const ENVIR &envir, const vec::fixed<3> &refPt) const override;
	virtual vec::fixed<6> hydroForce_2ndPot(const ENVIR &envir, const vec::fixed<3> &refPt) const override;
	virtual vec::fixed<6> hydroForce_convecAcc(const ENVIR &envir, const vec::fixed<3> &refPt) const override;
	virtual vec::fixed<6> hydroForce_axDiverg(const ENVIR &envir, const vec::fixed<3> &refPt) const override;
	virtual vec::fixed<6> hydroForce_accGradient(const ENVIR &envir, const vec::fixed<3> &refPt) const override;
	virtual vec::fixed<6> hydroForce_slendBodyRot(const ENVIR &envir, const vec::fixed<3> &refPt) const override;
	virtual vec::fixed<6> hydroForce_rem(const ENVIR &envir, const vec::fixed<3> &refPt) const override;

	/*****************************************************
		Printing
	*****************************************************/
	virtual std::string print() const override;

	/*****************************************************
		Clone for creating copies of the Morison Element
	*****************************************************/
	virtual MorisonRect* clone() const override;
};