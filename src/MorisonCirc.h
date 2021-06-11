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


public:
	virtual void evaluateQuantitiesAtBegin(const ENVIR &envir, const int hydroMode) override;


	/*****************************************************
		Constructors
	*****************************************************/
	MorisonCirc(const vec &node1Pos, const vec &node2Pos, const vec &cog, const int numIntPoints,
				const bool botPressFlag, const double axialCD_1, const double axialCa_1, const double axialCD_2, const double axialCa_2,
				const double diam, const double CD, double CM);

	/*****************************************************
		Forces acting on the Morison Element and added mass
	*****************************************************/
	virtual void make_local_base(arma::vec::fixed<3> &xvec, arma::vec::fixed<3> &yvec, arma::vec::fixed<3> &zvec, const arma::vec::fixed<3> &n1, const arma::vec::fixed<3> &n2) const override;
	virtual mat::fixed<6, 6> addedMass_perp(const double rho, const vec::fixed<3> &refPt, const int hydroMode) const override;
	virtual double A_perp(const int ii, const int jj, const vec::fixed<3> &x, const vec::fixed<3> &xG, const vec::fixed<3> &xvec, const vec::fixed<3> &yvec) const override;
	virtual double A_paral(const int ii, const int jj, const vec::fixed<3> &x, const vec::fixed<3> &xG, const vec::fixed<3> &zvec) const override;
	virtual mat::fixed<6, 6> addedMass_paral(const double rho, const vec::fixed<3> &refPt, const int hydroMode) const override;	

	// Forces up to second order - to be used in the evaluation of the total acceleration
	// Written in the global coordinate system.
	// Moments are given with respect to node1_sd.
	virtual vec::fixed<6> hydrostaticForce(const double rho, const double g, const vec::fixed<3> &refPt) const override;
	virtual vec::fixed<6> hydrostaticForce_sd(const double rho, const double g, const vec::fixed<3> &refPt) const override;
	virtual vec::fixed<6> hydroForce_1st(const ENVIR &envir, const vec::fixed<3> &refPt) const override;
	virtual vec::fixed<6> hydroForce_drag(const ENVIR &envir, const vec::fixed<3> &refPt) const override;
	virtual vec::fixed<6> hydroForce_relWaveElev(const ENVIR &envir, const vec::fixed<3> &refPt) const override;
	virtual vec::fixed<6> hydroForce_2ndPot(const ENVIR &envir, const vec::fixed<3> &refPt) const override;
	virtual vec::fixed<6> hydroForce_convecAcc(const ENVIR &envir, const vec::fixed<3> &refPt) const override;
	virtual vec::fixed<6> hydroForce_axDiverg(const ENVIR &envir, const vec::fixed<3> &refPt) const override;
	virtual vec::fixed<6> hydroForce_accGradient(const ENVIR &envir, const vec::fixed<3> &refPt) const override;
	virtual vec::fixed<6> hydroForce_slendBodyRot(const ENVIR &envir, const vec::fixed<3> &refPt) const override;
	virtual vec::fixed<6> hydroForce_rem(const ENVIR &envir, const vec::fixed<3> &refPt) const override;
	virtual void quantities4hydrostaticMatrix(double &zb, double &V, double &Awl, double &xwl, double &ywl, double &Ixx, double &Iyy, double &Ixy) const override;

	vec::fixed<6> hydrostaticForce_helper(const double rho, const double g, const vec::fixed<3> &refPt, const vec::fixed<3> &n1, const vec::fixed<3> &n2_in, const vec::fixed<3> &xvec, const vec::fixed<3> &yvec, const vec::fixed<3> &zvec) const;
	cx_vec::fixed<6> hydroForce_1st_coefs(const Wave &wave, double watDensity, double watDepth, double gravity) const;
	cx_vec::fixed<6> hydroForce_2ndPot_coefs(const Wave &wave_ii, const Wave &wave_jj, double watDensity, double watDepth, double gravity) const;

	/*****************************************************
		Printing
	*****************************************************/
	virtual std::string print() const override;

	/*****************************************************
		Clone for creating copies of the Morison Element
	*****************************************************/
	virtual MorisonCirc* clone() const override;
};