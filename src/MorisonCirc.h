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
	virtual void setPropertiesWithIFFT(const ENVIR &envir) override;


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
	virtual vec::fixed<6> hydrostaticForce(const double rho, const double g) const override;
	virtual vec::fixed<6> hydrodynamicForce(const ENVIR &envir, const int hydroMode, const vec::fixed<3> &refPt, const vec::fixed<3> &refPt_sd,
											vec::fixed<6> &force_drag, vec::fixed<6> &force_1, vec::fixed<6> &force_2,
											vec::fixed<6> &force_3, vec::fixed<6> &force_4, vec::fixed<6> &force_eta, vec::fixed<6> &force_rem) const override;

	// Analytical evaluation of the integration along the cylinder's length
	// of the term due to fluid acceleration in Morison's Equation.
	//
	// Written in the global coordinate system.
	// Moments are given with respect to node1.
	// TODO: Once things are finished, put these functions in MorisonElement
	vec::fixed<6> hydroForce_1st(const ENVIR &envir, const int hydroMode) const;
	cx_vec::fixed<6> hydroForce_1st_components(const Wave &wave, double watDensity, double watDepth, double gravity) const;
	vec::fixed<6> hydroForce_drag(const ENVIR &envir) const;
	vec::fixed<6> hydroForce_drag_fromIFFT(const ENVIR &envir) const;
	vec::fixed<6> hydroForce_drag_withoutIFFT(const ENVIR &envir) const;
	

	vec::fixed<6> morisonForce_inertia2nd(const ENVIR &envir) const;

	/*****************************************************
		Printing
	*****************************************************/
	virtual std::string print() const override;

	/*****************************************************
		Clone for creating copies of the Morison Element
	*****************************************************/
	virtual MorisonCirc* clone() const override;
};