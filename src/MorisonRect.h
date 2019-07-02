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
	MorisonRect(const vec &node1Pos, const vec &node2Pos, const vec &node3Pos, const vec &cog, const int numIntPoints, 
				const bool botPressFlag, const double axialCD, const double axialCa, const double diam_X, const double CD_X, const double CM_X,
				const double diam_Y, const double CD_Y, const double CM_Y,
				const double botArea, const double topArea);


	/*****************************************************
		Forces acting on the Morison Element and functions for node position/velocity/acceleration)
	*****************************************************/
	virtual void make_local_base(arma::vec::fixed<3> &xvec, arma::vec::fixed<3> &yvec, arma::vec::fixed<3> &zvec) const override;
	virtual mat::fixed<6, 6> addedMass_perp(const double rho) const override;
	virtual mat::fixed<6, 6> addedMass_paral(const double rho) const override;
	virtual vec::fixed<6> hydrostaticForce(const double rho, const double g, const double z_wl) const override;
	virtual vec::fixed<6> hydrodynamicForce(const ENVIR &envir, vec::fixed<6> &force_inertia, vec::fixed<6> &force_drag, vec::fixed<6> &force_froudeKrylov) const override;

	/*****************************************************
		Printing
	*****************************************************/
	virtual std::string print() const override;

	/*****************************************************
		Clone for creating copies of the Morison Element
	*****************************************************/
	virtual MorisonRect* clone() const override;
};