#pragma once

#include <memory> // For std::unique_ptr
#include <armadillo> // Linear algebra library with usage similar to MATLAB
#include "ENVIR.h"

using namespace arma;

class MorisonElement
{
protected:
	vec::fixed<3> m_cog2node1; // Position of node 1 with respect do the floater CoG (note that this is NOT the CoG of the FOWT), i.e. m_node1Pos(t) - CoG(t). Since the floater is treated as a rigid body, this is constant (it is written in the coordinate system attached to the floater).  
	vec::fixed<3> m_cog2node2;
	int m_numIntPoints;

	bool m_botPressFlag;
	double m_axialCD;
	double m_axialCa;		

	// At first, the position, velocity and acceleration of the nodes were not members of the class, and they were calculated directly in their respective getters.
	// However, this procedure involved the calculation of the same parameters several times, so I decided to keep them as members and calculate them with the function
	// updateNodesPosVelAcc(). This function MUST be called every time the FOWT position, velocity and acceleration are changed.
	vec::fixed<3> m_node1Pos;
	vec::fixed<3> m_node2Pos;
	vec::fixed<3> m_node1Vel;
	vec::fixed<3> m_node2Vel;
	vec::fixed<3> m_node1Acc;
	vec::fixed<3> m_node2Acc;
	vec::fixed<3> m_node1AccCentrip;
	vec::fixed<3> m_node2AccCentrip;

	// Nodes position at the beginning of the simulation
	vec::fixed<3> m_node1Pos_t0;
	vec::fixed<3> m_node2Pos_t0;

	// Nodes position and velocity considering only the mean and slow drift of the FOWT.
	// They are used in the calculation of the contribution of quadratic
	// terms (including the quadratic drag) to the hydrodynamic force.
	vec::fixed<3> m_node1Pos_sd;
	vec::fixed<3> m_node2Pos_sd;
	vec::fixed<3> m_node1Vel_sd;
	vec::fixed<3> m_node2Vel_sd;

public:
	MorisonElement(const vec &node1Pos, const vec &node2Pos, const vec &cog, const int numIntPoints, 
				   const bool botPressFlag, const double axialCD, const double axialCa);
	
	// Functions related to position, velocity and acceleration
	void updateNodesPosVelAcc(const vec::fixed<6> &floaterCoGpos, const vec::fixed<6> &floaterVel, const vec::fixed<6> &floaterAcc, const vec::fixed<6> &floaterCoGpos_SD, const vec::fixed<6> &floaterVel_SD);
	vec::fixed<3> node1Pos_t0() const;
	vec::fixed<3> node2Pos_t0() const;
	vec::fixed<3> node1Pos() const;
	vec::fixed<3> node2Pos() const;
	vec::fixed<3> node1Pos_sd() const;
	vec::fixed<3> node2Pos_sd() const;
	vec::fixed<3> node1Vel() const;
	vec::fixed<3> node2Vel() const;
	vec::fixed<3> node1Vel_sd() const;
	vec::fixed<3> node2Vel_sd() const;
	vec::fixed<3> node1Acc() const;
	vec::fixed<3> node2Acc() const;
	vec::fixed<3> node1AccCentrip() const;
	vec::fixed<3> node2AccCentrip() const;
	virtual void make_local_base(arma::vec::fixed<3> &xvec, arma::vec::fixed<3> &yvec, arma::vec::fixed<3> &zvec) const = 0;
	virtual void make_local_base_t0(arma::vec::fixed<3> &xvec, arma::vec::fixed<3> &yvec, arma::vec::fixed<3> &zvec) const = 0;
	virtual void make_local_base_sd(arma::vec::fixed<3> &xvec, arma::vec::fixed<3> &yvec, arma::vec::fixed<3> &zvec) const = 0;

	// Contribution to the added mass
	virtual mat::fixed<6, 6> addedMass_perp(const double rho, const vec::fixed<3> &refPt) const = 0;
	virtual double A_perp(const int ii, const int jj, const vec::fixed<3> &x, const vec::fixed<3> &xG, const vec::fixed<3> &xvec, const vec::fixed<3> &yvec) const = 0;
	virtual mat::fixed<6, 6> addedMass_paral(const double rho, const vec::fixed<3> &refPt) const = 0;

	virtual vec::fixed<6> hydrostaticForce(const double rho, const double g) const = 0;
	virtual vec::fixed<6> hydrodynamicForce(const ENVIR &envir, const int hydroMode, const mat::fixed<3, 3> &rotat, const vec::fixed<3> &refPt, const vec::fixed<3> &refPt_sd,
											vec::fixed<6> &force_inertia, vec::fixed<6> &force_drag, vec::fixed<6> &force_froudeKrylov,
											vec::fixed<6> &force_inertia_2nd_part1, vec::fixed<6> &force_inertia_2nd_part2, 
											vec::fixed<6> &force_inertia_2nd_part3, vec::fixed<6> &force_inertia_2nd_part4, 
											vec::fixed<6> &force_inertia_2nd_part5) const = 0;

	// Printers and getters
	virtual std::string print() const = 0;

	// Others
	virtual vec::fixed<3> findIntersectWL(const ENVIR &envir) const;

	/*****************************************************
		Clone for creating copies of the Morison Element
	*****************************************************/
	virtual MorisonElement* clone() const = 0;
};

