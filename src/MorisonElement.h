#pragma once

#include <memory> // For std::unique_ptr
#include <armadillo> // Linear algebra library with usage similar to MATLAB
#include "ENVIR.h"

using namespace arma;

class MorisonElement
{
protected:
	vec::fixed<3> m_cog2node1; // Position of node 1 with respect do the floater CoG (note that this is NOT the CoG of the FOWT), i.e. m_node1Pos(t) - CoG(t). Since the floater is treated as a rigid body, this is constant.  
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

	// Nodes position at the beginning of the simulation, which is useful for
	// calculating the hydrodynamic forces considering only first order terms
	vec::fixed<3> m_node1Pos_t0;
	vec::fixed<3> m_node2Pos_t0;

public:
	MorisonElement(const vec &node1Pos, const vec &node2Pos, const vec &cog, const int numIntPoints, 
				   const bool botPressFlag, const double axialCD, const double axialCa);
	
	// Functions related to position, velocity and acceleration
	void updateNodesPosVelAcc(const vec::fixed<6> &floaterCoGpos, const vec::fixed<6> &floaterVel, const vec::fixed<6> &floaterAcc);	
	vec::fixed<3> node1Pos() const;
	vec::fixed<3> node2Pos() const;
	vec::fixed<3> node1Vel() const;
	vec::fixed<3> node2Vel() const;
	vec::fixed<3> node1Acc() const;
	vec::fixed<3> node2Acc() const;
	vec::fixed<3> node1AccCentrip() const;
	vec::fixed<3> node2AccCentrip() const;
	virtual void make_local_base(arma::vec::fixed<3> &xvec, arma::vec::fixed<3> &yvec, arma::vec::fixed<3> &zvec) const = 0;

	// Contribution to the added mass
	virtual mat::fixed<6, 6> addedMass_perp(const double rho, const int hydroMode) const = 0;
	virtual mat::fixed<6, 6> addedMass_paral(const double rho, const int hydroMode) const = 0;

	// Forces and their auxiliaries
	virtual double findIntersectWL(const ENVIR &envir) const = 0;
	virtual vec::fixed<6> hydrostaticForce(const double rho, const double g, const double z_wl) const = 0;
	virtual vec::fixed<6> hydrodynamicForce(const ENVIR &envir, const int hydroMode, vec::fixed<6> &force_inertia, vec::fixed<6> &force_drag, vec::fixed<6> &force_froudeKrylov) const = 0;

	// Printers and getters
	virtual std::string print() const = 0;

	/*****************************************************
		Clone for creating copies of the Morison Element
	*****************************************************/
	virtual MorisonElement* clone() const = 0;
};

