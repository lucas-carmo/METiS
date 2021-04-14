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
	double m_axialCD_1;
	double m_axialCa_1;
	double m_axialCD_2;
	double m_axialCa_2;

	// At first, the position, velocity and acceleration of the nodes were not members of the class, and they were calculated directly in their respective getters.
	// However, this procedure involved the calculation of the same parameters several times, so I decided to keep them as members and calculate them with the function
	// updateNodesPosVelAcc(). This function MUST be called every time the FOWT position, velocity and acceleration are changed.
	vec::fixed<3> m_node1Pos;
	vec::fixed<3> m_node2Pos;
	vec::fixed<3> m_node1Vel;
	vec::fixed<3> m_node2Vel;
	vec::fixed<3> m_node1AccCentrip;
	vec::fixed<3> m_node2AccCentrip;

	// Nodes position at the beginning of the simulation
	vec::fixed<3> m_node1Pos_t0;
	vec::fixed<3> m_node2Pos_t0;

	// Nodes position and velocity considering only the mean and slow drift of the FOWT.
	vec::fixed<3> m_node1Pos_sd;
	vec::fixed<3> m_node2Pos_sd;
	vec::fixed<3> m_node1Vel_sd;
	vec::fixed<3> m_node2Vel_sd;
	
	// Vectors of the local base
	// They could be calculated every time from the nodes' position,
	// but as members they can be calculated only once when the cylinder is updated
	vec::fixed<3> m_xvec_t0; // Considering the initial position of the cylinder
	vec::fixed<3> m_yvec_t0;
	vec::fixed<3> m_zvec_t0;
	vec::fixed<3> m_xvec_sd; // Considering the slowly varying position of the cylinder
	vec::fixed<3> m_yvec_sd;
	vec::fixed<3> m_zvec_sd;
	vec::fixed<3> m_xvec; // Considering the instantaneous total position of the cylinder
	vec::fixed<3> m_yvec;
	vec::fixed<3> m_zvec;

	// Intersection with the instantaneous waterline
	vec::fixed<3> m_intersectWL;

	// Need the Z coordinate of the point that intersects the mean waterline at t=0
	// for the evaluation of the second-order wave forces considering the fixed body position
	vec::fixed<3> m_nodeWL;
	vec::fixed<3> m_cog2nodeWL;
	double m_Zwl{ 0 };

	// Rotation matrix from the body coordinate system to the global coordinate system
	mat::fixed<3, 3> m_RotatMatrix;
	mat::fixed<3, 3> m_RotatMatrix_sd;
	
	// Quantities calculated at the beginning of the simulaion using IFFT
	bool m_flagFixed{ false }; // Flag used to specify whether things will be evaluated about the mean position of the cylinder
	mat m_hydroForce_1st_Array;


public:
	MorisonElement(const vec &node1Pos, const vec &node2Pos, const vec &cog, const int numIntPoints, 
				   const bool botPressFlag, const double axialCD_1, const double axialCa_1, const double axialCD_2, const double axialCa_2);
	
	virtual void setPropertiesWithIFFT(const ENVIR &envir) = 0;

	// Functions related to position, velocity and acceleration
	void calcPosVel(const vec::fixed<6> &pos, const vec::fixed<6> &vel,
		            vec::fixed<3> &node1Pos, vec::fixed<3> &node2Pos, vec::fixed<3> &node1Vel, vec::fixed<3> &node2Vel,
		            vec::fixed<3> &xvec, vec::fixed<3> &yvec, vec::fixed<3> &zvec) ;	
	void updateMorisonElement(const ENVIR &envir, const vec::fixed<6> &floaterCoGpos, const vec::fixed<6> &floaterVel, 
							  const vec::fixed<6> &floaterCoGpos_SD, const vec::fixed<6> &floaterVel_SD);
	void updateMorisonElement(const vec::fixed<6> &floaterCoGpos, const vec::fixed<6> &floaterVel,
							  const vec::fixed<6> &floaterCoGpos_SD, const vec::fixed<6> &floaterVel_SD); // Version that does not find the intersection with the waterline

	vec::fixed<3> node1Pos_t0() const;
	vec::fixed<3> node2Pos_t0() const;
	vec::fixed<3> node1Pos_sd() const;
	vec::fixed<3> node2Pos_sd() const;
	vec::fixed<3> node1Pos() const;
	vec::fixed<3> node2Pos() const;

	vec::fixed<3> node1Vel_sd() const;
	vec::fixed<3> node2Vel_sd() const;
	vec::fixed<3> node1Vel() const;
	vec::fixed<3> node2Vel() const;	

	vec::fixed<3> node1AccCentrip() const;
	vec::fixed<3> node2AccCentrip() const;

	virtual void make_local_base(arma::vec::fixed<3> &xvec, arma::vec::fixed<3> &yvec, arma::vec::fixed<3> &zvec, const arma::vec::fixed<3> &n1, const arma::vec::fixed<3> &n2) const = 0;

	// Contribution to the added mass
	virtual mat::fixed<6, 6> addedMass_perp(const double rho, const vec::fixed<3> &refPt, const int hydroMode) const = 0;
	virtual double A_perp(const int ii, const int jj, const vec::fixed<3> &x, const vec::fixed<3> &xG, const vec::fixed<3> &xvec, const vec::fixed<3> &yvec) const = 0;
	virtual mat::fixed<6, 6> addedMass_paral(const double rho, const vec::fixed<3> &refPt, const int hydroMode) const = 0;
	virtual double A_paral(const int ii, const int jj, const vec::fixed<3> &x, const vec::fixed<3> &xG, const vec::fixed<3> &zvec) const = 0;

	// Forces up to second order - to be used in the evaluation of the total acceleration
	virtual vec::fixed<6> hydrostaticForce(const double rho, const double g) const = 0;
	virtual vec::fixed<6> hydrodynamicForce(const ENVIR &envir, const int hydroMode, const vec::fixed<3> &refPt, const vec::fixed<3> &refPt_sd,
											vec::fixed<6> &force_drag, vec::fixed<6> &force_1, vec::fixed<6> &force_2, 
											vec::fixed<6> &force_3, vec::fixed<6> &force_4, vec::fixed<6> &force_eta, vec::fixed<6> &force_rem) const = 0;
	

	// Printers and getters
	virtual std::string print() const = 0;

	// Others
	virtual vec::fixed<3> findIntersectWL(const ENVIR &envir) const;

	/*****************************************************
		Clone for creating copies of the Morison Element
	*****************************************************/
	virtual MorisonElement* clone() const = 0;
};

