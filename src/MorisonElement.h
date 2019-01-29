#pragma once

#include <memory> // For std::unique_ptr
#include <armadillo> // Linear algebra library with usage similar to MATLAB
#include "ENVIR.h"

using namespace arma;

class MorisonElement
{
protected:
	vec::fixed<3> m_cog2node1;
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

public:
	MorisonElement(vec cog2node1, vec cog2node2, int numIntPoints, 
				   bool botPressFlag, double axialCD, double axialCa);
	
	// Functions related to position, velocity and acceleration
	void updateNodesPosVelAcc(const vec::fixed<6> &floaterCoGpos, const vec::fixed<6> &floaterVel, const vec::fixed<6> &floaterAcc);
	mat::fixed<3, 3> RotatMatrix(const vec::fixed<3> &rotation) const;	
	vec::fixed<3> node1Pos() const;
	vec::fixed<3> node2Pos() const;
	vec::fixed<3> node1Vel() const;
	vec::fixed<3> node2Vel() const;
	vec::fixed<3> node1Acc() const;
	vec::fixed<3> node2Acc() const;

	// Forces
	virtual vec::fixed<6> hydrostaticForce(const ENVIR &envir) const = 0;
	
	virtual vec::fixed<6> hydrodynamicForce(const ENVIR &envir, vec::fixed<6> &force_inertia, vec::fixed<6> &force_drag, vec::fixed<6> &force_froudeKrylov) const = 0; 

	// Printers and getters
	virtual std::string print() const = 0;

	/*****************************************************
		Clone for creating copies of the Morison Element
	*****************************************************/
	virtual MorisonElement* clone() const = 0;
};

