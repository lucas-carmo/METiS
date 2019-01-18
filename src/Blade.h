#pragma once

#include <armadillo>
#include <vector>

using namespace arma; // For armadillo classes

class Blade
{
private:
    std::vector<double> m_span;
    std::vector<double> m_crvAC; // Local out-of-plane offset (when the blade-pitch angle is zero) of the aerodynamic center (reference point for the airfoil lift and drag forces), normal to the blade-pitch axis, as a result of blade curvature. It is positive downwind. Upwind turbines have negative crvAc for improved tower clearance.
    std::vector<double> m_swpAC; // Specifies the local in-plane offset (when the blade-pitch angle is zero) of the aerodynamic center (reference point for the airfoil lift and drag forces), normal to the blade-pitch axis, as a result of blade sweep; swpAC is positive in the direction opposite to the rotation.
    std::vector<double> m_crvAng; // Specifies the local angle (in degrees) from the blade-pitch axis of a vector normal to the plane of the airfoil, as a result of blade out-of-plane curvature (when the blade-pitch angle is zero). It is positive downwind. Upwind turbines have negative crvAng for improved tower clearance.
    std::vector<double> m_twist; // Specifies the local aerodynamic twist angle (in degrees) of the airfoil; it is the orientation of the local chord about the vector normal to the plane of the airfoil, positive to feather, leading edge upwind; the blade-pitch angle will be added to the local twist; 
    std::vector<double> m_chord; // Local chord length
    std::vector<int> m_airfoilID; // Airfoil data the local blade node is associated with

public:
	/*****************************************************
		Setters
	*****************************************************/
	void addBladeAeroLine(double span, double crvAC, double swpAC, double crvAng, double twist, double chord, int airfoilID);

	/*****************************************************
		Getters
	*****************************************************/
	unsigned int size() const;
	double span(unsigned int index) const;
	double crvAC(unsigned int index) const;
	double swpAC(unsigned int index) const;
	double crvAng(unsigned int index) const;
	double twist(unsigned int index) const;
	double chord(unsigned int index) const;
	int airoilID(unsigned int index) const;
};

