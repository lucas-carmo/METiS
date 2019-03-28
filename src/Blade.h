#pragma once

#include <armadillo>
#include <vector>

using namespace arma; // For armadillo classes

class Blade
{
private:
	// Properties that belong to each node
    std::vector<double> m_span;
    std::vector<double> m_crvAC; // Local out-of-plane offset (when the blade-pitch angle is zero) of the aerodynamic center (reference point for the airfoil lift and drag forces), normal to the blade-pitch axis, as a result of blade curvature. It is positive downwind. Upwind turbines have negative crvAc for improved tower clearance.
    std::vector<double> m_swpAC; // Specifies the local in-plane offset (when the blade-pitch angle is zero) of the aerodynamic center (reference point for the airfoil lift and drag forces), normal to the blade-pitch axis, as a result of blade sweep; swpAC is positive in the direction opposite to the rotation.
    std::vector<double> m_crvAng; // Specifies the local angle (in degrees) from the blade-pitch axis of a vector normal to the plane of the airfoil, as a result of blade out-of-plane curvature (when the blade-pitch angle is zero). It is positive downwind. Upwind turbines have negative crvAng for improved tower clearance.
    std::vector<double> m_twist; // Specifies the local aerodynamic twist angle (in degrees) of the airfoil; it is the orientation of the local chord about the vector normal to the plane of the airfoil, positive to feather, leading edge upwind; the blade-pitch angle will be added to the local twist; 
    std::vector<double> m_chord; // Local chord length
    std::vector<int> m_airfoilID; // Airfoil data the local blade node is associated with	
	std::vector< double > m_localSolidity;

	// Properties related to the BEMT method AND belong to each node
	// They are part of the Blade class to make the value calculated in the time step 'i' available as
	// first guess to the time step 'i+1'. So, to make it clear, the values stored in the members below
	// are the ones calculated in the last call of RNA::aeroForce, where the BEMT method is used.
	std::vector<double> m_a; // Axial induction factor (for each node)
	std::vector<double> m_ap; // Tangential induction factor (for each node)
	std::vector<double> m_phi; // Local inflow angle (for each node)

	// Radial distance of the node & coordinates of the nodes written in the hub coordinate system.
	// Since this needs to be calculated only in the first time step, I thought it better to keep them as members.
	std::vector<double> m_radius;
	std::vector<vec::fixed<3>> m_nodeCoord_hub;

	// Properties that belong to the whole blade
	double m_precone;
	double m_pitch;
	double m_initialAzimuth; // Azimuth position in t = 0

public:
	Blade();

	/*****************************************************
		Setters
	*****************************************************/
	void addBladeAeroNode(const double span, const double crvAC, const double swpAC, 
						  const double crvAng, const double twist, const double chord, 
						  const int airfoilID, const double hubRadius);

	void setPrecone(const double precone);
	void setPitch(const double pitch);
	void setInitialAzimuth(const double initialAzimuth);

	void setNodeRadius(const int index, const double hubRadius);

	/*****************************************************
		Getters
	*****************************************************/
	unsigned int size() const;
	double span(const unsigned int index) const;
	double crvAC(const unsigned int index) const;
	double swpAC(const unsigned int index) const;
	double crvAng(const unsigned int index) const;
	double twist(const unsigned int index) const;
	double chord(const unsigned int index) const;
	int airoilID(const unsigned int index) const;
	double radius(const unsigned int index) const;
	double axialIndFactor(const unsigned int index) const;
	double tangIndFactor(const unsigned int index) const;
	double phi(const unsigned int index) const; 
	double alpha(const unsigned int index) const;
	double localSolidity(const unsigned int index) const;

	double precone() const;
	double pitch() const;
	double initialAzimuth() const;

	/*****************************************************
		Calculate node position in different coordinate systems
	*****************************************************/	
	void setNodeCoord_hub(const int index);
	vec::fixed<3> nodeCoord_hub(const int index) const;	
	vec::fixed<3> nodeCoord_shaft(const int index, const double dAzimuth) const;	
	vec::fixed<3> nodeCoord_shaft(const vec::fixed<3> &nodeCoord_hub, const double dAzimuth) const;
	vec::fixed<3> nodeCoord_fowt(const vec::fixed<3> &nodeCoord_shaft, const double tilt, const double yaw, const double overhang, const double hubHeight2CoG) const;
	vec::fixed<3> nodeCoord_earth(const vec::fixed<6> &FOWTpos, const vec::fixed<3> &nodeCoord_tower) const;
};

