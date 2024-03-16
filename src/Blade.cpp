#include "Blade.h"
#include "auxFunctions.h" // For rotatMatrix()

using namespace arma;


Blade::Blade()
{
	m_precone = datum::nan;
	m_pitch = datum::nan;
	m_initialAzimuth = datum::nan;
}

/*****************************************************
	Setters
*****************************************************/
void Blade::addBladeAeroNode(const double span, const double crvAC, const double swpAC,
							 const double crvAng, const double twist, const double chord,
							 const int airfoilID, const double hubRadius)
{
	m_span.push_back(span);
	m_crvAC.push_back(crvAC);
	m_swpAC.push_back(swpAC);
	m_crvAng.push_back(crvAng);
	m_twist.push_back(twist);
	m_chord.push_back(chord);
	m_airfoilID.push_back(airfoilID);

	// Set node radial distance to the hub and its coordinates in the hub coordinate system
	m_radius.push_back(0);
	setNodeRadius(size() - 1, hubRadius); // update node radial distance

	m_nodeCoord_hub.push_back(vec::fixed<3> {0, 0, 0});
	setNodeCoord_hub(size() - 1); // update node position of the current blade node

	// Need the blades precone angle for calculating local solidity
	if (!is_finite(m_precone))
	{
		throw std::runtime_error("Need to specify blade precone before specifying blade aerodynamic properties.");
	}
	m_localSolidity.push_back(chord / (2 * datum::pi * m_radius.back() * cos(deg2rad(m_precone)) ));
}

void Blade::setPrecone(const double precone)
{
	m_precone = precone;
}

void Blade::setPitch(const double pitch)
{
	m_pitch = pitch;
}

void Blade::setInitialAzimuth(const double initialAzimuth)
{
	m_initialAzimuth = initialAzimuth;
}


/*****************************************************
	Getters
*****************************************************/
unsigned int Blade::size() const
{
	return static_cast<unsigned int>(m_airfoilID.size());
}

double Blade::span(const unsigned int index) const
{
	if (index < 0 || static_cast<unsigned int> (index) >= m_span.size())
	{
		throw std::runtime_error("Index out of range in Blade::span(const unsigned int index).");
	}

	return m_span[index];
}

double Blade::crvAC(const unsigned int index) const
{
	if (index < 0 || static_cast<unsigned int> (index) >= m_crvAC.size())
	{
		throw std::runtime_error("Index out of range in Blade::crvAC(const unsigned int index).");
	}

	return m_crvAC[index];
}

double Blade::swpAC(const unsigned int index) const
{
	if (index < 0 || static_cast<unsigned int> (index) >= m_swpAC.size())
	{
		throw std::runtime_error("Index out of range in Blade::swpAC(const unsigned int index).");
	}

	return m_swpAC[index];
}

double Blade::crvAng(const unsigned int index) const
{
	if (index < 0 || static_cast<unsigned int> (index) >= m_crvAng.size())
	{
		throw std::runtime_error("Index out of range in Blade::crvAng(const unsigned int index).");
	}

	return m_crvAng[index];
}

double Blade::twist(const unsigned int index) const
{
	if (index < 0 || static_cast<unsigned int> (index) >= m_twist.size())
	{
		throw std::runtime_error("Index out of range in Blade::twist(const unsigned int index).");
	}

	return m_twist[index];
}

double Blade::chord(const unsigned int index) const
{
	if (index < 0 || static_cast<unsigned int> (index) >= m_chord.size())
	{
		throw std::runtime_error("Index out of range in Blade::chord(const unsigned int index).");
	}

	return m_chord[index];
}

int Blade::airoilID(const unsigned int index) const
{
	if (index < 0 || static_cast<unsigned int> (index) >= m_airfoilID.size())
	{
		throw std::runtime_error("Index out of range in Blade::airfoilID(const unsigned int index).");
	}

	return m_airfoilID[index];
}

double Blade::radius(const unsigned int index) const
{
	if (index < 0 || static_cast<unsigned int> (index) >= m_radius.size())
	{
		throw std::runtime_error("Index out of range in Blade::radius(const unsigned int index).");
	}

	return m_radius[index];
}


double Blade::localSolidity(const unsigned int index) const
{
	if (index < 0 || static_cast<unsigned int> (index) >= m_localSolidity.size())
	{
		throw std::runtime_error("Index out of range in Blade::localSolidity(const unsigned int index).");
	}

	return m_localSolidity[index];
}

double Blade::precone() const
{
	return m_precone;
}

double Blade::pitch() const
{
	return m_pitch;
}

double Blade::initialAzimuth() const
{
	return m_initialAzimuth;
}

/*****************************************************
	Calculate node position in different coordinate systems
*****************************************************/
void Blade::setNodeRadius(const int index, const double hubRadius)
{
	if (index < 0 || static_cast<unsigned int> (index) >= m_radius.size())
	{
		throw std::runtime_error("Index out of range in Blade::setNodeCord_hub(const unsigned int index, const double hubRadius).");
	}

	m_radius[index] = hubRadius + span(index);
}


// Coordinates of a blade node written in the hub coordinate system.
// It requires the index of the node you are interested in and the hub radius.
// After the calculation, the value is stored in m_nodeCoord_hub(index) for future usage.
void Blade::setNodeCoord_hub(const int index)
{
	if (index < 0 || static_cast<unsigned int> (index) >= m_nodeCoord_hub.size())
	{
		throw std::runtime_error("Index out of range in Blade::setNodeCord_hub(const unsigned int index, const double hubRadius).");
	}

	vec::fixed<3> hubCoord;
	double r = radius(index);

	hubCoord[0] = r * sin(deg2rad(m_precone));
	hubCoord[1] = -r * sin(deg2rad(m_initialAzimuth)) * cos(deg2rad(m_precone));
	hubCoord[2] = r * cos(deg2rad(m_initialAzimuth)) * cos(deg2rad(m_precone));

	m_nodeCoord_hub.at(index) = hubCoord;
}

// Coordinates of a blade node written in the hub coordinate system.
vec::fixed<3> Blade::nodeCoord_hub(const int index) const
{
	if (index < 0 || static_cast<unsigned int> (index) >= size())
	{
		throw std::runtime_error("Index out of range in Blade::nodeCoord_hub(const unsigned int index).");
	}

	return m_nodeCoord_hub.at(index);
}

// Coordinates of a blade node written in the shaft coordinate system.
// dAzimuth must be given in degrees
vec::fixed<3> Blade::nodeCoord_shaft(const int index, const double dAzimuth) const
{
	// Do not need to check if index is out of range, since nodeCoord_hub already does that
	return nodeCoord_shaft(nodeCoord_hub(index), dAzimuth);
}

// Overload for when the input is the nodeCoord_hub of a certain node itself.
// dAzimuth must be given in degrees
vec::fixed<3> Blade::nodeCoord_shaft(const vec::fixed<3> &nodeCoord_hub, const double dAzimuth) const
{
	return ( rotatMatrix_deg(dAzimuth, 0, 0) * nodeCoord_hub );
}

// Coordinates of a blade node written in the fowt coordinate system.
// tilt and yaw must be given in degrees
vec::fixed<3> Blade::nodeCoord_fowt(const vec::fixed<3> &nodeCoord_shaft, const double tilt, const double yaw, const double overhang, const double hubHeight2CoG) const
{
	mat::fixed<3,3> rotat = rotatMatrix_deg(0,0,yaw) * rotatMatrix_deg(0, -tilt, 0);
	return (vec::fixed<3> {0,0,hubHeight2CoG} + rotat * (nodeCoord_shaft + vec::fixed<3> {overhang,0,0}) );
}

// Coordinates of a blade node written in the earth coordinate system.
vec::fixed<3> Blade::nodeCoord_earth(const vec::fixed<6> &FOWTpos, const vec::fixed<3> &nodeCoord_fowt) const
{
	return (FOWTpos.rows(0,2) + rotatMatrix(FOWTpos.rows(3,5)) * nodeCoord_fowt);
}

// Calculate local angle of attack for a given local inflow angle.
// All the angles are in degrees.
double Blade::alpha(const unsigned int index, const double phi) const
{
	if (index < 0 || static_cast<unsigned int> (index) >= size())
	{
		throw std::runtime_error("Index out of range in Blade::alpha(const unsigned int index, const double phi).");
	}

	return (phi - twist(index) - m_pitch);
}
