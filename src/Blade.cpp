#include "Blade.h"
#include "auxFunctions.h" // For rotatMatrix()

using namespace arma;

/*****************************************************
	Setters
*****************************************************/
// Blade(const double precone, const double azimuth, const double pitch)
// {
// 	m_precone = precone;
// 	m_azimuth = azimuth;
// 	m_pitch = pitch;
// Como eles são lidos em linhas diferentes, acho melhor deixar o constructor com NaN e 
// fazer uma função set pra cada um
// }


void Blade::addBladeAeroLine(const double span, const double crvAC, const double swpAC, 
							 const double crvAng, const double twist, const double chord, 
							 const int airfoilID)
{
	m_span.push_back(span);
	m_crvAC.push_back(crvAC);
	m_swpAC.push_back(swpAC);
	m_crvAng.push_back(crvAng);
	m_twist.push_back(twist);
	m_chord.push_back(chord);
	m_airfoilID.push_back(airfoilID);

	vec::fixed<3> nodeCoord_hub;
	nodeCoord_hub.fill(datum::nan);
	m_nodeCoord_hub.push_back(nodeCoord_hub);
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
	return m_span.at(index);
}

double Blade::crvAC(const unsigned int index) const
{
	return m_crvAC.at(index);
}

double Blade::swpAC(const unsigned int index) const
{
	return m_swpAC.at(index);
}

double Blade::crvAng(const unsigned int index) const
{
	return m_crvAng.at(index);
}

double Blade::twist(const unsigned int index) const
{
	return m_twist.at(index);
}

double Blade::chord(const unsigned int index) const
{
	return m_chord.at(index);
}

int Blade::airoilID(const unsigned int index) const
{
	return m_airfoilID.at(index);
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

// Coordinates of a blade node written in the hub coordinate system.
//
// It is only calculated once, requiring the index of the node you are interested in and the hub radius.
// After the calculation, the value is stored in m_nodeCoord_hub(index) for future usage.
vec::fixed<3> Blade::nodeCoord_hub(const unsigned int index, const double hubRadius)
{
	if (arma::is_finite(m_nodeCoord_hub.at(index)))
	{
		return m_nodeCoord_hub.at(index);
	}

	vec::fixed<3> hubCoord;
	double r = hubRadius + span(index);		

	hubCoord[0] = r * tan(precone());
	hubCoord[1] = -r * sin(initialAzimuth()) * cos(precone());
	hubCoord[2] = r * cos(initialAzimuth()) * cos(precone());

	m_nodeCoord_hub.at(index) = hubCoord;

	return hubCoord;
}

// An overload for when Blade::nodeCoord_hub(unsigned int index, double hubRadius) was already calculated
vec::fixed<3> Blade::nodeCoord_hub(const unsigned int index) const
{
	if (!arma::is_finite(m_nodeCoord_hub.at(index)))
	{
		throw std::runtime_error("Need to calculate Blade::m_nodeCoord_hub(unsigned int index, double hubRadius) at least once before calling Blade::nodeCoord_hub(unsigned int index).");	
	}
	return m_nodeCoord_hub.at(index);
}

// Coordinates of a blade node written in the shaft coordinate system.
vec::fixed<3> Blade::nodeCoord_shaft(const unsigned int index, const double dAzimuth) const
{	
	double angle = (initialAzimuth() + dAzimuth) * datum::pi / 180.;

	return ( rotatMatrix(vec::fixed<3> {angle, 0, 0}) * nodeCoord_hub(index) );
}

// Overload for when the input is the nodeCoord_hub of a certain node itself.
vec::fixed<3> Blade::nodeCoord_shaft(const vec::fixed<3> &nodeCoord_hub, const double dAzimuth) const
{
	double angle = (initialAzimuth() + dAzimuth) * datum::pi / 180.;
	return ( rotatMatrix(vec::fixed<3> {angle, 0, 0}) * nodeCoord_hub );
}

// Coordinates of a blade node written in the tower coordinate system.
vec::fixed<3> Blade::nodeCoord_tower(const vec::fixed<3> &nodeCoord_shaft, const double tilt, const double yaw, const double hubHeight) const
{
	mat::fixed<3,3> rotat = rotatMatrix(vec::fixed<3> {0, yaw*datum::pi/180., 0}) * rotatMatrix(vec::fixed<3> {0, -tilt*datum::pi/180., 0});
	return (rotat * nodeCoord_shaft + vec::fixed<3> {0,0,hubHeight} );
}

// Coordinates of a blade node written in the earth coordinate system.
vec::fixed<3> Blade::nodeCoord_earth(const vec::fixed<6> &FOWTpos, const vec::fixed<3> &nodeCoord_tower) const
{
	return (FOWTpos.rows(0,2) + rotatMatrix(FOWTpos.rows(3,5)) * nodeCoord_tower);
}