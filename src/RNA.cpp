#include "IO.h"
#include "RNA.h"
#include "auxFunctions.h"

RNA::RNA()
{
	m_hubRadius = arma::datum::nan; // Initialize with NaN in order to know whether it was already calculated or not, since it is needed for calling Blade::setNodeCoord_hub()
	m_hubHeight2CoG = arma::datum::nan; // Same thing
}


/*****************************************************
	Setters
*****************************************************/
void RNA::readRotorSpeed(const std::string &data)
{
	readDataFromString(data, m_rotorSpeed);
}

void RNA::readRotorTilt(const std::string &data)
{
	readDataFromString(data, m_rotorTilt);
}

void RNA::readRotorYaw(const std::string &data)
{
	readDataFromString(data, m_rotorYaw);
}

void RNA::readNumBlades(const std::string &data)
{
	unsigned int numBlades;
	readDataFromString(data, numBlades);
	m_blades.resize(numBlades);

	double dAzimuth = 360 / numBlades;
	for (unsigned int ii = 0; ii < numBlades; ++ii)
	{		
		m_blades.at(ii).setInitialAzimuth(ii * dAzimuth);
	}
}

void RNA::readBladePrecone(const std::string &data)
{
	double precone;
	readDataFromString(data, precone);

	if (numBlades() == 0)
	{
		throw std::runtime_error("Need to specify number of blades before their properties. Error in input line " + std::to_string(IO::getInLineNumber()) + ".");
	}

	for (unsigned int ii = 0; ii < numBlades(); ++ii)
	{		
		m_blades.at(ii).setPrecone(precone);
	}	
}

void RNA::readBladePitch(const std::string &data)
{
	double pitch;
	readDataFromString(data, pitch);

	if (numBlades() == 0)
	{
		throw std::runtime_error("Need to specify number of blades before their properties. Error in input line " + std::to_string(IO::getInLineNumber()) + ".");
	}

	for (unsigned int ii = 0; ii < numBlades(); ++ii)
	{		
		m_blades.at(ii).setPitch(pitch);
	}		
}

void RNA::readBladeAeroLine(const std::string &data)
{
	double span = 0;
	double crvAC = 0;
	double swpAC = 0;
	double crvAng = 0;
	double twist = 0;
	double chord = 0;
	int airfoilID = 0;

	// The seven properties provided by each line are separated by white spaces in the input string (whitespace or tab)
	std::vector<std::string> input = stringTokenize(data, " \t");

	// Check number of inputs
	if (input.size() != 7)
	{
		throw std::runtime_error("Unable to read the blade aerodynamic properties in input line " + std::to_string(IO::getInLineNumber()) + ". Wrong number of parameters.");
	}

	// Read data
	readDataFromString(input.at(0), span);
	readDataFromString(input.at(1), crvAC);
	readDataFromString(input.at(2), swpAC);
	readDataFromString(input.at(3), crvAng);
	readDataFromString(input.at(4), twist);
	readDataFromString(input.at(5), chord);
	readDataFromString(input.at(6), airfoilID);

	if (numBlades() == 0)
	{
		throw std::runtime_error("Need to specify number of blades before their properties. Error in input line " + std::to_string(IO::getInLineNumber()) + ".");
	}

	// Add the data read in the current line to each blade element
	for (unsigned int ii = 0; ii < numBlades(); ++ii)
	{		
		m_blades.at(ii).addBladeAeroNode(span, crvAC, swpAC, crvAng, twist, chord, airfoilID, hubRadius());
	}	
}

// Add a new empty airfoil to m_airfoils
void RNA::addAirfoil()
{
	m_airfoils.push_back(Airfoil());
}

// Read line of the input file whith the airfoil properties. These properties are appended
// to the last airfoil in m_airfoils
void RNA::readAirfoilLine(const std::string &data)
{
	double angle;
	double CL;
	double CD;
	double CM;

	// The four properties provided by each line are separated by white spaces in the input string (whitespace or tab)
	std::vector<std::string> input = stringTokenize(data, " \t");

	// Check number of inputs
	if (input.size() != 4)
	{
		throw std::runtime_error("Unable to read the airfoil properties in input line " + std::to_string(IO::getInLineNumber()) + ". Wrong number of parameters.");
	}

	// Read data
	readDataFromString(input.at(0), angle);
	readDataFromString(input.at(1), CL);
	readDataFromString(input.at(2), CD);
	readDataFromString(input.at(3), CM);

	// Check if there is any airfoil for appending the data read from the input file
	if (m_airfoils.empty())
	{
		addAirfoil();
	}

	// Add the data read in the current line to the last airfoil in m_airfoils
	try
	{
		m_airfoils.back().addAirfoilLine(angle, CL, CD, CM);
	}
	catch (std::exception &exception)
	{
		throw std::runtime_error(std::string(exception.what()) + " Error in input line " + std::to_string(IO::getInLineNumber()) + ".");
	}
}

void RNA::readHubRadius(const std::string &data)
{
	readDataFromString(data, m_hubRadius);
}

void RNA::readHubHeight(const std::string &data)
{
	readDataFromString(data, m_hubHeight);
}

void RNA::readOverhang(const std::string &data)
{
	readDataFromString(data, m_overhang);
}

void RNA::setHubHeight2CoG(const double zCoG)
{
	m_hubHeight2CoG = hubHeight() - zCoG;
}

/*****************************************************
	Getters
*****************************************************/
double RNA::rotorSpeed() const
{
	return m_rotorSpeed;
}

double RNA::rotorTilt() const
{
	return m_rotorTilt;
}

double RNA::rotorYaw() const
{
	return m_rotorYaw;
}

unsigned int RNA::numBlades() const
{
	return static_cast<unsigned int>(m_blades.size());
}

double RNA::bladePrecone(const unsigned int ii) const
{
	return m_blades.at(ii).precone();
}

double RNA::bladePitch(const unsigned int ii) const
{
	return m_blades.at(ii).pitch();
}

std::string RNA::printBladeAero() const
{
	std::string output = "";
	output = output +  "BlSpn (m)\t" + "BlCrvAC (m)\t" + "BlSwpAC (m)\t" + "BlCrvAng (deg)\t" + "BlTwist (deg)\t" + "BlChord (m)\t" + "BlAfID" + '\n';

	// For printing, print only the properties of the first blade
	for (unsigned int ii = 0; ii < m_blades.at(0).size(); ++ii)
	{
		output = output + std::to_string(m_blades.at(0).span(ii)) + "\t" + std::to_string(m_blades.at(0).crvAC(ii)) + "\t" + std::to_string(m_blades.at(0).swpAC(ii)) + "\t" + std::to_string(m_blades.at(0).crvAng(ii)) + "\t" + std::to_string(m_blades.at(0).twist(ii)) + "\t" + std::to_string(m_blades.at(0).chord(ii)) + "\t" + std::to_string(m_blades.at(0).airoilID(ii)) + "\n";
	}

	return output;
}

std::string RNA::printAirfoils() const
{
	std::string output = "";	

	for (unsigned int ii = 0; ii < m_airfoils.size(); ++ii)
	{
		output = output + "\nAirfoil #" + std::to_string(ii) + '\n';
		output = output + "\t\tAngle (deg)\t" + "CL\t" + "CD\t" + "CM\n";
		for (unsigned int jj = 0; jj < m_airfoils.at(ii).size(); ++jj) 
		{
			output = output + "\t\t" + std::to_string(m_airfoils.at(ii).angle(jj)) + "\t" + std::to_string(m_airfoils.at(ii).getCL(jj)) + "\t" + std::to_string(m_airfoils.at(ii).getCD(jj)) + "\t" + std::to_string(m_airfoils.at(ii).getCM(jj)) + "\n";
		}
		output = output + "\n";
	}

	return output;
}

//std::string printAirofoils() const;

double RNA::hubRadius() const
{
	return m_hubRadius;
}

double RNA::hubHeight() const
{
	return m_hubHeight;
}

double RNA::overhang() const
{
	return m_overhang;
}

double RNA::hubHeight2CoG() const
{
	if (!arma::is_finite(m_hubHeight2CoG))
	{
		throw std::runtime_error("Need to call Blade::setHubHeight2CoG(const double hubHeight) at least once before calling RNA::hubHeight2CoG().");	
	}
	return m_hubHeight2CoG;
}


/*****************************************************
	Caculation functions
*****************************************************/	
double RNA::dAzimuth(const double time) const
{
	//time*(param.rt.Spd / 60) * 360
	return (360 * time * rotorSpeed() / 60.);
}

// FOWTpos is the position of the center of the FOWT coordinate system with respect to the earth fixed system, written in the earth fixed system
// In other words, it is {m_floater.CoG(),0,0,0} + FOWT.m_disp, i.e. the initial position of the floater CoG + its displacement
vec::fixed<6> RNA::aeroForce(const ENVIR &envir, const vec::fixed<6> &FOWTpos, const vec::fixed<6> &FOWTvel) const
{
	// Node position written in the different coordinate systems
	vec::fixed<3> nodeCoord_hub{ 0,0,0 };
	vec::fixed<3> nodeCoord_shaft{ 0,0,0 };
	vec::fixed<3> nodeCoord_fowt{ 0,0,0 };
	vec::fixed<3> nodeCoord_earth{ 0,0,0 };
	
	// Initialize some variables that are needed for the wind calculation
	vec::fixed<3> windVel{ 0,0,0 }; // Wind velocity. In the end of the for below, it is written in the blade node coordinate system
	vec::fixed<3> nodeVel{ 0,0,0 }; // Blade node structural velocity, written in the blade node coordinate system
	vec::fixed<3> cog2node{ 0,0,0}; // Vector given by the difference between the position of the blade node and the origin of the fowt coordinate system (around which the body rotation is provided). Written in the earth coordinate system.
	double windRel_nVel{ 0 }; // Relative wind speed in the x direction of the node coordinate system (normal to the plane of rotation)
	double windRel_tVel{ 0 }; // Relative wind speed in the y direction of the node coordinate system (in the plane of rotation)
	double localTipSpeed{ 0 };
	double deltaAzimuth = dAzimuth(envir.time());
	double totalAzimuth{ 0 };

	// Initialize the aerodynamic coefficients: lift, drag, in the normal, and in the tangential directions
	double Cl{ 0 };
	double Cd{ 0 };
	double Cn{ 0 };
	double Ct{ 0 };

	// Initialize the flow angles and parameters needed for the BEMT method.
	// The values of a, ap and phi are actually stored as members of the Blade class, but calling 
	// m_blades.at(iiBlades).phi() all the time in the loop below is cumbersome and inefficient.
	// The value of alpha is calculated by a method of the Blade class.
	double a{ 0 }; // Axial induction factor
	double ap{ 0 }; // Tangential induction factor
	double phi{ 0 }; // Local inflow angle
	double alpha{ 0 }; // Local angle of attack

	// At first, the wind velocity is calculated in the earth coordinate system.
	// For the calculations, it needs to be written in the blade node coordinate system.
	// This is done using the matrices rigidBodyRotation, which passes the vector
	// from the earth cordinate system to the one attached to the body,
	// and then rotorRotation, which is different for each blade and depends on the
	// rotor configuration (i.e. tilt, yaw, current azimuth and blade precone)
	mat::fixed<3, 3> rigidBodyRotation = rotatMatrix(-FOWTpos.rows(0, 2));
	mat::fixed<3, 3> rotorRotation{ 0,0,0 };
	
	for (unsigned int iiBlades = 0; iiBlades < m_blades.size(); ++iiBlades)
	{
		totalAzimuth = deltaAzimuth + m_blades.at(iiBlades).initialAzimuth();
		rotorRotation = rotatMatrix_deg(vec::fixed<3> {0, -m_blades.at(iiBlades).precone(), 0}) * rotatMatrix_deg(vec::fixed<3> { -totalAzimuth, -rotorTilt(), -rotorYaw()});


		for (unsigned int iiNodes = 0; iiNodes < m_blades.at(iiBlades).size(); ++iiNodes)
		{
			nodeCoord_hub = m_blades.at(iiBlades).nodeCoord_hub(iiNodes);
			nodeCoord_shaft = m_blades.at(iiBlades).nodeCoord_shaft(iiNodes, deltaAzimuth);
			nodeCoord_fowt = m_blades.at(iiBlades).nodeCoord_fowt(nodeCoord_shaft, rotorTilt(), rotorYaw(), hubHeight2CoG(), overhang());
			nodeCoord_earth = m_blades.at(iiBlades).nodeCoord_earth(FOWTpos, nodeCoord_fowt);
			
			// windVel is first calculated in the global coordinate system. 
			// We need to convert it to the node coordinate system, in which:
			// - windVel[0] is the component that is normal to the rotation plan
			// - windVel[1] is the component that is in the rotation plan and in the tangential direction
			// - windVel[2] is the component that is in the rotation plan and in the radial direction
			windVel[0] = envir.windVel_X(nodeCoord_earth);			
			windVel = rotorRotation * (rigidBodyRotation * windVel);

			// Structural velocity of the nodes. Need to be written in the node coordinate system
			cog2node = rotatMatrix(FOWTpos.rows(0, 2)) * nodeCoord_fowt;
			nodeVel = FOWTvel.rows(0,2) + arma::cross(FOWTvel.rows(3,5), cog2node);  // nodeVel = linearVel + angVel ^ r ; this is written in the earth coordinate system			
			nodeVel = rotorRotation * (rigidBodyRotation * nodeVel); // Need to pass to the node coordinate system
			nodeVel += rotatMatrix_deg(vec::fixed<3> {-totalAzimuth, 0, 0}) * arma::cross(vec::fixed<3> {rotorSpeed()*2*datum::pi/60, 0, 0}, nodeCoord_shaft); // Need to add the velocity due to the rotor rotation, written in the node coordinate system

			// Normal and tangential componentes of the relative wind velocity
			windRel_nVel = windVel[0] - nodeVel[0];
			windRel_tVel = windVel[1] - nodeVel[1];

			// Local tip speed ratio
			localTipSpeed = std::abs(windRel_tVel/windRel_nVel);

			// Total local solidity


			/* 
			Começa o Brent
			*/

			// As a first guess, the values of phi, a and ap from the previous time step are used.
			phi = m_blades.at(iiBlades).phi(iiNodes);
			alpha = m_blades.at(iiBlades).alpha(iiNodes);
			a = m_blades.at(iiBlades).axialIndFactor(iiNodes);
			ap = m_blades.at(iiBlades).tangIndFactor(iiNodes);


			// Aerodynamic coefficients in the normal and tangential directions
			Cl = m_airfoils.at(m_blades.at(iiBlades).airoilID(iiNodes)).CL(alpha);
			Cd = m_airfoils.at(m_blades.at(iiBlades).airoilID(iiNodes)).CD(alpha);
			Cn = Cl * cos(deg2rad(phi)) + Cd * sin(deg2rad(phi));
			Ct = Cl * sin(deg2rad(phi)) - Cd * cos(deg2rad(phi));

			if (envir.time() == 0)
			{

			}


		}
	}
	return vec::fixed<6> {0,0,0,0,0,0};
}


//double k(const double phi,)