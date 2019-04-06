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

// Should be called when setting the RNA in FOWT::setRNA
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

double RNA::bladePrecone(const unsigned int indexBlade) const
{
	return m_blades.at(indexBlade).precone();
}

double RNA::bladePitch(const unsigned int indexBlade) const
{
	return m_blades.at(indexBlade).pitch();
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

// FOWTpos is the position of the center of the FOWT coordinate system with respect to the earth fixed system, written in the earth fixed system.
// In other words, it is {m_floater.CoG(),0,0,0} + FOWT.m_disp, i.e. the initial position of the floater CoG + its displacement
vec::fixed<6> RNA::aeroForce(const ENVIR &envir, const vec::fixed<6> &FOWTpos, const vec::fixed<6> &FOWTvel) const
{
	// Total aerodynamic force/moments acting on the hub
	vec::fixed<6> aeroForce{ 0,0,0,0,0,0 };

	// Forces at each node
	double nodeForce_n{ 0 };
	double nodeForce_t{ 0 };
	double nodeForce_m{ 0 };

	// Force acting on each blade (integration of the forces at each node)
	double dr{ 0 }; // Integration step between nodes
	double rm{0}; // Center of the integration step
	vec::fixed<6> bladeForce;

	// Initialize the aerodynamic coefficients: normal direction, tangential directions and moment
	double Cn{ 0 };
	double Ct{ 0 };
	double Cm{ 0 };

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
	double windRelVel{ 0 }; // Total wind relative velocity. It is the one that includes the correction due to the axial and tangential induction factors
	double localTipSpeed{ 0 };
	double deltaAzimuth = dAzimuth(envir.time());
	double totalAzimuth{ 0 };

	// Initialize the flow angles and parameters needed for the BEMT method.
	// The values of a, ap and phi are actually stored as members of the Blade class, but calling 
	// m_blades[iiBlades].phi() all the time in the loop below is cumbersome and inefficient.
	// The value of alpha is calculated by a method of the Blade class.
	double a{ 0 }; // Axial induction factor
	double ap{ 0 }; // Tangential induction factor
	double phi{ 0 }; // Local inflow angle
	double localSolidity{ 0 }; 
	double F{ 0 }; // Tip + hub loss factor	

	// At first, the wind velocity is calculated in the earth coordinate system.
	// For the calculations, it needs to be written in the blade node coordinate system.
	// This is done using the matrices rigidBodyRotation, which passes the vector
	// from the earth cordinate system to the one attached to the body,
	// and then rotorRotation, which is different for each blade and depends on the
	// rotor configuration (i.e. tilt, yaw, current azimuth and blade precone)
	mat::fixed<3, 3> rigidBodyRotation = rotatMatrix(-FOWTpos.rows(0, 2));
	mat::fixed<3, 3> rotorRotation{ 0,0,0 };
	
	// Variables that provide the brackets for Brent method
	double phi_min;
	double phi_max;

	for (unsigned int iiBlades = 0; iiBlades < m_blades.size(); ++iiBlades)
	{
		totalAzimuth = deltaAzimuth + m_blades[iiBlades].initialAzimuth();
		rotorRotation = rotatMatrix_deg(0, -m_blades[iiBlades].precone(), 0) * rotatMatrix_deg(-totalAzimuth, -rotorTilt(), -rotorYaw());

		bladeForce.zeros();
		for (unsigned int iiNodes = 0; iiNodes < m_blades[iiBlades].size(); ++iiNodes)
		{
			nodeCoord_hub = m_blades[iiBlades].nodeCoord_hub(iiNodes);
			nodeCoord_shaft = m_blades[iiBlades].nodeCoord_shaft(iiNodes, deltaAzimuth);
			nodeCoord_fowt = m_blades[iiBlades].nodeCoord_fowt(nodeCoord_shaft, rotorTilt(), rotorYaw(), hubHeight2CoG(), overhang());
			nodeCoord_earth = m_blades[iiBlades].nodeCoord_earth(FOWTpos, nodeCoord_fowt);
			
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
			localSolidity = numBlades() * m_blades[iiBlades].localSolidity(iiNodes);

			/************
				The solution of BEMT equations begins here.
				The solution method proposed by Ning, 2013 is used.
			************/
	
			// Bracket the solution in order to use Brent's method
			if (calcRes(90, iiBlades, iiNodes, localSolidity, localTipSpeed, envir.useTipLoss(), envir.useHubLoss()) >= 0)
			{
				phi_min = 1e-6;
				phi_max = 90;
			}
			else if (calcRes(-45, iiBlades, iiNodes, localSolidity, localTipSpeed, envir.useTipLoss(), envir.useHubLoss()) >= 0)
			{
				phi_min = -45;
				phi_max = -0.001;
			}
			else
			{
				phi_min = 90 + 0.001;
				phi_max = 180 - 0.001;
			}	

			phi = Brent(phi_min, phi_max, iiBlades, iiNodes, localSolidity, localTipSpeed, envir.useTipLoss(), envir.useHubLoss());

			Cm = RNA::Cm(phi, iiBlades, iiNodes);
			Cn = RNA::Cn(phi, iiBlades, iiNodes);
			Ct = RNA::Ct(phi, iiBlades, iiNodes);
			F = calcF(phi, iiNodes, envir.useTipLoss(), envir.useHubLoss());
			a = calcAxialIndFactor(calcK(phi, localSolidity, Cn, F), phi, F);				
			ap = calcTangIndFactor(calcKp(phi, localSolidity, Ct, F), F);

			// Wind relative velocity, including axial and tangential induction factors
			windRelVel = sqrt( pow(windRel_tVel * (1.0 + ap), 2) + pow(windRel_nVel * (1 - a), 2) );

			/****
			Calculate nodal forces
			****/
			nodeForce_m = 0.5 * Cm * pow(m_blades[iiBlades].chord(iiNodes), 2) * envir.airDensity() * pow(windRelVel, 2);
			nodeForce_n = 0.5 * Cn * m_blades[iiBlades].chord(iiNodes) * envir.airDensity() * pow(windRelVel, 2);
			nodeForce_t = 0.5 * Ct * m_blades[iiBlades].chord(iiNodes) * envir.airDensity() * pow(windRelVel, 2);

			/****
			Calculate blade forces	
			****/		
			if (iiNodes == 0)
			{
				dr = (m_blades[iiBlades].radius(iiNodes+1) - m_blades[iiBlades].radius(iiNodes)) / 2.0;
				rm = m_blades[iiBlades].radius(iiNodes) + dr/2;
			}
			else if (iiNodes == m_blades[iiBlades].size())
			{
				dr = (m_blades[iiBlades].radius(iiNodes) - m_blades[iiBlades].radius(iiNodes-1)) / 2.0;
				rm = m_blades[iiBlades].radius(iiNodes) - dr/2;
			}
			else
			{	
				// For nodes outside the extremities, we have rm[i] = rm[i-1] + dr[i-1]/2 + dr[i]/2, with i = iiNodes.
				// Since this is inside a loop where iiNodes is increasing, we use the current values of rm and dr
				// because they correspond to rm[i-1] and dr[i-1].
				rm += dr/2;			
				dr = (m_blades[iiBlades].radius(iiNodes+1) - m_blades[iiBlades].radius(iiNodes-1)) / 2.0;				 
				rm += dr/2;
			}
			
			bladeForce[0] += nodeForce_n * dr;
			bladeForce[1] += -nodeForce_t * dr;
			bladeForce[2] += 0;
			bladeForce.rows(3,5) += arma::cross(vec::fixed<3>{0,0,rm}, bladeForce.rows(0,2)); // Moment generated by the aerodynamic forces
			bladeForce[5] += nodeForce_m * dr; // Need to add the moment due to Cm 
		}

			/****
			Calculate hub forces - i.e. the sum of the forces acting on each blade written in the hub coord system
			****/
			aeroForce.rows(0,2) += rotatMatrix_deg(totalAzimuth, m_blades[iiBlades].precone(), 0) * bladeForce.rows(0,2);
			aeroForce.rows(3,5) += rotatMatrix_deg(totalAzimuth, m_blades[iiBlades].precone(), 0) * bladeForce.rows(3,5);
	}

	IO::print2outLine(IO::OUTFLAG_AD_HUB_FORCE, aeroForce);
	return aeroForce;
}

double RNA::Brent(const double phi_min, const double phi_max, const unsigned int bladeIndex, const unsigned int nodeIndex, const double localSolidity, const double localTipSpeed, const bool useTipLoss, const bool useHubLoss) const
{
    double FPP{pow(10,-11)};
    double nearzero{pow(10,-20)};
	double error{0};
    double tolerance{pow(10,-5)};

    double AA = phi_min;
    double BB = phi_max;
	double FA = calcRes(AA, bladeIndex, nodeIndex, localSolidity, localTipSpeed, useTipLoss, useHubLoss);
	double FB = calcRes(BB, bladeIndex, nodeIndex, localSolidity, localTipSpeed, useTipLoss, useHubLoss);

	bool not_done{true};
	int m{0};

	double CC, DD, EE, tol1, xm, phi, SS, PP, QQ, RR;

    if (FA*FB < 0)
    {
        double FC = FB;
        while (not_done)
        {
            if (FB*FC > 0)
            {
                CC = AA;
                FC = FA;
                DD = BB - AA;
                EE = DD;
            }
            if (abs(FC) < abs(FB))
            {
                AA = BB; BB = CC; CC = AA;
                FA = FB; FB = FC; FC = FA;
            }
            tol1 = 2.0*FPP*abs(BB) + 0.5*tolerance;
            xm = 0.5*(CC-BB);
            if ((abs(xm) <= tol1) || (abs(FA) < nearzero))
            {
                phi = BB;
                not_done = false;
            }
            else
            {
                if ((abs(EE) >= tol1) && (abs(FA) > abs(FB)))
                {
                    SS = FB/FA;
                    if (abs(AA-CC) < nearzero)
                    {
                        PP = 2.0*xm*SS;
                        QQ = 1.0 - SS;
                    }
                    else
                    {
                        QQ = FA/FC;
                        RR = FB/FC;
                        PP = SS*(2.0*xm*QQ*(QQ-RR) - (BB-AA)*(RR-1.0));
                        QQ = (QQ-1.0)*(RR-1.0)*(SS-1.0);
                    }
                    if (PP>nearzero)
                    {
                        QQ = - QQ;
                    }
                    PP = abs(PP);


                    if ((2.0*PP) < minimum(3.0*xm*QQ-abs(tol1*QQ),abs(EE*QQ)))
                    {EE = DD;   DD = PP/QQ;}
                    else
                    {
                        DD = xm;
                        EE = DD;
                    }
                }
                else
                {
                    DD = xm;
                    EE = DD;
                }
                AA = BB;
                FA = FB;
                if (abs(DD) > tol1)
                {
                    BB = BB + DD;
                }
                else
                {
                    if (xm>0)
                    {
                        BB = BB + abs(tol1);
                    }
                    else
                    {
                        BB = BB - abs(tol1);
                    }
                }
                FB = calcRes(BB, bladeIndex, nodeIndex, localSolidity, localTipSpeed, useTipLoss, useHubLoss);
                m = m+1;
            }
        }
    }
    else if (FB == 0)
    {
        phi = BB;
    }
    else
    {
		throw std::runtime_error("Interval without roots in RNA::Brent().");
    }

    return phi;
}



double RNA::calcRes(const double phi, const unsigned int bladeIndex, const unsigned int nodeIndex, const double localSolidity, const double localTipSpeed, const bool useTipLoss, const bool useHubLoss) const
{
    double F = calcF(phi, nodeIndex, useTipLoss, useHubLoss);

	if (F == 0)
	{
		return 0;
	}
    
	double k = calcK(phi, localSolidity, Cn(phi, bladeIndex, nodeIndex), F);
	double kp = calcKp(phi, localSolidity, Ct(phi, bladeIndex, nodeIndex), F);

	double a{0};
    if (phi > 0)
    {
        a = calcAxialIndFactor(k, phi, F);

        if (abs(a-1) < 1e-6)
            return -cos(deg2rad(phi))*(1-kp)/localTipSpeed;
        else
        {
            return sin(deg2rad(phi))/(1-a) - cos(deg2rad(phi))*(1-kp)/localTipSpeed;
        }
    }
    else // Propeler brake-region
        return sin(deg2rad(phi))*(1-k) - cos(deg2rad(phi))*(1-kp)/localTipSpeed;
}



double RNA::Cn(const double phi, const unsigned int bladeIndex, const unsigned int nodeIndex) const
{
	double alpha = m_blades[bladeIndex].alpha(nodeIndex, phi);
	double Cl = m_airfoils[m_blades[bladeIndex].airoilID(nodeIndex)].CL(alpha);
	double Cd = m_airfoils[m_blades[bladeIndex].airoilID(nodeIndex)].CD(alpha);	

	return Cl * cos(deg2rad(phi)) + Cd * sin(deg2rad(phi));
}

double RNA::Ct(const double phi, const unsigned int bladeIndex, const unsigned int nodeIndex) const
{
	double alpha = m_blades[bladeIndex].alpha(nodeIndex, phi);
	double Cl = m_airfoils[m_blades[bladeIndex].airoilID(nodeIndex)].CL(alpha);
	double Cd = m_airfoils[m_blades[bladeIndex].airoilID(nodeIndex)].CD(alpha);	

	return Cl * sin(deg2rad(phi)) - Cd * cos(deg2rad(phi));
}

double RNA::Cm(const double phi, const unsigned int bladeIndex, const unsigned int nodeIndex) const
{
	double alpha = m_blades[bladeIndex].alpha(nodeIndex, phi);
	return m_airfoils[m_blades[bladeIndex].airoilID(nodeIndex)].CM(alpha);
}

// But what if phi > 180? This leads to weird values of F
double RNA::calcF(const double phi, const int nodeIndex, const bool useTipLoss, const bool useHubLoss) const
{
	if (phi == 0)
	{
		return 0;
	}

	double Ftip = 1;
	double Fhub = 1;

	double r = m_blades[0].radius(nodeIndex);
	double R = m_blades[0].radius(m_blades[0].size() - 1); // Total radius of the blade is equal to the radius of the last node

	if (useTipLoss)
	{
		Ftip = (2.0 / datum::pi) * acos(exp (-(numBlades() / 2.0) * (R - r) / (r * sin(deg2rad(phi)))) );
	}

	if (useHubLoss)
	{
		Fhub = (2.0 / datum::pi) * acos(exp(-(numBlades() / 2.0) * (r - m_hubRadius) / (m_hubRadius * sin(deg2rad(phi)))) );
	}

	return Ftip * Fhub;
}


// Convenience parameter defined in the algorithm for solving the BEMT equations proposed by Ning, 2013. Related to the axial induction factor by a = k / (1+k)
double RNA::calcK(const double phi, const double localSolidity, const double Cn, const double F) const
{
	if (F == 0)
	{
		return 1e6;
	}
	else
	{
		return localSolidity * Cn / (4 * F * pow(sin(deg2rad(phi)), 2));
	}
}


// Convenience parameter defined in the algorithm for solving the BEMT equations proposed by Ning, 2013. Related to the tangential induction factor by ap = kp / (1-kp)
double RNA::calcKp(const double phi, const double localSolidity, const double Ct, const double F) const
{
	if (F == 0)
	{
		return 1e6;
	}
	else
	{
		return localSolidity * Ct / (4 * F * sin(deg2rad(phi)) * cos(deg2rad(phi)) );
	}
}



double RNA::calcAxialIndFactor(const double k, const double phi, const double F) const
{
	if (k <= -1)
	{
		throw std::runtime_error("The function RNA::calcAxialIndFactor(const double k, const double phi, const double F) requires k > -1.");
	}

	// The method proposed by Ning, 2013 uses Buhl's derivation for the region where momentum theory can not be applied. 
	// That region corresponds to 0.4 < a < 1.0, which is the same as k >= 2/3. The variables g1, g2, and g3 are used in that condition.
	double g1, g2, g3; 

	if (F == 0)
	{
		return 1;
	}

	if (phi > 0)
	{
		if (k < 2 / 3.) // This one corresponds to the momentum region
		{
			// Deal with the case where the denominator is zero.
			if (abs(k+1) < 1e-6) 
			{
				// Since we are in the momentum region, if we say that the singularity yields a value of 'a' larger than 1, which corresponds to the
				// propeller brake region, we do not introduce any artificial solutions to the BEMT equations.
				// We have chosen 10, but any value > 1 would work
				return 10;
			}

			return k / (k + 1);
		}

		else // This one corresponds to the empirical region. Buhl's correction with Beta = 0.4 is used. See Ning, 2013
		{
			g1 = 2.0 * F * k - (10 / 9. - F);
			g2 = 2.0 * F * k - F * (4 / 3. - F);
			g3 = 2 * F * k - (25 / 9. - 2 * F);

            if (abs(g3) < 1e-6)
			{
				return (1.0 - 1.0/(2.0*sqrt(g2)));
			}
			else
			{
				return (g1 - sqrt(g2))/g3;
			}
		}
	}
	else // If phi < 0, the flow is in the propeller-brake region
	{
        if (k > 1.0)
            return k/(k - 1.0);
        else
            return 0;
	}
}

double RNA::calcTangIndFactor(const double kp, const double F) const
{
    if (F == 0)
        return -1.0;
    else
        return kp/(1.0-kp);
}


