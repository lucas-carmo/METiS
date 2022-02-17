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
void RNA::setUseTipLoss(const bool useTipLoss)
{
	m_useTipLoss = useTipLoss;
}

void RNA::setUseHubLoss(const bool useHubLoss)
{
	m_useHubLoss = useHubLoss;
}

void RNA::setUseSkewCorr(const bool useSkewCorr)
{
	m_useSkewCorr = useSkewCorr;
}

void RNA::setRotorSpeed(const double rotorSpeed)
{
	m_rotorSpeed = rotorSpeed;
}

void RNA::setRotorTilt(const double rotorTilt)
{
	m_rotorTilt = rotorTilt;
}

void RNA::setRotorYaw(const double rotorYaw)
{
	m_rotorYaw = rotorYaw;
}

void RNA::setNumBlades(const unsigned int numBlades)
{
	m_blades.resize(numBlades);

	double dAzimuth = 360 / numBlades;
	for (unsigned int ii = 0; ii < numBlades; ++ii)
	{
		m_blades.at(ii).setInitialAzimuth(ii * dAzimuth);
	}
}

void RNA::setBladePrecone(const double precone)
{
	if (numBlades() == 0)
	{
		throw std::runtime_error("Need to specify number of blades before their properties. Error in input line " + std::to_string(IO::getInLineNumber()) + ".");
	}

	for (unsigned int ii = 0; ii < numBlades(); ++ii)
	{
		m_blades.at(ii).setPrecone(precone);
	}
}

void RNA::setBladePitch(const double pitch)
{
	if (numBlades() == 0)
	{
		throw std::runtime_error("Need to specify number of blades before their properties. Error in input line " + std::to_string(IO::getInLineNumber()) + ".");
	}

	for (unsigned int ii = 0; ii < numBlades(); ++ii)
	{
		m_blades.at(ii).setPitch(pitch);
	}
}

void RNA::addBladeAeroNode(const double span, const double crvAC, const double swpAC, const double crvAng,
	                         const double twist, const double chord, const int airfoilID)
{
	if (numBlades() == 0)
	{
		throw std::runtime_error("Need to specify number of blades before their properties. Error in RNA::setBladeAeroLine");
	}

	// Set the data read for each blade element
	for (unsigned int ii = 0; ii < numBlades(); ++ii)
	{
		m_blades.at(ii).addBladeAeroNode(span, crvAC, swpAC, crvAng, twist, chord, airfoilID, hubRadius());
	}
}

// Add a new empty airfoil to m_airfoils
void RNA::addEmptyAirfoil()
{
	m_airfoils.push_back(Airfoil());
}

// Read airfoil properties corresponding to a given angle of attack.
// These properties are appended to the last airfoil in m_airfoils.
void RNA::addAirfoilData(const double angle, const double CL, const double CD, const double CM)
{
	if (m_airfoils.empty())
	{
		throw std::runtime_error("m_airfoils can not be empty when calling RNA::addAirfoilLines(double angle, double CL, double CD, double CM).");
	}

	m_airfoils.back().addAirfoilData(angle, CL, CD, CM);
}

void RNA::setHubRadius(const double hubRadius)
{
	m_hubRadius = hubRadius;
}

void RNA::setHubHeight(const double hubHeight)
{
	m_hubHeight = hubHeight;
}

void RNA::setOverhang(const double overhang)
{
	m_overhang = overhang;
}

// Should be called when setting the RNA in FOWT::setRNA
void RNA::setHubHeight2CoG(const double zCoG)
{
	m_hubHeight2CoG = hubHeight() - zCoG;
}

/*****************************************************
	Getters
*****************************************************/
bool RNA::useTipLoss() const
{
	return m_useTipLoss;
}

bool RNA::useHubLoss() const
{
	return m_useHubLoss;
}

bool RNA::useSkewCorr() const
{
	return m_useSkewCorr;
}

double RNA::rotorSpeed() const
{
	// In rpm
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

double RNA::hubRadius() const
{
	if (!arma::is_finite(m_hubRadius))
	{
		throw std::runtime_error("Need to call RNA::setHubRadius(const double hubRadius) before calling RNA::hubRadius().");
	}
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
	return (360 * time * rotorSpeed() / 60.);
}

// FOWTpos is the position of the center of the FOWT coordinate system with respect to the earth fixed system, written in the earth fixed system.
// In other words, it is {m_FOWT.CoG(),0,0,0} + FOWT.m_disp, i.e. the initial position of the FOWT CoG + its displacement
arma::vec::fixed<6> RNA::aeroForce(const ENVIR &envir, const arma::vec::fixed<6> &FOWTpos, const arma::vec::fixed<6> &FOWTvel)
{
	// Total aerodynamic force/moments acting on the hub
	arma::vec::fixed<6> aeroForce{ 0,0,0,0,0,0 };

	// Forces at each node
	double nodeForce_n{ 0 };
	double nodeForce_t{ 0 };
	double nodeForce_m{ 0 };

	// Force acting on each blade (integration of the forces at each node)
	double dr{ 0 }; // Integration step between nodes
	double rm{0}; // Center of the integration step
	arma::vec::fixed<6> bladeForce;

	// Initialize the aerodynamic coefficients: normal direction, tangential directions and moment
	double Cn{ 0 };
	double Ct{ 0 };
	double Cm{ 0 };

	// Node position written in the different coordinate systems
	arma::vec::fixed<3> nodeCoord_hub{ 0,0,0 };
	arma::vec::fixed<3> nodeCoord_shaft{ 0,0,0 };
	arma::vec::fixed<3> nodeCoord_fowt{ 0,0,0 };
	arma::vec::fixed<3> nodeCoord_earth{ 0,0,0 };

	// Initialize some variables that are needed for the wind calculation
	arma::vec::fixed<3> windVel{ 0,0,0 }; // Wind velocity. In the end of the for below, it is written in the blade node coordinate system
	arma::vec::fixed<3> nodeVel{ 0,0,0 }; // Blade node structural velocity, written in the blade node coordinate system
	arma::vec::fixed<3> cog2node{ 0,0,0}; // Vector given by the difference between the position of the blade node and the origin of the fowt coordinate system (around which the body rotation is provided). Written in the earth coordinate system.
	double windRel_nVel{ 0 }; // Relative wind speed in the x direction of the node coordinate system (normal to the plane of rotation)
	double windRel_tVel{ 0 }; // Relative wind speed in the y direction of the node coordinate system (in the plane of rotation)
	double windRelVel{ 0 }; // Total wind relative velocity. It is the one that includes the correction due to the axial and tangential induction factors
	double localTipSpeed{ 0 };
	double deltaAzimuth = dAzimuth(envir.time());
	double totalAzimuth{ 0 };

	// Initialize the flow angles and parameters needed for the BEMT method.
	// The values of a, ap and phi are actually stored as members of the Blade class, but calling
	// m_blades.at(iiBlades).phi() all the time in the loop below is cumbersome and inefficient.
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
	arma::mat::fixed<3, 3> rigidBodyRotation = rotatMatrix(-FOWTpos.rows(3, 5));
	arma::mat::fixed<3, 3> rotorRotation{ 0,0,0 };

	// Variables that provide the brackets for Brent method
	double phi_min;
	double phi_max;

	for (unsigned int iiBlades = 0; iiBlades < m_blades.size(); ++iiBlades)
	{
		totalAzimuth = deltaAzimuth + m_blades.at(iiBlades).initialAzimuth();
		rotorRotation = rotatMatrix_deg(0, -m_blades.at(iiBlades).precone(), 0) * rotatMatrix_deg(-totalAzimuth, rotorTilt(), -rotorYaw());

		bladeForce.zeros();
		for (unsigned int iiNodes = 0; iiNodes < m_blades.at(iiBlades).size(); ++iiNodes)
		{
			nodeCoord_hub = m_blades.at(iiBlades).nodeCoord_hub(iiNodes);
			nodeCoord_shaft = m_blades.at(iiBlades).nodeCoord_shaft(iiNodes, deltaAzimuth);
			nodeCoord_fowt = m_blades.at(iiBlades).nodeCoord_fowt(nodeCoord_shaft, rotorTilt(), rotorYaw(), overhang(), hubHeight2CoG());
			nodeCoord_earth = m_blades.at(iiBlades).nodeCoord_earth(FOWTpos, nodeCoord_fowt);

			// windVel is first calculated in the global coordinate system.
			// We need to convert it to the node coordinate system, in which:
			// - windVel[0] is the component that is normal to the rotation plan
			// - windVel[1] is the component that is in the rotation plan and in the tangential direction
			// - windVel[2] is the component that is in the rotation plan and in the radial direction
			envir.windVel(windVel, nodeCoord_earth);
			windVel = rotorRotation * (rigidBodyRotation * windVel);

			// Structural velocity of the nodes. Need to be written in the node coordinate system
			cog2node = nodeCoord_fowt;
			nodeVel = FOWTvel.rows(0,2) + arma::cross(FOWTvel.rows(3,5), cog2node);  // nodeVel = linearVel + angVel ^ r ; this is written in the earth coordinate system
			nodeVel = rotorRotation * (rigidBodyRotation * nodeVel); // Need to pass to the node coordinate system
			nodeVel += rotatMatrix_deg(-totalAzimuth, 0, 0) * arma::cross(arma::vec::fixed<3> {rotorSpeed()*2*datum::pi/60, 0, 0}, nodeCoord_shaft); // Need to add the velocity due to the rotor rotation, written in the node coordinate system

			// Normal and tangential componentes of the relative wind velocity
			windRel_nVel = windVel[0] - nodeVel[0];
			windRel_tVel = windVel[1] - nodeVel[1];

			// Local tip speed ratio
			localTipSpeed = std::abs(windRel_tVel/windRel_nVel);

			// Total local solidity
			localSolidity = numBlades() * m_blades.at(iiBlades).localSolidity(iiNodes);

			/************
				The solution of BEMT equations begins here, using the solution method proposed by Ning, 2013.
			************/

			// Bracket the solution in order to use Brent's method
			if (calcRes(90, iiBlades, iiNodes, localSolidity, localTipSpeed) >= 0)
			{
				phi_min = 1e-10;
				phi_max = 90;
			}
			else if (calcRes(-45, iiBlades, iiNodes, localSolidity, localTipSpeed) >= 0)
			{
				phi_min = -45;
				phi_max = -0.001;
			}
			else
			{
				phi_min = 90 + 0.001;
				phi_max = 180 - 0.001;
			}

			phi = Brent(phi_min, phi_max, iiBlades, iiNodes, localSolidity, localTipSpeed);

			// Calculate force coefficients
			Cm = RNA::Cm(phi, iiBlades, iiNodes);
			Cn = RNA::Cn(phi, iiBlades, iiNodes);
			Ct = RNA::Ct(phi, iiBlades, iiNodes);

			F = calcF(phi, iiNodes);
			a = calcAxialIndFactor(calcK(phi, localSolidity, Cn, F), phi, F);
			if (F == 0)
			{
				phi = 0;
			}
			ap = calcTangIndFactor(calcKp(phi, localSolidity, Ct, F), F);


			// Wind relative velocity, including axial and tangential induction factors
			windRelVel = std::sqrt(std::pow(windRel_tVel * (1.0 + ap), 2) + std::pow(windRel_nVel * (1 - a), 2));

			/****
			Calculate nodal forces
			****/
			nodeForce_m = 0.5 * Cm * std::pow(m_blades.at(iiBlades).chord(iiNodes), 2) * envir.airDensity() * std::pow(windRelVel, 2);
			nodeForce_n = 0.5 * Cn * m_blades.at(iiBlades).chord(iiNodes) * envir.airDensity() * std::pow(windRelVel, 2);
			nodeForce_t = 0.5 * Ct * m_blades.at(iiBlades).chord(iiNodes) * envir.airDensity() * std::pow(windRelVel, 2);

			/****
			Calculate blade forces
			****/
			if (iiNodes == 0)
			{
				dr = (m_blades.at(iiBlades).radius(iiNodes + 1) - m_blades.at(iiBlades).radius(iiNodes)) / 2.0;
				rm = m_blades.at(iiBlades).radius(iiNodes) + dr / 2;
			}
			else if (iiNodes == m_blades.at(iiBlades).size() - 1)
			{
				dr = (m_blades.at(iiBlades).radius(iiNodes) - m_blades.at(iiBlades).radius(iiNodes - 1)) / 2.0;
				rm = m_blades.at(iiBlades).radius(iiNodes) - dr / 2;
			}
			else
			{
				// For nodes outside the extremities, we have rm[i] = rm[i-1] + dr[i-1]/2 + dr[i]/2, with i = iiNodes.
				// Since this is inside a loop where iiNodes is increasing, we use the current values of rm and dr
				// because they correspond to rm[i-1] and dr[i-1].
				rm += dr / 2;
				dr = (m_blades.at(iiBlades).radius(iiNodes + 1) - m_blades.at(iiBlades).radius(iiNodes - 1)) / 2.0;
				rm += dr / 2;
			}

			bladeForce[0] += nodeForce_n * dr;
			bladeForce[1] += -nodeForce_t * dr;
			bladeForce[2] += 0;
			bladeForce.rows(3, 5) += arma::cross(arma::vec::fixed<3>{0, 0, rm}, arma::vec::fixed<3>{nodeForce_n, -nodeForce_t, 0}) * dr; // Moment generated by the aerodynamic forces
			bladeForce[5] += nodeForce_m * dr; // Need to add the moment due to Cm
		}

		/****
		Calculate hub forces - i.e. the sum of the forces acting on each blade written in the hub coord system
		****/
		aeroForce.rows(0,2) += rotatMatrix_deg(m_blades.at(iiBlades).initialAzimuth(), m_blades.at(iiBlades).precone(), 0) * bladeForce.rows(0,2);
		aeroForce.rows(3,5) += rotatMatrix_deg(m_blades.at(iiBlades).initialAzimuth(), m_blades.at(iiBlades).precone(), 0) * bladeForce.rows(3,5);
	}	
	IO::print2outLine(IO::OUTFLAG_AD_HUB_FORCE, aeroForce);
	
	// Need to write aeroForce in the global reference plane + change the fulcrum to the FOWT CoG
	// 1) Write aeroForce in the shaft coordinate system
	aeroForce.rows(0, 2) = rotatMatrix_deg(deltaAzimuth, 0, 0) * aeroForce.rows(0, 2);
	aeroForce.rows(3, 5) = rotatMatrix_deg(deltaAzimuth, 0, 0) * aeroForce.rows(3, 5);

	// 2) Write aeroForce in the fowt coordinate system
	arma::mat::fixed<3, 3> rotatShaft2FOWT = rotatMatrix_deg(0, 0, -rotorYaw()) * rotatMatrix_deg(0, -rotorTilt(), 0);
	aeroForce.rows(0, 2) = rotatShaft2FOWT * aeroForce.rows(0, 2);
	aeroForce.rows(3, 5) = rotatShaft2FOWT * aeroForce.rows(3, 5);

	// 3) Change the fulcrum to the CoG of the FOWT. Distance is hubHeight2CoG() along the tower and overhang() along the shaft
	arma::vec::fixed<3> leverHub2CoG = arma::vec::fixed<3>{ 0,0,hubHeight2CoG() } + rotatShaft2FOWT * arma::vec::fixed<3> {overhang(), 0, 0};
	aeroForce.rows(3, 5) += arma::cross(leverHub2CoG, aeroForce.rows(0, 2));

	// 4) Write aeroForce in the global coordinate system
	aeroForce.rows(0, 2) = rotatMatrix(FOWTpos.rows(3, 5)) * aeroForce.rows(0, 2);
	aeroForce.rows(3, 5) = rotatMatrix(FOWTpos.rows(3, 5)) * aeroForce.rows(3, 5);	

	return aeroForce;
}

double RNA::Brent(const double phi_min, const double phi_max, const unsigned int bladeIndex, const unsigned int nodeIndex, const double localSolidity, const double localTipSpeed) const
{
	double FPP{std::pow(10,-11)};
	double nearzero{std::pow(10,-20)};
	double error{0};
	double tolerance{std::pow(10,-5)};

	double AA = phi_min;
	double BB = phi_max;
	double FA = calcRes(AA, bladeIndex, nodeIndex, localSolidity, localTipSpeed);
	double FB = calcRes(BB, bladeIndex, nodeIndex, localSolidity, localTipSpeed);

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
            if (std::abs(FC) < std::abs(FB))
            {
                AA = BB; BB = CC; CC = AA;
                FA = FB; FB = FC; FC = FA;
            }
            tol1 = 2.0*FPP*std::abs(BB) + 0.5*tolerance;
            xm = 0.5*(CC-BB);
            if ((std::abs(xm) <= tol1) || (std::abs(FA) < nearzero))
            {
                phi = BB;
                not_done = false;
            }
            else
            {
                if ((std::abs(EE) >= tol1) && (std::abs(FA) > std::abs(FB)))
                {
                    SS = FB/FA;
                    if (std::abs(AA-CC) < nearzero)
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
                    PP = std::abs(PP);


                    if ((2.0*PP) < minimum(3.0*xm*QQ-std::abs(tol1*QQ),std::abs(EE*QQ)))
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
                if (std::abs(DD) > tol1)
                {
                    BB = BB + DD;
                }
                else
                {
                    if (xm>0)
                    {
                        BB = BB + std::abs(tol1);
                    }
                    else
                    {
                        BB = BB - std::abs(tol1);
                    }
                }
                FB = calcRes(BB, bladeIndex, nodeIndex, localSolidity, localTipSpeed);
                m = m+1;
            }
        }
    }
    else if (almostEqual(FB, 0, tolerance))
    {
        phi = BB;
    }
	else if (almostEqual(FA, 0, tolerance))
	{
		phi = AA;
	}
    else
    {
		throw std::runtime_error("Interval without roots in RNA::Brent(). The interval is ]" + std::to_string(AA) + "," + std::to_string(BB) + "[" );
    }

    return phi;
}



double RNA::calcRes(const double phi, const unsigned int bladeIndex, const unsigned int nodeIndex, const double localSolidity, const double localTipSpeed) const
{
    double F = calcF(phi, nodeIndex);
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

        if (std::abs(a-1) < 1e-5)
            return -std::cos(deg2rad(phi))*(1-kp)/localTipSpeed;
        else
        {
            return std::sin(deg2rad(phi))/(1-a) - std::cos(deg2rad(phi))*(1-kp)/localTipSpeed;
        }
    }
    else // Propeler brake-region
        return std::sin(deg2rad(phi))*(1-k) - std::cos(deg2rad(phi))*(1-kp)/localTipSpeed;
}


double RNA::Cn(const double phi, const unsigned int bladeIndex, const unsigned int nodeIndex) const
{
	double alpha = m_blades.at(bladeIndex).alpha(nodeIndex, phi);

	// Airfoil IDs are counted starting at 1, that's why we need to add -1 to access the element.
	double Cl = m_airfoils.at(m_blades.at(bladeIndex).airoilID(nodeIndex)-1).CL(alpha);
	double Cd = m_airfoils.at(m_blades.at(bladeIndex).airoilID(nodeIndex)-1).CD(alpha);

	return Cl * std::cos(deg2rad(phi)) + Cd * std::sin(deg2rad(phi));
}


double RNA::Ct(const double phi, const unsigned int bladeIndex, const unsigned int nodeIndex) const
{
	double alpha = m_blades.at(bladeIndex).alpha(nodeIndex, phi);
	double Cl = m_airfoils.at(m_blades.at(bladeIndex).airoilID(nodeIndex)-1).CL(alpha);
	double Cd = m_airfoils.at(m_blades.at(bladeIndex).airoilID(nodeIndex)-1).CD(alpha);

	return Cl * std::sin(deg2rad(phi)) - Cd * std::cos(deg2rad(phi));
}


double RNA::Cm(const double phi, const unsigned int bladeIndex, const unsigned int nodeIndex) const
{
	double alpha = m_blades.at(bladeIndex).alpha(nodeIndex, phi);
	return m_airfoils.at(m_blades.at(bladeIndex).airoilID(nodeIndex)-1).CM(alpha);
}


// But what if phi > 180? This leads to weird values of F
double RNA::calcF(const double phi, const int nodeIndex) const
{
	if (phi == 0)
	{
		return 0;
	}

	double Ftip = 1;
	double Fhub = 1;

	double r = m_blades[0].radius(nodeIndex);
	double R = m_blades[0].radius(m_blades[0].size() - 1); // Total radius of the blade is equal to the radius of the last node

	if (m_useTipLoss)
	{
		Ftip = (2.0 / datum::pi) * std::acos(std::exp(-(numBlades() / 2.0) * (R - r) / (r * std::sin(deg2rad(phi)))) );
	}

	if (m_useHubLoss)
	{
		Fhub = (2.0 / datum::pi) * std::acos(std::exp(-(numBlades() / 2.0) * (r - m_hubRadius) / (m_hubRadius * std::abs(std::sin(deg2rad(phi))))) );
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
		return localSolidity * Cn / (4 * F * std::pow(std::sin(deg2rad(phi)), 2));
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
		return localSolidity * Ct / (4 * F * std::sin(deg2rad(phi)) * std::cos(deg2rad(phi)) );
	}
}


double RNA::calcAxialIndFactor(const double k, const double phi, const double F) const
{
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
			if (std::abs(k+1) < 1e-6)
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

            if (std::abs(g3) < 1e-6)
			{
				return (1.0 - 1.0/(2.0*std::sqrt(g2)));
			}
			else
			{
				return (g1 - std::sqrt(g2))/g3;
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
