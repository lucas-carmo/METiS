//#include "Floater.h"
#include "MorisonCirc.h"
#include "MorisonRect.h"
#include "IO.h"
#include "FOWT.h"
#include "Wave.h"

#include <iostream>
#include <vector>
#include <algorithm>    // std::binary_search
#include <utility> // For std::move
#include <stdexcept> // For std::exception
#include <armadillo>

using namespace arma;

/*****************************************************
	Constructors
*****************************************************/
Floater::Floater()
{
	// Assign nan to these variable in order to check later whether they were read from the input file or not
	m_mass = datum::nan;
	m_CoG.fill(datum::nan);	
	m_inertia.fill(datum::nan);
	m_addedMass_t0.fill(datum::nan);
}


/*****************************************************
	Overloaded operators
*****************************************************/
Floater& Floater::operator= (const Floater &floater)
{
	// Shallow copy for simple member variables
	m_mass = floater.m_mass;
	m_CoG = floater.m_CoG;
	m_inertia = floater.m_inertia;	
	m_addedMass_t0 = floater.m_addedMass_t0;
	m_disp = floater.m_disp;
	m_disp_sd = floater.m_disp_sd;

	// The member variables that need deep copying are:
	// - m_MorisonElements;

	// Check for self-assignment
    if (this == &floater)
		return *this;

	// In case there is data in m_MorisonElements
	m_MorisonElements.clear();

	// Resize m_MorisonElement to match the size of the one in the input floater
	m_MorisonElements.resize( floater.m_MorisonElements.size() );

	// m_MorisonElements.at(ii) is a std::unique_ptr<MorisonElements>
	for (int ii = 0; ii < floater.m_MorisonElements.size(); ++ii)
	{
		MorisonElement* rawPtr = floater.m_MorisonElements.at(ii)->clone();
		std::unique_ptr<MorisonElement> smartPtr(rawPtr);
		m_MorisonElements.at(ii) = std::move(smartPtr);
	}
		
    return *this;	
}


/*****************************************************
	Setters
*****************************************************/
void Floater::setMass(const double mass)
{
	m_mass = mass;
}

void Floater::setInertia(const vec::fixed<6> &inertia)
{
	m_inertia = inertia;
}

void Floater::setCoG(const vec::fixed<3> &cog)
{
	m_CoG = cog;
}

void Floater::setInstantSD(bool instantSD)
{
	m_instantSD = instantSD;
}

void Floater::evaluateQuantitiesAtBegin(const ENVIR &envir, const int hydroMode)
{
	for (int ii = 0; ii < m_MorisonElements.size(); ++ii)
	{
		// Evaluate fluid kinematics at each node
		m_MorisonElements.at(ii)->evaluateQuantitiesAtBegin(envir, hydroMode);
	}
}


/*
	Add circular cylinder Morison Element to m_MorisonElements
*/
void Floater::addMorisonCirc(vec::fixed<3> &node1_coord, vec::fixed<3> &node2_coord, const double diam, const double CD, const double CM, const unsigned int numIntPoints,
	const double axialCD_1, const double axialCa_1, const double axialCD_2, const double axialCa_2, const bool botPressFlag)
{
	// Since many times we need node1 to be below node2, it is better to swap them here in order to have less swaps in the future
	if (node1_coord[2] > node2_coord[2])
	{
		node1_coord.swap(node2_coord);
	}

	// Morison Elements need the position of the floater CoG
	// to define their relative position in the rigid body
	if (CoG().has_nan()) // If this is true, then the CoG of the floater wasn't read yet
	{
		throw std::runtime_error("Floater CoG must be specified before Morison Elements.");
	}

	// Create a circular cylinder Morison Element using the following constructor and add it to m_MorisonElements.
	m_MorisonElements.push_back(std::make_unique<MorisonCirc>(node1_coord, node2_coord, CoG(), numIntPoints,
		botPressFlag, axialCD_1, axialCa_1, axialCD_2, axialCa_2, diam, CD, CM));
}

/*
	Add rectangular cylinder Morison Element to m_MorisonElements
*/
void Floater::addMorisonRect(vec::fixed<3> &node1_coord, vec::fixed<3> &node2_coord, vec::fixed<3> &node3_coord, const double diam_X, const double diam_Y, 
	const double CD_X, const double CD_Y, const double CM_X, const double CM_Y, const unsigned int numIntPoints,
	const double axialCD_1, const double axialCa_1, const double axialCD_2, const double axialCa_2, const bool botPressFlag)
{
	// Morison Elements need the position of the floater CoG
	// to define their relative position in the rigid body
	if (CoG().has_nan()) // If this is true, then the CoG of the floater wasn't read yet
	{
		throw std::runtime_error("Floater CoG must be specified before Morison Elements.");
	}

	// Create a rectangular cylinder Morison Element using the following constructor and add it to m_MorisonElements.
	m_MorisonElements.push_back(std::make_unique<MorisonRect>(node1_coord, node2_coord, node3_coord, CoG(), numIntPoints,
		botPressFlag, axialCD_1, axialCa_1, axialCD_2, axialCa_2, diam_X, CD_X, CM_X, diam_Y, CD_Y, CM_Y));
}

mat::fixed<6, 6> Floater::addedMass_t0() const
{
	if (m_addedMass_t0.is_finite())
	{
		return m_addedMass_t0;
	}
	else
	{
		throw std::runtime_error("Tried to call Floater::addedMass_t0(), but it was not calculated yet. Try calling Floater::addedMass(const double density, const int hydroMode) first.");
	}
}

mat::fixed<6, 6> Floater::hydrostaticStiffness() const
{
	return m_Khs;
}

vec::fixed<6> Floater::CoGPos() const
{
	return join_cols(m_CoG, zeros(3,1)) + m_disp;
}

vec::fixed<6> Floater::CoGPos_1stOrd() const
{
	return join_cols(m_CoG, zeros(3, 1)) + m_disp_1stOrd;
}

vec::fixed<6> Floater::CoGPos_sd() const
{
	return join_cols(m_CoG, zeros(3, 1)) + m_disp_sd;
}

/*****************************************************
	Getters
*****************************************************/
vec::fixed<3> Floater::CoG() const
{
	return m_CoG;
}

double Floater::mass() const
{
	return m_mass;
}

mat::fixed<6,6> Floater::inertiaMatrix() const
{
	mat::fixed<6, 6> A(fill::zeros);
	A(0, 0) = A(1, 1) = A(2, 2) = m_mass;
	A(3, 3) = m_inertia(0);
	A(4, 4) = m_inertia(1);
	A(5, 5) = m_inertia(2);

	A(3, 4) = A(4, 3) = m_inertia(3);
	A(3, 5) = A(5, 3) = m_inertia(4);
	A(4, 5) = A(5, 4) = m_inertia(5);

	return A;
}

std::string Floater::printMass() const
{
	return std::to_string(m_mass);
}

std::string Floater::printInertia() const
{
	return "Ixx = " + std::to_string(m_inertia(0)) + " Iyy = " + std::to_string(m_inertia(1)) + " Izz = " + std::to_string(m_inertia(2)) +
		  " Ixy = " + std::to_string(m_inertia(3)) + " Ixz = " + std::to_string(m_inertia(4)) + " Iyy = " + std::to_string(m_inertia(5));
}

std::string Floater::printCoG() const
{
	std::string output = "(" + std::to_string( m_CoG(0) );
	for (int ii = 1; ii < m_CoG.n_elem; ++ii)
	{
		output = output + "," + std::to_string( m_CoG(ii) );
	}
	
	return output + ")";
}

std::string Floater::printMorisonElements() const
{
	std::string output = "";	
	output = output + "Number of Morison Elements: " + std::to_string(m_MorisonElements.size()) + "\n";
	for (int ii = 0; ii < m_MorisonElements.size(); ++ii)
	{
		output = output + "Morison Element #" + std::to_string(ii) + '\n';
		output = output + m_MorisonElements.at(ii)->print() + '\n';
	}
	return output;
}



/*****************************************************
	Forces, acceleration, displacement, etc
*****************************************************/
void Floater::update(const ENVIR &envir, const vec::fixed<6> &FOWTdisp, const vec::fixed<6> &FOWTvel, const vec::fixed<6> &FOWTdisp_1stOrd, const vec::fixed<6> &FOWTvel_1stOrd, const vec::fixed<6> &FOWTdisp_SD, const vec::fixed<6> &FOWTvel_SD)
{
	m_disp = FOWTdisp;
	m_disp_sd = FOWTdisp_SD;
	m_disp_1stOrd = FOWTdisp_1stOrd;

	for (int ii = 0; ii < m_MorisonElements.size(); ++ii)
	{
		m_MorisonElements.at(ii)->updateMorisonElement(envir, CoGPos(), FOWTvel, CoGPos_1stOrd(), FOWTvel_1stOrd, CoGPos_sd(), FOWTvel_SD);
	}
}

void Floater::setNode1stAcc(const vec::fixed<6> &FOWTacc_1stOrd)
{
	for (int ii = 0; ii < m_MorisonElements.size(); ++ii)
	{
		m_MorisonElements.at(ii)->updateAcc1stOrdNodeAtWL(FOWTacc_1stOrd);
	}
}

// Needed this function to evaluate the added mass at the beginning of the simulation
void Floater::setAddedMass_t0(const double density)
{
	m_addedMass_t0 = density*addedMass(1);
}

void Floater::setStiffnessMatrix(const double density, const double gravity)
{
	// Quantities required for evaluating the hydrostatic matrix
	double zb{ 0 }, V{ 0 }, Awl{ 0 }, xwl{ 0 }, ywl{ 0 }, Ixx{ 0 }, Iyy{ 0 }, Ixy{ 0 };

	for (int ii = 0; ii < m_MorisonElements.size(); ++ii)
	{
		double aux_zb{ 0 }, aux_V{ 0 }, aux_Awl{ 0 }, aux_xwl{ 0 }, aux_ywl{ 0 }, aux_Ixx{ 0 }, aux_Iyy{ 0 }, aux_Ixy{ 0 };
		m_MorisonElements.at(ii)->quantities4hydrostaticMatrix(aux_zb, aux_V, aux_Awl, aux_xwl, aux_ywl, aux_Ixx, aux_Iyy, aux_Ixy);

		Awl += aux_Awl;
		V += aux_V;
		xwl += aux_Awl * aux_xwl;
		zb += aux_V * aux_zb;
		Ixx += aux_Ixx + aux_Awl * aux_ywl * aux_ywl;
		Iyy += aux_Iyy + aux_Awl * aux_xwl*aux_xwl;
		Ixy += aux_Ixy + aux_Awl * aux_xwl*aux_ywl;
	}
	xwl *= 1 / Awl;
	ywl *= 1 / Awl;
	zb *= 1 / V;


	m_Khs.at(2, 2) = density * gravity * Awl;
	m_Khs.at(3, 3) = density * gravity * (Ixx + V * (zb - m_CoG.at(2)));
	m_Khs.at(4, 4) = density * gravity * (Iyy + V * (zb - m_CoG.at(2)));
	m_Khs.at(3, 4) = density * gravity * Ixy;
	m_Khs.at(4, 3) = m_Khs.at(3, 4);

	m_Vol = V;
}

mat::fixed<6, 6> Floater::addedMass(const double density, const int hydroMode) 
{
	mat::fixed<6, 6> A(fill::zeros);

	A = density*addedMass(hydroMode);
	return A;
}

vec::fixed<6> Floater::hydrodynamicForce_1stOrd(const ENVIR &envir) const
{
	vec::fixed<6> force_1stP(fill::zeros); // Due to first-order potential
	vec::fixed<3> refPt = CoGPos_sd().rows(0, 2);

	// Force acting on each element
	for (int ii = 0; ii < m_MorisonElements.size(); ++ii)
	{
		// If the cylinder is above the waterline, do not need to compute the hydrodynamic forces
		if (m_MorisonElements.at(ii)->node1Pos_sd().at(2) > 0)
			continue;

		force_1stP += m_MorisonElements.at(ii)->hydroForce_1st(envir, refPt);
	}

	IO::print2outLine(IO::OUTFLAG_HD_FORCE_1STP, force_1stP);
	return force_1stP;
}

vec::fixed<6> Floater::hydrodynamicForce_drag1stOrd(const ENVIR &envir) const
{
	vec::fixed<6> force_drag(fill::zeros); // Drag force
	vec::fixed<3> refPt = CoGPos_sd().rows(0, 2);

	// Force acting on each element
	for (int ii = 0; ii < m_MorisonElements.size(); ++ii)
	{
		// If the cylinder is above the waterline, do not need to compute the hydrodynamic forces
		if (m_MorisonElements.at(ii)->node1Pos_sd().at(2) > 0)
			continue;

		force_drag += m_MorisonElements.at(ii)->hydroForce_drag(envir, refPt, true);
	}

	return force_drag;
}

vec::fixed<6> Floater::hydrodynamicForce_dragTotal(const ENVIR &envir) const
{
	vec::fixed<6> force_drag(fill::zeros); // Drag force
	vec::fixed<3> refPt = CoGPos_sd().rows(0, 2);

	// Force acting on each element
	for (int ii = 0; ii < m_MorisonElements.size(); ++ii)
	{
		// If the cylinder is above the waterline, do not need to compute the hydrodynamic forces
		if (m_MorisonElements.at(ii)->node1Pos_sd().at(2) > 0)
			continue;

		force_drag += m_MorisonElements.at(ii)->hydroForce_drag(envir, refPt, false);
	}

	IO::print2outLine(IO::OUTFLAG_HD_FORCE_DRAG, force_drag);
	return force_drag;
}

// Moments are output with respect to the CoG of the floater
vec::fixed<6> Floater::hydrodynamicForce_2ndOrd_slenderbody(const ENVIR &envir, const vec::fixed<6> &F_1stOrd) const
{	
	// These forces below are used to output the different components to the output files and internally in MorisonElement.hydrodynamicForce().	
	vec::fixed<6> force_eta(fill::zeros);  // Due to relative wave elevation
	vec::fixed<6> force_conv(fill::zeros); // Due to convective acceleration
	vec::fixed<6> force_axdv(fill::zeros); // Due to axial-divergence acceleration
	vec::fixed<6> force_acgr(fill::zeros); // Due to gradient of fluid acceleration
	vec::fixed<6> force_rotn(fill::zeros); // Due to the rotation of the normal vector
	vec::fixed<6> force_2ndP(fill::zeros); // Component due to 2nd order incoming flow
	vec::fixed<6> force_rslb(fill::zeros); // Due to the rotation term from slender-body approximation
	vec::fixed<6> force_rem(fill::zeros);  // Remaining force components

	vec::fixed<3> refPt = CoGPos_sd().rows(0, 2);
	vec::fixed<3> Rvec = m_disp_1stOrd.rows(3, 5);

	// Force acting on each element
	for (int ii = 0; ii < m_MorisonElements.size(); ++ii)
	{
		// If the cylinder is above the waterline, do not need to compute the hydrodynamic forces
		if (m_MorisonElements.at(ii)->node1Pos_sd().at(2) > 0)
			continue;

		if (envir.waveStret() == 1)
		{
			force_eta += m_MorisonElements.at(ii)->hydroForce_relWaveElev_inertia(envir, refPt);
			force_eta += m_MorisonElements.at(ii)->hydroForce_relWaveElev_drag(envir, refPt);
		}

		force_acgr += m_MorisonElements.at(ii)->hydroForce_accGradient(envir, refPt);
		force_conv += m_MorisonElements.at(ii)->hydroForce_convecAcc(envir, refPt);
		force_axdv += m_MorisonElements.at(ii)->hydroForce_axDiverg(envir, refPt);
		force_2ndP += m_MorisonElements.at(ii)->hydroForce_2ndPot(envir, refPt);
		force_rslb += m_MorisonElements.at(ii)->hydroForce_slendBodyRot(envir, refPt);
		force_rem += m_MorisonElements.at(ii)->hydroForce_rem(envir, refPt);
	}

	force_rotn.rows(0, 2) = cross(Rvec, F_1stOrd.rows(0, 2));
	force_rotn.rows(3, 5) = cross(Rvec, F_1stOrd.rows(3, 5));

	vec::fixed<6> force = force_eta + force_conv + force_axdv + force_acgr + force_rotn + force_2ndP + force_rslb + force_rem;

	IO::print2outLine(IO::OUTFLAG_HD_FORCE_ETA, force_eta);
	IO::print2outLine(IO::OUTFLAG_HD_FORCE_CONV, force_conv);
	IO::print2outLine(IO::OUTFLAG_HD_FORCE_AXDV, force_axdv);
	IO::print2outLine(IO::OUTFLAG_HD_FORCE_ACGR, force_acgr);
	IO::print2outLine(IO::OUTFLAG_HD_FORCE_ROTN, force_rotn);
	IO::print2outLine(IO::OUTFLAG_HD_FORCE_2NDP, force_2ndP);
	IO::print2outLine(IO::OUTFLAG_HD_FORCE_RSLB, force_rslb);
	IO::print2outLine(IO::OUTFLAG_HD_FORCE_REM, force_rem);

	return force;
}

mat::fixed<6, 6> Floater::addedMass(const int hydroMode) const
{
	mat::fixed<6, 6> A(fill::zeros);
	vec::fixed<3> refPt(fill::zeros);

	for (int ii = 0; ii < m_MorisonElements.size(); ++ii)
	{
		// Evaluated at the fixed initial position
		if (hydroMode < 2)
		{
			refPt = CoGPos_sd().rows(0, 2);
		}
		else
		{
			refPt = CoGPos().rows(0,2);
		}

		A += m_MorisonElements.at(ii)->addedMass_perp(1, refPt, hydroMode) + m_MorisonElements.at(ii)->addedMass_paral(1, refPt, hydroMode);
	}
	return A;
}

vec::fixed<6> Floater::hydrostaticForce_stiffnessPart(bool flagUse1stOrd) const
{
	vec::fixed<6> disp{ m_disp };
	if (flagUse1stOrd)
	{
		disp = m_disp_1stOrd;
	}

	vec::fixed<6> force = -m_Khs * disp;

	// Print only the total hydrostatic force, not the intermediary evaluation with the first-order
	if (!flagUse1stOrd)
	{
		IO::print2outLine(IO::OUTFLAG_HS_FORCE, force);
	}
	return force;
}

vec::fixed<6> Floater::hydrostaticForce_staticBuoyancy(const double rho, const double g) const
{	
	return { 0,0, rho * g * m_Vol, 0, 0, 0 };
}

mat Floater::QTF_ToTime_AppN(const ENVIR& envir, vec& m_p12omega, vec& m_p12beta, const cx_mat& m_p12surge1, const cx_mat& m_p12sway1, const cx_mat& m_p12heave1, const cx_mat& m_p12roll1, const cx_mat& m_p12pitch1, const cx_mat& m_p12yaw1) // Generates time series from Newman's Approximation
{
	const Wave& waveQTF(envir.getWave(0));
		
	size_t nRows = m_p12omega.size();
    size_t nCols = m_p12omega.size();

	// Check the wave incidence requested in the input file
	double findwavedir = waveQTF.direction();

	uvec  wavedir = find(m_p12beta == findwavedir);
		
	// Assemble the QTF matrix for the requested wave incidence
	cx_mat m_p12surge = reshape(m_p12surge1.col(wavedir(0)), nRows, nCols);
	
	cx_mat m_p12sway  = reshape(m_p12sway1.col(wavedir(0)), nRows, nCols);
	cx_mat m_p12heave = reshape(m_p12heave1.col(wavedir(0)), nRows, nCols);
	cx_mat m_p12roll  = reshape(m_p12roll1.col(wavedir(0)), nRows, nCols);
	cx_mat m_p12pitch = reshape(m_p12pitch1.col(wavedir(0)), nRows, nCols);
	cx_mat m_p12yaw   = reshape(m_p12yaw1.col(wavedir(0)), nRows, nCols);
		
	

	mat m_p12surgeAmp = abs(m_p12surge);
	mat m_p12surgePha = arg(m_p12surge);
	mat m_p12surgeRe = real(m_p12surge);
	mat m_p12surgeIm = imag(m_p12surge);
		
	mat m_p12swayAmp = abs(m_p12sway);
	mat m_p12swayPha = arg(m_p12sway);
	mat m_p12swayRe = real(m_p12sway);
	mat m_p12swayIm = imag(m_p12sway);

	mat m_p12heaveAmp = abs(m_p12heave);
	mat m_p12heavePha = arg(m_p12heave);
	mat m_p12heaveRe = real(m_p12heave);
	mat m_p12heaveIm = imag(m_p12heave);

	mat m_p12rollAmp = abs(m_p12roll);
	mat m_p12rollPha = arg(m_p12roll);
	mat m_p12rollRe = real(m_p12roll);
	mat m_p12rollIm = imag(m_p12roll);

	mat m_p12pitchAmp = abs(m_p12pitch);
	mat m_p12pitchPha = arg(m_p12pitch);
	mat m_p12pitchRe = real(m_p12pitch);
	mat m_p12pitchIm = imag(m_p12pitch);

	mat m_p12yawAmp = abs(m_p12yaw);
	mat m_p12yawPha = arg(m_p12yaw);
	mat m_p12yawRe = real(m_p12yaw);
	mat m_p12yawIm = imag(m_p12yaw);

	// Get the mean drift force from QTF
	AppNauxiliar.surgeAmp = m_p12surgeAmp.diag();
    AppNauxiliar.surgePha = m_p12surgePha.diag();
    AppNauxiliar.surgeRe  = m_p12surgeRe.diag();
    AppNauxiliar.surgeIm  = m_p12surgeIm.diag();
	
    AppNauxiliar.swayAmp = m_p12swayAmp.diag();
    AppNauxiliar.swayPha = m_p12swayPha.diag();
    AppNauxiliar.swayRe  = m_p12swayRe.diag();
    AppNauxiliar.swayIm  = m_p12swayIm.diag();
    
    AppNauxiliar.heaveAmp = m_p12heaveAmp.diag();
    AppNauxiliar.heavePha = m_p12heavePha.diag();
    AppNauxiliar.heaveRe  = m_p12heaveRe.diag();
    AppNauxiliar.heaveIm  = m_p12heaveIm.diag();

    AppNauxiliar.rollAmp = m_p12rollAmp.diag();
    AppNauxiliar.rollPha = m_p12rollPha.diag();
    AppNauxiliar.rollRe  = m_p12rollRe.diag();
    AppNauxiliar.rollIm  = m_p12rollIm.diag();

    AppNauxiliar.pitchAmp = m_p12pitchAmp.diag();
    AppNauxiliar.pitchPha = m_p12pitchPha.diag();
    AppNauxiliar.pitchRe  = m_p12pitchRe.diag();
    AppNauxiliar.pitchIm  = m_p12pitchIm.diag();

    AppNauxiliar.yawAmp = m_p12yawAmp.diag();
    AppNauxiliar.yawPha = m_p12yawPha.diag();
    AppNauxiliar.yawRe  = m_p12yawRe.diag();
    AppNauxiliar.yawIm  = m_p12yawIm.diag();

	// Create a matrix with size nRows x nCols
    AppN.surgeAmp = zeros<mat>(nRows, nCols);
    AppN.swayAmp  = zeros<mat>(nRows, nCols);
    AppN.heaveAmp = zeros<mat>(nRows, nCols);
    AppN.rollAmp  = zeros<mat>(nRows, nCols);
    AppN.pitchAmp = zeros<mat>(nRows, nCols);
    AppN.yawAmp   = zeros<mat>(nRows, nCols);
    
    AppN.surgePha = zeros<mat>(nRows, nCols);
    AppN.swayPha  = zeros<mat>(nRows, nCols);
    AppN.heavePha = zeros<mat>(nRows, nCols);
    AppN.rollPha  = zeros<mat>(nRows, nCols);
    AppN.pitchPha = zeros<mat>(nRows, nCols);
    AppN.yawPha   = zeros<mat>(nRows, nCols);
    
    AppN.surgeRe  = zeros<mat>(nRows, nCols);
    AppN.swayRe   = zeros<mat>(nRows, nCols);
    AppN.heaveRe  = zeros<mat>(nRows, nCols);
    AppN.rollRe   = zeros<mat>(nRows, nCols);
    AppN.pitchRe  = zeros<mat>(nRows, nCols);
    AppN.yawRe    = zeros<mat>(nRows, nCols);
    
    AppN.surgeIm  = zeros<mat>(nRows, nCols);
    AppN.swayIm   = zeros<mat>(nRows, nCols);
    AppN.heaveIm  = zeros<mat>(nRows, nCols);
    AppN.rollIm   = zeros<mat>(nRows, nCols);
    AppN.pitchIm  = zeros<mat>(nRows, nCols);
    AppN.yawIm    = zeros<mat>(nRows, nCols);

	// Calculate Newman's Approximation
    for (size_t ii = 0; ii < nRows; ++ii) {
        for (size_t jj = 0; jj < nCols; ++jj) {

            AppN.surgeRe(ii,jj) = (AppNauxiliar.surgeRe(ii) + AppNauxiliar.surgeRe(jj))*0.5;
            AppN.surgeIm(ii,jj) = (AppNauxiliar.surgeIm(ii) + AppNauxiliar.surgeIm(jj))*0.5;

            AppN.swayRe(ii,jj) = (AppNauxiliar.swayRe(ii) + AppNauxiliar.swayRe(jj))*0.5;
            AppN.swayIm(ii,jj) = (AppNauxiliar.swayIm(ii) + AppNauxiliar.swayIm(jj))*0.5;

            AppN.heaveRe(ii,jj) = (AppNauxiliar.heaveRe(ii) + AppNauxiliar.heaveRe(jj))*0.5;
            AppN.heaveIm(ii,jj) = (AppNauxiliar.heaveIm(ii) + AppNauxiliar.heaveIm(jj))*0.5;

            AppN.rollRe(ii,jj) = (AppNauxiliar.rollRe(ii) + AppNauxiliar.rollRe(jj))*0.5;
            AppN.rollIm(ii,jj) = (AppNauxiliar.rollIm(ii) + AppNauxiliar.rollIm(jj))*0.5;

            AppN.pitchRe(ii,jj) = (AppNauxiliar.pitchRe(ii) + AppNauxiliar.pitchRe(jj))*0.5;
            AppN.pitchIm(ii,jj) = (AppNauxiliar.pitchIm(ii) + AppNauxiliar.pitchIm(jj))*0.5;

            AppN.yawRe(ii,jj) = (AppNauxiliar.yawRe(ii) + AppNauxiliar.yawRe(jj))*0.5;
            AppN.yawIm(ii,jj) = (AppNauxiliar.yawIm(ii) + AppNauxiliar.yawIm(jj))*0.5;
            
        }
    }

	// Store the matrix as a complex matrix
	AppN.surge = cx_mat(AppN.surgeRe, AppN.surgeIm);
	AppN.sway = cx_mat(AppN.swayRe, AppN.swayIm);
	AppN.heave = cx_mat(AppN.heaveRe, AppN.heaveIm);
	AppN.roll = cx_mat(AppN.rollRe, AppN.rollIm);
	AppN.pitch = cx_mat(AppN.pitchRe, AppN.pitchIm);
	AppN.yaw = cx_mat(AppN.yawRe, AppN.yawIm);
	AppN.omega  = m_p12omega;
	AppN.period = 2 * arma::datum::pi / m_p12omega;
	

	// Get start to perform QTF matrix interpolation
	vec xfWamit = AppN.omega;

	const vec &t = envir.getTimeArray();
	int npts = t.size();

	vec fFiltro1 = zeros<mat>(npts, 1);

	double wmax = 2* datum::pi/envir.timeStep();
	double dw = 2* datum::pi/envir.timeTotal();

	vec fEspectro = linspace<vec>(0, wmax, npts);
	int cont = 0;
	

    for (size_t ii = 0; ii < npts; ++ii) {
        if (fEspectro(ii) >= xfWamit.min() && fEspectro(ii) <= xfWamit.max()) {
            fFiltro1(cont) = fEspectro(ii);            
            cont++;
        }
    }

    
	vec xfFiltro = fFiltro1(span(0, cont-1));

	// Performs QTF matrix interpolation
   	// SURGE
    interp2(xfWamit, xfWamit, AppN.surgeRe, xfFiltro, xfFiltro, AppNInterp.surgeRe, "linear");
    interp2(xfWamit, xfWamit, AppN.surgeIm, xfFiltro, xfFiltro, AppNInterp.surgeIm, "linear");

    // SWAY
    interp2(xfWamit, xfWamit, AppN.swayRe, xfFiltro, xfFiltro, AppNInterp.swayRe, "linear");
    interp2(xfWamit, xfWamit, AppN.swayIm, xfFiltro, xfFiltro, AppNInterp.swayIm, "linear");

    // HEAVE
    interp2(xfWamit, xfWamit, AppN.heaveRe, xfFiltro, xfFiltro, AppNInterp.heaveRe, "linear");
    interp2(xfWamit, xfWamit, AppN.heaveIm, xfFiltro, xfFiltro, AppNInterp.heaveIm, "linear");

    // ROLL
    interp2(xfWamit, xfWamit, AppN.rollRe, xfFiltro, xfFiltro, AppNInterp.rollRe, "linear");
    interp2(xfWamit, xfWamit, AppN.rollIm, xfFiltro, xfFiltro, AppNInterp.rollIm, "linear");

    // PITCH
    interp2(xfWamit, xfWamit, AppN.pitchRe, xfFiltro, xfFiltro, AppNInterp.pitchRe, "linear");
    interp2(xfWamit, xfWamit, AppN.pitchIm, xfFiltro, xfFiltro, AppNInterp.pitchIm, "linear");

    // YAW
    interp2(xfWamit, xfWamit, AppN.yawRe, xfFiltro, xfFiltro, AppNInterp.yawRe, "linear");
    interp2(xfWamit, xfWamit, AppN.yawIm, xfFiltro, xfFiltro, AppNInterp.yawIm, "linear");

	// Store the matrix as a complex matrix
	AppNInterp.surge = cx_mat(AppNInterp.surgeRe, AppNInterp.surgeIm);
	AppNInterp.sway = cx_mat(AppNInterp.swayRe, AppNInterp.swayIm);
	AppNInterp.heave = cx_mat(AppNInterp.heaveRe, AppNInterp.heaveIm);
	AppNInterp.roll = cx_mat(AppNInterp.rollRe, AppNInterp.rollIm);
	AppNInterp.pitch = cx_mat(AppNInterp.pitchRe, AppNInterp.pitchIm);
	AppNInterp.yaw = cx_mat(AppNInterp.yawRe, AppNInterp.yawIm);

	AppNInterp.omega  = xfFiltro;
	AppNInterp.period = 2* datum::pi / xfFiltro;

	// Generates complex amplitudes
	int nWaves = envir.numberOfWaveComponents();
	cx_mat AmpC1 = zeros<cx_mat>(nWaves, 1);

	for (unsigned int iWave = 0; iWave < nWaves; ++iWave)
	{

		AmpC1(iWave,0) = envir.waveElev_coef(0, 0, iWave);
		//const Wave& waveQTF(envir.getWave(iWave));
		//AmpC(iWave, 0) = waveQTF.amp() * std::exp( cx_double(0,1) * waveQTF.phase() );

	}
	cx_mat AmpC = AmpC1(span(0, cont - 1),0);
	
	mat timeseriesAppN (npts, 6, fill::zeros);

	cx_mat AppNsurge = AppNInterp.surge;
	cx_mat AppNsway = AppNInterp.sway;
	cx_mat AppNheave = AppNInterp.heave;
	cx_mat AppNroll = AppNInterp.roll;
	cx_mat AppNpitch = AppNInterp.pitch;
	cx_mat AppNyaw = AppNInterp.yaw;

	// Generates time series from IFFT
	vec Tsurge = fourier_coefficients_qtf(envir, AppNsurge, AmpC); // Generates the time series
	vec Tsway  = fourier_coefficients_qtf(envir, AppNsway, AmpC); // Generates the time series
	vec Theave = fourier_coefficients_qtf(envir, AppNheave, AmpC); // Generates the time series
	vec Troll  = fourier_coefficients_qtf(envir, AppNroll, AmpC); // Generates the time series
	vec Tpitch = fourier_coefficients_qtf(envir, AppNpitch, AmpC); // Generates the time series
	vec Tyaw   = fourier_coefficients_qtf(envir, AppNyaw, AmpC); // Generates the time series

	timeseriesAppN.col(0) = Tsurge;
	timeseriesAppN.col(1) = Tsway;
	timeseriesAppN.col(2) = Theave;
	timeseriesAppN.col(3) = Troll;
	timeseriesAppN.col(4) = Tpitch;
	timeseriesAppN.col(5) = Tyaw;



	return timeseriesAppN;

}

vec::fixed<6> Floater::hydrodynamicForce_2ndOrd_AppN(const ENVIR& envir, vec& m_p12omega, vec& m_p12beta, const cx_mat& m_p12surge, const cx_mat& m_p12sway, const cx_mat& m_p12heave, const cx_mat& m_p12roll, const cx_mat& m_p12pitch, const cx_mat& m_p12yaw)
{

	vec::fixed<6> force(fill::zeros); // Vector to store the force at the requested time
	
    if (m_hydroForce_WAMIT_DiffLoads.is_empty()) // Check if it is necessary to generate the QTF matrix by Newman Approximation
	{
		m_hydroForce_WAMIT_DiffLoads = QTF_ToTime_AppN(envir, m_p12omega, m_p12beta, m_p12surge, m_p12sway, m_p12heave, m_p12roll, m_p12pitch, m_p12yaw); // Generates time series from Newman's Approximation
		
	}

	uword ind1 = envir.getInd4interp1(); // Get the index of the requested time	
	force = m_hydroForce_WAMIT_DiffLoads.row(ind1).t() * envir.ramp(); // Force at the requested time

	if (envir.shouldInterp()) // Check if it is necessary to interpolate in time
	{
		uword ind2 = envir.getInd4interp2(); // Get the index of the requested time + 1
		const vec &t = envir.getTimeArray(); // Get vector time
		force += (m_hydroForce_WAMIT_DiffLoads.row(ind2).t() - m_hydroForce_WAMIT_DiffLoads.row(ind1).t()) * (envir.time() - t(ind1)) / (t(ind2) - t(ind1)) * envir.ramp(); // Interpolated the force at the requested time, it's a column vector

	}

	return force;
}