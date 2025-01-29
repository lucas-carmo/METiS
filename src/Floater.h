#pragma once

#include "MorisonElement.h"
#include "ENVIR.h"
#include "FOWT.h"
#include "Wave.h"

#include <string>
#include <vector> // For std::vector
#include <memory> // For std::unique_ptr
#include <armadillo> // Linear algebra library with usage similar to MATLAB


using namespace arma; // For armadillo classes

class Floater
{
private:

	
	double m_mass;
	bool m_instantSD{ false }; // If FOWT.m_filterSD_omega < 0, the slow position is actually equal to the instantaneous position, which has an important effect in the functions of Floater class.
	vec::fixed<3> m_CoG; // Coordinates of center of gravity.
	vec::fixed<6> m_inertia; // Moments and products of inertia. It is a 6x1 array. Actually, it is a symmetric 3x3 matrix, hence 3 elements are simply repeated.
	mat::fixed<6, 6> m_addedMass_t0;
	mat::fixed<6, 6> m_Khs; // Hydrostatic stiffness matrix
	double m_Vol;
	std::vector<std::unique_ptr<MorisonElement>> m_MorisonElements;	

	vec::fixed<6> m_disp; // Floater displacement (i.e. [instantaneous position] - [initial position])
	vec::fixed<6> m_disp_sd; // Keep track of the slow position of the floater as well
	vec::fixed<6> m_disp_1stOrd; // Same thing, but first order

	

public:
	Floater();

	/*****************************************************
		Overloaded operators
	*****************************************************/
	Floater& operator=(const Floater &floater);

	/*****************************************************
		Setters
	*****************************************************/
	void setMass(const double mass);
	void setInertia(const vec::fixed<6> &inertia);
	void setCoG(const vec::fixed<3> &cog);
	void setInstantSD(bool instantSD);
	void evaluateQuantitiesAtBegin(const ENVIR &envir, const int hydroMode);

	void addMorisonCirc(vec::fixed<3> &node1_coord, vec::fixed<3> &node2_coord, const double diam, const double CD, const double CM, const unsigned int numIntPoints,
						const double botDiam, const double topDiam, const double axialCD, const double axialCa, const bool botPressFlag);
	void addMorisonRect(vec::fixed<3> &node1_coord, vec::fixed<3> &node2_coord, vec::fixed<3> &node3_coord, const double diam_X, const double diam_Y,
	const double CD_X, const double CD_Y, const double CM_X, const double CM_Y, const unsigned int numIntPoints,
	const double botArea, const double topArea, const double axialCD, const double axialCa, const bool botPressFlag);

	/*****************************************************
		Getters
	*****************************************************/
	vec::fixed<3> CoG() const;
	double mass() const;
	mat::fixed<6,6> inertiaMatrix() const;
	mat::fixed<6, 6> addedMass_t0() const;
	mat::fixed<6, 6> hydrostaticStiffness() const;

	// Floater position, including rotations, with respect to the CoG position
	vec::fixed<6> CoGPos() const;
	vec::fixed<6> CoGPos_1stOrd() const;
	vec::fixed<6> CoGPos_sd() const;

	std::string printMass() const;
	std::string printInertia() const;
	std::string printCoG() const;
	std::string printMorisonElements() const;


	/*****************************************************
		Forces, acceleration, displacement, etc
	*****************************************************/
	void update(const ENVIR &envir, const vec::fixed<6> &FOWTdisp, const vec::fixed<6> &FOWTvel, const vec::fixed<6> &FOWTdisp_1stOrd, const vec::fixed<6> &FOWTvel_1stOrd, const vec::fixed<6> &FOWTdisp_SD, const vec::fixed<6> &FOWTvel_SD);
	void setNode1stAcc(const vec::fixed<6> &FOWTacc_1stOrd);
	void setAddedMass_t0(const double density);
	void setStiffnessMatrix(const double density, const double gravity);
	mat::fixed<6, 6> addedMass(const int hydroMode) const; // Calculated considering unitary density
	mat::fixed<6, 6> addedMass(const double density, const int hydroMode);	
	vec::fixed<6> hydrodynamicForce_1stOrd(const ENVIR &envir) const;
	vec::fixed<6> hydrodynamicForce_drag1stOrd(const ENVIR &envir) const;
	vec::fixed<6> hydrodynamicForce_dragTotal(const ENVIR &envir) const;
	vec::fixed<6> hydrodynamicForce_2ndOrd_slenderbody(const ENVIR &envir, const vec::fixed<6> &F_1stOrd) const;
	vec::fixed<6> hydrostaticForce_stiffnessPart(bool flagUse1stOrd) const;
	vec::fixed<6> hydrostaticForce_staticBuoyancy(const double rho, const double g) const;	

	struct AppNstruct {

        mat surgeAmp, swayAmp, heaveAmp, rollAmp, pitchAmp, yawAmp;
        mat surgePha, swayPha, heavePha, rollPha, pitchPha, yawPha;
        mat surgeRe, swayRe, heaveRe, rollRe, pitchRe, yawRe;
        mat surgeIm, swayIm, heaveIm, rollIm, pitchIm, yawIm;
        cx_mat surge, sway, heave, roll, pitch, yaw;
        vec omega;
        vec period;
		vec beta;

    };
	
	struct AppNstructAux {

        vec surgeAmp, swayAmp, heaveAmp, rollAmp, pitchAmp, yawAmp;
        vec surgePha, swayPha, heavePha, rollPha, pitchPha, yawPha;
        vec surgeRe, swayRe, heaveRe, rollRe, pitchRe, yawRe;
        vec surgeIm, swayIm, heaveIm, rollIm, pitchIm, yawIm;
        cx_double surge, sway, heave, roll, pitch, yaw;
        vec omega;
        vec period;

    };

    AppNstructAux AppNauxiliar;
    AppNstruct AppN;
	AppNstruct AppNInterp;

	mat m_hydroForce_WAMIT_DiffLoads; // Auxiliary variable that stores the time series. Column 1-6: surge-yaw

	mat timeseriesAppN;

	vec::fixed<6> hydrodynamicForce_2ndOrd_AppN(const ENVIR &envir, vec& m_p12omega, vec& m_p12beta, const cx_mat& m_p12surge, const cx_mat& m_p12sway, const cx_mat& m_p12heave, const cx_mat& m_p12roll, const cx_mat& m_p12pitch, const cx_mat& m_p12yaw);
	mat QTF_ToTime_AppN(const ENVIR& envir, vec& m_p12omega, vec& m_p12beta, const cx_mat& m_p12surge, const cx_mat& m_p12sway, const cx_mat& m_p12heave, const cx_mat& m_p12roll, const cx_mat& m_p12pitch, const cx_mat& m_p12yaw);
};

