#pragma once

#include <string>
#include <armadillo>
#include "Floater.h"
#include "RNA.h"
#include "ENVIR.h"


using namespace arma; // For armadillo classes

class FOWT
{
private:
	/*
	Physical members of the FOWT (floater, rotor nacelle assembly (RNA), tower, mooring lines, ...)
	*/
	Floater m_floater;
	//Tower m_tower;
	RNA m_rna;
	vec::fixed<6> m_extConstForce;
	mat::fixed<6,6> m_extLinStiff;
	mat::fixed<6, 6> m_extLinDamp;

	/*
	Specification of the analysis
	*/
	int m_hydroMode{ 0 };
	int m_aeroMode{ 0 };
	int m_moorMode{ 0 };

	// Flags to specify the active degrees of freedom
	std::array<bool, 6> m_dofs = { 1, 1, 1, 1, 1, 1 };

	/* 
	FOWT properties derived from its subsystems
	*/
	double m_mass;
	vec::fixed<3> m_CoG;

	/*
	FOWT condition
	*/
	vec::fixed<6> m_disp; // m_disp(0:2) = Position with respect to the initial CoG (i.e. CoG(t) - CoG(0)) --- m_disp(3:5) = Rotation with respect to initial configuration. For now, we are considering small rotations
	vec::fixed<6> m_vel;
	vec::fixed<6> m_disp_1stOrd; // Same thing, but first order 
	vec::fixed<6> m_vel_1stOrd;

	// Axis system that follows the mean and slow drift.
	// They are evaluated by filtering the instantaneous position with the following parameters.
	double m_filterSD_omega;
	double m_filterSD_zeta;
	vec::fixed<6> m_disp_sd;
	vec::fixed<6> m_vel_sd;

public:
	FOWT();


	/*****************************************************
		Setters
	*****************************************************/
	void setHydroMode(const int hydroMode);
	void setAeroMode(const int aeroMode);
	void setMoorMode(const int moorMode);
	void setDoFs(std::array<bool, 6> &dofs);

	void setExtConstForce(const vec::fixed<6> &extConstForce);
	void setExtLinStiff(const mat::fixed<6,6> &extLinStiff);
	void setExtLinDamp(const mat::fixed<6, 6> &extLinDamp);

	void setFilderSD(const double omega, const double zeta);
	void setDispSD(const vec::fixed<6> &disp_sd);

	void setFloater(Floater &floater);
	void setRNA(RNA &rna);

	void setAddedMass_t0(const double density);
	void setStiffnessMatrix(const double density, const double gravity);
	void evaluateQuantitiesAtBegin(const ENVIR &envir); // If the body is fixed, most of the computationally expensive calculations can be performed only once at the beginning of the simulation

	/*****************************************************
		Getters
	*****************************************************/
	int hydroMode() const;
	int aeroMode() const;
	int moorMode() const;

	bool isDoFActive(const int index);

	double filterSD_omega() const;
	double filterSD_zeta() const;

	vec::fixed<3> CoG();
	double mass();

	vec::fixed<6> disp() const;
	vec::fixed<6> vel() const;
	vec::fixed<6> disp_1stOrd() const;
	vec::fixed<6> vel_1stOrd() const;
	vec::fixed<6> disp_sd() const;
 
	vec::fixed<6> constForce() const;
	std::string printLinStiff() const;
	std::string printFloater() const;
	std::string printRNA() const;
	std::string printHydroMode() const;
	std::string printAeroMode() const;
	std::string printMoorMode() const;
	std::string printDoF() const;

	/*****************************************************
		QTF and AppN
	*****************************************************/

	struct p12Struct {

        mat surgeAmp, swayAmp, heaveAmp, rollAmp, pitchAmp, yawAmp;
        mat surgePha, swayPha, heavePha, rollPha, pitchPha, yawPha;
        mat surgeRe, swayRe, heaveRe, rollRe, pitchRe, yawRe;
        mat surgeIm, swayIm, heaveIm, rollIm, pitchIm, yawIm;
        cx_mat surge, sway, heave, roll, pitch, yaw;
        vec omega;
        vec period;

    };

	p12Struct m_p12;

    p12Struct m_p12auxiliar;

	const vector<double> betaQTF;
	void readWAMIT_p12(const std::string &QTFPath, const vector<double> &betaQTF);
	/*****************************************************
		Forces, acceleration, displacement, etc
	*****************************************************/
	vec::fixed<12> calcAcceleration(const ENVIR &envir);
	void update(const ENVIR &envir, const vec::fixed<12> &disp, const vec::fixed<12> &vel);
	void update_sd(const vec::fixed<6> &disp, const double dt);
	void update_sd(const vec::fixed<6> &disp, const double dt, const double wf, const double zeta);

	vec::fixed<6> aeroForce(const ENVIR &envir);
	vec::fixed<6> mooringForce(bool flagUse1stOrd);
	vec::fixed<6> weightForce(const double gravity);
	vec::fixed<6> extLinearDamp(bool flagUse1stOrd);
};