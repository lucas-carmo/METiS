#include "FOWT.h"
#include "IO.h"
#include <iostream>
#include <string>
#include <iterator>
#include <armadillo>

#include <algorithm> // For std::find


/*****************************************************
	Constructors
*****************************************************/
FOWT::FOWT() : m_extLinStiff(fill::zeros), m_mass(datum::nan),
m_disp(fill::zeros), m_vel(fill::zeros),
m_disp_sd(fill::zeros), m_vel_sd(fill::zeros)
{
	m_CoG.fill(datum::nan);

	// Default values for the filter in case the user has not specified any.
	m_filterSD_omega = 0;
	m_filterSD_zeta = 0.20;
}

/*****************************************************
	Setters
*****************************************************/
void FOWT::setHydroMode(const int hydroMode)
{
	m_hydroMode = hydroMode;
}

void FOWT::setAeroMode(const int aeroMode)
{
	m_aeroMode = aeroMode;
}

void FOWT::setMoorMode(const int moorMode)
{
	m_moorMode = moorMode;
}

void FOWT::setDoFs(std::array<bool, 6> &dofs)
{
	m_dofs = dofs;
}


void FOWT::setExtLinStiff(const mat::fixed<6, 6> &extLinStiff)
{
	m_extLinStiff = extLinStiff;
}

void FOWT::setExtLinDamp(const mat::fixed<6, 6> &extLinDamp)
{
	m_extLinDamp = extLinDamp;
}

void FOWT::setExtConstForce(const vec::fixed<6> &extConstForce)
{
	m_extConstForce = extConstForce;
}

void FOWT::setFilderSD(const double omega, const double zeta)
{
	m_filterSD_omega = omega;
	m_filterSD_zeta = zeta;
}

void FOWT::setDispSD(const vec::fixed<6> &disp_sd)
{
	m_disp_sd = disp_sd;
}


void FOWT::setFloater(Floater &floater)
{
	m_floater = floater;

	if (m_filterSD_omega < 0)
	{
		m_floater.setInstantSD(true);
	}
}

void FOWT::setRNA(RNA &rna)
{
	m_rna = rna;

	// Need to set the vertical distance between the hub and the CoG
	if (!arma::is_finite(m_floater.CoG()))
	{
		throw std::runtime_error("Need to set the floater CoG before calling FOWT::setRNA(RNA &rna)");
	}

	m_rna.setHubHeight2CoG(CoG().at(2));
}

void FOWT::setAddedMass_t0(const double density)
{
	m_floater.setAddedMass_t0(density);
}

void FOWT::setStiffnessMatrix(const double density, const double gravity)
{
	m_floater.setStiffnessMatrix(density, gravity);
}

void FOWT::evaluateQuantitiesAtBegin(const ENVIR &envir)
{
	if (envir.getTimeArray().size() == 0)
		throw std::runtime_error("Time array should be set before calling FOWT::evaluateQuantitiesAtBegin.");
	if (hydroMode() != 0)
		m_floater.evaluateQuantitiesAtBegin(envir, m_hydroMode);
}




/*****************************************************
	Getters
*****************************************************/
int FOWT::hydroMode() const
{
	return m_hydroMode;
}

int FOWT::aeroMode() const
{
	return m_aeroMode;
}

int FOWT::moorMode() const
{
	return m_moorMode;
}

bool FOWT::isDoFActive(const int index)
{
	if (index >= 0 && index < 6)
		return m_dofs.at(index);
	else
		throw std::runtime_error("Invalid index in FOWT::isDoFActive(const int index).");
}

double FOWT::filterSD_omega() const
{
	return m_filterSD_omega;
}

double FOWT::filterSD_zeta() const
{
	return m_filterSD_zeta;
}

vec::fixed<3> FOWT::CoG()
{
	// If CoG was not calculated yet, calculate it
	if (!arma::is_finite(m_CoG))
	{
		m_CoG = m_floater.CoG();
	}

	return m_CoG;
}

double FOWT::mass()
{
	// If m_mass was not calculated yet, calculate it
	if (!arma::is_finite(m_mass))
	{
		m_mass = m_floater.mass();
	}

	return m_mass;
}

vec::fixed<6> FOWT::disp() const
{
	return m_disp;
}

vec::fixed<6> FOWT::vel() const
{
	return m_vel;
}

vec::fixed<6> FOWT::disp_1stOrd() const
{
	return m_disp_1stOrd;
}

vec::fixed<6> FOWT::vel_1stOrd() const
{
	return m_vel_1stOrd;
}

vec::fixed<6> FOWT::disp_sd() const
{
	return m_disp_sd;
}

vec::fixed<6> FOWT::constForce() const
{
	return m_extConstForce;
}

std::string FOWT::printLinStiff() const
{
	std::string output = "\n\t";
	for (int ii = 0; ii < m_extLinStiff.n_rows; ++ii)
	{
		for (int jj = 0; jj < m_extLinStiff.n_cols; ++jj)
		{
			output = output + std::to_string(m_extLinStiff.at(ii, jj));
			(jj == 5) ? (output = output + "\n\t") : (output = output + " \t ; \t");
		}
	}

	return output + "\n";
}


std::string FOWT::printFloater() const
{
	std::string output = "";

	output = output + "\tMass:\t" + m_floater.printMass() + "\n";
	output = output + "\tCoG:\t" + m_floater.printCoG() + "\n";
	output = output + "\tInertia Matrix:\t" + m_floater.printInertia() + "\n";
	output = output + "\tMorison Elements:\n" + m_floater.printMorisonElements() + "\n";

	mat::fixed<6, 6> A = m_floater.addedMass_t0();
	output = output + "\tDimensional added Mass at initial position (same units as input):\n";
	for (int ii = 0; ii < 6; ++ii)
	{
		output = output + "\t\t";
		for (int jj = 0; jj < 6; ++jj)
		{
			output = output + std::to_string(A(ii, jj));

			(jj == 5) ? (output = output + '\n') : (output = output + " \t\t ; \t\t");
		}
	}

	mat::fixed<6, 6> K = m_floater.hydrostaticStiffness();
	output = output + "\tHydrostatic stiffness (same units as input):\n";
	for (int ii = 0; ii < 6; ++ii)
	{
		output = output + "\t\t";
		for (int jj = 0; jj < 6; ++jj)
		{
			output = output + std::to_string(K(ii, jj));

			(jj == 5) ? (output = output + '\n') : (output = output + " \t\t ; \t\t");
		}
	}

	output = output + '\n';

	return output;
}


std::string FOWT::printRNA() const
{
	std::string output = "";

	output = output + "Use Tip Loss:\t" + std::to_string(m_rna.useTipLoss()) + "\n";
	output = output + "Use Hub Loss:\t" + std::to_string(m_rna.useHubLoss()) + "\n";
	output = output + "Use Skew Correction:\t" + std::to_string(m_rna.useSkewCorr()) + "\n";
	output = output + "\tRotor Speed:\t" + std::to_string(m_rna.rotorSpeed()) + "\n";
	output = output + "\tRotor Tilt :\t" + std::to_string(m_rna.rotorTilt()) + "\n";
	output = output + "\tRotor Yaw:\t" + std::to_string(m_rna.rotorYaw()) + "\n";
	output = output + "\tNumBlades:\t" + std::to_string(m_rna.numBlades()) + "\n";
	output = output + "\tBlades Precone:\t" + std::to_string(m_rna.bladePrecone(0)) + "\n";
	output = output + "\tBlades Pitch:\t" + std::to_string(m_rna.bladePitch(0)) + "\t" + std::to_string(m_rna.bladePitch(1)) + "\t" + std::to_string(m_rna.bladePitch(2)) + "\n";
	output = output + "\tHub Radius:\t" + std::to_string(m_rna.hubRadius()) + "\n";
	output = output + "\tHub Height:\t" + std::to_string(m_rna.hubHeight()) + "\n";
	output = output + "\tOverhang:\t" + std::to_string(m_rna.overhang()) + "\n";
	output = output + "\tBlade aerodynamic properties:\n" + m_rna.printBladeAero();
	output = output + "\tAirfoils properties:\n" + m_rna.printAirfoils();

	return output;
}

std::string FOWT::printHydroMode() const
{
	return std::to_string(m_hydroMode);
}

std::string FOWT::printAeroMode() const
{
	return std::to_string(m_aeroMode);
}

std::string FOWT::printMoorMode() const
{
	return std::to_string(m_moorMode);
}


std::string FOWT::printDoF() const
{
	std::string output = "";
	for (int ii = 0; ii < 6; ++ii)
	{
		output += std::to_string(m_dofs[ii]) + ' ';
	}

	return output;
}

/*****************************************************
	QTF and AppN
*****************************************************/


// Read .12 WAMIT files (Total second-order forces by direct method)
void FOWT::readWAMIT_p12(const std::string &QTFPath)
{

	const bool flagFreq = false;

	// Read data
	std::ifstream file(QTFPath);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file");
    }

    std::vector<std::vector<double>> dataIn;
    std::string line;

    while (std::getline(file, line)) {

        std::istringstream iss(line);
        std::vector<double> row((std::istream_iterator<double>(iss)), std::istream_iterator<double>());
        dataIn.push_back(row);

    }
    file.close();
    

    if (flagFreq) {
        // Convert period (in s) to angular frequency (in rad/s)
        for (auto& row : dataIn) {

            row[0] = 2 * datum::pi / row[0];
            row[1] = 2 * datum::pi / row[1];

        }
    }



    // Sort data
    std::sort(dataIn.begin(), dataIn.end(), [](const std::vector<double>& a, const std::vector<double>& b) {

        return std::tie(a[4], a[0], a[1]) < std::tie(b[4], b[0], b[1]);

    });

	// Assign the different data to individual column vectors
    std::vector<double> wOrT1, wOrT2, beta1, beta2, I, mod, pha, reP, imP;
    for (const auto& row : dataIn) {

        wOrT1.push_back(row[0]); // angular frequency(w) or period(T) (row)
        wOrT2.push_back(row[1]); // angular frequency (w) or period (T) (column)
        beta1.push_back(row[2]); // row wave heading
        beta2.push_back(row[3]); // columnwave heading
        I.push_back(row[4]);     // DoF index
        mod.push_back(row[5]);   // Absolute value
        pha.push_back(row[6]);   // Phase
        reP.push_back(row[7]);	 // Real part
        imP.push_back(row[8]);	 // Imaginary part

    }
    

    // wOrT_unique will be used to indicate the row and column of each value in the QTF matrices
    // Unique frequencies
    std::vector<double> wOrT_unique = wOrT1;    
    wOrT_unique.insert(wOrT_unique.end(), wOrT2.begin(), wOrT2.end());
    std::sort(wOrT_unique.begin(), wOrT_unique.end());
    wOrT_unique.erase(std::unique(wOrT_unique.begin(), wOrT_unique.end()), wOrT_unique.end());

	// Unique incidences
    std::vector<double> betaQTF = beta1;    
    betaQTF.insert(betaQTF.end(), beta2.begin(), beta2.end());
    std::sort(betaQTF.begin(), betaQTF.end());
    betaQTF.erase(std::unique(betaQTF.begin(), betaQTF.end()), betaQTF.end());



    size_t nRows    = wOrT_unique.size();
    size_t nCols    = wOrT_unique.size();
	size_t nBetaQTF = betaQTF.size();	

	m_p12.surge = zeros<cx_mat>(nRows * nCols, nBetaQTF);
	m_p12.sway = zeros<cx_mat>(nRows * nCols, nBetaQTF);
	m_p12.heave = zeros<cx_mat>(nRows * nCols, nBetaQTF);
	m_p12.roll = zeros<cx_mat>(nRows * nCols, nBetaQTF);
	m_p12.pitch = zeros<cx_mat>(nRows * nCols, nBetaQTF);
	m_p12.yaw = zeros<cx_mat>(nRows * nCols, nBetaQTF);
	
	// loop to store each QTF matrix in an other matrix. Each column of this matrix represents a wave incidence
	for (size_t jj = 0; jj < nBetaQTF; ++jj) {

		// Create a matrix with size nRows x nCols for each field of W12
		m_p12auxiliar.surgeAmp = zeros<mat>(nRows, nCols);
		m_p12auxiliar.swayAmp  = zeros<mat>(nRows, nCols);
		m_p12auxiliar.heaveAmp = zeros<mat>(nRows, nCols);
		m_p12auxiliar.rollAmp  = zeros<mat>(nRows, nCols);
		m_p12auxiliar.pitchAmp = zeros<mat>(nRows, nCols);
		m_p12auxiliar.yawAmp   = zeros<mat>(nRows, nCols);
		
		m_p12auxiliar.surgePha = zeros<mat>(nRows, nCols);
		m_p12auxiliar.swayPha  = zeros<mat>(nRows, nCols);
		m_p12auxiliar.heavePha = zeros<mat>(nRows, nCols);
		m_p12auxiliar.rollPha  = zeros<mat>(nRows, nCols);
		m_p12auxiliar.pitchPha = zeros<mat>(nRows, nCols);
		m_p12auxiliar.yawPha   = zeros<mat>(nRows, nCols);
		
		m_p12auxiliar.surgeRe  = zeros<mat>(nRows, nCols);
		m_p12auxiliar.swayRe   = zeros<mat>(nRows, nCols);
		m_p12auxiliar.heaveRe  = zeros<mat>(nRows, nCols);
		m_p12auxiliar.rollRe   = zeros<mat>(nRows, nCols);
		m_p12auxiliar.pitchRe  = zeros<mat>(nRows, nCols);
		m_p12auxiliar.yawRe    = zeros<mat>(nRows, nCols);
		
		m_p12auxiliar.surgeIm  = zeros<mat>(nRows, nCols);
		m_p12auxiliar.swayIm   = zeros<mat>(nRows, nCols);
		m_p12auxiliar.heaveIm  = zeros<mat>(nRows, nCols);
		m_p12auxiliar.rollIm   = zeros<mat>(nRows, nCols);
		m_p12auxiliar.pitchIm  = zeros<mat>(nRows, nCols);
		m_p12auxiliar.yawIm    = zeros<mat>(nRows, nCols);

		m_p12auxiliar.omega    = zeros<mat>(nRows, 1);
		m_p12auxiliar.period   = zeros<mat>(nRows, 1);

    
    
		// if flagFreq == 1: the lower triangle is written
		// otherwise, it is the upper triangle
		for (size_t ii = 0; ii < wOrT1.size(); ++ii) {
			if (beta1[ii] == betaQTF[jj] && beta2[ii] == betaQTF[jj]) {
				size_t rr = std::distance(wOrT_unique.begin(), std::find(wOrT_unique.begin(), wOrT_unique.end(), wOrT1[ii])); // Row index
				size_t cc = std::distance(wOrT_unique.begin(), std::find(wOrT_unique.begin(), wOrT_unique.end(), wOrT2[ii])); // Column index

				if (I[ii] == 1){
					m_p12auxiliar.surgeAmp(rr, cc) = mod[ii];
					m_p12auxiliar.surgePha(rr, cc) = pha[ii];
					m_p12auxiliar.surgeRe(rr, cc) = reP[ii];
					m_p12auxiliar.surgeIm(rr, cc) = imP[ii];
				}
				
				if (I[ii] == 2){
					m_p12auxiliar.swayAmp(rr, cc) = mod[ii];
					m_p12auxiliar.swayPha(rr, cc) = pha[ii];
					m_p12auxiliar.swayRe(rr, cc) = reP[ii];
					m_p12auxiliar.swayIm(rr, cc) = imP[ii];
				}
				if (I[ii] == 3){
					m_p12auxiliar.heaveAmp(rr, cc) = mod[ii];
					m_p12auxiliar.heavePha(rr, cc) = pha[ii];
					m_p12auxiliar.heaveRe(rr, cc) = reP[ii];
					m_p12auxiliar.heaveIm(rr, cc) = imP[ii];
				}

				if (I[ii] == 4){
					m_p12auxiliar.rollAmp(rr, cc) = mod[ii];
					m_p12auxiliar.rollPha(rr, cc) = pha[ii];
					m_p12auxiliar.rollRe(rr, cc) = reP[ii];
					m_p12auxiliar.rollIm(rr, cc) = imP[ii];
				}

				if (I[ii] == 5){
					m_p12auxiliar.pitchAmp(rr, cc) = mod[ii];
					m_p12auxiliar.pitchPha(rr, cc) = pha[ii];
					m_p12auxiliar.pitchRe(rr, cc) = reP[ii];
					m_p12auxiliar.pitchIm(rr, cc) = imP[ii];
				}

				if (I[ii] == 6){
					m_p12auxiliar.yawAmp(rr, cc) = mod[ii];
					m_p12auxiliar.yawPha(rr, cc) = pha[ii];
					m_p12auxiliar.yawRe(rr, cc) = reP[ii];
					m_p12auxiliar.yawIm(rr, cc) = imP[ii];
				}

			}
			
		}

		m_p12auxiliar.period = wOrT_unique;


		
		// Symmetrize matrices
		if (flagFreq) {
			
			// Add the symmetric upper triangle to the QTF matrices.
			// It could be donne with a loop using 'eval', but it is too ugly and hard
			// to debug.Although much longer, the way used below is far simpler.
			m_p12auxiliar.surgeAmp = m_p12auxiliar.surgeAmp + trimatl(m_p12auxiliar.surgeAmp, -1).t();
			m_p12auxiliar.surgePha = m_p12auxiliar.surgePha - trimatl(m_p12auxiliar.surgePha, -1).t();
			m_p12auxiliar.surgeRe  = m_p12auxiliar.surgeRe  + trimatl(m_p12auxiliar.surgeRe, -1).t();
			m_p12auxiliar.surgeIm  = m_p12auxiliar.surgeIm  - trimatl(m_p12auxiliar.surgeIm, -1).t();
			
			m_p12auxiliar.swayAmp = m_p12auxiliar.swayAmp + trimatl(m_p12auxiliar.swayAmp, -1).t();
			m_p12auxiliar.swayPha = m_p12auxiliar.swayPha - trimatl(m_p12auxiliar.swayPha, -1).t();
			m_p12auxiliar.swayRe  = m_p12auxiliar.swayRe  + trimatl(m_p12auxiliar.swayRe, -1).t();
			m_p12auxiliar.swayIm  = m_p12auxiliar.swayIm  - trimatl(m_p12auxiliar.swayIm, -1).t();

			m_p12auxiliar.heaveAmp = m_p12auxiliar.heaveAmp + trimatl(m_p12auxiliar.heaveAmp, -1).t();
			m_p12auxiliar.heavePha = m_p12auxiliar.heavePha - trimatl(m_p12auxiliar.heavePha, -1).t();
			m_p12auxiliar.heaveRe  = m_p12auxiliar.heaveRe  + trimatl(m_p12auxiliar.heaveRe, -1).t();
			m_p12auxiliar.heaveIm  = m_p12auxiliar.heaveIm  - trimatl(m_p12auxiliar.heaveIm, -1).t();

			m_p12auxiliar.rollAmp = m_p12auxiliar.rollAmp + trimatl(m_p12auxiliar.rollAmp, -1).t();
			m_p12auxiliar.rollPha = m_p12auxiliar.rollPha - trimatl(m_p12auxiliar.rollPha, -1).t();
			m_p12auxiliar.rollRe  = m_p12auxiliar.rollRe  + trimatl(m_p12auxiliar.rollRe, -1).t();
			m_p12auxiliar.rollIm  = m_p12auxiliar.rollIm  - trimatl(m_p12auxiliar.rollIm, -1).t();

			m_p12auxiliar.pitchAmp = m_p12auxiliar.pitchAmp + trimatl(m_p12auxiliar.pitchAmp, -1).t();
			m_p12auxiliar.pitchPha = m_p12auxiliar.pitchPha - trimatl(m_p12auxiliar.pitchPha, -1).t();
			m_p12auxiliar.pitchRe  = m_p12auxiliar.pitchRe  + trimatl(m_p12auxiliar.pitchRe, -1).t();
			m_p12auxiliar.pitchIm  = m_p12auxiliar.pitchIm  - trimatl(m_p12auxiliar.pitchIm, -1).t();

			m_p12auxiliar.yawAmp = m_p12auxiliar.yawAmp + trimatl(m_p12auxiliar.yawAmp, -1).t();
			m_p12auxiliar.yawPha = m_p12auxiliar.yawPha - trimatl(m_p12auxiliar.yawPha, -1).t();
			m_p12auxiliar.yawRe  = m_p12auxiliar.yawRe  + trimatl(m_p12auxiliar.yawRe, -1).t();
			m_p12auxiliar.yawIm  = m_p12auxiliar.yawIm  - trimatl(m_p12auxiliar.yawIm, -1).t();

		} else {

			// Add the symmetric lower triangle to the QTF matrices.
			// It could be donne with a loop using 'eval', but it is too ugly and hard
			// to debug.Although much longer, the way used below is far simpler.

			m_p12auxiliar.surgeAmp = trimatu(m_p12auxiliar.surgeAmp, 1).t()    + m_p12auxiliar.surgeAmp;
			m_p12auxiliar.surgePha = -1*trimatu(m_p12auxiliar.surgePha, 1).t() + m_p12auxiliar.surgePha;
			m_p12auxiliar.surgeRe  = trimatu(m_p12auxiliar.surgeRe, 1).t()     + m_p12auxiliar.surgeRe;
			m_p12auxiliar.surgeIm  = -1*trimatu(m_p12auxiliar.surgeIm, 1).t()  + m_p12auxiliar.surgeIm;

			m_p12auxiliar.swayAmp = trimatu(m_p12auxiliar.swayAmp, 1).t()      + m_p12auxiliar.swayAmp;
			m_p12auxiliar.swayPha = -1 * trimatu(m_p12auxiliar.swayPha, 1).t() + m_p12auxiliar.swayPha;
			m_p12auxiliar.swayRe  = trimatu(m_p12auxiliar.swayRe, 1).t()       + m_p12auxiliar.swayRe;
			m_p12auxiliar.swayIm  = -1 * trimatu(m_p12auxiliar.swayIm, 1).t()  + m_p12auxiliar.swayIm;

			m_p12auxiliar.heaveAmp = trimatu(m_p12auxiliar.heaveAmp, 1).t()      + m_p12auxiliar.heaveAmp;
			m_p12auxiliar.heavePha = -1 * trimatu(m_p12auxiliar.heavePha, 1).t() + m_p12auxiliar.heavePha;
			m_p12auxiliar.heaveRe  = trimatu(m_p12auxiliar.heaveRe, 1).t()       + m_p12auxiliar.heaveRe;
			m_p12auxiliar.heaveIm  = -1 * trimatu(m_p12auxiliar.heaveIm, 1).t()  + m_p12auxiliar.heaveIm;

			m_p12auxiliar.rollAmp = trimatu(m_p12auxiliar.rollAmp, 1).t()      + m_p12auxiliar.rollAmp;
			m_p12auxiliar.rollPha = -1 * trimatu(m_p12auxiliar.rollPha, 1).t() + m_p12auxiliar.rollPha;
			m_p12auxiliar.rollRe  = trimatu(m_p12auxiliar.rollRe, 1).t()       + m_p12auxiliar.rollRe;
			m_p12auxiliar.rollIm  = -1 * trimatu(m_p12auxiliar.rollIm, 1).t()  + m_p12auxiliar.rollIm;

			m_p12auxiliar.pitchAmp = trimatu(m_p12auxiliar.pitchAmp, 1).t()      + m_p12auxiliar.pitchAmp;
			m_p12auxiliar.pitchPha = -1 * trimatu(m_p12auxiliar.pitchPha, 1).t() + m_p12auxiliar.pitchPha;
			m_p12auxiliar.pitchRe  = trimatu(m_p12auxiliar.pitchRe, 1).t()       + m_p12auxiliar.pitchRe;
			m_p12auxiliar.pitchIm  = -1 * trimatu(m_p12auxiliar.pitchIm, 1).t()  + m_p12auxiliar.pitchIm;

			m_p12auxiliar.yawAmp = trimatu(m_p12auxiliar.yawAmp, 1).t()      + m_p12auxiliar.yawAmp;
			m_p12auxiliar.yawPha = -1 * trimatu(m_p12auxiliar.yawPha, 1).t() + m_p12auxiliar.yawPha;
			m_p12auxiliar.yawRe  = trimatu(m_p12auxiliar.yawRe, 1).t()       + m_p12auxiliar.yawRe;
			m_p12auxiliar.yawIm  = -1 * trimatu(m_p12auxiliar.yawIm, 1).t()  + m_p12auxiliar.yawIm;

		}

		// Create a matrix with size nRows x nCols for each field of W12
		m_p12aux.surgeAmp = zeros<mat>(nRows, nCols);
		m_p12aux.swayAmp  = zeros<mat>(nRows, nCols);
		m_p12aux.heaveAmp = zeros<mat>(nRows, nCols);
		m_p12aux.rollAmp  = zeros<mat>(nRows, nCols);
		m_p12aux.pitchAmp = zeros<mat>(nRows, nCols);
		m_p12aux.yawAmp   = zeros<mat>(nRows, nCols);
		
		m_p12aux.surgePha = zeros<mat>(nRows, nCols);
		m_p12aux.swayPha  = zeros<mat>(nRows, nCols);
		m_p12aux.heavePha = zeros<mat>(nRows, nCols);
		m_p12aux.rollPha  = zeros<mat>(nRows, nCols);
		m_p12aux.pitchPha = zeros<mat>(nRows, nCols);
		m_p12aux.yawPha   = zeros<mat>(nRows, nCols);
		
		m_p12aux.surgeRe  = zeros<mat>(nRows, nCols);
		m_p12aux.swayRe   = zeros<mat>(nRows, nCols);
		m_p12aux.heaveRe  = zeros<mat>(nRows, nCols);
		m_p12aux.rollRe   = zeros<mat>(nRows, nCols);
		m_p12aux.pitchRe  = zeros<mat>(nRows, nCols);
		m_p12aux.yawRe    = zeros<mat>(nRows, nCols);
		
		m_p12aux.surgeIm  = zeros<mat>(nRows, nCols);
		m_p12aux.swayIm   = zeros<mat>(nRows, nCols);
		m_p12aux.heaveIm  = zeros<mat>(nRows, nCols);
		m_p12aux.rollIm   = zeros<mat>(nRows, nCols);
		m_p12aux.pitchIm  = zeros<mat>(nRows, nCols);
		m_p12aux.yawIm    = zeros<mat>(nRows, nCols);

		// It is necessary to "flip" the matrix
		for (size_t ii = 0; ii < nRows; ++ii) {
			for (size_t jj = 0; jj < nCols; ++jj) {

				m_p12aux.surgeAmp(ii,jj) = m_p12auxiliar.surgeAmp(nRows - 1 - ii, nCols - 1 - jj);
				m_p12aux.surgePha(ii,jj) = m_p12auxiliar.surgePha(nRows - 1 - ii, nCols - 1 - jj);
				m_p12aux.surgeRe(ii,jj) = m_p12auxiliar.surgeRe(nRows - 1 - ii, nCols - 1 - jj);
				m_p12aux.surgeIm(ii,jj) = m_p12auxiliar.surgeIm(nRows - 1 - ii, nCols - 1 - jj);

				m_p12aux.swayAmp(ii,jj) = m_p12auxiliar.swayAmp(nRows - 1 - ii, nCols - 1 - jj);
				m_p12aux.swayPha(ii,jj) = m_p12auxiliar.swayPha(nRows - 1 - ii, nCols - 1 - jj);
				m_p12aux.swayRe(ii,jj) = m_p12auxiliar.swayRe(nRows - 1 - ii, nCols - 1 - jj);
				m_p12aux.swayIm(ii,jj) = m_p12auxiliar.swayIm(nRows - 1 - ii, nCols - 1 - jj);

				m_p12aux.heaveAmp(ii,jj) = m_p12auxiliar.heaveAmp(nRows - 1 - ii, nCols - 1 - jj);
				m_p12aux.heavePha(ii,jj) = m_p12auxiliar.heavePha(nRows - 1 - ii, nCols - 1 - jj);
				m_p12aux.heaveRe(ii,jj) = m_p12auxiliar.heaveRe(nRows - 1 - ii, nCols - 1 - jj);
				m_p12aux.heaveIm(ii,jj) = m_p12auxiliar.heaveIm(nRows - 1 - ii, nCols - 1 - jj);

				m_p12aux.rollAmp(ii,jj) = m_p12auxiliar.rollAmp(nRows - 1 - ii, nCols - 1 - jj);
				m_p12aux.rollPha(ii,jj) = m_p12auxiliar.rollPha(nRows - 1 - ii, nCols - 1 - jj);
				m_p12aux.rollRe(ii,jj) = m_p12auxiliar.rollRe(nRows - 1 - ii, nCols - 1 - jj);
				m_p12aux.rollIm(ii,jj) = m_p12auxiliar.rollIm(nRows - 1 - ii, nCols - 1 - jj);

				m_p12aux.pitchAmp(ii,jj) = m_p12auxiliar.pitchAmp(nRows - 1 - ii, nCols - 1 - jj);
				m_p12aux.pitchPha(ii,jj) = m_p12auxiliar.pitchPha(nRows - 1 - ii, nCols - 1 - jj);
				m_p12aux.pitchRe(ii,jj) = m_p12auxiliar.pitchRe(nRows - 1 - ii, nCols - 1 - jj);
				m_p12aux.pitchIm(ii,jj) = m_p12auxiliar.pitchIm(nRows - 1 - ii, nCols - 1 - jj);

				m_p12aux.yawAmp(ii,jj) = m_p12auxiliar.yawAmp(nRows - 1 - ii, nCols - 1 - jj);
				m_p12aux.yawPha(ii,jj) = m_p12auxiliar.yawPha(nRows - 1 - ii, nCols - 1 - jj);
				m_p12aux.yawRe(ii,jj) = m_p12auxiliar.yawRe(nRows - 1 - ii, nCols - 1 - jj);
				m_p12aux.yawIm(ii,jj) = m_p12auxiliar.yawIm(nRows - 1 - ii, nCols - 1 - jj);

			}
		}
	

		// Store the matrix as a complex matrix
		m_p12aux.surge  = cx_mat(m_p12aux.surgeRe, m_p12aux.surgeIm);
		m_p12aux.sway   = cx_mat(m_p12aux.swayRe, m_p12aux.swayIm);
		m_p12aux.heave  = cx_mat(m_p12aux.heaveRe, m_p12aux.heaveIm);
		m_p12aux.roll   = cx_mat(m_p12aux.rollRe, m_p12aux.rollIm);
		m_p12aux.pitch  = cx_mat(m_p12aux.pitchRe, m_p12aux.pitchIm);
		m_p12aux.yaw    = cx_mat(m_p12aux.yawRe, m_p12aux.yawIm);


		


		// Stores each QTF matrix in an other matrix. Each column of this matrix represents a wave incidence

		cx_vec Colsurge = vectorise(m_p12aux.surge);
		cx_vec Colsway = vectorise(m_p12aux.sway);
		cx_vec Colheave = vectorise(m_p12aux.heave);
		cx_vec Colroll = vectorise(m_p12aux.roll);
		cx_vec Colpitch = vectorise(m_p12aux.pitch);
		cx_vec Colyaw = vectorise(m_p12aux.yaw);
		

		m_p12.surge.col(jj) = Colsurge;
		m_p12.sway.col(jj)  = Colsway;
		m_p12.heave.col(jj) = Colheave;
		m_p12.roll.col(jj)  = Colroll;
		m_p12.pitch.col(jj) = Colpitch;
		m_p12.yaw.col(jj)   = Colyaw;



	}

	m_p12.omega  = 2 * arma::datum::pi / reverse(m_p12auxiliar.period);
	m_p12.period = reverse(m_p12auxiliar.period);
	m_p12.beta   = betaQTF;


}



/*****************************************************
	Forces, acceleration, displacement, etc
*****************************************************/
// Update FOWT displacement, velocity, acceleration and any other necessary state.
void FOWT::update(const ENVIR &envir, const vec::fixed<12> &disp, const vec::fixed<12> &vel)
{
	m_disp_1stOrd = disp.rows(0, 5);
	m_vel_1stOrd = vel.rows(0, 5);
	m_disp = disp.rows(6, 11);
	m_vel = vel.rows(6, 11);

	if (m_filterSD_omega < 0)
	{
		m_disp_sd = disp;
		m_vel_sd = vel;
	}

	// Aqui tem que passar os deslocamentos com relacao ao CoG do floater. Calcular aqui mesmo baseado na posicao do centro de referencia de movimento 
	m_floater.update(envir, m_disp, m_vel, m_disp_1stOrd, m_vel_1stOrd, m_disp_sd, m_vel_sd);
}


// Evaluate (and update) the axis system that follows the slow position.
// Done so by applying a second-order low-pass filter to the displacement.
void FOWT::update_sd(const vec::fixed<6> &disp, const double dt)
{
	double wf = m_filterSD_omega;
	double zeta = m_filterSD_zeta;

	// If the filtering frequency is zero, there is nothing to be done here
	if (wf == 0)
	{
		return;
	}

	// If the filtering frequency is below zero, the slow drift position is actually 
	// equal to the instantaneous position, and it is updated in FOWT::update()
	if (wf > 0)
	{
		// Integrated using RK4. However, since the simulation time step must be small enough to 
		// capture wave frequency motions, which are way faster than slow drift motions, even an
		// explicit Euler method would work
		vec::fixed<6> acc_sd_k1 = -2 * wf*zeta*m_vel_sd - wf * wf * (m_disp_sd - m_disp);
		vec::fixed<6> vel_sd_k1 = acc_sd_k1 * dt;
		vec::fixed<6> disp_sd_k1 = m_vel_sd * dt;

		vec::fixed<6> acc_sd_k2 = -2 * wf*zeta*(m_vel_sd + vel_sd_k1 / 2) - wf * wf * (m_disp_sd + disp_sd_k1 / 2 - (m_disp + (disp - m_disp) / 2));
		vec::fixed<6> vel_sd_k2 = acc_sd_k2 * dt;
		vec::fixed<6> disp_sd_k2 = (m_vel_sd + vel_sd_k1 / 2) * dt;

		vec::fixed<6> acc_sd_k3 = -2 * wf*zeta*(m_vel_sd + vel_sd_k2 / 2) - wf * wf * (m_disp_sd + disp_sd_k2 / 2 - (m_disp + (disp - m_disp) / 2));
		vec::fixed<6> vel_sd_k3 = acc_sd_k3 * dt;
		vec::fixed<6> disp_sd_k3 = (m_vel_sd + vel_sd_k2 / 2) * dt;

		vec::fixed<6> acc_sd_k4 = -2 * wf*zeta*(m_vel_sd + vel_sd_k3) - wf * wf * (m_disp_sd + disp_sd_k3 - disp);
		vec::fixed<6> vel_sd_k4 = acc_sd_k3 * dt;
		vec::fixed<6> disp_sd_k4 = (m_vel_sd + vel_sd_k3) * dt;

		m_vel_sd += (vel_sd_k1 + 2 * vel_sd_k2 + 2 * vel_sd_k3 + vel_sd_k4) / 6;
		m_disp_sd += (disp_sd_k1 + 2 * disp_sd_k2 + 2 * disp_sd_k3 + disp_sd_k4) / 6;
	}
}



vec::fixed<12> FOWT::calcAcceleration(const ENVIR &envir)
{	
	mat::fixed<6, 6> addedMass(m_floater.addedMass_t0());
	IO::print2outLine(IO::OUTFLAG_ADDED_MASS_DIAG, addedMass.diag());

	mat::fixed<6, 6> inertiaMatrix = addedMass + m_floater.inertiaMatrix();	

	/*
		Calculate the force and acceleration due to first-order wave forces
	*/
	// The first order hydrodynamic forces (including hydrostatic and added mass) may be needed ahead
	// if second-order forces are required. Hence, assign the following two to avoid recalculating later.
	// The others are not necessary, as their evaluation is very fast.
	vec::fixed<6> hydrodynamicForce_1stOrd = m_floater.hydrodynamicForce_1stOrd(envir);
	vec::fixed<6> hydrostaticForce_1stOrd = m_floater.hydrostaticForce_stiffnessPart(true);
	vec::fixed<6> extLinDamp_1stOrd = extLinearDamp(true);
	vec::fixed<6> force_1stOrd = hydrodynamicForce_1stOrd + hydrostaticForce_1stOrd
							   + m_floater.hydrodynamicForce_drag1stOrd(envir) + mooringForce(true) + extLinDamp_1stOrd
							   + m_floater.hydrostaticForce_staticBuoyancy(envir.watDensity(), envir.gravity()) + weightForce(envir.gravity());

	// Calculate the acceleration only if at least one dof is activated 
	// (i.e. if at least one element of m_dofs is equal to 'true') 
	vec::fixed<6> acc_1stOrd(fill::zeros);
	if (std::find(m_dofs.begin(), m_dofs.end(), true) != m_dofs.end())
	{
		// Avoid coupling effects when a DoF is disabled and the others are not. 
		// Please note that the force was printed BEFORE this is done, in such a way 
		// that the full 6 component force vector is printed, even if the DoF is not active. 
		for (int ii = 0; ii < 6; ++ii)
		{
			if (!m_dofs[ii])
			{
				force_1stOrd.at(ii) = 0;

				// Even when the dof is disabled, it is necessary to provide a non-zero inertia 
				// entry at the main diagonal to avoid a singular matrix 
				inertiaMatrix.row(ii).zeros();
				inertiaMatrix.col(ii).zeros();
				inertiaMatrix.at(ii, ii) = 1;
			}
		}

		// Solve inertiaMatrix * acc = force 
		// Armadillo will throw its own exception if this computation fails. 
		acc_1stOrd = arma::solve(inertiaMatrix, force_1stOrd);
	}

	// TODO: N�o precisa calcular de novo se for s� 1a ordem, pq o resultado � igual ao de cima.

	// Calculate second order forces, if necessary
	vec::fixed<6> hydroForce_2ndOrd(fill::zeros);

	// Calculate second order forces with slender body approximation
	if (m_hydroMode == 2)
	{		
		m_floater.setNode1stAcc(acc_1stOrd);
		hydroForce_2ndOrd = m_floater.hydrodynamicForce_2ndOrd_slenderbody(envir, hydrodynamicForce_1stOrd + hydrostaticForce_1stOrd + extLinDamp_1stOrd);	
	}

	// Calculate second order forcer with Newman's approximation
	else if (m_hydroMode == 3)
	{		
		hydroForce_2ndOrd = m_floater.hydrodynamicForce_2ndOrd_AppN(envir, m_p12.omega, m_p12.beta, m_p12.surge, m_p12.sway, m_p12.heave, m_p12.roll, m_p12.pitch, m_p12.yaw);
	}
	
	// Calculate second order forcer with Full QTF
	else if (m_hydroMode == 4)
	{		
		// chamar funcao para ler .12d proveniente do arquivo IO
		//hydroForce_2ndOrd = hydrodynamicForce_2ndOrd_FullQTF(parametros);	
	}

	// Calculate the total force acting on the FOWT
	// Remember that the weight is already included in hydroForce_1stOrd, but the drag is not!
	vec::fixed<6> FaddMass = -(m_floater.addedMass(envir.watDensity(), m_hydroMode) - m_floater.addedMass_t0())* acc_1stOrd;
	vec::fixed<6> force = hydrodynamicForce_1stOrd + hydroForce_2ndOrd + m_floater.hydrostaticForce_stiffnessPart(false) 
						+ m_floater.hydrodynamicForce_dragTotal(envir) + mooringForce(false) + extLinearDamp(false)
						+ m_floater.hydrostaticForce_staticBuoyancy(envir.watDensity(), envir.gravity()) + weightForce(envir.gravity()) + aeroForce(envir) + FaddMass;
	IO::print2outLine(IO::OUTFLAG_TOTAL_FORCE, force);
	IO::print2outLine(IO::OUTFLAG_HD_FORCE, hydrodynamicForce_1stOrd + hydroForce_2ndOrd);

	// Calculate the acceleration only if at least one dof is activated
	// (i.e. if at least one element of m_dofs is equal to 'true')
	vec::fixed<6> acc_total(fill::zeros);
	if (std::find(m_dofs.begin(), m_dofs.end(), true) != m_dofs.end())
	{
		// Avoid coupling effects when a DoF is disabled and the others are not.
		// For doing so, set the calculated force to zero if the dof is deactivated.
		// Please note that the force was printed BEFORE this is done, in such a way
		// that the full 6 component force vector is printed, even if the DoF is not active.
		for (int ii = 0; ii < 6; ++ii)
		{
			if (!m_dofs[ii])
			{
				force.at(ii) = 0;

				// Even when the dof is disabled, it is necessary to provide a non-zero inertia
				// entry at the main diagonal to avoid a singular matrix
				inertiaMatrix.row(ii).zeros();
				inertiaMatrix.col(ii).zeros();
				inertiaMatrix.at(ii, ii) = 1;
			}
		}

		// Solve inertiaMatrix * acc = force
		// Armadillo will throw its own exception if this computation fails.
		acc_total = arma::solve(inertiaMatrix, force);		
	}
	IO::print2outLine(IO::OUTFLAG_HD_ADD_MASS_FORCE, FaddMass);
	return join_cols(acc_1stOrd, acc_total);
}

vec::fixed<6> FOWT::aeroForce(const ENVIR &envir)
{
	if (m_aeroMode == 1)
	{		
		return m_rna.aeroForce(envir, m_disp + join_cols(CoG(), vec::fixed<3> {0, 0, 0}), m_vel);
	}

	return vec::fixed<6> {0, 0, 0, 0, 0, 0};
}

vec::fixed<6> FOWT::mooringForce(bool flagUse1stOrd)
{
	vec::fixed<6> force{ 0, 0, 0, 0, 0, 0 };
	if (m_moorMode == 1)
	{
		if (flagUse1stOrd)
		{
			force = (-m_extLinStiff * m_disp_1stOrd + m_extConstForce);
		}
		else
		{
			force = (-m_extLinStiff * m_disp + m_extConstForce);
		}
	}

	IO::print2outLine(IO::OUTFLAG_MOOR_FORCE, force);

	return force;
}

vec::fixed<6> FOWT::weightForce(const double gravity)
{
	return vec::fixed<6> {0, 0, -gravity * mass(), 0, 0, 0};
}

vec::fixed<6> FOWT::extLinearDamp(bool flagUse1stOrd)
{
	vec::fixed<6> force{ 0, 0, 0, 0, 0, 0 };
	if (flagUse1stOrd)
	{
		force = -m_extLinDamp * m_vel_1stOrd;
	}
	else
	{
		force = -m_extLinDamp * m_vel;
	}

	//IO::print2outLine(IO::OUTFLAG_LINDAMP_FORCE, force);

	return force;
}
