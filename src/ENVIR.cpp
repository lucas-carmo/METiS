#include "ENVIR.h"
#include "IO.h"
#include "auxFunctions.h"

#include <iostream>
#include <vector>
#include <algorithm>    // std::binary_search
#include <utility> // For std::move
#include <complex> 

using namespace arma;

/*****************************************************
	Constructors
*****************************************************/
ENVIR::ENVIR()
{
	// Initialize with NaN so we can check whether they were defined later
	m_gravity = arma::datum::nan;
	m_watDepth = arma::datum::nan;
	m_timeStep = arma::datum::nan;
	m_timeTotal = arma::datum::nan;
	m_timeRamp = arma::datum::nan;	
}

void ENVIR::setCurrentTime(const double time)
{
	m_time = time;
}

/*****************************************************
	Setters
*****************************************************/
void ENVIR::setTimeStep(const double timeStep)
{
	m_timeStep = timeStep;
}

void ENVIR::setPrintStep(const double printStep)
{
	m_printStep = printStep;
}

void ENVIR::setTimeTotal(const double timeTotal)
{
	m_timeTotal = timeTotal;
}

void ENVIR::setTimeRamp(const double timeRamp)
{
	m_timeRamp = timeRamp;
}

void ENVIR::setGravity(const double gravity)
{
	m_gravity = gravity;
}

void ENVIR::setWatDens(const double watDens)
{
	m_watDens = watDens;
}

void ENVIR::setWatDepth(const double watDepth)
{
	m_watDepth = watDepth;
}

void ENVIR::setWaveStret(const unsigned int waveStret)
{
	m_waveStret = waveStret;
}

void ENVIR::setAirDens(const double airDens)
{
	m_airDens = airDens;
}

void ENVIR::setWindRefVel(const double windRefVel)
{
	m_windRefVel = windRefVel;
}

void ENVIR::setWindDir(const double windDir)
{
	m_windDir = windDir;
}

void ENVIR::setWindRefHeight(const double windRefHeight)
{
	m_windRefHeight = windRefHeight;
}

void ENVIR::setWindExp(const double windExp)
{
	m_windExp = windExp;
}

void ENVIR::addNode(const unsigned int nodeID, const double nodeCoordX, const double nodeCoordY, const double nodeCoordZ)
{
	if (m_nodesID.size() != 0) // If this is not the first node that will be added to m_nodesID
	{
		if (nodeID <= m_nodesID.back()) // Then verify if its ID is larger than the previous one, thus garanteeing that m_nodesID is in ascending order (this is needed to use binary search to find nodes IDs)
		{
			throw std::runtime_error("Nodes must be provided in ascending order, but node " + std::to_string(nodeID) + " was specified after node " + std::to_string(m_nodesID.back()) + ". In ENVIR::addNode");
		}
	}

	m_nodesID.push_back(nodeID);
	m_nodesCoord.push_back(vec::fixed<3> {nodeCoordX, nodeCoordY, nodeCoordZ});
}


void ENVIR::addRegularWave(const std::string &waveType, const double height, const double freqORperiod, const double direction, const double phase)
{
	// Check whether the water depth was defined
	if (!is_finite(m_watDepth))
	{
		throw std::runtime_error("The water depth must be specified before the waves. In ENVIR::addRegularWave.");
	}

	// Check whether the acceleration of gravity was defined
	if (!is_finite(m_gravity))
	{
		throw std::runtime_error("The acceleration of gravity must be specified before the waves. In ENVIR::addRegularWave.");
	}

	m_wave.push_back(Wave(waveType, height, freqORperiod, direction, phase, m_watDepth, m_gravity));
}


void ENVIR::addJonswap(const double Hs, const double Tp, const double gamma, const double direction, const double wlow, const double whigh, const int numberOfRegularWaves, const double dwMax)
{
	arma::vec w{ 0 };
	arma::vec dw{ 0 };

	// If the number of frequencies is not specified, the frequency resolution is determined by the total simulation time.
	// A Harmonic Deterministic Amplitude Scheme (HDAS) is used.
	// See Merigaud, A. and Ringwood, John V. - Free-Surface Time-Series Generation for Wave Energy Applications - IEEE J. of Oceanic Eng. - 2018
	if (!is_finite(m_timeTotal) || !is_finite(m_timeStep) || !is_finite(m_timeRamp))
	{
		throw std::runtime_error("The total simulation time and time step must be specified before the irregular waves. In ENVIR::addJonswap.");
	}

	double dw0 = 2 * arma::datum::pi / m_timeTotal;

	if (numberOfRegularWaves <= 0)
	{
		if (m_wave.size() != 0)
		{
			throw std::runtime_error("Wave options that use IFFT to evaluate wave kinematics, such as JONSWAP without specifying the number of components or externally generated wave elevation, can not be specified with other waves. In ENVIR::addJonswap.");
		}

		w = arma::regspace(0, dw0, 2 * arma::datum::pi /m_timeStep);
		dw = dw0 * arma::ones<arma::vec>(w.size());
		m_flagIFFT = true;
	}

	// Otherwise, use the Equal Area Deterministic Amplitude Scheme (EADAS)
	else
	{
		if (getFlagIFFT())
		{
			throw std::runtime_error("Wave options that use IFFT to evaluate wave kinematics, such as JONSWAP without specifying the number of components or externally generated wave elevation, can not be specified with other waves. In ENVIR::addJonswap.");
		}

		double dE = Hs * Hs / 16. / numberOfRegularWaves; // m0 = (Hs/4)^2 divided by the number of strips, which is equal to the number of wave components		

		// The final size of w (and dw) is probably different from numberOfRegularWaves due to the limitations given by dwMax
		// and the high cut-off frequency (whigh). There can be both fewer or more frequencies than what was originally asked by the used.
		// In any case, this first guess of vector size is used to reduce the calls of insert_rows.
		w = arma::zeros<arma::vec>(numberOfRegularWaves);
		dw = w;

		// Find the frequencies that bound each strip		
		//
		// The objective is to find dw such that the area between w[ii] and w[ii-1] is equal to dE.
		// By considering a trapezoidal integration, this is equivalent to solving
		// f(x) = x * (S(w_a+x) + S(w_a)) - 2*dE = 0
		// For doing so, the bissection method is used (it would be nice to implement a general function for solving this kind o problem instead of implementing it directly here)
		// However, this leads to a very coarse discretization in the parts of the spectrum where there is little energy.
		// Hence, if the x found is too big, i.e. larger than dwMax, it is break down in steps equal do dwMax
		double w_a = wlow;
		int ii = 0;
		double E_sum = 0;
		while (true)
		{
			// Need to bracket the solution first
			double a = 0; // It is clear that f(0) has a negative value
			double b = dw0; // First guess for a positive value is the dw obtained from a harmonic scheme

			// Increase b until f(b) > 0			
			while ((b * (JONSWAP(w_a + b, Tp, Hs, gamma) + JONSWAP(w_a, Tp, Hs, gamma)) - 2 * dE) < 0)
			{
				b = b + dw0; // Use steps of dw0
			}

			double x_j = (a + b) / 2;
			double eps4dw = 1e-6*dE;
			while (std::abs(x_j * (JONSWAP(w_a + x_j, Tp, Hs, gamma) + JONSWAP(w_a, Tp, Hs, gamma)) - 2 * dE) > eps4dw)
			{
				// Test with limit a
				if ((a * (JONSWAP(w_a + a, Tp, Hs, gamma) + JONSWAP(w_a, Tp, Hs, gamma)) - 2 * dE)
					* (x_j * (JONSWAP(w_a + x_j, Tp, Hs, gamma) + JONSWAP(w_a, Tp, Hs, gamma)) - 2 * dE) < 0)
				{
					b = x_j;
				}

				// Test with limit b
				else if ((b * (JONSWAP(w_a + b, Tp, Hs, gamma) + JONSWAP(w_a, Tp, Hs, gamma)) - 2 * dE)
					* (x_j * (JONSWAP(w_a + x_j, Tp, Hs, gamma) + JONSWAP(w_a, Tp, Hs, gamma)) - 2 * dE) < 0)
				{
					a = x_j;
				}

				// It is very unlikely that this product will be exactly zero, which means that x_j is the exact solution, but this will be covered anywawy
				else
				{
					break;
				}

				// New guess for next step
				x_j = (a + b) / 2;
			}

			if (x_j == 0)
			{
				throw std::runtime_error("Something went wrong in the frequency discretization in ENVIR::addJonswap for frequency w = " + std::to_string(w_a) + "rad/s");
			}

			// x_j is the frequency discretization that leads to an energy dE between w_a and w_a+x_j		
			if (x_j > dwMax)
			{
				x_j = dwMax;
			}
			if (x_j < 0)
			{
				break;
			}

			dw(ii) = x_j;
			w(ii) = w_a + 0.5*dw(ii); // The frequency is taken at the midpoint between the boundaries
			E_sum += dw(ii) * JONSWAP(w(ii), Tp, Hs, gamma);

			if ((w(ii) >= whigh) || (E_sum >= Hs * Hs / 16))
			{
				// There may be cases that ii is smaller than numberOfRegularWaves.
				// These extra elements are excluded from the vector.
				if (ii < w.size() - 1)
				{
					w.shed_rows(ii + 1, w.size() - 1);
				}
				break;
			}

			w_a += dw(ii);
			ii = ii + 1;

			if (ii > w.size() - 1)
			{
				w.insert_rows(ii, 1);
				dw.insert_rows(ii, 1);
			}

		}
	}

	// Wave parameters that will be calculated using the JONSWAP spectrum
	double height(0);
	double phase(0);

	// Loop to create the waves
	int ii = 0;
	double Sw = 0;
	for (int ii = 0; ii < w.size(); ++ii)
	{	
		Sw = 0;
		if (w.at(ii) > wlow && w.at(ii) < whigh)
		{
			Sw = JONSWAP(w.at(ii), Tp, Hs, gamma);
		}

		height = 2 * std::sqrt(2 * Sw * dw(ii));
		phase = 360 * arma::randu(1, 1).at(0, 0);
		addRegularWave("TRWave", height, 2 * arma::datum::pi / w.at(ii), direction, phase);
	}
}

void ENVIR::addWaveElevSeries(const std::string &elevFlPath, const double direction)
{
	// Open input file
	std::ifstream elevFl;
	elevFl.open(elevFlPath);
	if (!elevFl)
	{
		throw std::runtime_error("Unable to open file " + elevFlPath + " for reading.");
	}

	// Check if other waves were already specified 
	if (m_wave.size() != 0)
	{
		throw std::runtime_error("Wave options that use IFFT to evaluate wave kinematics, such as JONSWAP without specifying the number of components or externally generated wave elevation, can not be specified with other waves. In ENVIR::addWaveElevSeries.");
	}

	// Count number of lines in file
	std::string line;
	int nlines{ 0 };
	while (std::getline(elevFl, line))
	{
		++nlines;
	}

	elevFl.clear(); // Clear and return to the beginning
	elevFl.seekg(0);

	// Create vectors to store the time and wave elevation, which must be arranged as two columns in the input file
	std::vector<double> time(nlines), elev(nlines);
	int ii = 0;
	while (std::getline(elevFl, line))
	{
		std::vector<std::string> input = stringTokenize(line, " \t");
		if (input.size() != 2)
		{
			throw std::runtime_error("Unable to read time series of wave elevation from file " + elevFlPath + ". Wrong number of parameters in line " + std::to_string(ii+1) + ".");
		}
		
		time.at(ii) = string2num<double>(input.at(0));
		elev.at(ii) = string2num<double>(input.at(1));

		++ii;
	}
	
	if (time.at(0) != 0)
	{
		throw std::runtime_error("Time series of wave elevation provided in " + elevFlPath + " should start at t = 0.");
	}

	// Frequency vector 
	double dw = 2 * arma::datum::pi / time.at(time.size()-1);
	vec w = arma::regspace(0, dw, (time.size()-1) * dw);
	m_flagIFFT = true;

	// FFT and assignment of the computed values to the vector of waves
	cx_vec amp(mkl_fft_real(elev));
	amp *= 1. / time.size(); 
	amp.rows(std::floor(amp.size()/2), amp.size()-1).zeros();
	amp.rows(1, amp.size() - 1) *= 2;
	for (int ii = 0; ii < w.size(); ++ii)
	{
		addRegularWave("TRWave", 2*std::abs(amp.at(ii)), 2 * arma::datum::pi / w.at(ii), direction, -std::arg(amp.at(ii)) * 180. / arma::datum::pi);
	}

	// Make sure that the time step and total simulation time are equal to the wave elevation file
	double newTimeStep = (time.at(time.size() - 1) - time.at(0)) / (time.size()-1);
	IO::print2log("Setting time step to " + std::to_string(newTimeStep) + "s to match the one from the the wave elevation file. Previous value was " + std::to_string(m_timeStep) + ".");
	m_timeStep = newTimeStep;

	IO::print2log("Setting total simulation time to " + std::to_string(time.at(time.size() - 1)) + "s to match the one from the the wave elevation file. Previous value was " + std::to_string(m_timeTotal) + ".");
	m_timeTotal = time.at(time.size() - 1);
}


void ENVIR::addWaveProbe(const unsigned int ID)
{
	// Check whether nodes were specified
	if (this->isNodeEmpty())
	{
		throw std::runtime_error("Nodes should be specified before adding wave locations. In: ENVIR::addWaveProbe");
	}

	std::vector<unsigned int>::const_iterator iter = std::find(m_waveProbeID.begin(), m_waveProbeID.end(), ID); // Check if this wave probe was already included

	if (iter == m_waveProbeID.end())
	{
		m_waveProbe.push_back(this->getNode(ID)); // Get the node coordinate and add it to m_waveProbe
		m_waveProbeID.push_back(ID);
	}
}

void ENVIR::evaluateWaveKinematics()
{
	typedef std::vector<cx_double> cx_stdvec;

	m_timeArray = arma::linspace(0, m_timeTotal, std::ceil(m_timeTotal / m_timeStep) + 1);
	m_timeRampArray = zeros(m_timeArray.size(), 1);
	for (int it = 0; it < m_timeArray.size(); ++it)
	{
		m_timeRampArray.at(it) = ramp(m_timeArray.at(it));
	}

	vec w{ zeros(m_wave.size()) };
	for (int iWave = 0; iWave < m_wave.size(); ++iWave)
	{
		w.at(iWave) = m_wave.at(iWave).angFreq();
	}
	
	if (IO::isOutputActive(IO::OUTFLAG_WAVE_ELEV))
	{
		m_waveElevArray = zeros(m_timeArray.size(), m_waveProbeID.size());
		for (int iProbe = 0; iProbe < m_waveProbeID.size(); ++iProbe)
		{
			cx_mat amp(m_wave.size(), 1); // Complex amplitude			
			for (int iWave = 0; iWave < m_wave.size(); ++iWave)
			{
				amp.at(iWave) = waveElev_coef(m_waveProbe.at(iProbe).at(0), m_waveProbe.at(iProbe).at(1), iWave);
			}

			if (getFlagIFFT())
			{
				m_waveElevArray.col(iProbe) = m_wave.size() * mkl_ifft_real(amp);
			}
			else
			{
				for (int it = 0; it < m_timeArray.size(); ++it)
				{
					cx_vec sinCos{ cos(w * m_timeArray.at(it)), sin(w * m_timeArray.at(it)) };
					m_waveElevArray.at(it, iProbe) = arma::accu(real(amp % sinCos));
				}
			}

			m_waveElevArray.col(iProbe) %= m_timeRampArray;
		}
	}

	if (IO::isOutputActive(IO::OUTFLAG_WAVE_VEL))
	{
		m_waveVel1stArray_x = zeros(m_timeArray.size(), m_waveProbeID.size());
		m_waveVel1stArray_z = zeros(m_timeArray.size(), m_waveProbeID.size());
		m_waveVel1stArray_y = zeros(m_timeArray.size(), m_waveProbeID.size());

		for (int iProbe = 0; iProbe < m_waveProbeID.size(); ++iProbe)
		{
			cx_mat amp(m_wave.size(), 3); // Complex amplitude
			for (int iWave = 0; iWave < m_wave.size(); ++iWave)
			{
				amp.row(iWave) = u1_coef(m_waveProbe.at(iProbe).at(0), m_waveProbe.at(iProbe).at(1), m_waveProbe.at(iProbe).at(2), iWave).st();
			}
			
			if (getFlagIFFT())
			{
				for (int idof = 0; idof < 3; ++idof)
				{
					cx_stdvec aux = conv_to<cx_stdvec>::from(amp.col(0));
					m_waveVel1stArray_x.col(iProbe) = m_wave.size() * mat(mkl_ifft_real(aux));

					aux = conv_to<cx_stdvec>::from(amp.col(1));
					m_waveVel1stArray_y.col(iProbe) = m_wave.size() * mat(mkl_ifft_real(aux));

					aux = conv_to<cx_stdvec>::from(amp.col(2));
					m_waveVel1stArray_z.col(iProbe) = m_wave.size() * mat(mkl_ifft_real(aux));
				}
			}
			else
			{
				for (int it = 0; it < m_timeArray.size(); ++it)
				{
					cx_vec sinCos{ cos(w * m_timeArray.at(it)), sin(w * m_timeArray.at(it)) };
					m_waveVel1stArray_x.at(it, iProbe) = arma::accu(real(amp.col(0) % sinCos));
					m_waveVel1stArray_y.at(it, iProbe) = arma::accu(real(amp.col(1) % sinCos));
					m_waveVel1stArray_z.at(it, iProbe) = arma::accu(real(amp.col(2) % sinCos));
				}
			}

			m_waveVel1stArray_x.col(iProbe) %= m_timeRampArray;
			m_waveVel1stArray_y.col(iProbe) %= m_timeRampArray;
			m_waveVel1stArray_z.col(iProbe) %= m_timeRampArray;
		}
	}
}

/*****************************************************
	Getters
*****************************************************/
double ENVIR::timeStep() const
{
	return m_timeStep;
}

double ENVIR::printStep() const
{
	return m_printStep;
}

double ENVIR::timeTotal() const
{
	return m_timeTotal;
}

double ENVIR::time() const
{
	return m_time;
}

double ENVIR::gravity() const
{
	return m_gravity;
}

double ENVIR::watDensity() const
{
	return m_watDens;
}

double ENVIR::watDepth() const
{
	return m_watDepth;
}

unsigned int ENVIR::waveStret() const
{
	return m_waveStret;
}

double ENVIR::airDensity() const
{
	return m_airDens;
}

double ENVIR::windRefVel() const
{
	return m_windRefVel;
}

double ENVIR::windRefHeight() const
{
	return m_windRefHeight;
}

double ENVIR::windDir() const
{
	return m_windDir;
}

double ENVIR::windExp() const
{
	return m_windExp;
}

unsigned int ENVIR::numberOfWaveComponents() const
{
	return m_wave.size();
}

const Wave& ENVIR::getWave(unsigned int waveIndex) const
{
	return m_wave.at(waveIndex);
}

double ENVIR::waveElevAtProbe(const unsigned int ID) const
{
	// This should never occur in production code, but may be useful for debugging
	if (m_timeArray.size() == 0)
	{
		return waveElev(m_waveProbe[ID][0], m_waveProbe[ID][1]);
	}
	else
	{
		uword ind1 = m_ind4interp1;
		double eta = m_waveElevArray.at(ind1, ID);
		if (shouldInterp())
		{
			uword ind2 = getInd4interp2();
			eta += (m_waveElevArray.at(ind2, ID) - m_waveElevArray.at(ind1, ID)) * (m_time - m_timeArray.at(ind1)) / (m_timeArray.at(ind2) - m_timeArray.at(ind1));
		}
		return eta;
	}
}

vec::fixed<3> ENVIR::waveVelAtProbe(const unsigned int ID) const
{
	// This should never occur in production code, but may be useful for debugging
	if (m_timeArray.size() == 0)
	{		
		return u1(m_waveProbe.at(ID), 0);
	}
	else
	{
		uword ind1 = m_ind4interp1;
		double vel_x = m_waveVel1stArray_x.at(ind1, ID);
		double vel_y = m_waveVel1stArray_y.at(ind1, ID);
		double vel_z = m_waveVel1stArray_z.at(ind1, ID);
		if (shouldInterp())
		{
			uword ind2 = getInd4interp2();
			vel_x += (m_waveVel1stArray_x.at(ind2, ID) - m_waveVel1stArray_x.at(ind1, ID)) * (m_time - m_timeArray.at(ind1)) / (m_timeArray.at(ind2) - m_timeArray.at(ind1));
			vel_y += (m_waveVel1stArray_y.at(ind2, ID) - m_waveVel1stArray_y.at(ind1, ID)) * (m_time - m_timeArray.at(ind1)) / (m_timeArray.at(ind2) - m_timeArray.at(ind1));
			vel_z += (m_waveVel1stArray_z.at(ind2, ID) - m_waveVel1stArray_z.at(ind1, ID)) * (m_time - m_timeArray.at(ind1)) / (m_timeArray.at(ind2) - m_timeArray.at(ind1));
		}
		return { vel_x, vel_y, vel_z };
	}
}

bool ENVIR::getFlagIFFT() const
{
	return m_flagIFFT;
}

bool ENVIR::shouldInterp() const
{
	return m_shouldInterp;
}

uword ENVIR::getInd4interp1() const
{
	return m_ind4interp1;
}

uword ENVIR::getInd4interp2() const
{
	return m_ind4interp2;
}


const vec& ENVIR::getTimeArray() const
{
	return m_timeArray;
}

const vec& ENVIR::getRampArray() const
{
	return m_timeRampArray;
}

/*****************************************************
	Printing
*****************************************************/
std::string ENVIR::printTimeStep() const
{
	return std::to_string(m_timeStep);
}

std::string ENVIR::printPrintStep() const
{
	return std::to_string(m_printStep);
}

std::string ENVIR::printTimeTotal() const
{
	return std::to_string(m_timeTotal);
}

std::string ENVIR::printTimeRamp() const
{
	return std::to_string(m_timeRamp);
}

std::string ENVIR::printNodes() const
{
	std::string output = "";
	for (int ii = 0; ii < m_nodesID.size(); ++ii)
	{
		output = output + "( " + std::to_string(m_nodesID.at(ii)) +
			", " + std::to_string(m_nodesCoord.at(ii)(0)) +
			", " + std::to_string(m_nodesCoord.at(ii)(1)) +
			", " + std::to_string(m_nodesCoord.at(ii)(2)) + " )\n";
	}
	return output;
}

std::string ENVIR::printGrav() const
{
	return std::to_string(m_gravity);
}

std::string ENVIR::printWatDens() const
{
	return std::to_string(m_watDens);
}

std::string ENVIR::printWatDepth() const
{
	return std::to_string(m_watDepth);
}

std::string ENVIR::printWave() const
{
	std::string output = "";
	output = output + "Number of waves: " + std::to_string(m_wave.size()) + "\n";
	output = output + "Wave characteristics (printing up to 100 components): " + std::to_string(m_wave.size()) + "\n";
	for (int ii = 0; ii < m_wave.size(); ++ii)
	{
		if (ii == 100)
			break;
		output = output + "Wave #" + std::to_string(ii) + "\n";
		output = output + "Height: " + std::to_string(m_wave.at(ii).height()) + "\n";
		output = output + "Period: " + std::to_string(m_wave.at(ii).period()) + "\n";
		output = output + "Wave number: " + std::to_string(m_wave.at(ii).waveNumber()) + "\n";
		output = output + "Length: " + std::to_string(m_wave.at(ii).length()) + "\n";
		output = output + "Direction: " + std::to_string(m_wave.at(ii).direction()) + "\n";
		output = output + "Phase: " + std::to_string(m_wave.at(ii).phase()) + "\n\n";
	}
	return output;
}

std::string ENVIR::printWaveProbe() const
{
	std::string output = "";
	for (int ii = 0; ii < m_waveProbe.size(); ++ii)
	{
		output = output + "Location #" + std::to_string(ii) + ": (" + std::to_string(m_waveProbe.at(ii).at(0))
			+ "," + std::to_string(m_waveProbe.at(ii).at(1)) + "," + std::to_string(m_waveProbe.at(ii).at(2)) + ")\n";
	}
	return output;
}


void ENVIR::printWaveCharact() const
{
	for (int ii = 0; ii < m_waveProbe.size(); ++ii)
	{
			if (IO::isOutputActive(IO::OUTFLAG_WAVE_ELEV))
				IO::print2outLine(IO::OUTFLAG_WAVE_ELEV, m_waveProbeID[ii], waveElevAtProbe(ii));

			if (IO::isOutputActive(IO::OUTFLAG_WAVE_VEL))
				IO::print2outLine(IO::OUTFLAG_WAVE_VEL, m_waveProbeID[ii], waveVelAtProbe(ii));

			if (IO::isOutputActive(IO::OUTFLAG_WAVE_ACC))
				IO::print2outLine(IO::OUTFLAG_WAVE_ACC, m_waveProbeID[ii], du1dt(m_waveProbe[ii], 0));

			if (IO::isOutputActive(IO::OUTFLAG_WAVE_ACC_2ND))
				IO::print2outLine(IO::OUTFLAG_WAVE_ACC_2ND, m_waveProbeID[ii], du2dt(m_waveProbe[ii]));

			if (IO::isOutputActive(IO::OUTFLAG_WAVE_PRES))
				IO::print2outLine(IO::OUTFLAG_WAVE_PRES, m_waveProbeID[ii], wavePressure(m_waveProbe[ii]));

			if (IO::isOutputActive(IO::OUTFLAG_WAVE_PRES_2ND))
				IO::print2outLine(IO::OUTFLAG_WAVE_PRES_2ND, m_waveProbeID[ii], wavePressure_2ndOrd(m_waveProbe[ii]));
	}
}


/*****************************************************
	Other functions
*****************************************************/
// Return coordinates of a node based on its ID
// Throws a std::runtime_error if the node could not be found
arma::vec::fixed<3> ENVIR::getNode(unsigned int ID) const
{
	std::vector<unsigned int>::const_iterator iter = std::find(m_nodesID.begin(), m_nodesID.end(), ID); // Find node by its ID.
	vec::fixed<3> node_coord(fill::zeros);
	if (iter != m_nodesID.end())
	{
		auto index = std::distance(m_nodesID.begin(), iter); // Get index by the distance between the iterator found above and m_nodes.begin()
		node_coord = m_nodesCoord.at(index);
	}
	else
	{
		throw std::runtime_error("Unable to find node with ID " + std::to_string(ID) + ". Error in input line " + std::to_string(IO::getInLineNumber()) + ".");
	}

	return node_coord;
}

bool ENVIR::isNodeEmpty() const
{
	return m_nodesID.empty();
}

bool ENVIR::isWaveProbeEmpty() const
{
	return m_waveProbe.empty();
}

void ENVIR::stepTime()
{
	stepTime(m_timeStep);
}

void ENVIR::stepTime(double const step)
{
	m_time += step;

	// Find indices for interpolation of the time vector
	m_ind4interp1 = index_min(abs(m_timeArray - m_time));	

	m_shouldInterp = true ;
	double t1 = m_timeArray.at(m_ind4interp1);
	if (abs(t1 - m_time) < arma::datum::eps)
	{
		m_shouldInterp = false;
	}
	else
	{
		(t1 < m_time) ? (m_ind4interp2 = m_ind4interp1 + 1) : (m_ind4interp2 = m_ind4interp1, m_ind4interp1 = m_ind4interp2 - 1);
	}

	// Sometimes at the end of the simulation ind2 may be larger than the size of the array
	// due to a time step that is slightly larger than it should. It is ok to set both ind1 to the last index
	// and avoid interpolations in this case
	if (m_ind4interp1 > m_timeArray.size() - 1 || m_ind4interp2 > m_timeArray.size() - 1)
	{
		m_ind4interp1 = m_timeArray.size() - 1;
		m_shouldInterp = false;
	}
}


double ENVIR::ramp() const
{
	return ramp(m_time);
}

double ENVIR::ramp(double time) const
{
	double ramp{ 1 };

	if (time < m_timeRamp)
	{
		ramp = 0.5 * (1 - cos(datum::pi * time / m_timeRamp));
	}

	return ramp;
}


// We consider linear Airy waves, with velocity potential:
// phi = g*A/w * cosh(k(z+h))/cosh(k*h) * sin(k*x - w*t)
cx_double ENVIR::waveElev_coef(const double x, const double y, const unsigned int waveIndex) const
{
	double A = m_wave.at(waveIndex).amp();
	double k = m_wave.at(waveIndex).waveNumber();
	double beta = m_wave.at(waveIndex).direction() * arma::datum::pi / 180.;
	double phase = m_wave.at(waveIndex).phase() * arma::datum::pi / 180.;

	return { A*cos(k*cos(beta)*x + k * sin(beta)*y + phase) , -A*sin(k*cos(beta)*x + k * sin(beta)*y + phase) };
}

double ENVIR::waveElev(const double x, const double y, const double time) const
{
	double elev{ 0 };

	for (int ii = 0; ii < m_wave.size(); ++ii)
	{
		if (m_wave.at(ii).amp() != 0)
		{
			double w = m_wave.at(ii).angFreq();
			elev += (waveElev_coef(x, y, ii) * cx_double { cos(w * time), sin(w * time) }).real();
		}
	}

	return elev;
}

double ENVIR::waveElev(const double x, const double y) const
{	
	return ramp() * waveElev(x, y, m_time);
}


double ENVIR::wavePressure(const vec::fixed<3> &coord, const unsigned int waveIndex) const
{
	double p(0);

	double x = coord[0];
	double y = coord[1];
	double z = coord[2];
	double h = m_watDepth;
	double t = m_time;
	double rho = m_watDens;
	double g = m_gravity;
	double w = m_wave.at(waveIndex).angFreq();
	double A = m_wave.at(waveIndex).amp();
	double k = m_wave.at(waveIndex).waveNumber();
	double cosBeta = m_wave.at(waveIndex).cosBeta();
	double sinBeta = m_wave.at(waveIndex).sinBeta();
	double phase = m_wave.at(waveIndex).phase() * arma::datum::pi / 180.;

	// This formulation is valid only below the mean water level, i.e. z <= 0
	if (z <= 0)
	{
		// When k*h is too high, which happens for deep water/short waves, sinh(k*h) and cosh(k*h) become too large and are considered "inf".
		// Hence, we chose a threshold of 10, above which the deep water approximation is employed.
		double khz(0);
		if (k*h >= 10)
		{
			khz = exp(k*z);
		}
		else
		{
			khz = cosh(k * (z + h)) / cosh(k*h);
		}

		return ramp() * rho * g * A * khz * cos(k*cosBeta*x + k * sinBeta*y - w * t + phase);
	}

	else
	{
		return 0;
	}
}

double ENVIR::wavePressure(const vec::fixed<3> &coord) const
{
	double p(0);

	for (int ii = 0; ii < m_wave.size(); ++ii)
	{
		p += ENVIR::wavePressure(coord, ii);
	}

	return p;
}

cx_vec::fixed<3> ENVIR::u1_coef(const double x, const double y, const double z, const unsigned int waveIndex) const
{
	cx_vec::fixed<3> vel(fill::zeros);
	double h = m_watDepth;
	double t = m_time;
	double w = m_wave.at(waveIndex).angFreq();
	double A = m_wave.at(waveIndex).amp();
	double k = m_wave.at(waveIndex).waveNumber();
	double cosBeta = m_wave.at(waveIndex).cosBeta();
	double sinBeta = m_wave.at(waveIndex).sinBeta();
	double phase = m_wave.at(waveIndex).phase() * arma::datum::pi / 180.;
	double khz_xy(0), khz_z(0);

	if (z <= 0 && k > 0)
	{
		// When k*h is too high, which happens for deep water/short waves, sinh(k*h) and cosh(k*h) become too large and are considered "inf".
		// Hence, we chose a threshold of 10, above which the deep water approximation is employed.
		if (k*h >= 10)
		{
			khz_xy = exp(k*z);
			khz_z = khz_xy;
		}
		else
		{
			khz_xy = cosh(k * (z + h)) / sinh(k*h);
			khz_z = sinh(k * (z + h)) / sinh(k*h);
		}

		vel.row(0) = { cosBeta * khz_xy * cos(k*cosBeta*x + k * sinBeta*y + phase), -cosBeta * khz_xy * sin(k*cosBeta*x + k * sinBeta*y + phase) };
		vel.row(1) = { sinBeta * khz_xy * cos(k*cosBeta*x + k * sinBeta*y + phase), -sinBeta * khz_xy * sin(k*cosBeta*x + k * sinBeta*y + phase) };
		vel.row(2) = { khz_z * sin(k*cosBeta*x + k * sinBeta*y + phase), khz_z * cos(k*cosBeta*x + k * sinBeta*y + phase) };
	}

	return w * A * vel;
}

vec::fixed<3> ENVIR::u1(const vec::fixed<3> &coord, const double zwl) const
{
	arma::vec::fixed<3> vel = { 0,0,0 };

	if (m_waveStret <= 1 && zwl != 0)
	{
		throw std::runtime_error("zwl = " + std::to_string(zwl) + "is incompatible with wave streching mode " + std::to_string(m_waveStret));
	}

	double z = coord.at(2);
	if (z > zwl)
	{
		return vel;
	}

	if (m_waveStret == 2)
	{
		z = m_watDepth * (m_watDepth + z) / (m_watDepth + zwl) - m_watDepth;
	}

	for (int ii = 0; ii < m_wave.size(); ++ii)
	{
		if (m_wave.at(ii).amp() != 0)
		{
			double w = m_wave.at(ii).angFreq();
			vel += real(u1_coef(coord.at(0), coord.at(1), z, ii) * cx_double { cos(w * m_time), sin(w * m_time) });
		}
	}
	return vel * ramp();
}

vec::fixed<3> ENVIR::du1dt_eachWave(const vec::fixed<3> &coord, const unsigned int waveIndex) const
{
	arma::vec::fixed<3> acc = { 0,0,0 };

	// More friendly notation
	double x = coord[0];
	double y = coord[1];
	double z = coord[2];
	double h = m_watDepth;
	double t = m_time;
	double w = m_wave.at(waveIndex).angFreq();
	double A = m_wave.at(waveIndex).amp();
	double k = m_wave.at(waveIndex).waveNumber();
	double cosBeta = m_wave.at(waveIndex).cosBeta();
	double sinBeta = m_wave.at(waveIndex).sinBeta();
	double phase = m_wave.at(waveIndex).phase() * arma::datum::pi / 180.;
	double khz_xy(0), khz_z(0);

	if (z <= 0 && k > 0)
	{
		// When k*h is too high, which happens for deep water/short waves, sinh(k*h) and cosh(k*h) become too large and are considered "inf".
		// Hence, we chose a threshold of 10, above which the deep water approximation is employed.
		if (k*h >= 10)
		{
			khz_xy = exp(k*z);
			khz_z = khz_xy;
		}
		else
		{
			khz_xy = cosh(k * (z + h)) / sinh(k*h);
			khz_z = sinh(k * (z + h)) / sinh(k*h);
		}

		acc[0] = w * w * A * khz_xy * cosBeta * sin(k*cosBeta*x + k * sinBeta*y - w * t + phase);
		acc[1] = w * w * A * khz_xy * sinBeta * sin(k*cosBeta*x + k * sinBeta*y - w * t + phase);
		acc[2] = -w * w * A * khz_z * cos(k*cosBeta*x + k * sinBeta*y - w * t + phase);
	}

	return acc * ramp();
}


vec::fixed<3> ENVIR::du1dt(const vec::fixed<3> &coord, const double zwl) const
{
	arma::vec::fixed<3> acc = { 0,0,0 };

	if (m_waveStret <= 1 && zwl != 0)
	{
		throw std::runtime_error("zwl = " + std::to_string(zwl) + "is incompatible with wave streching mode " + std::to_string(m_waveStret));
	}

	double z = coord.at(2);
	if (z > zwl)
	{
		return acc;
	}

	if (m_waveStret == 2)
	{
		z = m_watDepth * (m_watDepth + z) / (m_watDepth + zwl) - m_watDepth;
	}

	for (int ii = 0; ii < m_wave.size(); ++ii)
	{
		acc += ENVIR::du1dt_eachWave(vec::fixed<3>({ coord.at(0), coord.at(1), z }), ii);
	}

	return acc;
}

vec::fixed<3> ENVIR::du1dx_eachWave(const vec::fixed<3> &coord, const unsigned int waveIndex) const
{
	arma::vec::fixed<3> acc = { 0,0,0 };

	// More friendly notation
	double x = coord[0];
	double y = coord[1];
	double z = coord[2];
	double h = m_watDepth;
	double t = m_time;
	double w = m_wave.at(waveIndex).angFreq();
	double A = m_wave.at(waveIndex).amp();
	double k = m_wave.at(waveIndex).waveNumber();
	double cosBeta = m_wave.at(waveIndex).cosBeta();
	double sinBeta = m_wave.at(waveIndex).sinBeta();
	double phase = m_wave.at(waveIndex).phase() * arma::datum::pi / 180.;
	double khz_xy(0), khz_z(0);

	if (z <= 0 && k > 0)
	{
		// When k*h is too high, which happens for deep water/short waves, sinh(k*h) and cosh(k*h) become too large and are considered "inf".
		// Hence, we chose a threshold of 10, above which the deep water approximation is employed.
		if (k*h >= 10)
		{
			khz_xy = exp(k*z);
			khz_z = khz_xy;
		}
		else
		{
			khz_xy = cosh(k * (z + h)) / sinh(k*h);
			khz_z = sinh(k * (z + h)) / sinh(k*h);
		}

		acc[0] = -w * A * k * khz_xy * cosBeta * cosBeta * sin(k*cosBeta*x + k * sinBeta*y - w * t + phase);
		acc[1] = -w * A * k * khz_xy * sinBeta * cosBeta * sin(k*cosBeta*x + k * sinBeta*y - w * t + phase);
		acc[2] = w * A * k * khz_z * cosBeta * cos(k*cosBeta*x + k * sinBeta*y - w * t + phase);
	}

	return acc * ramp();
}

vec::fixed<3> ENVIR::du1dx(const vec::fixed<3> &coord, const double zwl) const
{
	arma::vec::fixed<3> acc = { 0,0,0 };

	if (m_waveStret <= 1 && zwl != 0)
	{
		throw std::runtime_error("zwl = " + std::to_string(zwl) + "is incompatible with wave streching mode " + std::to_string(m_waveStret));
	}

	double z = coord.at(2);
	if (z > zwl)
	{
		return acc;
	}

	if (m_waveStret == 2)
	{
		z = m_watDepth * (m_watDepth + z) / (m_watDepth + zwl) - m_watDepth;
	}

	for (int ii = 0; ii < m_wave.size(); ++ii)
	{
		acc += ENVIR::du1dx_eachWave(vec::fixed<3>({ coord.at(0), coord.at(1), z }), ii);
	}
	return acc;
}

vec::fixed<3> ENVIR::du1dy_eachWave(const vec::fixed<3> &coord, const unsigned int waveIndex) const
{
	arma::vec::fixed<3> acc = { 0,0,0 };

	// More friendly notation
	double x = coord[0];
	double y = coord[1];
	double z = coord[2];
	double h = m_watDepth;
	double t = m_time;
	double w = m_wave.at(waveIndex).angFreq();
	double A = m_wave.at(waveIndex).amp();
	double k = m_wave.at(waveIndex).waveNumber();
	double cosBeta = m_wave.at(waveIndex).cosBeta();
	double sinBeta = m_wave.at(waveIndex).sinBeta();
	double phase = m_wave.at(waveIndex).phase() * arma::datum::pi / 180.;
	double khz_xy(0), khz_z(0);

	if (z <= 0 && k > 0)
	{
		// When k*h is too high, which happens for deep water/short waves, sinh(k*h) and cosh(k*h) become too large and are considered "inf".
		// Hence, we chose a threshold of 10, above which the deep water approximation is employed.
		if (k*h >= 10)
		{
			khz_xy = exp(k*z);
			khz_z = khz_xy;
		}
		else
		{
			khz_xy = cosh(k * (z + h)) / sinh(k*h);
			khz_z = sinh(k * (z + h)) / sinh(k*h);
		}

		acc[0] = -w * A * k * khz_xy * cosBeta * sinBeta * sin(k*cosBeta*x + k * sinBeta*y - w * t + phase);
		acc[1] = -w * A * k * khz_xy * sinBeta * sinBeta * sin(k*cosBeta*x + k * sinBeta*y - w * t + phase);
		acc[2] = w * A * k * khz_z * sinBeta * cos(k*cosBeta*x + k * sinBeta*y - w * t + phase);
	}

	return acc * ramp();
}

vec::fixed<3> ENVIR::du1dy(const vec::fixed<3> &coord, const double zwl) const
{
	arma::vec::fixed<3> acc = { 0,0,0 };

	if (m_waveStret <= 1 && zwl != 0)
	{
		throw std::runtime_error("zwl = " + std::to_string(zwl) + "is incompatible with wave streching mode " + std::to_string(m_waveStret));
	}

	double z = coord.at(2);
	if (z > zwl)
	{
		return acc;
	}

	if (m_waveStret == 2)
	{
		z = m_watDepth * (m_watDepth + z) / (m_watDepth + zwl) - m_watDepth;
	}

	for (int ii = 0; ii < m_wave.size(); ++ii)
	{
		acc += ENVIR::du1dy_eachWave(vec::fixed<3>({ coord.at(0), coord.at(1), z }), ii);
	}
	return acc;
}

vec::fixed<3> ENVIR::du1dz_eachWave(const vec::fixed<3> &coord, const unsigned int waveIndex) const
{
	arma::vec::fixed<3> acc = { 0,0,0 };

	// More friendly notation
	double x = coord[0];
	double y = coord[1];
	double z = coord[2];
	double h = m_watDepth;
	double t = m_time;
	double w = m_wave.at(waveIndex).angFreq();
	double A = m_wave.at(waveIndex).amp();
	double k = m_wave.at(waveIndex).waveNumber();
	double cosBeta = m_wave.at(waveIndex).cosBeta();
	double sinBeta = m_wave.at(waveIndex).sinBeta();
	double phase = m_wave.at(waveIndex).phase() * arma::datum::pi / 180.;
	double khz_xy(0), khz_z(0);

	if (z <= 0 && k > 0)
	{
		// When k*h is too high, which happens for deep water/short waves, sinh(k*h) and cosh(k*h) become too large and are considered "inf".
		// Hence, we chose a threshold of 10, above which the deep water approximation is employed.
		if (k*h >= 10)
		{
			khz_xy = exp(k*z);
			khz_z = khz_xy;
		}
		else
		{
			khz_xy = cosh(k * (z + h)) / sinh(k*h);
			khz_z = sinh(k * (z + h)) / sinh(k*h);
		}

		acc[0] = w * A * k * khz_z * cosBeta * cos(k*cosBeta*x + k * sinBeta*y - w * t + phase);
		acc[1] = w * A * k * khz_z * sinBeta * cos(k*cosBeta*x + k * sinBeta*y - w * t + phase);
		acc[2] = w * A * k * khz_xy * sin(k*cosBeta*x + k * sinBeta*y - w * t + phase);
	}

	return acc * ramp();
}

vec::fixed<3> ENVIR::du1dz(const vec::fixed<3> &coord, const double zwl) const
{
	arma::vec::fixed<3> acc = { 0,0,0 };

	if (m_waveStret <= 1 && zwl != 0)
	{
		throw std::runtime_error("zwl = " + std::to_string(zwl) + "is incompatible with wave streching mode " + std::to_string(m_waveStret));
	}

	double z = coord.at(2);
	if (z > zwl)
	{
		return acc;
	}

	if (m_waveStret == 2)
	{
		z = m_watDepth * (m_watDepth + z) / (m_watDepth + zwl) - m_watDepth;
	}

	for (int ii = 0; ii < m_wave.size(); ++ii)
	{
		acc += ENVIR::du1dz_eachWave(vec::fixed<3>({ coord.at(0), coord.at(1), z }), ii);
	}
	return acc;
}

vec::fixed<3> ENVIR::du2dt(const vec::fixed<3> &coord, const unsigned int waveIndex1, const unsigned int waveIndex2) const
{
	arma::vec::fixed<3> acc = { 0,0,0 };

	// In this case, the calculation below would lead to 0/0, but the limit is actually 0. This is fine, as the second order potential should not contribute to the mean drift.
	if (waveIndex1 == waveIndex2)
	{
		return acc;
	}

	// More friendly notation
	double z = coord[2];
	double h = m_watDepth;
	double t = m_time;
	double g = m_gravity;

	double w1 = m_wave.at(waveIndex1).angFreq();
	double A1 = m_wave.at(waveIndex1).amp();
	double k1 = m_wave.at(waveIndex1).waveNumber();
	double b1 = m_wave.at(waveIndex1).direction() * arma::datum::pi / 180.;
	double p1 = m_wave.at(waveIndex1).phase() * arma::datum::pi / 180.;
	double cosB1 = m_wave.at(waveIndex1).cosBeta();
	double sinB1 = m_wave.at(waveIndex1).sinBeta();

	double A2 = m_wave.at(waveIndex2).amp();
	double w2 = m_wave.at(waveIndex2).angFreq();
	double k2 = m_wave.at(waveIndex2).waveNumber();
	double b2 = m_wave.at(waveIndex2).direction() * arma::datum::pi / 180.;
	double p2 = m_wave.at(waveIndex2).phase() * arma::datum::pi / 180.;
	double cosB2 = m_wave.at(waveIndex2).cosBeta();
	double sinB2 = m_wave.at(waveIndex2).sinBeta();

	// This formulation is valid only below the mean water level, i.e. z <= 0
	if (z <= 0 && k1 > 0 && k2 > 0)
	{
		arma::vec::fixed<3> k1_k2 = { k1 * cosB1 - k2 * cosB2, k1 * sinB1 - k2 * sinB2, 0 };
		double norm_k1_k2 = arma::norm(k1_k2);

		// Isso aqui soh depende das propriedades das ondas e da profundidades. Da pra calcular previamente uma vez soh.
		double aux = ((w2 - w1) / (w1 * w2)) * k1 * k2 * (std::cos(b1 - b2) + std::tanh(k1*h) * std::tanh(k2*h))
			- 0.5 * (k1*k1 / (w1 * pow(std::cosh(k1*h), 2)) - k2 * k2 / (w2 * pow(std::cosh(k2*h), 2)));
		aux = aux / (g * norm_k1_k2 * std::tanh(norm_k1_k2 * h) - (w1 - w2)*(w1 - w2));

		double khz_xy = std::cosh(norm_k1_k2 * (z + h)) / std::cosh(norm_k1_k2 * h);
		double khz_z = std::sinh(norm_k1_k2 * (z + h)) / std::cosh(norm_k1_k2 * h);

		acc[0] = 0.5 * A1 * A2 * (w1 - w2) * (k1 * cosB1 - k2 * cosB2) * g*g * aux * khz_xy
			* std::sin(dot(k1_k2, coord) - (w1 - w2) * t + p1 - p2);

		acc[1] = 0.5 * A1 * A2 * (w1 - w2) * (k1 * sinB1 - k2 * sinB2) * g*g * aux * khz_xy
			* std::sin(dot(k1_k2, coord) - (w1 - w2) * t + p1 - p2);

		acc[2] = -0.5 * A1 * A2 * (w1 - w2) * g*g * aux * norm_k1_k2 * khz_z
			* std::cos(dot(k1_k2, coord) - (w1 - w2) * t + p1 - p2);
	}

	return acc * ramp();
}

vec::fixed<3> ENVIR::du2dt(const vec::fixed<3> &coord) const
{
	arma::vec::fixed<3> acc = { 0,0,0 };

	// When i == j, du2dt = {0,0,0}, so it is safe to skip this part of the loop.
	// Besides, as only the real part of the second-order difference-frequency potential is used,
	// the acceleration due to a pair ij is equal to ji.
	for (int ii = 0; ii < m_wave.size(); ++ii)
	{
		for (int jj = ii + 1; jj < m_wave.size(); ++jj)
		{
			acc += 2 * ENVIR::du2dt(coord, ii, jj);
		}
	}
	return acc;
}

double ENVIR::wavePressure_2ndOrd(const vec::fixed<3> &coord, const unsigned int waveIndex1, const unsigned int waveIndex2) const
{
	double p{ 0 };

	// In this case, the calculation below would lead to 0/0, but the limit is actually 0. This is fine, as the second order potential should not contribute to the mean drift.
	if (waveIndex1 == waveIndex2)
	{
		return p;
	}

	// More friendly notation
	double z = coord[2];
	double h = m_watDepth;
	double t = m_time;
	double g = m_gravity;

	double w1 = m_wave.at(waveIndex1).angFreq();
	double A1 = m_wave.at(waveIndex1).amp();
	double k1 = m_wave.at(waveIndex1).waveNumber();
	double b1 = m_wave.at(waveIndex1).direction() * arma::datum::pi / 180.;
	double phase1 = m_wave.at(waveIndex1).phase() * arma::datum::pi / 180.;
	double cosB1 = m_wave.at(waveIndex1).cosBeta();
	double sinB1 = m_wave.at(waveIndex1).sinBeta();

	double A2 = m_wave.at(waveIndex2).amp();
	double w2 = m_wave.at(waveIndex2).angFreq();
	double k2 = m_wave.at(waveIndex2).waveNumber();
	double b2 = m_wave.at(waveIndex2).direction() * arma::datum::pi / 180.;
	double phase2 = m_wave.at(waveIndex2).phase() * arma::datum::pi / 180.;
	double cosB2 = m_wave.at(waveIndex2).cosBeta();
	double sinB2 = m_wave.at(waveIndex2).sinBeta();

	// This formulation is valid only below the mean water level, i.e. z <= 0
	if (z <= 0 && k1 > 0 && k2 > 0)
	{
		arma::vec::fixed<3> k1_k2 = { k1 * cosB1 - k2 * cosB2, k1 * sinB1 - k2 * sinB2, 0 };
		double norm_k1_k2 = arma::norm(k1_k2);

		// Isso aqui soh depende das propriedades das ondas e da profundidades. Da pra calcular previamente uma vez soh.
		double aux = ((w2 - w1) / (w1 * w2)) * k1 * k2 * (std::cos(b1 - b2) + std::tanh(k1*h) * std::tanh(k2*h))
			- 0.5 * (k1*k1 / (w1 * pow(std::cosh(k1*h), 2)) - k2 * k2 / (w2 * pow(std::cosh(k2*h), 2)));
		aux = aux / (g * norm_k1_k2 * std::tanh(norm_k1_k2 * h) - (w1 - w2)*(w1 - w2));

		p = 0.5 * m_watDens * A1 * A2 * (w1 - w2) * g*g * aux * std::cosh(norm_k1_k2 * (z + h)) / std::cosh(norm_k1_k2 * h)
			* std::cos(dot(k1_k2, coord) - (w1 - w2) * t + phase1 - phase2);
	}

	return p * ramp();
}

double ENVIR::wavePressure_2ndOrd(const vec::fixed<3> &coord) const
{
	double p{ 0 };

	// When i == j, p = 0, so it is safe to skip this part of the loop.
	// Besides, as only the real part of the second-order difference-frequency potential is used,
	// the acceleration due to a pair ij is equal to ji.
	for (int ii = 0; ii < m_wave.size(); ++ii)
	{
		for (int jj = ii + 1; jj < m_wave.size(); ++jj)
		{
			p += 2 * ENVIR::wavePressure_2ndOrd(coord, ii, jj);
		}
	}
	return p;
}

vec::fixed<3> ENVIR::dadx_eachWave(const vec::fixed<3> &coord, const unsigned int waveIndex) const
{
	arma::vec::fixed<3> acc = { 0,0,0 };

	// More friendly notation
	double x = coord[0];
	double y = coord[1];
	double z = coord[2];
	double h = m_watDepth;
	double t = m_time;
	double w = m_wave.at(waveIndex).angFreq();
	double A = m_wave.at(waveIndex).amp();
	double k = m_wave.at(waveIndex).waveNumber();
	double cosBeta = m_wave.at(waveIndex).cosBeta();
	double sinBeta = m_wave.at(waveIndex).sinBeta();
	double phase = m_wave.at(waveIndex).phase() * arma::datum::pi / 180.;
	double khz_xy(0), khz_z(0);

	if (z <= 0 && k > 0)
	{
		// When k*h is too high, which happens for deep water/short waves, sinh(k*h) and cosh(k*h) become too large and are considered "inf".
		// Hence, we chose a threshold of 10, above which the deep water approximation is employed.
		if (k*h >= 10)
		{
			khz_xy = exp(k*z);
			khz_z = khz_xy;
		}
		else
		{
			khz_xy = cosh(k * (z + h)) / sinh(k*h);
			khz_z = sinh(k * (z + h)) / sinh(k*h);
		}

		acc[0] = w * w * A * khz_xy * k * cosBeta * cosBeta * cos(k*cosBeta*x + k * sinBeta*y - w * t + phase);
		acc[1] = w * w * A * khz_xy * k * cosBeta * sinBeta * cos(k*cosBeta*x + k * sinBeta*y - w * t + phase);
		acc[2] = -w * w * A * khz_z * k * cosBeta * cos(k*cosBeta*x + k * sinBeta*y - w * t + phase);
	}

	return acc * ramp();
}

vec::fixed<3> ENVIR::dadx(const vec::fixed<3> &coord, const double zwl) const
{
	arma::vec::fixed<3> acc = { 0,0,0 };

	if (m_waveStret <= 1 && zwl != 0)
	{
		throw std::runtime_error("zwl = " + std::to_string(zwl) + "is incompatible with wave streching mode " + std::to_string(m_waveStret));
	}

	double z = coord.at(2);
	if (z > zwl)
	{
		return acc;
	}

	if (m_waveStret == 2)
	{
		z = m_watDepth * (m_watDepth + z) / (m_watDepth + zwl) - m_watDepth;
	}

	for (int ii = 0; ii < m_wave.size(); ++ii)
	{
		acc += ENVIR::dadx_eachWave(vec::fixed<3>({ coord.at(0), coord.at(1), z }), ii);
	}
	return acc;
}

vec::fixed<3> ENVIR::dady_eachWave(const vec::fixed<3> &coord, const unsigned int waveIndex) const
{
	arma::vec::fixed<3> acc = { 0,0,0 };

	// More friendly notation
	double x = coord[0];
	double y = coord[1];
	double z = coord[2];
	double h = m_watDepth;
	double t = m_time;
	double w = m_wave.at(waveIndex).angFreq();
	double A = m_wave.at(waveIndex).amp();
	double k = m_wave.at(waveIndex).waveNumber();
	double cosBeta = m_wave.at(waveIndex).cosBeta();
	double sinBeta = m_wave.at(waveIndex).sinBeta();
	double phase = m_wave.at(waveIndex).phase() * arma::datum::pi / 180.;
	double khz_xy(0), khz_z(0);

	if (z <= 0 && k > 0)
	{
		// When k*h is too high, which happens for deep water/short waves, sinh(k*h) and cosh(k*h) become too large and are considered "inf".
		// Hence, we chose a threshold of 10, above which the deep water approximation is employed.
		if (k*h >= 10)
		{
			khz_xy = exp(k*z);
			khz_z = khz_xy;
		}
		else
		{
			khz_xy = cosh(k * (z + h)) / sinh(k*h);
			khz_z = sinh(k * (z + h)) / sinh(k*h);
		}

		acc[0] = w * w * A * khz_xy * k * sinBeta * cosBeta * cos(k*cosBeta*x + k * sinBeta*y - w * t + phase);
		acc[1] = w * w * A * khz_xy * k * sinBeta * sinBeta * cos(k*cosBeta*x + k * sinBeta*y - w * t + phase);
		acc[2] = -w * w * A * khz_z * k * sinBeta * cos(k*cosBeta*x + k * sinBeta*y - w * t + phase);
	}


	return acc * ramp();
}

vec::fixed<3> ENVIR::dady(const vec::fixed<3> &coord, const double zwl) const
{
	arma::vec::fixed<3> acc = { 0,0,0 };

	if (m_waveStret <= 1 && zwl != 0)
	{
		throw std::runtime_error("zwl = " + std::to_string(zwl) + "is incompatible with wave streching mode " + std::to_string(m_waveStret));
	}

	double z = coord.at(2);
	if (z > zwl)
	{
		return acc;
	}

	if (m_waveStret == 2)
	{
		z = m_watDepth * (m_watDepth + z) / (m_watDepth + zwl) - m_watDepth;
	}

	for (int ii = 0; ii < m_wave.size(); ++ii)
	{
		acc += ENVIR::dady_eachWave(vec::fixed<3>({ coord.at(0), coord.at(1), z }), ii);
	}
	return acc;
}

vec::fixed<3> ENVIR::dadz_eachWave(const vec::fixed<3> &coord, const unsigned int waveIndex) const
{
	arma::vec::fixed<3> acc = { 0,0,0 };

	// More friendly notation
	double x = coord[0];
	double y = coord[1];
	double z = coord[2];
	double h = m_watDepth;
	double t = m_time;
	double w = m_wave.at(waveIndex).angFreq();
	double A = m_wave.at(waveIndex).amp();
	double k = m_wave.at(waveIndex).waveNumber();
	double cosBeta = m_wave.at(waveIndex).cosBeta();
	double sinBeta = m_wave.at(waveIndex).sinBeta();
	double phase = m_wave.at(waveIndex).phase() * arma::datum::pi / 180.;
	double khz_xy(0), khz_z(0);

	if (z <= 0 && k > 0)
	{
		// When k*h is too high, which happens for deep water/short waves, sinh(k*h) and cosh(k*h) become too large and are considered "inf".
		// Hence, we chose a threshold of 10, above which the deep water approximation is employed.
		if (k*h >= 10)
		{
			khz_xy = exp(k*z);
			khz_z = khz_xy;
		}
		else
		{
			khz_xy = sinh(k * (z + h)) / sinh(k*h);
			khz_z = cosh(k * (z + h)) / sinh(k*h);
		}

		acc[0] = w * w * k * A * khz_xy * cosBeta * sin(k*cosBeta*x + k * sinBeta*y - w * t + phase);
		acc[1] = w * w * k * A * khz_xy * sinBeta * sin(k*cosBeta*x + k * sinBeta*y - w * t + phase);
		acc[2] = -w * w * k * A * khz_z * cos(k*cosBeta*x + k * sinBeta*y - w * t + phase);
	}

	return acc * ramp();
}

vec::fixed<3> ENVIR::dadz(const vec::fixed<3> &coord, const double zwl) const
{
	arma::vec::fixed<3> acc = { 0,0,0 };

	if (m_waveStret <= 1 && zwl != 0)
	{
		throw std::runtime_error("zwl = " + std::to_string(zwl) + "is incompatible with wave streching mode " + std::to_string(m_waveStret));
	}

	double z = coord.at(2);
	if (z > zwl)
	{
		return acc;
	}

	if (m_waveStret == 2)
	{
		z = m_watDepth * (m_watDepth + z) / (m_watDepth + zwl) - m_watDepth;
	}

	for (int ii = 0; ii < m_wave.size(); ++ii)
	{
		acc += ENVIR::dadz_eachWave(vec::fixed<3>({ coord.at(0), coord.at(1), z }), ii);
	}
	return acc;
}



double ENVIR::windVel_X(const vec::fixed<3> &coord) const
{
	return (ramp() * windRefVel() * std::cos(m_windDir * arma::datum::pi / 180.) * pow(coord[2] / windRefHeight(), windExp()));
}

double ENVIR::windVel_Y(const vec::fixed<3> &coord) const
{
	return (-ramp() * windRefVel() * std::sin(m_windDir * arma::datum::pi / 180.) * pow(coord[2] / windRefHeight(), windExp()));
}








// JONSWAP spectrum considering angular frequency (rad/s)
double JONSWAP(const double w, const double Tp, const double Hs, const double gamma)
{
	if (w == 0)
	{
		return 0;
	}

	double wp = 2 * arma::datum::pi / Tp;
	double sigma(0);
	double A(0);
	double Sw(0);

	if (w <= wp)
	{
		sigma = 0.07;
	}
	else
	{
		sigma = 0.09;
	}

	A = std::exp(-0.5 * std::pow((w / wp - 1) / sigma, 2));
	return 320 * Hs*Hs * pow(Tp, -4) * std::pow(w, -5) * std::exp(-1950 * pow(Tp, -4) * std::pow(w, -4)) * std::pow(gamma, A);
}