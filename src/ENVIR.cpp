#include "ENVIR.h"
#include "IO.h"
#include "auxFunctions.h"

#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>    // std::binary_search
#include <utility> // For std::move
#include <complex> 

using namespace arma;

// TODO: Wheeler stretching is not considered when IFFT is used. Need to fix that later.

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

	m_windRefVel = arma::datum::nan;
	m_windRefHeight = arma::datum::nan;
	m_windExp = arma::datum::nan;
	m_windDir = arma::datum::nan;
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

void ENVIR::setWindFromTurbFile(const std::string &fileName)
{
	std::ifstream is(fileName, std::ifstream::binary);
	if (!is)
	{
		throw std::runtime_error("Unable to read turbulence input file " + fileName + ".");
	}
	m_turbFileName = fileName;
	m_flagTurbWind = true;
	int nffc = 3; // I dont know what this means exactly, but it seems to be the number of spatial dimensions (X, Y and Z)

	// Need wind direction in some of the calculations below
	if (!arma::is_finite(m_windDir))
	{
		throw std::runtime_error("Need to specify WindDir before turbulence file.");
	}

	//----------------------------
	// get the header information
	//----------------------------
	// TurbSim format identifier
	// Non periodic: 7
	// Periodic: 8	
	int16_t tmp;
	is.read(reinterpret_cast<char*>(&tmp), sizeof(tmp));
	if (tmp == 8) m_flagPeriodicTurb = true;
		
	int32_t nz{ 0 }; // the number of grid points vertically, INT(4)
	is.read(reinterpret_cast<char*>(&nz), sizeof(nz));

	int32_t ny{ 0 }; // the number of grid points laterally, INT(4)
	is.read(reinterpret_cast<char*>(&ny), sizeof(ny));

	int32_t ntwr{ 0 }; // the number of tower points, INT(4)
	is.read(reinterpret_cast<char*>(&ntwr), sizeof(ntwr));

	int32_t nt{ 0 }; // the number of time steps, INT(4)
	is.read(reinterpret_cast<char*>(&nt), sizeof(nt));

	float dz{ 0 }; // grid spacing in vertical direction, REAL(4), in m
	is.read(reinterpret_cast<char*>(&dz), sizeof(float));

	float dy{ 0 }; // grid spacing in lateral direction, REAL(4), in m
	is.read(reinterpret_cast<char*>(&dy), sizeof(float));

	is.read(reinterpret_cast<char*>(&m_windDt), sizeof(float)); // grid spacing in delta time, REAL(4), in m / s
	m_windTimeTotal = (nt - 1) * m_windDt;

	is.read(reinterpret_cast<char*>(&m_windRefVel), sizeof(float)); // the mean wind speed at hub height, REAL(4), in m / s

	is.read(reinterpret_cast<char*>(&m_windRefHeight), sizeof(float)); // height of the hub, REAL(4), in m

	float z1{ 0 }; // height of the bottom of the grid, REAL(4), in m
	is.read(reinterpret_cast<char*>(&z1), sizeof(float));

	std::array<float, 3> Vslope;
	std::array<float, 3> Voffset;
	is.read(reinterpret_cast<char*>(&Vslope[0]), sizeof(float)); // the U - component slope for scaling, REAL(4)
	is.read(reinterpret_cast<char*>(&Voffset[0]), sizeof(float)); // the U - component offset for scaling, REAL(4)
	is.read(reinterpret_cast<char*>(&Vslope[1]), sizeof(float)); // the V - component slope for scaling, REAL(4)
	is.read(reinterpret_cast<char*>(&Voffset[1]), sizeof(float)); // the V - component offset for scaling, REAL(4)
	is.read(reinterpret_cast<char*>(&Vslope[2]), sizeof(float)); // the W - component slope for scaling, REAL(4)
	is.read(reinterpret_cast<char*>(&Voffset[2]), sizeof(float)); // the W - component offset for scaling, REAL(4)

	int32_t nchar{ 0 }; // the number of characters in the description string, max 200, INT(4)
	is.read(reinterpret_cast<char*>(&nchar), sizeof(nchar));

	char *asciiSTR = new char[nchar];
	is.read(asciiSTR, nchar); // the ASCII integer representation of the character string			
	std::cout << asciiSTR << "\n";
	delete[] asciiSTR;

	//-------------------------
	// get the grid information
	//-------------------------
	int32_t nPts = ny * nz;
	int32_t nv = nffc * nPts;       // the size of one time step
	int32_t nvTwr = nffc * ntwr;

	m_windVelocity.set_size(nt);
	m_windVelocity.fill(arma::zeros<arma::fcube>(nffc, ny, nz));
	arma::fcube twrVelocity(nt, nffc, ntwr);

	// Loop the time steps
	for (int it = 0; it < nt; ++it)
	{
		// --------------------
		// get the grid points
		// --------------------
		std::vector<int16_t> v(nv);
		for (int jj = 0; jj < nv; ++jj)
		{
			is.read(reinterpret_cast<char*>(&v[jj]), sizeof(int16_t));
		}

		int ip = 0;
		for (int iz = 0; iz < nz; ++iz)
		{
			for (int iy = 0; iy < ny; ++iy)
			{
				for (int k = 0; k < nffc; ++k)
				{
					m_windVelocity.at(it).at(k, iy, iz) = (v[ip] - Voffset[k]) / Vslope[k];
					ip = ip + 1;
				}

				// Rotate based on wind direction
				float auxU = m_windVelocity.at(it).at(0, iy, iz);
				float auxV = m_windVelocity.at(it).at(1, iy, iz);
				m_windVelocity.at(it).at(0, iy, iz) = auxU * m_windDirCos + auxV * m_windDirSin;
				m_windVelocity.at(it).at(1, iy, iz) = -auxU * m_windDirSin + auxV * m_windDirCos;								
			}

		}		
		//---------------------
		//get the tower points
		//---------------------
		// NOTE:: Tower part wasnt tested yet
		if (nvTwr > 0)
		{
			std::vector<int16_t> v(nvTwr);
			for (int jj = 0; jj < nvTwr; ++jj)
			{
				is.read(reinterpret_cast<char*>(&v[jj]), sizeof(int16_t));
			}

			// scale the data
			for (int k = 0; k < nffc; ++k)
			{
				for (int jj = 0; jj < nvTwr; ++jj)
				{
					twrVelocity(it, k, jj) = (v[k + jj * nffc] - Voffset[k]) / Vslope[k];
				}
			}
		}
	}

	m_windGrid_y = arma::regspace<fvec>(0, 1, ny - 1) * dy - (ny - 1) / 2.0 * dy;
	m_windGrid_z = arma::regspace<fvec>(0, 1, nz - 1) * dz + z1;

	// Just like in NREL's InflowWind, the grid is located at a distance gridwidth/2 downwind 
	// of the tower for non periodic wind, "so that no aerodynamic analysis nodes are outside
	// of the wind domain even if the rotor is initialized with a +/-90 degree yaw".
	double m_windGrid_x0{ 0 };
	if (!m_flagPeriodicTurb)
	{
		m_windGrid_x0 = m_windGrid_y.back();
	}	

	// Rotate grid considering wind dir (horizontal only)	
	m_windGrid_x =  m_windGrid_x0 * m_windDirCos + m_windGrid_y * m_windDirSin;
	m_windGrid_y = -m_windGrid_x0 * m_windDirSin + m_windGrid_y * m_windDirCos;

	// m_windGrid_x and m_windGrid_y may be in ascending or desceding order depending on the rotation angle,
	// but it is better to have m_windGrid_y in ascending order to make it easier to determine its 
	// lower and upper bound. Hence we sort it and rearrange m_windGrid_x accordingly.
	if (m_windGrid_y[1] < m_windGrid_y[0])
	{
		m_windGrid_x = arma::reverse(m_windGrid_x);
		m_windGrid_y = arma::reverse(m_windGrid_y);
	}
}

void ENVIR::setWindFromUnifFile(const std::string & fileName)
{
	// Need wind reference height, not exactly in this function, but when
	// the wind velocity is calculated. This check should be moved to
	// functions specifically concerned with consistency check of inputs,
	// but it is located here for the moment.
	if (!arma::is_finite(m_windRefHeight))
	{
		throw std::runtime_error("Need to specify 'windHeight' before wind file.");
	}

	std::ifstream is;
	is.open(fileName);
	if (!is)
	{
		throw std::runtime_error("Unable to read uniform wind file " + fileName + ".");
	}
	m_wnd_fileName = fileName;
	m_flagUnifWind = true;
	
	// Count number of lines in file
	std::string line;
	int nlines{ 0 }; // Number of lines with actual data (excluding header and comments)
	while (std::getline(is, line))
	{
		if (line.at(0) != '!') ++nlines; // Do not count comments, which are identified by "!"
	}

	// Clear and return to the beginning
	is.clear();
	is.seekg(0);

	// Initialize vectors to the correct size
	m_wnd_time.zeros(nlines);
	m_wnd_windSpeed.zeros(nlines);
	m_wnd_windDir.zeros(nlines);
	m_wnd_vertSpeed.zeros(nlines);
	m_wnd_horizSheer.zeros(nlines);
	m_wnd_windExp.zeros(nlines);
	m_wnd_linVertSheer.zeros(nlines);
	m_wnd_gustSpeed.zeros(nlines);

	// Read data to the respective vectors
	int ii = 0;
	nlines = 0;
	while (std::getline(is, line))
	{		
		++nlines;
		if (line.at(0) == '!') continue;		

		std::vector<std::string> input = stringTokenize(line, " \t");
		if (input.size() != 8)
		{
			throw std::runtime_error("Unable to read wind file '" + m_wnd_fileName + "'. Wrong number of parameters in line " + std::to_string(nlines) + ".");
		}

		m_wnd_time.at(ii)         = string2num<float>(input.at(0));
		m_wnd_windSpeed.at(ii)    = string2num<float>(input.at(1));
		m_wnd_windDir.at(ii)      = string2num<float>(input.at(2));
		m_wnd_vertSpeed.at(ii)    = string2num<float>(input.at(3));
		m_wnd_horizSheer.at(ii)   = string2num<float>(input.at(4));
		m_wnd_windExp.at(ii)      = string2num<float>(input.at(5));
		m_wnd_linVertSheer.at(ii) = string2num<float>(input.at(6));
		m_wnd_gustSpeed.at(ii)    = string2num<float>(input.at(7));
		++ii;
	}
}

void ENVIR::setWindRefLength(const double windRefLength)
{
	m_wnd_refLength = windRefLength;
}

void ENVIR::setWindRefVel(const double windRefVel)
{
	m_windRefVel = windRefVel;
}

void ENVIR::setWindDir(const double windDir)
{
	m_windDir = windDir;
	m_windDirCos = std::cos(m_windDir * arma::datum::pi / 180.);
	m_windDirSin = std::sin(m_windDir * arma::datum::pi / 180.);
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

		w = arma::regspace(0, dw0, 2 * arma::datum::pi / m_timeStep);
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

void ENVIR::addWaveElevSeries(const std::string &elevFlPath, const double direction, const double wlow, const double whigh)
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
			throw std::runtime_error("Unable to read time series of wave elevation from file " + elevFlPath + ". Wrong number of parameters in line " + std::to_string(ii + 1) + ".");
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
	double dw = 2 * arma::datum::pi / time.at(time.size() - 1);
	vec w = arma::regspace(0, dw, (time.size() - 1) * dw);
	m_flagIFFT = true;

	// FFT and assignment of the computed values to the vector of waves
	cx_vec A(mkl_fft_real(elev));
	A *= 1. / time.size();
	A.rows(std::floor(A.size() / 2), A.size() - 1).zeros();
	A.rows(1, A.size() - 1) *= 2;
	for (int ii = 0; ii < w.size(); ++ii)
	{
		if (w.at(ii) < wlow || (w.at(ii) > whigh && whigh > 0))
		{
			A.at(ii) = 0;
		}
		addRegularWave("TRWave", 2 * std::abs(A.at(ii)), 2 * arma::datum::pi / w.at(ii), direction, -std::arg(A.at(ii)) * 180. / arma::datum::pi);
	}

	// Make sure that the total simulation time are equal to the wave elevation file
	IO::print2log("Setting total simulation time to " + std::to_string(time.at(time.size() - 1)) + "s to match the one from the the wave elevation file. Previous value was " + std::to_string(m_timeTotal) + ".");
	m_timeTotal = time.at(time.size() - 1);
}


void ENVIR::addWindProbe(const unsigned int ID)
{
	// Check whether nodes were specified
	if (this->isNodeEmpty())
	{
		throw std::runtime_error("Nodes should be specified before adding wind probes. In: ENVIR::addWindProbe");
	}

	// Check if this wind probe was already included
	std::vector<unsigned int>::const_iterator iter = std::find(m_windProbeID.begin(), m_windProbeID.end(), ID);

	if (iter == m_windProbeID.end())
	{
		// Get the node coordinate and add it to m_windProbe
		m_windProbe.push_back(this->getNode(ID));
		m_windProbeID.push_back(ID);
	}
}

void ENVIR::addWaveProbe(const unsigned int ID)
{
	// Check whether nodes were specified
	if (this->isNodeEmpty())
	{
		throw std::runtime_error("Nodes should be specified before adding wave probes. In: ENVIR::addWaveProbe");
	}

	// Check if this wave probe was already included
	std::vector<unsigned int>::const_iterator iter = std::find(m_waveProbeID.begin(), m_waveProbeID.end(), ID);

	if (iter == m_waveProbeID.end())
	{
		// Get the node coordinate and add it to m_waveProbe
		m_waveProbe.push_back(this->getNode(ID));
		m_waveProbeID.push_back(ID);
	}
}

void ENVIR::evaluateWaveKinematics()
{
	typedef std::vector<cx_double> cx_stdvec;

	m_timeArray = arma::regspace(0, m_timeStep, m_timeTotal);

	// Due to roundoff errors, sometimes the size of the vector of wave components
	// is different from the size of the time simulation array in cases where IFFT
	// should be used. This can not happen, then the size of the time vector is fixed below.
	if (getFlagIFFT() && m_timeArray.size() != m_wave.size())
	{
		m_timeArray = arma::linspace(0, m_timeTotal, m_wave.size());
		double newTimeStep = m_timeArray.at(1) - m_timeArray.at(0);
		IO::print2log("Setting time step to " + std::to_string(newTimeStep) + "s to avoid incompatible sizes in IFFT calculation. Previous value was " + std::to_string(m_timeStep) + ".");
		m_timeStep = newTimeStep;
	}

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


	// 2nd Order quantities follow a slightly different logic.
	// If IFFT is going to be used, the frequencies are equally spaced and hence
	// the elements of the matrix of 2nd order coefficients can be grouped in a single
	// amplitude vector that is used as an input for the IFFT algorithm. The vector of
	// frequencies is the vector of the difference frequencies.
	// See 2011, Incorporating Irregular Nonlinear Waves in Coupled Simulation of Offshore Wind Turbines, Agarwal and Manuel, App. Ocean Research 
	//
	// If, however, IFFT is not going to be used, the wave frequencies are not necessarily
	// equally spaced, hence we need to perform a computationally expensive double sum
	if (IO::isOutputActive(IO::OUTFLAG_WAVE_PRES_2ND))
	{
		m_wavePress2ndArray = zeros(m_timeArray.size(), 1);
		for (int iProbe = 0; iProbe < m_waveProbeID.size(); ++iProbe)
		{
			if (getFlagIFFT())
			{
				cx_mat amp_dw(m_wave.size(), 1, fill::zeros); // Vector with the complex amplitude at the difference frequency
				for (int iWave = 0; iWave < m_wave.size(); ++iWave)
				{
					if (m_wave.at(iWave).amp() == 0) continue;
					for (int jWave = 0; jWave <= iWave; ++jWave)
					{
						if (m_wave.at(jWave).amp() == 0) continue;
						cx_double aux = wavePressure_2ndOrd_coef(m_waveProbe.at(iProbe), iWave, jWave);
						if (iWave != jWave)
						{
							aux *= 2; // Because we are summing only the upper part of the matrix
						}
						amp_dw.at(iWave - jWave) += aux;
					}
				}

				m_wavePress2ndArray.col(iProbe) = m_wave.size() * mkl_ifft_real(amp_dw);
			}
			else
			{
				for (int it = 0; it < m_timeArray.size(); ++it)
				{
					for (int iWave = 0; iWave < m_wave.size(); ++iWave)
					{
						for (int jWave = 0; jWave <= iWave; ++jWave)
						{
							double w_ii{ m_wave.at(iWave).angFreq() }, w_jj{ m_wave.at(jWave).angFreq() };
							cx_double sinCos{ cos((w_ii - w_jj) * m_timeArray.at(it)), sin((w_ii - w_jj) * m_timeArray.at(it)) };
							cx_double amp_dw = wavePressure_2ndOrd_coef(m_waveProbe.at(iProbe), iWave, jWave);
							if (iWave != jWave)
							{
								amp_dw *= 2;
							}

							m_wavePress2ndArray.at(it, iProbe) += real(amp_dw * sinCos);
						}
					}
				}
			}

			m_wavePress2ndArray.col(iProbe) %= m_timeRampArray;
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

bool ENVIR::getFlagWindTurb() const
{
	return m_flagTurbWind;
}

bool ENVIR::getFlagWindUnif() const
{
	return m_flagUnifWind;
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

double ENVIR::wavePres2ndAtProbe(const unsigned int ID) const
{
	// This should never occur in production code, but may be useful for debugging
	if (m_timeArray.size() == 0)
	{
		return wavePressure_2ndOrd(m_waveProbe[ID]);
	}
	else
	{
		uword ind1 = m_ind4interp1;
		double p = m_wavePress2ndArray.at(ind1, ID);
		if (shouldInterp())
		{
			uword ind2 = getInd4interp2();
			p += (m_wavePress2ndArray.at(ind2, ID) - m_wavePress2ndArray.at(ind1, ID)) * (m_time - m_timeArray.at(ind1)) / (m_timeArray.at(ind2) - m_timeArray.at(ind1));
		}
		return p;
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

std::string ENVIR::getTurbFileName() const
{
	return m_turbFileName;
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
			IO::print2outLine(IO::OUTFLAG_WAVE_PRES_2ND, m_waveProbeID[ii], wavePres2ndAtProbe(ii));
	}
}

void ENVIR::printWindVelocity() const
{
	for (int ii = 0; ii < m_windProbe.size(); ++ii)
	{	
		vec::fixed<3> windVel(fill::zeros);
		this->windVel(windVel, m_windProbe.at(ii));
		IO::print2outLine(IO::OUTFLAG_WIND_VEL, m_windProbeID[ii], windVel);
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

	m_shouldInterp = true;
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

// Prof. Pesce:
// I would prefer a hyperbolic tangent modulation... It goes smoothly and asymptotically to a constant. 
// In the case of the cosine modulation its second derivative with respect to time is non-null at t=Tr 
// and this is a source of a localized impulse.
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
double ENVIR::waveElev(const double x, const double y) const
{
	return ramp() * waveElev(x, y, m_time);
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

cx_double ENVIR::waveElev_coef(const double x, const double y, const unsigned int waveIndex) const
{
	double A = m_wave.at(waveIndex).amp();
	double k = m_wave.at(waveIndex).waveNumber();
	double beta = m_wave.at(waveIndex).direction() * arma::datum::pi / 180.;
	double phase = m_wave.at(waveIndex).phase() * arma::datum::pi / 180.;

	return { A*cos(k*cos(beta)*x + k * sin(beta)*y + phase) , -A * sin(k*cos(beta)*x + k * sin(beta)*y + phase) };
}

double ENVIR::wavePressure(const vec::fixed<3> &coord) const
{
	return ramp() * wavePressure(coord, m_time);
}

double ENVIR::wavePressure(const vec::fixed<3>& coord, const double time) const
{
	double p(0);

	for (int ii = 0; ii < m_wave.size(); ++ii)
	{
		if (m_wave.at(ii).amp() != 0)
		{
			double w = m_wave.at(ii).angFreq();
			p += (wavePressure_coef(coord.at(0), coord.at(1), coord.at(2), ii) * cx_double { cos(w * time), sin(w * time) }).real();
		}
	}

	return p;
}

cx_double ENVIR::wavePressure_coef(const double x, const double y, const double z, const unsigned int waveIndex) const
{
	cx_double p(0);

	double h = m_watDepth;
	double rho = m_watDens;
	double g = m_gravity;
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

		p = cx_double(cos(k*cosBeta*x + k * sinBeta*y + phase), -sin(k*cosBeta*x + k * sinBeta*y + phase));
		p *= rho * g * A * khz;
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

		vel.row(0) = cx_double(cosBeta * khz_xy * cos(k*cosBeta*x + k * sinBeta*y + phase), -cosBeta * khz_xy * sin(k*cosBeta*x + k * sinBeta*y + phase));
		vel.row(1) = cx_double(sinBeta * khz_xy * cos(k*cosBeta*x + k * sinBeta*y + phase), -sinBeta * khz_xy * sin(k*cosBeta*x + k * sinBeta*y + phase));
		vel.row(2) = cx_double(khz_z * sin(k*cosBeta*x + k * sinBeta*y + phase), khz_z * cos(k*cosBeta*x + k * sinBeta*y + phase));
	}

	return w * A * vel;
}

double ENVIR::wavePressure_2ndOrd(const vec::fixed<3> &coord) const
{
	return ramp() * wavePressure_2ndOrd(coord, m_time);
}

double ENVIR::wavePressure_2ndOrd(const vec::fixed<3> &coord, const double time) const
{
	double p{ 0 };

	// When i == j, p = 0, so it is safe to skip this part of the loop.
	// Besides, as only the real part of the second-order difference-frequency potential is used,
	// the acceleration due to a pair ij is equal to ji.
	for (unsigned int ii = 0; ii < m_wave.size(); ++ii)
	{
		for (unsigned int jj = ii + 1; jj < m_wave.size(); ++jj)
		{
			double w_ii{ m_wave.at(ii).angFreq() }, w_jj{ m_wave.at(jj).angFreq() };
			cx_double sinCos({ cos((w_ii - w_jj) * time), sin((w_ii - w_jj) *time) });

			p += 2 * real(wavePressure_2ndOrd_coef(coord, ii, jj) * sinCos);
		}
	}

	return p;
}

cx_double ENVIR::wavePressure_2ndOrd_coef(const vec::fixed<3>& coord, const unsigned int waveIndex1, const unsigned int waveIndex2) const
{
	cx_double p{ 0 };

	// In this case, the calculation below would lead to 0/0, but the limit is actually 0. This is fine, as the second order potential should not contribute to the mean drift.
	if (waveIndex1 == waveIndex2)
	{
		return p;
	}

	// More friendly notation
	double z = coord[2];
	double h = m_watDepth;
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

		double kh{ 0 };
		if (norm_k1_k2 * h >= 10)
		{
			kh = exp(norm_k1_k2 * z);

		}
		else
		{
			kh = std::cosh(norm_k1_k2 * (z + h)) / std::cosh(norm_k1_k2 * h);
		}

		p = cx_double(cos(dot(k1_k2, coord) + phase1 - phase2), -sin(dot(k1_k2, coord) + phase1 - phase2));
		p *= 0.5 * m_watDens * A1 * A2 * (w1 - w2) * g*g * aux * kh;
	}

	return p;
}

vec::fixed<3> ENVIR::u1(const vec::fixed<3> &coord, const double zwl) const
{
	return ramp() * u1(coord, zwl, m_time);
}

vec::fixed<3> ENVIR::u1(const vec::fixed<3>& coord, const double zwl, const double time) const
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
	return vel;
}

vec::fixed<3> ENVIR::du1dt(const vec::fixed<3> &coord, const double zwl) const
{
	return ramp() * du1dt(coord, zwl, m_time);
}

vec::fixed<3> ENVIR::du1dt(const vec::fixed<3>& coord, const double zwl, const double time) const
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
		if (m_wave.at(ii).amp() != 0)
		{
			double w = m_wave.at(ii).angFreq();
			acc += real(du1dt_coef(coord.at(0), coord.at(1), z, ii) * cx_double { cos(w * m_time), sin(w * m_time) });
		}
	}

	return acc;
}

cx_vec::fixed<3> ENVIR::du1dt_coef(const double x, const double y, const double z, const unsigned int waveIndex) const
{
	cx_vec::fixed<3> acc(fill::zeros);

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

		acc[0] = cx_double(sin(k*cosBeta*x + k * sinBeta*y + phase), cos(k*cosBeta*x + k * sinBeta*y + phase));
		acc[0] *= w * w * A * khz_xy * cosBeta;

		acc[1] = cx_double(sin(k*cosBeta*x + k * sinBeta*y + phase), cos(k*cosBeta*x + k * sinBeta*y + phase));
		acc[1] *= w * w * A * khz_xy * sinBeta;

		acc[2] = cx_double(cos(k*cosBeta*x + k * sinBeta*y + phase), -sin(k*cosBeta*x + k * sinBeta*y + phase));
		acc[2] *= -w * w * A * khz_z;
	}

	return acc;
}

vec::fixed<3> ENVIR::du2dt(const vec::fixed<3> &coord) const
{
	return ramp() * du2dt(coord, m_time);
}

vec::fixed<3> ENVIR::du2dt(const vec::fixed<3>& coord, const double time) const
{
	arma::vec::fixed<3> acc = { 0,0,0 };

	// When i == j, p = 0, so it is safe to skip this part of the loop.
	// Besides, as only the real part of the second-order difference-frequency potential is used,
	// the acceleration due to a pair ij is equal to ji.
	for (unsigned int ii = 0; ii < m_wave.size(); ++ii)
	{
		for (unsigned int jj = ii + 1; jj < m_wave.size(); ++jj)
		{
			double w_ii{ m_wave.at(ii).angFreq() }, w_jj{ m_wave.at(jj).angFreq() };
			cx_double sinCos({ cos((w_ii - w_jj) * time), sin((w_ii - w_jj) *time) });

			acc += 2 * real(du2dt_coef(coord, ii, jj) * sinCos);
		}
	}

	return acc;
}

cx_vec::fixed<3> ENVIR::du2dt_coef(const vec::fixed<3>& coord, const unsigned int waveIndex1, const unsigned int waveIndex2) const
{
	cx_vec::fixed<3> acc(fill::zeros);

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

		acc.at(0) = cx_double(sin(dot(k1_k2, coord) + p1 - p2), cos(dot(k1_k2, coord) + p1 - p2));
		acc.at(0) *= 0.5 * A1 * A2 * (w1 - w2) * (k1 * cosB1 - k2 * cosB2) * g*g * aux * khz_xy;

		acc.at(1) = cx_double(sin(dot(k1_k2, coord) + p1 - p2), cos(dot(k1_k2, coord) + p1 - p2));
		acc.at(1) *= 0.5 * A1 * A2 * (w1 - w2) * (k1 * sinB1 - k2 * sinB2) * g*g * aux * khz_xy;

		acc.at(2) = cx_double(cos(dot(k1_k2, coord) + p1 - p2), -sin(dot(k1_k2, coord) + p1 - p2));
		acc.at(2) *= -0.5 * A1 * A2 * (w1 - w2) * g*g * aux * norm_k1_k2 * khz_z;
	}

	return acc;
}

vec::fixed<3> ENVIR::du1dx(const vec::fixed<3> &coord, const double zwl) const
{
	return ramp() * du1dx(coord, zwl, m_time);
}

vec::fixed<3> ENVIR::du1dx(const vec::fixed<3>& coord, const double zwl, const double time) const
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
		if (m_wave.at(ii).amp() != 0)
		{
			double w = m_wave.at(ii).angFreq();
			acc += real(du1dx_coef(coord.at(0), coord.at(1), z, ii) * cx_double { cos(w * m_time), sin(w * m_time) });
		}
	}

	return acc;
}

cx_vec::fixed<3> ENVIR::du1dx_coef(const double x, const double y, const double z, const unsigned int waveIndex) const
{
	cx_vec::fixed<3> acc(fill::zeros);

	// More friendly notation
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

		cx_double i(0, 1);
		cx_double aux = w * A * k * cosBeta * cx_double(sin(k*cosBeta*x + k * sinBeta*y + phase), cos(k*cosBeta*x + k * sinBeta*y + phase));
		acc.at(0) = -aux * khz_xy * cosBeta;
		acc.at(1) = -aux * khz_xy * sinBeta;
		acc.at(2) = -aux * i * khz_z;
	}

	return acc;
}

vec::fixed<3> ENVIR::du1dy(const vec::fixed<3> &coord, const double zwl) const
{
	return ramp() * du1dy(coord, zwl, m_time);
}

vec::fixed<3> ENVIR::du1dy(const vec::fixed<3>& coord, const double zwl, const double time) const
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
		if (m_wave.at(ii).amp() != 0)
		{
			double w = m_wave.at(ii).angFreq();
			acc += real(du1dy_coef(coord.at(0), coord.at(1), z, ii) * cx_double { cos(w * m_time), sin(w * m_time) });
		}
	}
	return acc;
}

cx_vec::fixed<3> ENVIR::du1dy_coef(const double x, const double y, const double z, const unsigned int waveIndex) const
{
	cx_vec::fixed<3> acc(fill::zeros);

	// More friendly notation
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

		cx_double i(0, 1);
		cx_double aux = w * A * k * sinBeta * cx_double(sin(k*cosBeta*x + k * sinBeta*y + phase), cos(k*cosBeta*x + k * sinBeta*y + phase));
		acc.at(0) = -aux * khz_xy * cosBeta;
		acc.at(1) = -aux * khz_xy * sinBeta;
		acc.at(2) = -aux * i * khz_z;
	}

	return acc;
}

vec::fixed<3> ENVIR::du1dz(const vec::fixed<3> &coord, const double zwl) const
{
	return ramp() * du1dz(coord, zwl, m_time);
}

vec::fixed<3> ENVIR::du1dz(const vec::fixed<3>& coord, const double zwl, const double time) const
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
		if (m_wave.at(ii).amp() != 0)
		{
			double w = m_wave.at(ii).angFreq();
			acc += real(du1dz_coef(coord.at(0), coord.at(1), z, ii) * cx_double { cos(w * m_time), sin(w * m_time) });
		}
	}
	return acc;
}

cx_vec::fixed<3> ENVIR::du1dz_coef(const double x, const double y, const double z, const unsigned int waveIndex) const
{
	cx_vec::fixed<3> acc(fill::zeros);

	// More friendly notation
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

		cx_double i(0, 1);
		cx_double aux = w * A * k * cx_double(cos(k*cosBeta*x + k * sinBeta*y + phase), -sin(k*cosBeta*x + k * sinBeta*y + phase));
		acc.at(0) = aux * khz_z * cosBeta;
		acc.at(1) = aux * khz_z * sinBeta;
		acc.at(2) = aux * i * khz_xy;
	}

	return acc;
}

vec::fixed<3> ENVIR::da1dx(const vec::fixed<3> &coord, const double zwl) const
{
	return ramp() * da1dx(coord, zwl, m_time);
}

vec::fixed<3> ENVIR::da1dx(const vec::fixed<3> &coord, const double zwl, const double time) const
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
		if (m_wave.at(ii).amp() != 0)
		{
			double w = m_wave.at(ii).angFreq();
			acc += real(da1dx_coef(coord.at(0), coord.at(1), z, ii) * cx_double { cos(w * m_time), sin(w * m_time) });
		}
	}
	return acc;
}

cx_vec::fixed<3> ENVIR::da1dx_coef(const double x, const double y, const double z, const unsigned int waveIndex) const
{
	cx_vec::fixed<3> acc(fill::zeros);

	// More friendly notation
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

		acc.at(0) = cx_double(cos(k*cosBeta*x + k * sinBeta*y + phase), -sin(k*cosBeta*x + k * sinBeta*y + phase));
		acc.at(0) *= w * w * A * k * khz_xy * cosBeta * cosBeta;

		acc.at(1) = cx_double(cos(k*cosBeta*x + k * sinBeta*y + phase), -sin(k*cosBeta*x + k * sinBeta*y + phase));
		acc.at(1) *= w * w * A * k * khz_xy * cosBeta * sinBeta;

		acc.at(2) = cx_double(sin(k*cosBeta*x + k * sinBeta*y + phase), cos(k*cosBeta*x + k * sinBeta*y + phase));
		acc.at(2) *= w * w * A * k * khz_z * cosBeta;
	}

	return acc;
}

vec::fixed<3> ENVIR::da1dy(const vec::fixed<3> &coord, const double zwl) const
{
	return ramp() * da1dy(coord, zwl, m_time);
}

vec::fixed<3> ENVIR::da1dy(const vec::fixed<3> &coord, const double zwl, const double time) const
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
		if (m_wave.at(ii).amp() != 0)
		{
			double w = m_wave.at(ii).angFreq();
			acc += real(da1dy_coef(coord.at(0), coord.at(1), z, ii) * cx_double { cos(w * m_time), sin(w * m_time) });
		}
	}
	return acc;
}

cx_vec::fixed<3> ENVIR::da1dy_coef(const double x, const double y, const double z, const unsigned int waveIndex) const
{
	cx_vec::fixed<3> acc(fill::zeros);

	// More friendly notation
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

		acc.at(0) = cx_double(cos(k*cosBeta*x + k * sinBeta*y + phase), -sin(k*cosBeta*x + k * sinBeta*y + phase));
		acc.at(0) *= w * w * A * k * khz_xy * cosBeta * sinBeta;

		acc.at(1) = cx_double(cos(k*cosBeta*x + k * sinBeta*y + phase), -sin(k*cosBeta*x + k * sinBeta*y + phase));
		acc.at(1) *= w * w * A * k * khz_xy * sinBeta * sinBeta;

		acc.at(2) = cx_double(sin(k*cosBeta*x + k * sinBeta*y + phase), cos(k*cosBeta*x + k * sinBeta*y + phase));
		acc.at(2) *= w * w * A * k * khz_z * sinBeta;
	}

	return acc;
}

vec::fixed<3> ENVIR::da1dz(const vec::fixed<3> &coord, const double zwl) const
{
	return ramp() * da1dz(coord, zwl, m_time);
}

vec::fixed<3> ENVIR::da1dz(const vec::fixed<3> &coord, const double zwl, const double time) const
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
		if (m_wave.at(ii).amp() != 0)
		{
			double w = m_wave.at(ii).angFreq();
			acc += real(da1dz_coef(coord.at(0), coord.at(1), z, ii) * cx_double { cos(w * m_time), sin(w * m_time) });
		}
	}
	return acc;
}

cx_vec::fixed<3> ENVIR::da1dz_coef(const double x, const double y, const double z, const unsigned int waveIndex) const
{
	cx_vec::fixed<3> acc(fill::zeros);

	// More friendly notation
	double h = m_watDepth;
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

		acc.at(0) = cx_double(sin(k*cosBeta*x + k * sinBeta*y + phase), cos(k*cosBeta*x + k * sinBeta*y + phase));
		acc.at(0) *= w * w * A * k * khz_z * cosBeta;

		acc.at(1) = cx_double(sin(k*cosBeta*x + k * sinBeta*y + phase), cos(k*cosBeta*x + k * sinBeta*y + phase));
		acc.at(1) *= w * w * A * k * khz_z * sinBeta;

		acc.at(2) = cx_double(cos(k*cosBeta*x + k * sinBeta*y + phase), -sin(k*cosBeta*x + k * sinBeta*y + phase));
		acc.at(2) *= -w * w * A * k * khz_xy;
	}

	return acc;
}

vec::fixed<3> ENVIR::gradP1(const vec::fixed<3>& coord, const double zwl) const
{
	return ramp() * gradP1(coord, zwl, m_time);
}

vec::fixed<3> ENVIR::gradP1(const vec::fixed<3>& coord, const double zwl, const double time) const
{
	vec::fixed<3> gradP1(fill::zeros);

	if (m_waveStret <= 1 && zwl != 0)
	{
		throw std::runtime_error("zwl = " + std::to_string(zwl) + "is incompatible with wave streching mode " + std::to_string(m_waveStret));
	}

	double z = coord.at(2);
	if (z > zwl)
	{
		return gradP1;
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
			gradP1 += real(gradP1_coef(coord.at(0), coord.at(1), z, ii) * cx_double { cos(w * m_time), sin(w * m_time) });
		}
	}
	return gradP1;
}

cx_vec::fixed<3> ENVIR::gradP1_coef(const double x, const double y, const double z, const unsigned int waveIndex) const
{
	cx_vec::fixed<3> gradP(fill::zeros);

	// More friendly notation
	double h = m_watDepth;
	double rho = m_watDens;
	double g = m_gravity;
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

		gradP.at(0) = k*cosBeta*khz_xy*cx_double(-sin(k*cosBeta*x + k * sinBeta*y + phase), -cos(k*cosBeta*x + k * sinBeta*y + phase));
		gradP.at(1) = k*sinBeta*khz_xy*cx_double(-sin(k*cosBeta*x + k * sinBeta*y + phase), -cos(k*cosBeta*x + k * sinBeta*y + phase));
		gradP.at(2) = k*khz_z*cx_double(cos(k*cosBeta*x + k * sinBeta*y + phase), -sin(k*cosBeta*x + k * sinBeta*y + phase));
		gradP *= rho * g * A;
	}

	return gradP;
}

void ENVIR::windVel(vec::fixed<3> &windVel, const vec::fixed<3> &coord) const
{
	if (m_flagTurbWind)
	{				
		// Distance from a point P to a line whose vertices are v1 and v2 can be calculated by
		// ds = ||cross(P-v1, v2-v1)|| / norm(v2-v1)
		// This is what is done below, but the math was developed for the special case at hand
		// in order to save computational time
		double dx0 = coord.at(0) - m_windGrid_x.at(0);
		double dy0 = coord.at(1) - m_windGrid_y.at(0);
		double ds = dx0 * m_windDirCos - dy0 * m_windDirSin;

		// Find the point at the grid that corresponds to "coord".
		// It is P = (x,y) = coord - ds * V/||V||, with V/||V|| simply the unitary vector
		// pointing at the direction of the wind propagation.
		double x = coord.at(0) - ds * m_windDirCos;
		double y = coord.at(1) + ds * m_windDirSin;
		
		if (y < m_windGrid_y.front() || y > m_windGrid_y.back() ||
			coord.at(2) < m_windGrid_z.front() || coord.at(2) > m_windGrid_z.back() )
		{
			throw std::runtime_error("Point (" + std::to_string(x) + "," + std::to_string(y) + "," + std::to_string(coord.at(2)) + "), obtained from (" +
				std::to_string(coord.at(0)) + "," + std::to_string(coord.at(1)) + "," + std::to_string(coord.at(2)) +
				"), is outside of turbulence grid domain (function ENVIR::windVel).");
		}

		// This point is probably not one of the grid nodes, hence, it is necessary to interpolate.
		// To do so, find the indices of the lower bounds that will be used in the interpolation
		// for each of the grid vectors.
		// 
		// Since we are dealing with a plane given by the Z coordinate + the line described by the vectors
		// m_windGrid_x and m_windGrid_y, only one of the horizontal coordinates is actually needed, as the
		// index for the other one would be the same.
		arma::uword ind1_y = arma::index_min(abs(m_windGrid_y - y));
		if (m_windGrid_y(ind1_y) > y)
		{
			ind1_y -= 1;
		}

		arma::uword ind1_z = arma::index_min(abs(m_windGrid_z - coord.at(2)));
		if (m_windGrid_z(ind1_z) > coord.at(2))
		{
			ind1_z -= 1;
		}
		
		// The turbulence at "coord" is given by the one at the point found above, but not at the current time step.
		// Based on Taylor frozen turbulence hypothesis, the time taken for the turbulent wind to go from
		// the grid to the desired point "coord" is simply dt = ds/||V||. If this is a positive value,
		// it means that the wind at "coord" is the same that was at the grid in t-dt. Otherwise, it means that the wind
		// at coord will be at P in t+dt. Since the grid is located downwind of the turbine, the latter is expected to occur.
		double dt = ds / m_windRefVel;
		double t = time() - dt;
		

		if (m_flagPeriodicTurb)
		{
			double intPart{ 0 };
			double fracPart = modf(t / m_windTimeTotal, &intPart);
			t = fracPart * m_windTimeTotal;
			t += ((t < 0) ? m_windTimeTotal : 0); // If t is negative, we are walking backwards from the end of the time array
		}
		else
		{	
			if (t > m_windTimeTotal)
			{
				throw std::runtime_error("Required simulation time " + std::to_string(t) +
					"s is outside range of turbulence field (" + std::to_string(m_windTimeTotal) + "s).");
			}
			else if (t < 0)
			{
				throw std::runtime_error("Required simulation time for turbulence " + std::to_string(t) +
					"s is negative, meaning that a larger turbulence grid is required.");
			}
		}

		// Time also needs to be interpolated
		int ind1_t = floor( t / m_windDt);		
		double t1 = ind1_t * m_windDt;
		double t2 = t1 + m_windDt;		

		double auxVelocity_t1{ 0 };
		double auxVelocity_t2{ 0 };
		for (int ii = 0; ii < 3; ++ii)
		{			
			// Spatial interpolation at t1
			auxVelocity_t1 = bilinearInterpolation(m_windGrid_y(ind1_y), m_windGrid_y(ind1_y + 1), m_windGrid_z(ind1_z), m_windGrid_z(ind1_z + 1),
				                                  m_windVelocity(ind1_t)(ii, ind1_y, ind1_z), m_windVelocity.at(ind1_t).at(ii, ind1_y, ind1_z + 1),
				                                  m_windVelocity(ind1_t)(ii, ind1_y + 1, ind1_z), m_windVelocity.at(ind1_t).at(ii, ind1_y + 1, ind1_z + 1),
				                                  y, coord.at(2));

			// Spatial interpolation at t2
			auxVelocity_t2 = bilinearInterpolation(m_windGrid_y(ind1_y), m_windGrid_y(ind1_y + 1), m_windGrid_z(ind1_z), m_windGrid_z(ind1_z + 1),
                                                   m_windVelocity.at(ind1_t + 1).at(ii, ind1_y, ind1_z), m_windVelocity.at(ind1_t + 1).at(ii, ind1_y, ind1_z + 1),
                                                   m_windVelocity.at(ind1_t + 1).at(ii, ind1_y + 1, ind1_z), m_windVelocity.at(ind1_t + 1).at(ii, ind1_y + 1, ind1_z + 1),
                                                   y, coord.at(2));
			// Temporal interpolation
			windVel.at(ii) = auxVelocity_t1 + (auxVelocity_t2 - auxVelocity_t1) * (t - t1) / (t2 - t1);
		}
	}
	else if (m_flagUnifWind)
	{
		// Declare the variables that will be used to compute the wind velocities at the current time step
		float U{ 0 };
		float theta{ 0 };
		float W{ 0 };
		float hSheer{ 0 };
		float windExp{ 0 };
		float vSheer{ 0 };
		float ugust{ 0 };

		// If we are below the first record in the file, quantities are constant and equal to the first record
		if (m_time < m_wnd_time[0])
		{
			U = m_wnd_windSpeed[0];
			theta = m_wnd_windDir[0];
			W = m_wnd_vertSpeed[0];
			hSheer = m_wnd_horizSheer[0];
			windExp = m_wnd_windExp[0];
			vSheer = m_wnd_linVertSheer[0];
			ugust = m_wnd_gustSpeed[0];
		}

		// Same thing if it is above the last record
		else if (m_time > m_wnd_time.back())
		{
			U = m_wnd_windSpeed.back();
			theta = m_wnd_windDir.back();
			W = m_wnd_vertSpeed.back();
			hSheer = m_wnd_horizSheer.back();
			windExp = m_wnd_windExp.back();
			vSheer = m_wnd_linVertSheer.back();
			ugust = m_wnd_gustSpeed.back();
		}

		// Otherwise, quantities are interpolated between the values provided in m_wnd_time
		else
		{
			// The first step is to get this index
			arma::uword i1 = index_min(abs(m_wnd_time - m_time));
			arma::uword i2;

			// Make sure i1 is the lower index and i2 the upper index
			if (m_wnd_time(i1) <= m_time)
			{
				i2 = i1 + 1;
			}
			else
			{
				i2 = i1;
				i1 = i2 - 1;
			}
			float t1 = m_wnd_time(i1);
			float t2 = m_wnd_time(i2);

			// Interpolate values
			U     = m_wnd_windSpeed.at(i1) + (m_wnd_windSpeed.at(i2) - m_wnd_windSpeed.at(i1)) * (m_time - t1) / (t2 - t1);
			theta = m_wnd_windDir.at(i1) + (m_wnd_windDir.at(i2) - m_wnd_windDir.at(i1)) * (m_time - t1) / (t2 - t1);
			W = m_wnd_vertSpeed.at(i1) + (m_wnd_vertSpeed.at(i2) - m_wnd_vertSpeed.at(i1)) * (m_time - t1) / (t2 - t1);
			hSheer = m_wnd_horizSheer.at(i1) + (m_wnd_horizSheer.at(i2) - m_wnd_horizSheer.at(i1)) * (m_time - t1) / (t2 - t1);
			windExp = m_wnd_windExp.at(i1) + (m_wnd_windExp.at(i2) - m_wnd_windExp.at(i1)) * (m_time - t1) / (t2 - t1);
			vSheer = m_wnd_linVertSheer.at(i1) + (m_wnd_linVertSheer.at(i2) - m_wnd_linVertSheer.at(i1)) * (m_time - t1) / (t2 - t1);
			ugust = m_wnd_gustSpeed.at(i1) + (m_wnd_gustSpeed.at(i2) - m_wnd_gustSpeed.at(i1)) * (m_time - t1) / (t2 - t1);
		}
		float sinTheta = std::sin(theta*arma::datum::pi / 180.);
		float cosTheta = std::cos(theta*arma::datum::pi / 180.);

		// Compute horizontal velocity using the formula from InflowWind User's Guide (2016)
		U = U * pow(coord[2] / m_windRefHeight, windExp) +
			U * (hSheer / m_wnd_refLength) * (coord[0] * sinTheta + coord[1] * cosTheta) +
			U * (vSheer / m_wnd_refLength) * (coord[2] - m_windRefHeight) +
			ugust;

		windVel.at(0) =  U * cosTheta;
		windVel.at(1) = -U * sinTheta;
		windVel.at(2) =  W;		
	}
	else
	{		
		double U = ramp() * windRefVel()  * pow(coord[2] / m_windRefHeight, m_windExp);
		windVel.at(0) = U * m_windDirCos;
		windVel.at(1) = -U * m_windDirSin;
		windVel.at(2) = 0;
	}
}

mat ENVIR::timeSeriesFromAmp(cx_mat &inAmp, const vec &w) const
{
	mat out(m_timeArray.size() , inAmp.n_cols, fill::zeros);
	uword ncols = inAmp.n_cols;
	if (getFlagIFFT())
	{
		out = m_wave.size() * mkl_ifft_real(inAmp) % repmat(getRampArray(), 1, inAmp.n_cols);
	}
	else
	{
		for (unsigned int it = 0; it < m_timeArray.size(); ++it)
		{
			cx_vec sinCos{ cos(w * m_timeArray.at(it)), sin(w * m_timeArray.at(it)) };

			out.row(it) = sum(real(inAmp % repmat(sinCos, 1, inAmp.n_cols)), 0);
		}

		out %= repmat(getRampArray(), 1, inAmp.n_cols);
	}

	return out;
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