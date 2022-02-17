// This is a complete mess. I need to take some time and organize this class.
// - The function readInputFile has a lot of repetition. I should make a template that reads the input line based on the aux variable used in each block. Or overload a function that does it,

#include "IO.h"

#include <fstream> // Include file input/output classes
#include <iostream>
#include <sstream>
#include <stdexcept> // For std::exception
#include <iomanip> // For input/output manipulators
#include <sys/types.h> // Using this and stat.h to check if directories and files exist
#include <sys/stat.h>
#include <cstdlib> // for exit()
#include <array>


/*****************************************************
Defining and initializing static member variables
*****************************************************/
std::string IO::m_inFilePath = "";
std::ifstream IO::m_inFl;
unsigned int IO::m_inLineNumber = 0;

std::string IO::m_logFilePath = "";
std::ofstream IO::m_logFl;

std::string IO::m_sumFilePath = "";
std::ofstream IO::m_sumFl;

std::string IO::m_outFilePath = "";
std::ofstream IO::m_outFl;
const unsigned int IO::m_outColumnWidth = 25;
const unsigned int IO::m_outNumPrecision = 4;
std::array<bool, IO::OUTFLAG_SIZE> IO::m_whichResult2Output;

std::stringstream IO::m_outLineHeader;
std::stringstream IO::m_outLine;
bool IO::m_shouldWriteOutLineHeader = false;
bool IO::m_shouldWriteOutLine = false;



/*****************************************************
Class functions related to Input
*****************************************************/

//	This is the main input function.
//	It is responsible for reading the input file line by line and assigning what is read to the FOWT and ENVIR classes
void IO::readInputFile(FOWT &fowt, ENVIR &envir)
{
	// Classes that are members of FOWT and ENVIR
	Floater floater;
	RNA rna;

	// Mean displacement for hydrodynamic calculations
	vec::fixed<6> meanDisp(fill::zeros);

	// Initial displacement and velocity
	vec::fixed<6> disp0(fill::zeros);
	vec::fixed<6> vel0(fill::zeros);

	// Read file line by line
	std::string strInput{""};
	while (m_inFl)
	{
		IO::readLineInputFile(strInput);

		/**************************
		Read data based on keywords
		**************************/
		/*
			Read data to envir
		*/
		if (caseInsCompare(getKeyword(strInput), "TimeStep"))
		{
			envir.setTimeStep(string2num<double>(getData(strInput)));
		}

		else if (caseInsCompare(getKeyword(strInput), "PrintStep"))
		{
			envir.setPrintStep(string2num<double>(getData(strInput)));
		}

		else if (caseInsCompare(getKeyword(strInput), "TimeTotal"))
		{
			envir.setTimeTotal(string2num<double>(getData(strInput)));
		}

		else if (caseInsCompare(getKeyword(strInput), "TimeRamp"))
		{
			envir.setTimeRamp(string2num<double>(getData(strInput)));
		}

		else if (caseInsCompare(getKeyword(strInput), "Grav"))
		{
			envir.setGravity(string2num<double>(getData(strInput)));
		}

		else if (caseInsCompare(getKeyword(strInput), "WatDens"))
		{
			envir.setWatDens(string2num<double>(getData(strInput)));
		}

		else if (caseInsCompare(getKeyword(strInput), "WatDepth"))
		{
			envir.setWatDepth(string2num<double>(getData(strInput)));
		}

		else if (caseInsCompare(getKeyword(strInput), "WaveStret"))
		{
			envir.setWaveStret(string2num<unsigned int>(getData(strInput)));
		}

		else if (caseInsCompare(getKeyword(strInput), "AirDens"))
		{
			envir.setAirDens(string2num<double>(getData(strInput)));
		}

		else if (caseInsCompare(getKeyword(strInput), "WindTurbFile"))
		{
			if (arma::is_finite(envir.windRefVel()) || arma::is_finite(envir.windRefHeight()) || arma::is_finite(envir.windExp()))
			{
				throw std::runtime_error("Cannot read turbulence file given in input line " + std::to_string(IO::getInLineNumber()) + " because either WindHeight, WindVel or WindExp were also specified.");
			}
			envir.setWindFromTurbFile(getData(strInput));
		}

		else if (caseInsCompare(getKeyword(strInput), "WindVel"))
		{
			if (envir.getFlagWindTurb())
			{
				throw std::runtime_error("Cannot set WindVel given in input line " + std::to_string(IO::getInLineNumber()) + " because a turbulence file was already provided.");
			}
			envir.setWindRefVel(string2num<double>(getData(strInput)));
		}

		else if (caseInsCompare(getKeyword(strInput), "WindDir"))
		{
			envir.setWindDir(string2num<double>(getData(strInput)));
		}

		else if (caseInsCompare(getKeyword(strInput), "WindHeight"))
		{
			if (envir.getFlagWindTurb())
			{
				throw std::runtime_error("Cannot set WindHeight given in input line " + std::to_string(IO::getInLineNumber()) + " because a turbulence file was already provided.");
			}
			envir.setWindRefHeight(string2num<double>(getData(strInput)));
		}

		else if (caseInsCompare(getKeyword(strInput), "WindExp"))
		{
			if (envir.getFlagWindTurb())
			{
				throw std::runtime_error("Cannot set WindExp given in input line " + std::to_string(IO::getInLineNumber()) + " because a turbulence file was already provided.");
			}
			envir.setWindExp(string2num<double>(getData(strInput)));
		}

		// Waves, nodes and other inputs, are special because they have
		// several lines to specify their characteristics.
		// In theses cases, we need to loop all the lines and read the inputs
		// line by line.
		else if (caseInsCompare(getKeyword(strInput), "Wave"))
		{
			// Read next line, since current line is just the main keyword
			IO::readLineInputFile(strInput);

			// Loop the lines in the Wave section of the input file and add the waves to envir
			while (!caseInsCompare(getKeyword(strInput), "END")) // The END keyword indicates the end of the list of waves
			{
				if (!m_inFl) // Signal if the end of file is reached before the end keyword
				{
					throw std::runtime_error("End of file reached before END keyword in WAVE specification.");
					return;
				}

				// Check if this line is a regular wave
				if (caseInsCompare(getKeyword(strInput), "TRWave") || caseInsCompare(getKeyword(strInput), "FRWave") || caseInsCompare(getKeyword(strInput), "WRWave"))
				{
					// Need 4 inputs separated by a space or a tab:
					// wave height, period or frequency or angular frequency, direction of propagation, and phase
					std::vector<std::string> input = stringTokenize(getData(strInput), " \t");
					if (input.size() != 4)
					{
						throw std::runtime_error("Unable to read the wave in input line " + std::to_string(IO::getInLineNumber()) + ". Wrong number of parameters.");
					}

					// Can not specify individual wave components if an option that uses IFFT to calculate wave kinematics was already specified
					if (envir.getFlagIFFT())
					{
						throw std::runtime_error("Wave options that use IFFT to evaluate wave kinematics, such as JONSWAP without specifying the number of components or externally generated wave elevation, can not be specified with other waves. In ENVIR::addJonswap.");
					}

					double aux_height = string2num<double>(input.at(0));
					double aux_freqORperiod = string2num<double>(input.at(1));
					double aux_direction = string2num<double>(input.at(2));
					double aux_phase = string2num<double>(input.at(3));

					envir.addRegularWave(getKeyword(strInput), aux_height, aux_freqORperiod, aux_direction, aux_phase);
				}

				// Check if it is a JONSWAP spectrum
				else if (caseInsCompare(getKeyword(strInput), "JONSW"))
				{
					// Need at least 7 inputs, with an additional optional, separated by a space or a tab for defining the JONSWAP spectrum:
					// significant height, peak period, gamma, direction of propagation (at present, unidirectional seas are implemented)
					// lowest frequency limit (rad/s), highest frequency limit (rad/s) (beyond this limits, the spectrum is set to zero), 
					// seed for random number generator, and, optionally, the number of wave components.
					std::vector<std::string> input = stringTokenize(getData(strInput), " \t");
					if (input.size() != 7 && input.size() != 9)
					{
						throw std::runtime_error("Unable to read JONSWAP spectrum in input line " + std::to_string(IO::getInLineNumber()) + ". Wrong number of parameters.");
					}

					double aux_Hs = string2num<double>(input.at(0));
					double aux_Tp = string2num<double>(input.at(1));
					double aux_gamma = string2num<double>(input.at(2));
					double aux_dir = string2num<double>(input.at(3));
					double aux_wlow = string2num<double>(input.at(4));
					double aux_whigh = string2num<double>(input.at(5));
					
					// Use the seed for the RNG
					if (caseInsCompare(input.at(6), "?"))
					{
						arma::arma_rng::set_seed_random();
					}
					else
					{
						double seed = string2num<double>(input.at(6));
						arma::arma_rng::set_seed(seed);
					}

					// If the number of components is specified, the Equal Area method is used to discretize the wave spectrum
					// Otherwise, the number of components is calculated from the total simulation time
					int aux_nComponents = -1;
					double aux_dwMax = 0;
					if (input.size() == 9)
					{
						aux_nComponents = string2num<int>(input.at(7));
						aux_dwMax = string2num<double>(input.at(8));
						if (aux_dwMax <= 0)
						{
							throw std::runtime_error("Maximum delta omega specified in JONSWAP spectrum in input line " + std::to_string(IO::getInLineNumber()) + " must be > 0.");
						}
					}					
					envir.addJonswap(aux_Hs, aux_Tp, aux_gamma, aux_dir, aux_wlow, aux_whigh, aux_nComponents, aux_dwMax);
				}

				// Check if the input is a file with the wave elevation series at z = 0
				else if (caseInsCompare(getKeyword(strInput), "ELEV"))
				{
					// Need 4 inputs separated by a comma
					// Path to the file with the wave elevation series, wave direction, lowest frequency and high frequency limits 
					// (which determine the region outside of which the wave amplitude is set to zero)
					std::vector<std::string> input = stringTokenize(getData(strInput), ",");
					if (input.size() != 4)
					{
						throw std::runtime_error("Unable to read the wave in input line " + std::to_string(IO::getInLineNumber()) + 
							". Wrong number of parameters. Since the file can have whitespaces in its path, make sure you are using commas to separete the file path and the wave incidence.");
					}
					std::string elevFlPath = input.at(0);
					double waveDir = string2num<double>(input.at(1));
					double wlow = string2num<double>(input.at(2));
					double whigh = string2num<double>(input.at(3));
					

					envir.addWaveElevSeries(elevFlPath, waveDir, wlow, whigh);
				}

				// Otherwise, there could be a typo or something of the kind.
				else
				{
					throw std::runtime_error("Unknown wave type '" + getKeyword(strInput) + "' in line " + std::to_string(IO::getInLineNumber()) + ".");
				}

				IO::readLineInputFile(strInput);
			}
		}

		// Just like for the waves, a list of nodes is supposed to follow the "Nodes" keyword
		else if (caseInsCompare(getKeyword(strInput), "Nodes"))
		{
			// Read next line, since current line is just the main keyword
			IO::readLineInputFile(strInput);

			while (!caseInsCompare(getKeyword(strInput), "END"))
			{
				if (!m_inFl) // Signal if the end of file is reached before the end keyword
				{
					throw std::runtime_error("End of file reached before END keyword in NODES specification.");
					return;
				}

				// Nodes are specified by a vec with four components: ID, X coord, Y coord, and Z coord.
				// They are separated by commas in the input string.
				std::vector<std::string> input = stringTokenize(strInput, ",");
				if (input.size() != 4)
				{
					throw std::runtime_error("Unable to read the node in input line " + std::to_string(IO::getInLineNumber()) + ". Wrong number of parameters.");
					return;
				}

				unsigned int aux_nodeID = string2num<unsigned int>(input.at(0));
				double aux_nodeCoordX = string2num<double>(input.at(1));
				double aux_nodeCoordY = string2num<double>(input.at(2));
				double aux_nodeCoordZ = string2num<double>(input.at(3));
				envir.addNode(aux_nodeID, aux_nodeCoordX, aux_nodeCoordY, aux_nodeCoordZ);

				// Done with this line. Read next one.
				IO::readLineInputFile(strInput);
			}
		}

		/*
			Read data to fowt
		*/
		else if (caseInsCompare(getKeyword(strInput), "Hydro"))
		{
			fowt.setHydroMode(string2num<int>(getData(strInput)));
		}

		else if (caseInsCompare(getKeyword(strInput), "Aero"))
		{
			fowt.setAeroMode(string2num<int>(getData(strInput)));
		}

		else if (caseInsCompare(getKeyword(strInput), "Moor"))
		{
			fowt.setMoorMode(string2num<int>(getData(strInput)));
		}

		else if (caseInsCompare(getKeyword(strInput), "DOFS"))
		{
			// The flags for each of the six degrees of freedom are separated by white spaces in the input string (whitespace or tab)
			std::vector<std::string> input = stringTokenize(getData(strInput), " \t");
			if (input.size() != 6)
			{
				throw std::runtime_error("Unable to read the DoFs in input line " + std::to_string(IO::getInLineNumber()) + ". Wrong number of parameters.");
			}

			// Read to an auxiliar array before passing to fowt
			std::array<bool, 6> aux_activeDoFs;
			for (int ii = 0; ii < aux_activeDoFs.size(); ++ii)
			{
				aux_activeDoFs.at(ii) = string2num<bool>(input.at(ii));
			}

			fowt.setDoFs(aux_activeDoFs);
		}
		
		// The mooring line stiffness is a 6x6 matrix with columns separated by whitespaces or tabs.
		// Different rows correspond to different lines in the input file.
		else if (caseInsCompare(getKeyword(strInput), "ExtLinStiff"))
		{
		arma::mat::fixed<6, 6> aux_extStiff(fill::zeros);
		for (unsigned int countRows = 0; countRows < 6; ++countRows)
		{
			IO::readLineInputFile(strInput);

			std::vector<std::string> input = stringTokenize(strInput, " \t");
			if (input.size() != 6)
			{
				throw std::runtime_error("Unable to read row of linear stiffness matrix in input line " + std::to_string(IO::getInLineNumber()) + ". Wrong number of parameters.");
			}

			// Read to an auxiliar matrix before passing to fowt
			for (int ii = 0; ii < aux_extStiff.n_rows; ++ii)
			{
				aux_extStiff.at(countRows, ii) = string2num<double>(input.at(ii));
			}
		}
		fowt.setExtLinStiff(aux_extStiff);
		}

		// The external linear damping is a 6x6 matrix with columns separated by whitespaces or tabs.
		// Different rows correspond to different lines in the input file.
		else if (caseInsCompare(getKeyword(strInput), "extLinDamp"))
		{
		arma::mat::fixed<6, 6> aux_extDamp(fill::zeros);
		for (unsigned int countRows = 0; countRows < 6; ++countRows)
		{
			IO::readLineInputFile(strInput);

			std::vector<std::string> input = stringTokenize(strInput, " \t");
			if (input.size() != 6)
			{
				throw std::runtime_error("Unable to read row of linear damping matrix in input line " + std::to_string(IO::getInLineNumber()) + ". Wrong number of parameters.");
			}

			// Read to an auxiliar matrix before passing to fowt
			for (int ii = 0; ii < aux_extDamp.n_rows; ++ii)
			{
				aux_extDamp.at(countRows, ii) = string2num<double>(input.at(ii));
			}
		}
		fowt.setExtLinDamp(aux_extDamp);
		}


		else if (caseInsCompare(getKeyword(strInput), "ExtConstForce"))
		{
			// The 6 components of the external constant force are separated by commas in the input string
			std::vector<std::string> input = stringTokenize(getData(strInput), ",");
			if (input.size() != 6)
			{
				throw std::runtime_error("Unable to read external constant force in input line " + std::to_string(IO::getInLineNumber()) + ". Wrong number of parameters.");
			}

			// Read data to an auxiliary temporary variable
			vec::fixed<6> aux;
			for (int ii = 0; ii < aux.size(); ++ii)
			{
				aux.at(ii) = string2num<double>(input.at(ii));
			}

			fowt.setExtConstForce(aux);
		}

		else if (caseInsCompare(getKeyword(strInput), "FiltSlowDrift"))
		{
			// Need 2 inputs separated by a space or a tab:
			// angular frequency and damping levels (in % of critical damping)
			std::vector<std::string> input = stringTokenize(getData(strInput), " \t");
			if (input.size() != 2)
			{
				throw std::runtime_error("Unable to read slow drift filter in input line " + std::to_string(IO::getInLineNumber()) + ". Wrong number of parameters.");
			}
			fowt.setFilderSD(string2num<double>(input.at(0)), string2num<double>(input.at(1)));
		}

		// Read data to floater - which is a part of FOWT
		else if (caseInsCompare(getKeyword(strInput), "FloaterMass"))
		{
			floater.setMass(string2num<double>(getData(strInput)));
		}

		else if (caseInsCompare(getKeyword(strInput), "FloaterInertia"))
		{
			// The different components of the inertia matrix are separated by commas in the input string
			std::vector<std::string> input = stringTokenize(getData(strInput), ",");
			if (input.size() != 6)
			{
				throw std::runtime_error("Unable to read the floater inertia matrix in input line " + std::to_string(IO::getInLineNumber()) + ". Wrong number of parameters.");
			}

			vec::fixed<6> aux;
			for (int ii = 0; ii < aux.size(); ++ii)
			{
				aux(ii) = string2num<double>(input.at(ii));
			}

			floater.setInertia(aux);
		}

		else if (caseInsCompare(getKeyword(strInput), "FloaterCoG"))
		{
			// The coordinates of the CoG are separated by commas in the input string
			std::vector<std::string> input = stringTokenize(getData(strInput), ",");
			if (input.size() != 3)
			{
				throw std::runtime_error("Unable to read the CoG in input line " + std::to_string(IO::getInLineNumber()) + ". Wrong number of parameters.");
			}

			vec::fixed<3> aux;
			for (int ii = 0; ii < aux.size(); ++ii)
			{
				aux(ii) = string2num<double>(input.at(ii));
			}

			floater.setCoG(aux);
		}

		else if (caseInsCompare(getKeyword(strInput), "meanDisp"))
		{
		// The components of the mean displacement of the FOWT for hydrodynamic calculation are separated by commas in the input string
		std::vector<std::string> input = stringTokenize(getData(strInput), ",");
		if (input.size() != 6)
		{
			throw std::runtime_error("Unable to read meanDisp in input line " + std::to_string(IO::getInLineNumber()) + ". Wrong number of parameters.");
		}

		for (int ii = 0; ii < 6; ++ii)
		{
			meanDisp(ii) = string2num<double>(input.at(ii));
		}
		}

		else if (caseInsCompare(getKeyword(strInput), "Disp0"))
		{
			// The components of the initial displacement of the FOWT are separated by commas in the input string
			std::vector<std::string> input = stringTokenize(getData(strInput), ",");
			if (input.size() != 6)
			{
				throw std::runtime_error("Unable to read Disp0 in input line " + std::to_string(IO::getInLineNumber()) + ". Wrong number of parameters.");
			}

			for (int ii = 0; ii < 6; ++ii)
			{
				disp0(ii) = string2num<double>(input.at(ii));
			}
		}

		else if (caseInsCompare(getKeyword(strInput), "Vel0"))
		{
			// The components of the initial displacement of the FOWT are separated by commas in the input string
			std::vector<std::string> input = stringTokenize(getData(strInput), ",");
			if (input.size() != 6)
			{
				throw std::runtime_error("Unable to read Vel0 in input line " + std::to_string(IO::getInLineNumber()) + ". Wrong number of parameters.");
			}

			for (int ii = 0; ii < 6; ++ii)
			{
				vel0(ii) = string2num<double>(input.at(ii));
			}	
		}

		else if (caseInsCompare(getKeyword(strInput), "Morison_circ")) // A list of circular cylinder Morison Elements is supposed to follow the "Morison_circ" keyword
		{
			IO::readLineInputFile(strInput); // Read next line, since current line is just the main keyword

			while (!caseInsCompare(getKeyword(strInput), "END"))
			{
				if (!m_inFl) // Signal if the end of file is reached before the end keyword
				{
					throw std::runtime_error("End of file reached before END keyword in MORISON_CIRC specification.");
					return;
				}

				// The eleven properties of a circular cylinder Morison's Element are separated by white spaces in the input string.
				std::vector<std::string> input = stringTokenize(strInput, " \t");
				if (input.size() != 11)
				{
					throw std::runtime_error("Unable to read the circular cylinder in input line " + std::to_string(IO::getInLineNumber()) + ". Wrong number of parameters.");
				}

				// In order to specify the nodes corresponding to the extremities of the cylinder, it is
				// necessary to have read the nodes from the input file.
				if (envir.isNodeEmpty())
				{
					throw std::runtime_error("Nodes should be specified before Morison Elements. Error in input line " + std::to_string(IO::getInLineNumber()));
				}

				// Aux variables to handle the data read from the input file
				vec::fixed<3> aux_node1_coord = envir.getNode(string2num<unsigned int>(input.at(0)));
				vec::fixed<3> aux_node2_coord = envir.getNode(string2num<unsigned int>(input.at(1)));
				double aux_diam = string2num<double>(input.at(2));
				double aux_CD = string2num<double>(input.at(3));
				double aux_CM = string2num<double>(input.at(4));
				unsigned int aux_numIntPoints = string2num<unsigned int>(input.at(5));
				double aux_axialCD_1 = string2num<double>(input.at(6));
				double aux_axialCa_1 = string2num<double>(input.at(7));
				double aux_axialCD_2 = string2num<double>(input.at(8));
				double aux_axialCa_2 = string2num<double>(input.at(9));
				bool aux_botPressFlag = string2num<bool>(input.at(10));

				floater.addMorisonCirc(aux_node1_coord, aux_node2_coord, aux_diam, aux_CD, aux_CM, aux_numIntPoints, aux_axialCD_1, aux_axialCa_1, aux_axialCD_2, aux_axialCa_2, aux_botPressFlag);

				// Go to next line
				IO::readLineInputFile(strInput);
			}
		}

		else if (caseInsCompare(getKeyword(strInput), "Morison_rect")) // A list of rectangular cylinder Morison Elements is supposed to follow the "Morison_rect" keyword
		{
			IO::readLineInputFile(strInput); // Read next line, since current line is just the main keyword

			while (!caseInsCompare(getKeyword(strInput), "END"))
			{
				if (!m_inFl) // Signal if the end of file is reached before the end keyword
				{
					throw std::runtime_error("End of file reached before END keyword in MORISON_RECT specification.");
					return;
				}

				// The fifteen properties of a circular cylinder Morison's Element are separated by white spaces in the input string.
				std::vector<std::string> input = stringTokenize(strInput, " \t");
				if (input.size() != 15)
				{
					throw std::runtime_error("Unable to read the rectangular cylinder in input line " + std::to_string(IO::getInLineNumber()) + ". Wrong number of parameters.");
				}

				// Check whether nodes were specified
				if (envir.isNodeEmpty())
				{
					throw std::runtime_error("Nodes should be specified before Morison Elements. Error in input line " + std::to_string(IO::getInLineNumber()));
				}

				// Aux variables to handle the data read from the input file
				vec::fixed<3> aux_node1_coord = envir.getNode(string2num<unsigned int>(input.at(0)));
				vec::fixed<3> aux_node2_coord = envir.getNode(string2num<unsigned int>(input.at(1)));
				vec::fixed<3> aux_node3_coord = envir.getNode(string2num<unsigned int>(input.at(2)));
				double aux_diam_X = string2num<double>(input.at(3));
				double aux_CD_X = string2num<double>(input.at(4));
				double aux_CM_X = string2num<double>(input.at(5));
				double aux_diam_Y = string2num<double>(input.at(6));
				double aux_CD_Y = string2num<double>(input.at(7));
				double aux_CM_Y = string2num<double>(input.at(8));
				unsigned int aux_numIntPoints = string2num<unsigned int>(input.at(9));
				double aux_axialCD_1 = string2num<double>(input.at(10));
				double aux_axialCa_1 = string2num<double>(input.at(11));
				double aux_axialCD_2 = string2num<double>(input.at(12));
				double aux_axialCa_2 = string2num<double>(input.at(13));
				bool aux_botPressFlag = string2num<bool>(input.at(14));

				floater.addMorisonRect(aux_node1_coord, aux_node2_coord, aux_node3_coord, aux_diam_X, aux_diam_Y, aux_CD_X, aux_CD_Y,
					aux_CM_X, aux_CM_Y, aux_numIntPoints, aux_axialCD_1, aux_axialCa_1, aux_axialCD_2, aux_axialCa_2, aux_botPressFlag);

				IO::readLineInputFile(strInput);
			}
		}

		// Read data to rna - which is a part of FOWT
		else if (caseInsCompare(getKeyword(strInput), "UseTipLoss"))
		{
			rna.setUseTipLoss(string2num<bool>(getData(strInput)));
		}

		else if (caseInsCompare(getKeyword(strInput), "UseHubLoss"))
		{
			rna.setUseHubLoss(string2num<bool>(getData(strInput)));
		}

		else if (caseInsCompare(getKeyword(strInput), "UseSkewCorr"))
		{
			rna.setUseSkewCorr(string2num<bool>(getData(strInput)));
		}

		else if (caseInsCompare(getKeyword(strInput), "RotSpeed"))
		{
			rna.setRotorSpeed(string2num<double>(getData(strInput)));
		}

		else if (caseInsCompare(getKeyword(strInput), "RotTilt"))
		{
			rna.setRotorTilt(string2num<double>(getData(strInput)));
		}

		else if (caseInsCompare(getKeyword(strInput), "RotYaw"))
		{
			rna.setRotorYaw(string2num<double>(getData(strInput)));
		}

		else if (caseInsCompare(getKeyword(strInput), "BldPitch"))
		{
			rna.setBladePitch(string2num<double>(getData(strInput)));
		}

		else if (caseInsCompare(getKeyword(strInput), "BldPrecone"))
		{
			rna.setBladePrecone(string2num<double>(getData(strInput)));
		}

		else if (caseInsCompare(getKeyword(strInput), "NumBlades"))
		{
			rna.setNumBlades(string2num<unsigned int>(getData(strInput)));
		}

		else if (caseInsCompare(getKeyword(strInput), "HubRadius"))
		{
			rna.setHubRadius(string2num<double>(getData(strInput)));
		}

		else if (caseInsCompare(getKeyword(strInput), "HubHeight"))
		{
			rna.setHubHeight(string2num<double>(getData(strInput)));
		}

		else if (caseInsCompare(getKeyword(strInput), "Overhang"))
		{
			rna.setOverhang(string2num<double>(getData(strInput)));
		}

		else if (caseInsCompare(getKeyword(strInput), "Blades_aero"))
		{
			IO::readLineInputFile(strInput); // Read next line, since current line is just the main keyword

			while (!caseInsCompare(getKeyword(strInput), "END"))
			{
				if (!m_inFl) // Signal if the end of file is reached before the end keyword
				{
					throw std::runtime_error("End of file reached before END keyword in BLADES_AERO specification.");
					return;
				}

				// The seven properties provided by each line are separated by white spaces in the input string (whitespace or tab)
				std::vector<std::string> input = stringTokenize(strInput, " \t");
				if (input.size() != 7)
				{
					throw std::runtime_error("Unable to read the blade aerodynamic properties in input line " + std::to_string(IO::getInLineNumber()) + ". Wrong number of parameters.");
				}

				// Aux variables to handle the data read from the input file
				double aux_span = string2num<double>(input.at(0));
				double aux_crvAC = string2num<double>(input.at(1));
				double aux_swpAC = string2num<double>(input.at(2));
				double aux_crvAng = string2num<double>(input.at(3));
				double aux_twist = string2num<double>(input.at(4));
				double aux_chord = string2num<double>(input.at(5));
				int aux_airfoilID = string2num<int>(input.at(6));

				rna.addBladeAeroNode(aux_span, aux_crvAC, aux_swpAC, aux_crvAng, aux_twist, aux_chord, aux_airfoilID);

				IO::readLineInputFile(strInput);
			}
		}


		else if (caseInsCompare(getKeyword(strInput), "Airfoil_data"))
		{
			IO::readLineInputFile(strInput); // Read next line, since current line is just the main keyword

			rna.addEmptyAirfoil();
			while (!caseInsCompare(getKeyword(strInput), "END"))
			{
				if (!m_inFl) // Signal if the end of file is reached before the end keyword
				{
					throw std::runtime_error("End of file reached before END keyword in AIRFOIL_DATA specification.");
					return;
				}

				// The four properties provided by each line are separated by white spaces in the input string (whitespace or tab)
				std::vector<std::string> input = stringTokenize(strInput, " \t");
				if (input.size() != 4)
				{
					throw std::runtime_error("Unable to read the airfoil properties in input line " + std::to_string(IO::getInLineNumber()) + ". Wrong number of parameters.");
				}

				// Aux variables to handle the data read from the input file
				double aux_angle = string2num<double>(input.at(0));
				double aux_CL = string2num<double>(input.at(1));
				double aux_CD = string2num<double>(input.at(2));
				double aux_CM = string2num<double>(input.at(3));

				rna.addAirfoilData(aux_angle, aux_CL, aux_CD, aux_CM);

				IO::readLineInputFile(strInput);
			}
		}


		else if (caseInsCompare(getKeyword(strInput), "Blades_elasto"))
		{
			IO::readLineInputFile(strInput); // Read next line, since current line is just the main keyword

			while (!caseInsCompare(getKeyword(strInput), "END"))
			{
		 		if (!m_inFl) // Signal if the end of file is reached before the end keyword
		 		{
		 			throw std::runtime_error("End of file reached before END keyword in BLADES_ELASTO specification.");
		 			return;
		 		}

		 		// Implement in the future
		 		IO::readLineInputFile(strInput);
			}
		}

		else if (caseInsCompare(getKeyword(strInput), "Tower_Aero"))
		{
			IO::readLineInputFile(strInput); // Read next line, since current line is just the main keyword

			while (!caseInsCompare(getKeyword(strInput), "END"))
			{
		 		if (!m_inFl) // Signal if the end of file is reached before the end keyword
		 		{
		 			throw std::runtime_error("End of file reached before END keyword in TOWER_AERO specification.");
		 			return;
		 		}

				// Implement in the future

		 		IO::readLineInputFile(strInput);
			}
		}

		else if (caseInsCompare(getKeyword(strInput), "Tower_Elasto"))
		{
			IO::readLineInputFile(strInput); // Read next line, since current line is just the main keyword

			while (!caseInsCompare(getKeyword(strInput), "END"))
			{
		 		if (!m_inFl) // Signal if the end of file is reached before the end keyword
		 		{
		 			throw std::runtime_error("End of file reached before END keyword in TOWER_ELASTO specification.");
		 			return;
		 		}

				// Implement in the future

		 		IO::readLineInputFile(strInput);
			}
		}

		else if (caseInsCompare(getKeyword(strInput), "Output")) // List of parameters that will be output
		{
			IO::readLineInputFile(strInput); // Read next line, since current line is just the keyword

			while (!caseInsCompare(getKeyword(strInput), "END"))
			{
				if (!m_inFl) // Signal if the end of file is reached before the end keyword
				{
					throw std::runtime_error("End of file reached before END keyword in OUTPUT specification.");
					return;
				}

				IO::setResults2Output(strInput, envir);

				IO::readLineInputFile(strInput);
			}
		}

		else if (!caseInsCompare(getKeyword(strInput), "END_OF_INPUT_FILE"))
		{
			print2log("WARNING: Unknown keyword '" + getKeyword(strInput) + "' in line " + std::to_string(IO::getInLineNumber()) + ".");
		}
	}
	
	fowt.setFloater(floater);
	fowt.setRNA(rna);	

	// Set initial state.
	//
	// The slow drift position is always necessary, even when the filter frequency is taken as 0,
	// in which case it is equal to the initial fixed position though the whole simulation. 
	// However, in this condition the filter is skipepd in fowt_updata_sd, hence we change the frequency 
	// before calling it and them restore it back to its original value.	
	double wf = fowt.filterSD_omega();
	fowt.setFilderSD(1, fowt.filterSD_zeta());
	fowt.setDispSD(meanDisp);
	fowt.update(envir, fill::zeros, fill::zeros);
	fowt.setAddedMass_t0(envir.watDensity()); // The added mass matrix hydrostatics is evaluated considering the mean position
	fowt.setStiffnessMatrix(envir.watDensity(), envir.gravity()); // The hydrostatics is evaluated considering the instantaneous position, thus it DOES NOT take into account the mean displacement
	fowt.update(envir, join_cols(vec::fixed<6>(fill::zeros), disp0), join_cols(vec::fixed<6>(fill::zeros), vel0));
	fowt.setFilderSD(wf, fowt.filterSD_zeta());
}


/*
	Set output files path and open them. There are three output files:
	- Formatted output file, with the time series of selected outputs;
	- Summary file, with a summary of the simulation;
	- Log file, with errors and warnings.
*/
void IO::setFiles(const std::string &inFlPath)
{
	// Check whether we are not trying to reset the input file
	if (!m_inFilePath.empty())
	{
		throw std::runtime_error("You can not reset the input file.");
	}
	m_inFilePath = inFlPath;

	// Open input file
	m_inFl.open(m_inFilePath);
	if (!m_inFl)
	{
		throw std::runtime_error("Unable to open file " + m_inFilePath + " for reading.");
	}

	// Get path of the folder where the input file is located
	std::string folderPath = getFileFolder(m_inFilePath);

	// Name of the input file, without extension
	std::string inFlName = getFileName(m_inFilePath);

	// The output files are created in the same folder as the input file.
	//
	// If a file with the same name of a given output file already exists in folderPath, then a _1 is appended to its name. If the latter exists as well, append_2 instead. Keep this process until output_n does not exist and is then created.
	// For consistency, the same _n is used for all the output files, even if a previous number does not exist for one of the files. This number _n is the largest available number starting from 1.
	// So, if the files file_out, file_sum_2, file_log_3  already exist, the files file_log_1, file_sum_1, and file_out_1 are created.
	struct stat info;

	// Set 'log', 'summary' and 'formatted output' files
	m_logFilePath = folderPath + filesep + inFlName + "_log.txt";
	m_sumFilePath = folderPath + filesep + inFlName + "_sum.txt";
	m_outFilePath = folderPath + filesep + inFlName + "_out.txt";

	int index = 1; // This is the part where we verify if the file exists and append numbers as needed
	while (stat(m_logFilePath.c_str(), &info) == 0 || stat(m_sumFilePath.c_str(), &info) == 0 || stat(m_outFilePath.c_str(), &info) == 0)
	{
		m_logFilePath = folderPath + filesep + inFlName + "_log_" + std::to_string(index) + ".txt";
		m_sumFilePath = folderPath + filesep + inFlName + "_sum_" + std::to_string(index) + ".txt";
		m_outFilePath = folderPath + filesep + inFlName + "_out_" + std::to_string(index) + ".txt";
		++index;
	}

	// Open the output files
	m_logFl.open(m_logFilePath); // Check if we can open it
	if (!m_logFl)
	{
		throw std::runtime_error("Unable to open file " + m_logFilePath + " for writting.");
	}

	m_sumFl.open(m_sumFilePath);
	if (!m_sumFl)
	{
		throw std::runtime_error("Unable to open file " + m_sumFilePath + " for writting.");
	}

	m_outFl.open(m_outFilePath);
	if (!m_outFl)
	{
		throw std::runtime_error("Unable to open file " + m_outFilePath + " for writting.");
	}

	// Print output files path to the console
	std::cout << "Printing log file to: '"              << m_logFilePath << "'\n";
	std::cout << "Printing summary file to: '"          << m_sumFilePath << "'\n";
	std::cout << "Printing formatted output file to: '" << m_outFilePath << "'\n";
}


// Read line from input file to string "strInput".
// The function deals with empty lines and comments using functions "hasContent" and
// "thereIsCommentInString". Besides, it updates the line number counter inLineNumber
void IO::readLineInputFile(std::string &strInput)
{
	strInput = "";

	// Read next file line to string strInput and update line number counter.
    // Repeat this process until the line has some content or end of file is achieved.
	while (!hasContent(strInput) && m_inFl)
	{
		std::getline(m_inFl, strInput);
		++m_inLineNumber;

		// Remove comments from line
		if (thereIsCommentInString(strInput))
		{
			removeComments(strInput);
		}
	}

	if (!m_inFl)
	{
		strInput = "END_OF_INPUT_FILE";
	}
}


unsigned int IO::getInLineNumber()
{
	return m_inLineNumber;
}


// Set the flags specifying which variables must be output
// and other necessary information (like the IDs of the points where the wave elevation will be output)
void IO::setResults2Output(std::string strInput, ENVIR &envir)
{
	const std::string keyword = getKeyword(strInput);
	bool isOutput = false;

	if (caseInsCompare(keyword, "fowt_disp"))
	{
		m_whichResult2Output.at(IO::OUTFLAG_FOWT_DISP) = true;
		m_whichResult2Output.at(IO::OUTFLAG_FOWT_DISP_1ST) = true;
		isOutput = true;
	}

	if (caseInsCompare(keyword, "fowt_vel"))
	{
		m_whichResult2Output.at(IO::OUTFLAG_FOWT_VEL) = true;
		m_whichResult2Output.at(IO::OUTFLAG_FOWT_VEL_1ST) = true;
		isOutput = true;
	}

	if (caseInsCompare(keyword, "fowt_acc"))
	{
		m_whichResult2Output.at(IO::OUTFLAG_FOWT_ACC) = true;
		m_whichResult2Output.at(IO::OUTFLAG_FOWT_ACC_1ST) = true;
		isOutput = true;
	}

	if (caseInsCompare(keyword, "fowt_disp_sd"))
	{
		m_whichResult2Output.at(IO::OUTFLAG_FOWT_DISP_SD) = true;
		isOutput = true;
	}

	if (caseInsCompare(keyword, "hd_drag_force"))
	{
		m_whichResult2Output.at(IO::OUTFLAG_HD_FORCE_DRAG) = true;
		isOutput = true;
	}

	if (caseInsCompare(keyword, "hd_force_1stP"))
	{
		m_whichResult2Output.at(IO::OUTFLAG_HD_FORCE_1STP) = true;
		isOutput = true;
	}

	if (caseInsCompare(keyword, "hd_force_eta"))
	{
		m_whichResult2Output.at(IO::OUTFLAG_HD_FORCE_ETA) = true;
		isOutput = true;
	}

	if (caseInsCompare(keyword, "hd_force_conv"))
	{
		m_whichResult2Output.at(IO::OUTFLAG_HD_FORCE_CONV) = true;
		isOutput = true;
	}

	if (caseInsCompare(keyword, "hd_force_axDv"))
	{
		m_whichResult2Output.at(IO::OUTFLAG_HD_FORCE_AXDV) = true;
		isOutput = true;
	}

	if (caseInsCompare(keyword, "hd_force_2ndP"))
	{
		m_whichResult2Output.at(IO::OUTFLAG_HD_FORCE_2NDP) = true;
		isOutput = true;
	}

	if (caseInsCompare(keyword, "hd_force_acgr"))
	{
		m_whichResult2Output.at(IO::OUTFLAG_HD_FORCE_ACGR) = true;
		isOutput = true;
	}

	if (caseInsCompare(keyword, "hd_force_rotN"))
	{
		m_whichResult2Output.at(IO::OUTFLAG_HD_FORCE_ROTN) = true;
		isOutput = true;
	}

	if (caseInsCompare(keyword, "hd_force_RSLB"))
	{
		m_whichResult2Output.at(IO::OUTFLAG_HD_FORCE_RSLB) = true;
		isOutput = true;
	}

	if (caseInsCompare(keyword, "hd_force_Rem"))
	{
		m_whichResult2Output.at(IO::OUTFLAG_HD_FORCE_REM) = true;
		isOutput = true;
	}

	if (caseInsCompare(keyword, "hd_add_mass_force"))
	{
		m_whichResult2Output.at(IO::OUTFLAG_HD_ADD_MASS_FORCE) = true;
		isOutput = true;
	}

	if (caseInsCompare(keyword, "hd_force"))
	{
		m_whichResult2Output.at(IO::OUTFLAG_HD_FORCE) = true;
		isOutput = true;
	}

	if (caseInsCompare(keyword, "hs_force"))
	{
		m_whichResult2Output.at(IO::OUTFLAG_HS_FORCE) = true;
		isOutput = true;
	}

	if (caseInsCompare(keyword, "moor_force"))
	{
		m_whichResult2Output.at(IO::OUTFLAG_MOOR_FORCE) = true;
		isOutput = true;
	}

	if (caseInsCompare(keyword, "ad_hub_force"))
	{
		m_whichResult2Output.at(IO::OUTFLAG_AD_HUB_FORCE) = true;
		isOutput = true;
	}

	if (caseInsCompare(keyword, "added_mass"))
	{
		m_whichResult2Output.at(IO::OUTFLAG_ADDED_MASS_DIAG) = true;
		isOutput = true;
	}

	if (caseInsCompare(keyword, "total_force"))
	{
		m_whichResult2Output.at(IO::OUTFLAG_TOTAL_FORCE) = true;
		isOutput = true;
	}

	// Add wave probes for wave elev, velocity, acceleration or pressure
	bool addWaveProbe{ false };
	if (caseInsCompare(keyword, "wave_elev"))
	{
		m_whichResult2Output.at(IO::OUTFLAG_WAVE_ELEV) = true;
		addWaveProbe = true;
		isOutput = true;
	}

	if (caseInsCompare(keyword, "wave_vel"))
	{
		m_whichResult2Output.at(IO::OUTFLAG_WAVE_VEL) = true;
		addWaveProbe = true;
		isOutput = true;
	}

	if (caseInsCompare(keyword, "wave_acc"))
	{
		m_whichResult2Output.at(IO::OUTFLAG_WAVE_ACC) = true;
		addWaveProbe = true;
		isOutput = true;
	}

	if (caseInsCompare(keyword, "wave_pres"))
	{
		m_whichResult2Output.at(IO::OUTFLAG_WAVE_PRES) = true;
		addWaveProbe = true;
		isOutput = true;
	}

	if (caseInsCompare(keyword, "wave_acc_2nd"))
	{
		m_whichResult2Output.at(IO::OUTFLAG_WAVE_ACC_2ND) = true;
		addWaveProbe = true;
		isOutput = true;
	}

	if (caseInsCompare(keyword, "wave_pres_2nd"))
	{
		m_whichResult2Output.at(IO::OUTFLAG_WAVE_PRES_2ND) = true;
		addWaveProbe = true;
		isOutput = true;
	}

	if (addWaveProbe)
	{
		if (!getData(strInput).empty())
		{
			// The wave locations are specified by node IDs separated by tabs or white-spaces
			std::vector<std::string> input = stringTokenize(getData(strInput), " \t");

			if (input.empty())
			{
				throw std::runtime_error("You should specify at least one node ID for defining a wave location. Error in input line " + std::to_string(IO::getInLineNumber()) + ".");
			}

			for (int ii = 0; ii < input.size(); ++ii)
			{
				envir.addWaveProbe(string2num<unsigned int>(input.at(ii)));
			}
		}
	}

	if (!isOutput)
	{
		print2log("WARNING: Unknown output option '" + keyword + "' in line " + std::to_string(IO::getInLineNumber()) + ".");
	}
}



void IO::checkInputs(const FOWT &fowt, const ENVIR &envir)
{

}






/*****************************************************
Class functions related to output
*****************************************************/

// Returns the header that is later printed to the console and to the summary file (perhaps to other files in the future).
// I know that making it a function is not optimal in terms of performance, but since it is irrelevant
// compared to the time spent with the rest of the code, I chose this solutions because it is easier to maintain.
std::string IO::METiS_Header()
{
	std::string header("");
	header += ">--------<>--------<>--------<+>--------<>--------<>--------<\n";
	header += ">                                                           <\n";
	header += ">                        METiS - USP                        <\n";
	header += ">          Morison Equation Time Domain Simulation          <\n";
	header += ">              University of Sao Paulo - Brazil             <\n";
	header += ">                                                           <\n";
	header += ">                                                   v. " + g_METIS_VERSION + "<\n";
	header += ">--------<>--------<>--------<+>--------<>--------<>--------<\n";

	return header;
}

void IO::print2log(const std::string &str)
{
	if (m_logFl) // If we are able to write to the log file, do it
	{
		m_logFl << str << std::endl;
	}

	// Write to the console as well
	std::cout << "\n" << str << std::endl;
}


// Set m_shouldWriteOutLine to true or false. This is the variable
// that tells when the program should write to the formatted output file.
// Same thing applies to m_shouldWriteOutLineHeader
void IO::print2outLine_turnOn()
{
	m_shouldWriteOutLine = true;
}

void IO::print2outLine_turnOff()
{
	m_shouldWriteOutLine = false;
}

void IO::print2outLineHeader_turnOn()
{
	m_shouldWriteOutLineHeader = true;
}

void IO::print2outLineHeader_turnOff()
{
	m_shouldWriteOutLineHeader = false;
}

// Different functions to print to the formatted output file string streams depending on the output flag. This one is for
// outputs represented by vectors with six components (forces with moments, FOWT displacements, etc)
void IO::print2outLine(const OutFlag &flag, const arma::vec::fixed<6> &vector_6)
{
	// Check whether the specified flag is indeed one that requires a vector with six components
	if ((flag != IO::OUTFLAG_FOWT_DISP) && (flag != IO::OUTFLAG_FOWT_VEL) && (flag != IO::OUTFLAG_FOWT_ACC) && 
		(flag != IO::OUTFLAG_FOWT_DISP_1ST) && (flag != IO::OUTFLAG_FOWT_VEL_1ST) && (flag != IO::OUTFLAG_FOWT_ACC_1ST) && (flag != IO::OUTFLAG_FOWT_DISP_SD) &&
		(flag != IO::OUTFLAG_TOTAL_FORCE) && (flag != IO::OUTFLAG_HD_FORCE) && (flag != IO::OUTFLAG_HS_FORCE) && (flag != IO::OUTFLAG_MOOR_FORCE) &&
		(flag != IO::OUTFLAG_HD_FORCE_DRAG) && (flag != IO::OUTFLAG_HD_FORCE_1STP) && (flag != IO::OUTFLAG_HD_FORCE_ETA) &&
		(flag != IO::OUTFLAG_HD_FORCE_CONV) && (flag != IO::OUTFLAG_HD_FORCE_AXDV) && (flag != IO::OUTFLAG_HD_FORCE_ACGR) &&
		(flag != IO::OUTFLAG_HD_FORCE_ROTN) && (flag != IO::OUTFLAG_HD_FORCE_2NDP) && (flag != IO::OUTFLAG_HD_FORCE_RSLB) && (flag != IO::OUTFLAG_HD_FORCE_REM) &&
		(flag != IO::OUTFLAG_HD_ADD_MASS_FORCE) && (flag != IO::OUTFLAG_AD_HUB_FORCE) && (flag != IO::OUTFLAG_ADDED_MASS_DIAG) && (flag != OUTFLAG_DEBUG_VEC_6)
	   )
	{
		throw std::runtime_error("Unknown output flag in function IO::print2outLine(const OutFlag &flag, const arma::vec::fixed<6> &force).");
	}

	// If the print header flag is true and if this is one of the requested output variables,
	// then print the header based on the output flag
	if (m_shouldWriteOutLineHeader && (m_whichResult2Output.at(flag) || flag == OUTFLAG_DEBUG_VEC_6))
	{
		if (flag == OUTFLAG_DEBUG_VEC_6)
		{
			for (int ii = 1; ii <= 6; ++ii)
			{
				print2outLineHeader("DEBUG_VEC_" + std::to_string(ii));
			}
		}

		if (flag == OUTFLAG_HD_FORCE_DRAG)
		{
			for (int ii = 1; ii <= 6; ++ii)
			{
				print2outLineHeader("hd_drag_force_" + std::to_string(ii));
			}
		}

		if (flag == OUTFLAG_HD_FORCE_1STP)
		{
			for (int ii = 1; ii <= 6; ++ii)
			{
				print2outLineHeader("hd_force_1stP_" + std::to_string(ii));
			}
		}

		if (flag == OUTFLAG_HD_FORCE_ETA)
		{
			for (int ii = 1; ii <= 6; ++ii)
			{
				print2outLineHeader("hd_force_eta_" + std::to_string(ii));
			}
		}

		if (flag == OUTFLAG_HD_FORCE_CONV)
		{
			for (int ii = 1; ii <= 6; ++ii)
			{
				print2outLineHeader("hd_force_conv_" + std::to_string(ii));
			}
		}

		if (flag == OUTFLAG_HD_FORCE_AXDV)
		{
			for (int ii = 1; ii <= 6; ++ii)
			{
				print2outLineHeader("hd_force_axdv_" + std::to_string(ii));
			}
		}

		if (flag == OUTFLAG_HD_FORCE_ACGR)
		{
			for (int ii = 1; ii <= 6; ++ii)
			{
				print2outLineHeader("hd_force_acgr_" + std::to_string(ii));
			}
		}

		if (flag == OUTFLAG_HD_FORCE_ROTN)
		{
			for (int ii = 1; ii <= 6; ++ii)
			{
				print2outLineHeader("hd_force_rotn_" + std::to_string(ii));
			}
		}

		if (flag == OUTFLAG_HD_FORCE_2NDP)
		{
			for (int ii = 1; ii <= 6; ++ii)
			{
				print2outLineHeader("hd_force_2ndP_" + std::to_string(ii));
			}
		}

		if (flag == OUTFLAG_HD_FORCE_RSLB)
		{
			for (int ii = 1; ii <= 6; ++ii)
			{
				print2outLineHeader("hd_force_rslb_" + std::to_string(ii));
			}
		}

		if (flag == OUTFLAG_HD_FORCE_REM)
		{
			for (int ii = 1; ii <= 6; ++ii)
			{
				print2outLineHeader("hd_force_rem_" + std::to_string(ii));
			}
		}

		if (flag == OUTFLAG_HD_ADD_MASS_FORCE)
		{
			for (int ii = 1; ii <= 6; ++ii)
			{
				print2outLineHeader("hd_add_mass_force_" + std::to_string(ii));
			}
		}

		if (flag == OUTFLAG_HD_FORCE)
		{
			for (int ii = 1; ii <= 6; ++ii)
			{
				print2outLineHeader("hd_force_" + std::to_string(ii));
			}
		}

		if (flag == OUTFLAG_HS_FORCE)
		{
			for (int ii = 1; ii <= 6; ++ii)
			{
				print2outLineHeader("hs_force_" + std::to_string(ii));
			}
		}

		if (flag == OUTFLAG_AD_HUB_FORCE)
		{
			for (int ii = 1; ii <= 6; ++ii)
			{
				print2outLineHeader("ad_hub_force_" + std::to_string(ii));
			}
		}

		if (flag == OUTFLAG_MOOR_FORCE)
		{
			for (int ii = 1; ii <= 6; ++ii)
			{
				print2outLineHeader("moor_force_" + std::to_string(ii));
			}
		}

		if (flag == OUTFLAG_ADDED_MASS_DIAG)
		{
			for (int ii = 1; ii <= 6; ++ii)
			{
				print2outLineHeader("A_" + std::to_string(ii) + std::to_string(ii));
			}
		}

		if (flag == OUTFLAG_TOTAL_FORCE)
		{
			for (int ii = 1; ii <= 6; ++ii)
			{
				print2outLineHeader("total_force_" + std::to_string(ii));
			}
		}

		if (flag == OUTFLAG_FOWT_DISP)
		{
			print2outLineHeader("surge");
			print2outLineHeader("sway");
			print2outLineHeader("heave");
			print2outLineHeader("roll");
			print2outLineHeader("pitch");
			print2outLineHeader("yaw");
		}

		if (flag == OUTFLAG_FOWT_DISP_1ST)
		{
			print2outLineHeader("surge_1st");
			print2outLineHeader("sway_1st");
			print2outLineHeader("heave_1st");
			print2outLineHeader("roll_1st");
			print2outLineHeader("pitch_1st");
			print2outLineHeader("yaw_1st");
		}

		if (flag == OUTFLAG_FOWT_DISP_SD)
		{
			print2outLineHeader("surge_sd");
			print2outLineHeader("sway_sd");
			print2outLineHeader("heave_sd");
			print2outLineHeader("roll_sd");
			print2outLineHeader("pitch_sd");
			print2outLineHeader("yaw_sd");
		}

		if (flag == OUTFLAG_FOWT_VEL)
		{
			print2outLineHeader("surge_vel");
			print2outLineHeader("sway_vel");
			print2outLineHeader("heave_vel");
			print2outLineHeader("roll_vel");
			print2outLineHeader("pitch_vel");
			print2outLineHeader("yaw_vel");
		}

		if (flag == OUTFLAG_FOWT_VEL_1ST)
		{
			print2outLineHeader("surge_vel_1st");
			print2outLineHeader("sway_vel_1st");
			print2outLineHeader("heave_vel_1st");
			print2outLineHeader("roll_vel_1st");
			print2outLineHeader("pitch_vel_1st");
			print2outLineHeader("yaw_vel_1st");
		}

		if (flag == OUTFLAG_FOWT_ACC)
		{
			print2outLineHeader("surge_acc");
			print2outLineHeader("sway_acc");
			print2outLineHeader("heave_acc");
			print2outLineHeader("roll_acc");
			print2outLineHeader("pitch_acc");
			print2outLineHeader("yaw_acc");
		}

		if (flag == OUTFLAG_FOWT_ACC_1ST)
		{
			print2outLineHeader("surge_acc_1st");
			print2outLineHeader("sway_acc_1st");
			print2outLineHeader("heave_acc_1st");
			print2outLineHeader("roll_acc_1st");
			print2outLineHeader("pitch_acc_1st");
			print2outLineHeader("yaw_acc_1st");
		}
	}

	// If the printing flag is true and if this is one of the requested output variables,
	// then print it to the output line	
	if (m_shouldWriteOutLine && (m_whichResult2Output.at(flag) || flag == OUTFLAG_DEBUG_VEC_6))
	{
		for (int ii = 0; ii < 6; ++ii)
		{
			print2outLine(vector_6.at(ii));
		}
	}
}


// This one is useful for wave velocity and acceleration, which need the ID of the node where they are calculated
// and the value itself, which three component vector. Other future outputs may profit from this function as well.
void IO::print2outLine(const OutFlag &flag, const int ID, const arma::vec::fixed<3> &vector_3)
{
	if ((flag != OUTFLAG_WAVE_VEL) && (flag != OUTFLAG_WAVE_ACC) && (flag != OUTFLAG_WAVE_ACC_2ND) && (flag != OUTFLAG_DEBUG_VEC_3))
	{
		throw std::runtime_error("Unknown output flag in function IO::print2outLine(const OutFlag &flag, const int ID, const arma::vec::fixed<3> &vector_3).");
	}

	// If the print header flag is true and if this is one of the requested output variables,
	// then print the header based on the output flag
	if (m_shouldWriteOutLineHeader && (m_whichResult2Output.at(flag) || flag == OUTFLAG_DEBUG_VEC_3))
	{
		if (flag == OUTFLAG_DEBUG_VEC_3)
		{
			print2outLineHeader("DEBUG_VEC_" + std::to_string(ID) + "_x");
			print2outLineHeader("DEBUG_VEC_" + std::to_string(ID) + "_y");
			print2outLineHeader("DEBUG_VEC_" + std::to_string(ID) + "_z");
		}

		if (flag == OUTFLAG_WAVE_VEL)
		{
			print2outLineHeader("wave_vel_" + std::to_string(ID) + "_x");
			print2outLineHeader("wave_vel_" + std::to_string(ID) + "_y");
			print2outLineHeader("wave_vel_" + std::to_string(ID) + "_z");
		}

		if (flag == OUTFLAG_WAVE_ACC)
		{
			print2outLineHeader("wave_acc_" + std::to_string(ID) + "_x");
			print2outLineHeader("wave_acc_" + std::to_string(ID) + "_y");
			print2outLineHeader("wave_acc_" + std::to_string(ID) + "_z");
		}

		if (flag == OUTFLAG_WAVE_ACC_2ND)
		{
			print2outLineHeader("wave_acc_2nd_" + std::to_string(ID) + "_x");
			print2outLineHeader("wave_acc_2nd_" + std::to_string(ID) + "_y");
			print2outLineHeader("wave_acc_2nd_" + std::to_string(ID) + "_z");
		}
	}

	// If the printing flag is true and if this is one of the requested output variables,
	// then print it to the output line
	if (m_shouldWriteOutLine && (m_whichResult2Output.at(flag) || flag == OUTFLAG_DEBUG_VEC_3))
	{
		for (int ii = 0; ii < 3; ++ii)
		{
			print2outLine(vector_3.at(ii));
		}
	}
}


// This one is useful for wave elevation and wave pressure, which need the ID of the node where it is calculated
// and the value itself, which is a double. Other future outputs may profit from this function as well.
void IO::print2outLine(const OutFlag &flag, const int ID, const double num)
{
	if ( (flag != OUTFLAG_WAVE_ELEV) && (flag != OUTFLAG_WAVE_PRES) && (flag != OUTFLAG_WAVE_PRES_2ND) && (flag != OUTFLAG_DEBUG_NUM))
	{
		throw std::runtime_error("Unknown output flag in function IO::print2outLine(const OutFlag &flag, const int ID, const double num).");
	}

	// If the print header flag is true and if this is one of the requested output variables,
	// then print the header based on the output flag
	if (m_shouldWriteOutLine && (m_whichResult2Output.at(flag) || flag == OUTFLAG_DEBUG_NUM))
	{
		if (flag == OUTFLAG_DEBUG_NUM)
		{
			print2outLineHeader("DEBUG_NUM_" + std::to_string(ID));
		}

		if (flag == OUTFLAG_WAVE_ELEV)
		{
			print2outLineHeader("wave_elev_" + std::to_string(ID));
		}

		if (flag == OUTFLAG_WAVE_PRES)
		{
			print2outLineHeader("wave_pres_" + std::to_string(ID));
		}

		if (flag == OUTFLAG_WAVE_PRES_2ND)
		{
			print2outLineHeader("wave_pres_2nd_" + std::to_string(ID));
		}
	}

	// If the printing flag is true and if this is one of the requested output variables,
	// then print it to the output line
	if (m_shouldWriteOutLine && (m_whichResult2Output.at(flag) || flag == OUTFLAG_DEBUG_NUM))
	{
		print2outLine(num);
	}
}


// Different basic functions to print to the formatted output file string streams depending on the variable type.
// The format parameters (column width and precision) are member variables of the IO class and
// can be adjusted changing the values of m_outColumnWidth and m_outNumPrecision
void IO::print2outLine(const std::string &str)
{
	m_outLine << std::setw(IO::m_outColumnWidth) << str;
}

void IO::print2outLine(const double num)
{
	m_outLine << std::setw(IO::m_outColumnWidth) << std::scientific << std::setprecision(IO::m_outNumPrecision) << num;
}

void IO::print2outLine_decimal(const double num)
{
	m_outLine << std::setw(IO::m_outColumnWidth) << std::fixed << std::setprecision(IO::m_outNumPrecision) << num;
}

void IO::print2outLine(const int num)
{
	m_outLine << std::setw(IO::m_outColumnWidth) << num;
}

void IO::print2outLineHeader(const std::string &str)
{
	m_outLineHeader << std::setw(IO::m_outColumnWidth) << str;
}


// Since the data is output in columns, it is necessary to add a new line at each new print step
void IO::newLineOutFile()
{
	if (!m_outFl)
	{
		throw std::runtime_error("Unable to write to formatted output file.");
	}

	m_outFl << '\n';
}

void IO::printOutLineHeader2outFile()
{
	if (!m_outFl)
	{
		throw std::runtime_error("Unable to write to formatted output file.");
	}

	if (m_shouldWriteOutLineHeader)
	{
		m_outFl << std::setw(IO::m_outColumnWidth) << "Time"; // First column must always be the current time
		m_outFl << m_outLineHeader.str() << '\n'; // Then comes the results that are stored in the header stringstream
		m_outLineHeader.str(""); // Need to clear the stream for the next time step
	}
}

void IO::printOutLine2outFile()
{
	if (!m_outFl)
	{
		throw std::runtime_error("Unable to write to formatted output file.");
	}

	if (m_shouldWriteOutLine)
	{
		m_outFl << m_outLine.str() << '\n';
		m_outLine.str(""); // Need to clear the stream for the next time step
	}
}


// Print the members of fowt and envir. Useful for debugging.
void IO::printSumFile(const FOWT &fowt, const ENVIR &envir)
{
	if (!m_sumFl)
	{
		throw std::runtime_error("Unable to open file " + m_sumFilePath + " for writting.");
	}

	m_sumFl << IO::METiS_Header();
	m_sumFl << "\n\n";

	m_sumFl << "ENVIR:\n";
	m_sumFl << "Base time Step:\t" << envir.timeStep() << '\n';
	m_sumFl << "Print Step:\t" << envir.printStep() << '\n';
	m_sumFl << "Total Time:\t" << envir.timeTotal() << '\n';
	m_sumFl << "Time Ramp:\t" << envir.printTimeRamp() << '\n';
	m_sumFl << "Gravity:\t" << envir.gravity() << '\n';
	m_sumFl << "Water Density:\t" << envir.watDensity() << '\n';
	m_sumFl << "Water Depth:\t" << envir.watDepth() << '\n';
	m_sumFl << "Wave Stretching:\t" << envir.waveStret() << '\n';
	m_sumFl << "Air density:\t" << envir.airDensity() << '\n';
	m_sumFl << "Wind Ref velocity:\t" << envir.windRefVel() << '\n';
	m_sumFl << "Wind direction:\t" << envir.windDir() << '\n';
	m_sumFl << "Wind Height:\t" << envir.windRefHeight() << '\n';
	m_sumFl << "Wind exp:\t" << envir.windExp() << '\n';
	m_sumFl << "Nodes: \n" << envir.printNodes() << '\n';
	m_sumFl << "Wave Locations: " << envir.printWaveProbe() << '\n';
	m_sumFl << "\n" << envir.printWave() << '\n';

	m_sumFl << "\n\n";
	m_sumFl << "FOWT:\n";
	m_sumFl << "Hydro Mode:\t" << fowt.printHydroMode() << "\n";
	m_sumFl << "Aero Mode:\t" << fowt.printAeroMode() << "\n";
	m_sumFl << "Moor Mode:\t" << fowt.printMoorMode() << "\n";
	m_sumFl << "Filter SD:\t" << fowt.filterSD_omega() << "\t" << fowt.filterSD_zeta() << "\n";
	m_sumFl << "DOFs:\t" << fowt.printDoF() << '\n';
	m_sumFl << "MeanDisp:\t" << fowt.disp_sd() << '\n';

	if (fowt.moorMode() == 0)
	{
		m_sumFl << "No stiffness considered because Moor Mode = 0\n";
	}
	else
	{
		m_sumFl << "Linear Stiffness:\t" << fowt.printLinStiff() << '\n';
		m_sumFl << "Const force:\t" << fowt.constForce() << '\n';
	}

	if (fowt.hydroMode() == 0)
	{
		m_sumFl << "No floater considered because Hydro Mode = 0\n";
	}
	else
	{
		m_sumFl << "Floater:\n" << fowt.printFloater();
	}

	if (fowt.aeroMode() == 0)
	{
		m_sumFl << "No RNA considered because Aero Mode = 0\n";
	}
	else
	{
		m_sumFl << "RNA:\n" << fowt.printRNA();
	}

	m_sumFl << "\n\n";
	m_sumFl << "Output Variables:\n" << IO::printOutVar();
}


// Some printing functions
std::string IO::printOutVar()
{
	std::string output = "";
	for (int ii = 0; ii < IO::OUTFLAG_SIZE; ++ii)
	{
		bool printFlag(true);
		switch (ii)
		{
		case IO::OUTFLAG_FOWT_DISP:
			output += "FOWT rigid displacement: ";
			break;

		case IO::OUTFLAG_FOWT_VEL:
			output += "FOWT rigid velocity: ";
			break;

		case IO::OUTFLAG_FOWT_ACC:
			output += "FOWT rigid acceleration: ";
			break;

		case IO::OUTFLAG_FOWT_DISP_SD:
			output += "FOWT rigid displacement - Slow: ";
			break;

		case IO::OUTFLAG_WAVE_ELEV:
			output += "Wave Elevation: ";
			break;

		case IO::OUTFLAG_WAVE_VEL:
			output += "Wave Velocity: ";
			break;

		case IO::OUTFLAG_WAVE_ACC:
			output += "Wave Acceleration: ";
			break;

		case IO::OUTFLAG_WAVE_ACC_2ND:
			output += "Wave Acceleration - 2nd order: ";
			break;

		case IO::OUTFLAG_WAVE_PRES:
			output += "Wave Pressure: ";
			break;

		case IO::OUTFLAG_WAVE_PRES_2ND:
			output += "Wave Pressure - 2nd order: ";
			break;		

		case IO::OUTFLAG_HD_FORCE_DRAG:
			output += "Hydrodynamic force - Drag: ";
			break;

		case IO::OUTFLAG_HD_FORCE_1STP:
			output += "Hydrodynamic force - 1stp - 1st order pot: ";
			break;

		case IO::OUTFLAG_HD_FORCE_ETA:
			output += "Hydrodynamic force - Eta - Wave elevation: ";
			break;

		case IO::OUTFLAG_HD_FORCE_CONV:
			output += "Hydrodynamic force - CONV - Conv acc: ";
			break;

		case IO::OUTFLAG_HD_FORCE_AXDV:
			output += "Hydrodynamic force - AXDIV - Ax-diverg acc: ";
			break;

		case IO::OUTFLAG_HD_FORCE_ACGR:
			output += "Hydrodynamic force - ACGR - Acceleration gradient: ";
			break;

		case IO::OUTFLAG_HD_FORCE_ROTN:
			output += "Hydrodynamic force - ROTN - Rotation of normal vector: ";
			break;

		case IO::OUTFLAG_HD_FORCE_2NDP:
			output += "Hydrodynamic force - 2NDP - 2nd order pot: ";
			break;

		case IO::OUTFLAG_HD_FORCE_RSLB:
			output += "Hydrodynamic force - RLSB - Rotation from slender-body app: ";
			break;

		case IO::OUTFLAG_HD_FORCE_REM:
			output += "Hydrodynamic force - Remaining: ";
			break;		

		case IO::OUTFLAG_HD_ADD_MASS_FORCE:
			output += "Force due to added mass: ";
			break;

		case IO::OUTFLAG_HD_FORCE:
			output += "Hydrodynamic force: ";
			break;

		case IO::OUTFLAG_HS_FORCE:
			output += "Hydrostatic force: ";
			break;

		case IO::OUTFLAG_MOOR_FORCE:
			output += "Mooring force: ";
			break;

		case IO::OUTFLAG_AD_HUB_FORCE:
			output += "Aerodynamic forces (Hub): ";
			break;

		case IO::OUTFLAG_ADDED_MASS_DIAG:
			output += "Main diag. added mass: ";
			break;

		case IO::OUTFLAG_TOTAL_FORCE:
			output += "Total force: ";
			break;

		// Options that are not printted to the sum file:
		// - Debug options, as they are for development usage
		// - 1st order quantities, because they are activated with their respective total flags
		case IO::OUTFLAG_DEBUG_NUM:
			printFlag = false;
			break;

		case OUTFLAG_DEBUG_VEC_3:
			printFlag = false;
			break;

		case OUTFLAG_DEBUG_VEC_6:
			printFlag = false;
			break;

		case IO::OUTFLAG_FOWT_DISP_1ST:
			printFlag = false;
			break;

		case IO::OUTFLAG_FOWT_VEL_1ST:
			printFlag = false;
			break;

		case IO::OUTFLAG_FOWT_ACC_1ST:
			printFlag = false;
			break;

		default:
			output += "Unknown specifier in output flags.";
			break;
		}

		if (printFlag)
			output += std::to_string(m_whichResult2Output.at(ii)) + "\n";
	}
	return output;
}

bool IO::isOutputActive(const OutFlag &flag)
{
	return m_whichResult2Output.at(flag);
}




/*****************************************************
	Additional functions related to input/output
*****************************************************/
std::string getKeyword(const std::string &str)
{
	// Get the first part of the string, the one before the first '\t' or white-space
	return str.substr(0, str.find_first_of(" \t"));
}

// Get the part of the string after the keyword, excluding the '\t' or white-space
// (i.e. the part of the string after the first '\t' or white-space).
std::string getData(const std::string &str)
{
	// Check if input string is empty
	if (str.empty())
	{
		throw std::runtime_error("Empty string passed to getData(). Error in input line " + std::to_string(IO::getInLineNumber()) + ".");
	}

	std::vector<std::string> str_tokenized = stringTokenize(str, " \t");

	// If str_tokenized has only one element, i.e. only the keyword, then return an empty string
	if (str_tokenized.size() == 1)
	{
		return "";
	}

	std::string str_out = "";
	for (int ii = 1; ii < str_tokenized.size(); ++ii) // Start at one to skip keyword
	{
		str_out += str_tokenized.at(ii) + " ";
	}
	return str_out;
}
