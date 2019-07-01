#include "IO.h"

#include <fstream> // Include file input/output classes
#include <iostream>
#include <sstream>
#include <stdexcept> // For std::exception
#include <iomanip> // For input/output manipulators
#include <sys/types.h> // Using this and stat.h to check if directories and files exist
#include <sys/stat.h>
#include <cstdlib> // for exit()


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
const unsigned int IO::m_outColumnWidth = 18;
const unsigned int IO::m_outNumPrecision = 4;
std::array<bool, IO::OUTFLAG_SIZE> IO::m_whichResult2Output;

std::stringstream IO::m_outLineHeader; 
std::stringstream IO::m_outLine;
bool IO::m_shouldWriteOutLineHeader = true;
bool IO::m_shouldWriteOutLine = true;



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

	// Read file line by line
	while (m_inFl)
	{
		std::string strInput;
		IO::readLineInputFile(strInput);

		/**************************
		Read data based on keywords
		**************************/
		/*
			Read data to envir
		*/
		if (caseInsCompare(getKeyword(strInput), "Hydro"))
		{
			fowt.readHydroMode(getData(strInput));
			continue;
		}

		if (caseInsCompare(getKeyword(strInput), "Aero"))
		{
			fowt.readAeroMode(getData(strInput));
			continue;
		}

		if (caseInsCompare(getKeyword(strInput), "Moor"))
		{
			fowt.readMoorMode(getData(strInput));
			continue;
		}

		if (caseInsCompare(getKeyword(strInput), "DOFS"))
		{
			fowt.readDOFs(getData(strInput));
			continue;
		}		

		if (caseInsCompare(getKeyword(strInput), "TimeStep"))
		{
			envir.readTimeStep(getData(strInput));
			continue;
		}

		else if (caseInsCompare(getKeyword(strInput), "TimeTotal"))
		{
			envir.readTimeTotal(getData(strInput));
			continue;
		}

		else if (caseInsCompare(getKeyword(strInput), "TimeRamp"))
		{
			envir.readTimeRamp(getData(strInput));
			continue;
		}

		else if (caseInsCompare(getKeyword(strInput), "UseTipLoss"))
		{
			envir.readUseTipLoss(getData(strInput));
			continue;
		}

		else if (caseInsCompare(getKeyword(strInput), "UseHubLoss"))
		{
			envir.readUseHubLoss(getData(strInput));
			continue;
		}

		else if (caseInsCompare(getKeyword(strInput), "UseSkewCorr"))
		{
			envir.readUseSkewCorr(getData(strInput));
			continue;
		}

		else if (caseInsCompare(getKeyword(strInput), "Grav"))
		{
			envir.readGrav(getData(strInput));
			continue;
		}

		else if (caseInsCompare(getKeyword(strInput), "WatDens"))
		{
			envir.readWatDens(getData(strInput));
			continue;
		}

		else if (caseInsCompare(getKeyword(strInput), "WatDepth"))
		{
			envir.readWatDepth(getData(strInput));
			continue;
		}

		else if (caseInsCompare(getKeyword(strInput), "AirDens"))
		{
			envir.readAirDens(getData(strInput));
			continue;
		}

		else if (caseInsCompare(getKeyword(strInput), "WindVel"))
		{
			envir.readWindRefVel(getData(strInput));
			continue;
		}

		else if (caseInsCompare(getKeyword(strInput), "WindHeight"))
		{
			envir.readWindRefHeight(getData(strInput));
			continue;
		}

		else if (caseInsCompare(getKeyword(strInput), "WindExp"))
		{
			envir.readWindExp(getData(strInput));
			continue;
		}

		else if (caseInsCompare(getKeyword(strInput), "Wave")) // A list of Waves is supposed to follow the "Wave keyword"
		{
			IO::readLineInputFile(strInput); // Read next line, since current line is just the main keyword

			while (!caseInsCompare(getKeyword(strInput), "END")) // The END keyword indicates the end of the list of waves
			{
				if (!m_inFl) // Signal if the end of file is reached before the end keyword
				{
					throw std::runtime_error("End of file reached before END keyword in WAVE specification.");
					return;
				}

				envir.addWave(Wave(strInput)); // Add this wave to the environment

				IO::readLineInputFile(strInput);
			}
			continue;
		}

		else if (caseInsCompare(getKeyword(strInput), "Nodes")) // A list of Nodes is supposed to follow the "Wave keyword"
		{
			IO::readLineInputFile(strInput); // Read next line, since current line is just the main keyword

			while (!caseInsCompare(getKeyword(strInput), "END"))
			{
				if (!m_inFl) // Signal if the end of file is reached before the end keyword
				{
					throw std::runtime_error("End of file reached before END keyword in NODES specification.");
					return;
				}

				envir.addNode(strInput); // Add this node to the floater

				IO::readLineInputFile(strInput);
			}
			continue;
		}


		/*
			Read data to fowt
		*/
		else if (caseInsCompare(getKeyword(strInput), "LinStiff"))
		{
			fowt.readLinStiff(getData(strInput));
			continue;
		}

		else if (caseInsCompare(getKeyword(strInput), "FloaterMass"))
		{
			floater.readMass(getData(strInput));
			continue;
		}

		else if (caseInsCompare(getKeyword(strInput), "FloaterInertia"))
		{
			floater.readInertia(getData(strInput));
			continue;
		}

		else if (caseInsCompare(getKeyword(strInput), "FloaterCoG"))
		{
			floater.readCoG(getData(strInput));
			continue;
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

				floater.addMorisonCirc(strInput, envir); // Add this Morison Element to the floater. Need envir for the nodes location.

				IO::readLineInputFile(strInput);
			}
			continue;
		}


		else if (caseInsCompare(getKeyword(strInput), "Morison_rect")) // A list of rectangular cylinder Morison Elements is supposed to follow the "Morison_circ" keyword
		{
			IO::readLineInputFile(strInput); // Read next line, since current line is just the main keyword

			while (!caseInsCompare(getKeyword(strInput), "END"))
			{
				if (!m_inFl) // Signal if the end of file is reached before the end keyword
				{
					throw std::runtime_error("End of file reached before END keyword in MORISON_RECT specification.");
					return;
				}

				floater.addMorisonRect(strInput, envir); // Add this Morison Element to the floater. Need envir for the nodes location.

				IO::readLineInputFile(strInput);
			}
			continue;
		}


		else if (caseInsCompare(getKeyword(strInput), "RotSpeed"))
		{
			rna.readRotorSpeed(getData(strInput));
			continue;
		}


		else if (caseInsCompare(getKeyword(strInput), "RotTilt"))
		{
			rna.readRotorTilt(getData(strInput));
			continue;
		}


		else if (caseInsCompare(getKeyword(strInput), "RotYaw"))
		{
			rna.readRotorYaw(getData(strInput));
			continue;
		}


		else if (caseInsCompare(getKeyword(strInput), "BldPitch"))
		{
			rna.readBladePitch(getData(strInput));
			continue;
		}


		else if (caseInsCompare(getKeyword(strInput), "BldPrecone"))
		{
			rna.readBladePrecone(getData(strInput));
			continue;
		}


		else if (caseInsCompare(getKeyword(strInput), "NumBlades"))
		{
			rna.readNumBlades(getData(strInput));
			continue;
		}


		else if (caseInsCompare(getKeyword(strInput), "HubRadius"))
		{
			rna.readHubRadius(getData(strInput));
			continue;
		}


		else if (caseInsCompare(getKeyword(strInput), "HubHeight"))
		{
			rna.readHubHeight(getData(strInput));
			continue;
		}


		else if (caseInsCompare(getKeyword(strInput), "Overhang"))
		{
			rna.readOverhang(getData(strInput));
			continue;
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

				rna.readBladeAeroLine(strInput);

				IO::readLineInputFile(strInput);
			}
			continue;
		}


		else if (caseInsCompare(getKeyword(strInput), "Airfoil_data"))
		{
			rna.addAirfoil(); // Add new airfoil to rna
			IO::readLineInputFile(strInput); // Read next line, since current line is just the main keyword

			while (!caseInsCompare(getKeyword(strInput), "END"))
			{
				if (!m_inFl) // Signal if the end of file is reached before the end keyword
				{
					throw std::runtime_error("End of file reached before END keyword in AIRFOIL_DATA specification.");
					return;
				}

				rna.readAirfoilLine(strInput);

				IO::readLineInputFile(strInput);
			}
			continue;
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
			continue;
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
			continue;
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
			continue;
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
			continue;
		}

		else if (!caseInsCompare(getKeyword(strInput), "END_OF_INPUT_FILE"))
		{
			print2log("WARNING: Unknown keyword '" + getKeyword(strInput) + "' in line " + std::to_string(IO::getInLineNumber()) + ".");
		}
	}

	fowt.setFloater(floater);
	fowt.setRNA(rna);
}



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
	std::getline(m_inFl, strInput); // Read next file line to string strInput
	++m_inLineNumber; // Update line number counter

					  // Repeat this process until the line has some content or end of file is achieved
	while (!hasContent(strInput) && m_inFl)
	{
		std::getline(m_inFl, strInput);
		++m_inLineNumber;
	}

	// Remove comments from line
	if (thereIsCommentInString(strInput))
	{
		removeComments(strInput);
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
	if (caseInsCompare(getKeyword(strInput), "fowt_disp"))
	{
		m_whichResult2Output.at(IO::OUTFLAG_FOWT_DISP) = true;
	}

	if (caseInsCompare(getKeyword(strInput), "fowt_vel"))
	{
		m_whichResult2Output.at(IO::OUTFLAG_FOWT_VEL) = true;
	}

	if (caseInsCompare(getKeyword(strInput), "fowt_acc"))
	{
		m_whichResult2Output.at(IO::OUTFLAG_FOWT_ACC) = true;
	}

	if (caseInsCompare(getKeyword(strInput), "wave_elev"))
	{
		m_whichResult2Output.at(IO::OUTFLAG_WAVE_ELEV) = true;

		if (!getData(strInput).empty())
		{
			envir.addWaveLocation(getData(strInput)); // Add the wave locations to envir. Pass only the part of the string after the keyword
		}
	}

	if (caseInsCompare(getKeyword(strInput), "wave_vel")) // Similar to wave_elev
	{
		m_whichResult2Output.at(IO::OUTFLAG_WAVE_VEL) = true;

		if (!getData(strInput).empty())
		{
			envir.addWaveLocation(getData(strInput));
		}
	}

	if (caseInsCompare(getKeyword(strInput), "wave_acc")) // Similar to wave_elev
	{
		m_whichResult2Output.at(IO::OUTFLAG_WAVE_ACC) = true;

		if (!getData(strInput).empty())
		{
			envir.addWaveLocation(getData(strInput));
		}
	}

	if (caseInsCompare(getKeyword(strInput), "hd_force"))
	{
		m_whichResult2Output.at(IO::OUTFLAG_HD_FORCE) = true;
	}

	if (caseInsCompare(getKeyword(strInput), "hs_force"))
	{
		m_whichResult2Output.at(IO::OUTFLAG_HS_FORCE) = true;
	}

	if (caseInsCompare(getKeyword(strInput), "total_force"))
	{
		m_whichResult2Output.at(IO::OUTFLAG_TOTAL_FORCE) = true;
	}	

	if (caseInsCompare(getKeyword(strInput), "hd_inertia_force"))
	{
		m_whichResult2Output.at(IO::OUTFLAG_HD_INERTIA_FORCE) = true;
	}

	if (caseInsCompare(getKeyword(strInput), "hd_drag_force"))
	{
		m_whichResult2Output.at(IO::OUTFLAG_HD_DRAG_FORCE) = true;
	}

	if (caseInsCompare(getKeyword(strInput), "hd_fk_force"))
	{
		m_whichResult2Output.at(IO::OUTFLAG_HD_FK_FORCE) = true;
	}

	if (caseInsCompare(getKeyword(strInput), "ad_hub_force"))
	{
		m_whichResult2Output.at(IO::OUTFLAG_AD_HUB_FORCE) = true;
	}
	
}



static void checkInputs(const FOWT &fowt, const ENVIR &envir)
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
	std::cout << "\n\n" << str << std::endl;
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

// Different functions to print to the formatted output file string streams depending on the output flag
void IO::print2outLine(const OutFlag &flag, const arma::vec::fixed<6> &vector_6)
{
	// Check whether the specified flag is indeed one that requires a vector with six components
	if ((flag != OUTFLAG_FOWT_DISP) && (flag != OUTFLAG_FOWT_VEL) && (flag != OUTFLAG_FOWT_ACC) && 
		(flag != OUTFLAG_HD_FORCE) && (flag != OUTFLAG_HS_FORCE) && (flag != OUTFLAG_TOTAL_FORCE) &&
		(flag != OUTFLAG_HD_INERTIA_FORCE) && (flag != OUTFLAG_HD_DRAG_FORCE) && (flag != OUTFLAG_HD_FK_FORCE) &&
		(flag != OUTFLAG_AD_HUB_FORCE)
	   )
	{
		throw std::runtime_error("Unknown output flag in function IO::print2outLine(const OutFlag &flag, const arma::vec::fixed<6> &force).");
	}

	// If the print header flag is true and if this is one of the requested output variables,
	// then print the header based on the output flag	
	if (m_shouldWriteOutLineHeader && m_whichResult2Output.at(flag))
	{
		if (flag == OUTFLAG_HD_FORCE)
		{
			for (int ii = 1; ii <= 6; ++ii)
			{
				print2outLineHeader("HD_Force_" + std::to_string(ii));
			}
		}

		if (flag == OUTFLAG_HS_FORCE)
		{
			for (int ii = 1; ii <= 6; ++ii)
			{
				print2outLineHeader("HS_Force_" + std::to_string(ii));
			}
		}

		if (flag == OUTFLAG_TOTAL_FORCE)
		{
			for (int ii = 1; ii <= 6; ++ii)
			{
				print2outLineHeader("TOTAL_Force_" + std::to_string(ii));
			}
		}		

		if (flag == OUTFLAG_HD_INERTIA_FORCE)
		{
			for (int ii = 1; ii <= 6; ++ii)
			{
				print2outLineHeader("HD_INRT_Force_" + std::to_string(ii));
			}
		}

		if (flag == OUTFLAG_HD_DRAG_FORCE)
		{
			for (int ii = 1; ii <= 6; ++ii)
			{
				print2outLineHeader("HD_DRAG_Force_" + std::to_string(ii));
			}
		}

		if (flag == OUTFLAG_HD_FK_FORCE)
		{
			for (int ii = 1; ii <= 6; ++ii)
			{
				print2outLineHeader("HD_FK_Force_" + std::to_string(ii));
			}
		}

		if (flag == OUTFLAG_AD_HUB_FORCE)
		{
			for (int ii = 1; ii <= 6; ++ii)
			{
				print2outLineHeader("AD_HUB_Force_" + std::to_string(ii));
			}
		}		

		if (flag == OUTFLAG_FOWT_DISP)
		{
			print2outLineHeader("Surge");
			print2outLineHeader("Sway");
			print2outLineHeader("Heave");
			print2outLineHeader("Roll");
			print2outLineHeader("Pitch");
			print2outLineHeader("Yaw");
		}

		if (flag == OUTFLAG_FOWT_VEL)
		{
			print2outLineHeader("Surge_Vel");
			print2outLineHeader("Sway_Vel");
			print2outLineHeader("Heave_Vel");
			print2outLineHeader("Roll_Vel");
			print2outLineHeader("Pitch_Vel");
			print2outLineHeader("Yaw_Vel");
		}

		if (flag == OUTFLAG_FOWT_ACC)
		{
			print2outLineHeader("Surge_Acc");
			print2outLineHeader("Sway_Acc");
			print2outLineHeader("Heave_Acc");
			print2outLineHeader("Roll_Acc");
			print2outLineHeader("Pitch_Acc");
			print2outLineHeader("Yaw_Acc");
		}
	}


	// If the printing flag is true and if this is one of the requested output variables,
	// then print it to the output line
	if (m_shouldWriteOutLine && m_whichResult2Output.at(flag))
	{
		for (int ii = 0; ii < 6; ++ii)
		{
			print2outLine(vector_6.at(ii));
		}
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

	m_outFl << std::setw(IO::m_outColumnWidth) << "Time"; // First column must always be the current time
	m_outFl << m_outLineHeader.str() << '\n'; // Then comes the results that are stored in the header stringstream
	m_outLineHeader.str(""); // Need to clear the stream for the next time step
}

void IO::printOutLine2outFile()
{
	if (!m_outFl)
	{
		throw std::runtime_error("Unable to write to formatted output file.");
	}

	m_outFl << m_outLine.str() << '\n';
	m_outLine.str(""); // Need to clear the stream for the next time step
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
	m_sumFl << "Time Step:\t" << envir.timeStep() << '\n';
	m_sumFl << "Total Time:\t" << envir.timeTotal() << '\n';
	m_sumFl << "Time Ramp:\t" << envir.printTimeRamp() << '\n';
	m_sumFl << "Use Tip Loss:\t" << envir.useTipLoss() << '\n';
	m_sumFl << "Use Hub Loss:\t" << envir.useHubLoss() << '\n';
	m_sumFl << "Use Skew Correction:\t" << envir.useSkewCorr() << '\n';
	m_sumFl << "Gravity:\t" << envir.gravity() << '\n';
	m_sumFl << "Water Density:\t" << envir.watDensity() << '\n';
	m_sumFl << "Water Depth:\t" << envir.watDepth() << '\n';
	m_sumFl << "Air density:\t" << envir.airDensity() << '\n';
	m_sumFl << "Wind X velocity:\t" << envir.windRefVel() << '\n';
	m_sumFl << "Wind Height:\t" << envir.windRefHeight() << '\n';
	m_sumFl << "Wind exp:\t" << envir.windExp() << '\n';
	m_sumFl << "Nodes: \n" << envir.printNodes() << '\n';
	m_sumFl << "Wave Locations: " << envir.printWaveLocation() << '\n';
	m_sumFl << "\n" << envir.printWave() << '\n';

	m_sumFl << "\n\n";
	m_sumFl << "FOWT:\n";
	m_sumFl << "Hydro Mode:\t" << fowt.printHydroMode() << "\n";
	m_sumFl << "Aero Mode:\t" << fowt.printAeroMode() << "\n";
	m_sumFl << "Moor Mode:\t" << fowt.printMoorMode() << "\n";
	m_sumFl << "DOFs:\t" << fowt.printDoF() << '\n';
	m_sumFl << "Linear Stiffness:\t" << fowt.printLinStiff() << '\n';
	m_sumFl << "Floater:\n" << fowt.printFloater();
	m_sumFl << "RNA:\n" << fowt.printRNA();

	m_sumFl << "\n\n";
	m_sumFl << "Output Variables:\n" << IO::printOutVar();
}


// Some printing functions
std::string IO::printOutVar()
{
	std::string output = "";
	for (int ii = 0; ii < IO::OUTFLAG_SIZE; ++ii)
	{
		switch (ii)
		{
		case IO::OUTFLAG_FOWT_DISP:
			output += "FOWT rigid motion: ";
			break;

		case IO::OUTFLAG_FOWT_VEL:
			output += "FOWT rigid motion: ";
			break;

		case IO::OUTFLAG_FOWT_ACC:
			output += "FOWT rigid motion: ";
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

		case IO::OUTFLAG_HD_FORCE:
		    output += "Hydrodynamic force: ";
			break;

		case IO::OUTFLAG_HS_FORCE:
		    output += "Hydrostatic force: ";
			break;		

		case IO::OUTFLAG_HD_INERTIA_FORCE:
			output += "Hydrodynamic inertial force: ";
			break;

		case IO::OUTFLAG_HD_DRAG_FORCE:
			output += "Hydrodynamic drag force: ";
			break;

		case IO::OUTFLAG_HD_FK_FORCE:
			output += "Hydrodynamic Froude-Krylov force: ";
			break;

		case IO::OUTFLAG_AD_HUB_FORCE:
			output += "Aerodynamic forces (Hub): ";
			break;


		case IO::OUTFLAG_TOTAL_FORCE:
			output += "Total force: ";
			break;

		default:
			output += "Unknown specifier in output flags.";
			break;
		}
        output += std::to_string(m_whichResult2Output.at(ii)) + "\n";
	}
	return output;
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
		str_out += str_tokenized.at(ii) + "\t";
	}
	return str_out;
}