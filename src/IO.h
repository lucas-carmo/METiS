#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <sstream>
#include <algorithm> // Defines a collection of functions especially designed to be used on ranges of elements.
#include "FOWT.h"
#include "ENVIR.h"
#include "auxFunctions.h" // Include functions that are useful for string manipulation (among others) everytime IO.h is included

// Forward declaration of METiS Version
extern const std::string g_METIS_VERSION;

class IO
{
public:

	// Enum for the indices of the variables that can be output
	// They are used in the functions:
	//		- setResults2Output(std::string strInput, FOWT &fowt, ENVIR &envir);
	//		- in the different print2outLine(const OutFlag &flag, 'some output variable');
	//		- printOutVar()
	enum OutFlag
	{
		OUTFLAG_FOWT_DISP,
		OUTFLAG_FOWT_VEL,
		OUTFLAG_FOWT_ACC,
		OUTFLAG_FOWT_DISP_SD,
//
		OUTFLAG_WAVE_ELEV,
		OUTFLAG_WAVE_VEL,
		OUTFLAG_WAVE_ACC,
		OUTFLAG_WAVE_PRES,
		OUTFLAG_WAVE_ACC_2ND,
//
		OUTFLAG_HD_INERTIA_FORCE,
		OUTFLAG_HD_DRAG_FORCE,
		OUTFLAG_HD_FK_FORCE,
		OUTFLAG_HD_2ND_FORCE_PART1,
		OUTFLAG_HD_FORCE,
//
		OUTFLAG_HS_FORCE,
//
		OUTFLAG_AD_HUB_FORCE,
//
		OUTFLAG_TOTAL_FORCE,
//
		OUTFLAG_SIZE
	};

private:
	// Members related to input
	static std::string m_inFilePath;
	static std::ifstream m_inFl;
	static unsigned int m_inLineNumber; // Stores the current line number while reading the input file

	// Members related to log file
	static std::string m_logFilePath;
	static std::ofstream m_logFl;

	// Members related to summary file
	static std::string m_sumFilePath;
	static std::ofstream m_sumFl;

	// Members related to formatted output file
	static std::string m_outFilePath;
	static std::ofstream m_outFl;
	static const unsigned int m_outColumnWidth;
	static const unsigned int m_outNumPrecision;
	static std::array<bool, IO::OUTFLAG_SIZE> m_whichResult2Output;

	static std::stringstream m_outLineHeader; // String stream with the header identifying each column of the formatted output file
	static std::stringstream m_outLine; // String stream with the data that is output at each time step (FOWT displacement, hydro force components, anything that is a function of time)
	static bool m_shouldWriteOutLineHeader;
	static bool m_shouldWriteOutLine;


public:
	// Set input file and output files based on the input file path
	static void setFiles(const std::string &inFlPath);

	// Functions related to Input
	static void readLineInputFile(std::string &strInput);
	static unsigned int getInLineNumber();
	static void readInputFile(FOWT &fowt, ENVIR &envir);
	static void setResults2Output(std::string strInput, ENVIR &envir);
	static void checkInputs(const FOWT &fowt, const ENVIR &envir);

	/*
	Functions related to output
	*/
	static std::string METiS_Header();
	static void print2log(const std::string &str);

	// To summary file
	static void printSumFile(const FOWT &fowt, const ENVIR &envir);

	// To formatted output file
	static void print2outLineHeader_turnOn();
	static void print2outLineHeader_turnOff();
	static void print2outLine_turnOn();
	static void print2outLine_turnOff();

	// Functions that write to the stringstreams m_outLineHeader and m_outLine
	static void print2outLine(const OutFlag &flag, const arma::vec::fixed<6> &vector_6);
	static void print2outLine(const OutFlag &flag, const int ID, const double num);
	static void print2outLine(const OutFlag &flag, const int ID, const arma::vec::fixed<3> &vector_3);

	static void print2outLine(const std::string &str);
	static void print2outLine(const double num);
	static void print2outLine_decimal(const double num);
	static void print2outLine(const int num);
	static void print2outLineHeader(const std::string &str);

	// Functions that actually write to the output file
	static void printOutLineHeader2outFile();
	static void printOutLine2outFile();
	static void newLineOutFile();

	// Other printing functions
	static std::string printOutVar();
};



/*****************************************************
    Additional functions and templates related to input/output
*****************************************************/
std::string getKeyword(const std::string &str);

// Get the part of the string after the keyword, excluding the '\t' or white-space
std::string getData(const std::string &str);
