#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <sstream>
#include <algorithm> // Defines a collection of functions especially designed to be used on ranges of elements.
#include <cctype> // This header declares a set of functions to classify and transform individual characters, like toupper
#include <cwctype> // Same thing for wide characters
#include "FOWT.h"
#include "ENVIR.h"

// Forward declaration of METiS Version
extern const std::string g_METIS_VERSION;

// Define file separator for current platform
const std::string filesep =
#ifdef _WIN32
	"\\";
#else
	"/";
#endif


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
		OUTFLAG_FOWT_POS,
//		
		OUTFLAG_WAVE_ELEV,
		OUTFLAG_WAVE_VEL,
		OUTFLAG_WAVE_ACC,
//		
		OUTFLAG_HD_FORCE,
		OUTFLAG_HS_FORCE,		
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
	static std::stringstream m_outLine; // String stream with the data that is output at each time step (FOWT position, hydro force components, anything that is a function of time)
	static bool m_shouldWriteOutLineHeader;
	static bool m_shouldWriteOutLine;	


public:
	// Set input file and output files based on the input file path
	static void setFiles(const std::string &inFlPath);

	// Functions related to Input	
	static void readLineInputFile(std::string &strInput);
	static unsigned int getInLineNumber();
	static void readInputFile(FOWT &fowt, ENVIR &envir);	
	static void setResults2Output(std::string strInput, FOWT &fowt, ENVIR &envir);

	/*
	Functions related to output
	*/
	static std::string METiS_Header();
	static void writeErrorMessage(const std::string &str);
	static void writeWarningMessage(const std::string &str);	

	// To summary file
	static void printSumFile(const FOWT &fowt, const ENVIR &envir);
	
	// To formatted output file
	static void print2outLineHeader_turnOn(); 
	static void print2outLineHeader_turnOff();		
	static void print2outLine_turnOn(); 
	static void print2outLine_turnOff();	

	// Functions that write to the stringstreams m_outLineHeader and m_outLine
	static void print2outLine(const OutFlag &flag, const arma::vec::fixed<6> &vector_6);
	static void print2outLine(const std::string &str);
	static void print2outLine(const double num);
	static void print2outLine(const int num);
	static void print2outLineHeader(const std::string &str);

		// Functions that actually write to the output file
	static void printOutLineHeader2outFile();
	static void printOutLine2outFile();
	static void newLineOutFile();

	// Some printing functions
	static std::string printOutVar();
};



/*****************************************************
    Additional functions related to input/output 
*****************************************************/

// Verify whether a string contains a comment, marked by a '%'
bool thereIsCommentInString(const std::string& str);

// Verify whether a string has content, i.e. if it is not empty, it is not just
// white spaces or tabs, and it does not start with a comment mark ('%')
bool hasContent(const std::string &str);

// Remove comments from string, marked by a '%'
void removeComments(std::string &str);

std::string getKeyword(const std::string &str);

// Get the part of the string after the keyword, excluding the '\t' or white-space
std::string getData(const std::string &str);

// Tokenize a string using a given delimiter
std::vector<std::string> stringTokenize(const std::string &str, const std::string &delim = " \t");


// Case-insensitive string comparison
// Found at https://www.safaribooksonline.com/library/view/c-cookbook/0596007612/ch04s14.html

// Convert lowercase letter to uppercase and compare if they are equal
inline bool caseInsCharCompareN(char a, char b);

// Same thing for wchar
inline bool caseInsCharCompareW(wchar_t a, wchar_t b);

// Compare each character of the strings to see if they match
bool caseInsCompare(const std::string& s1, const std::string& s2);

// Same thing for wstring
bool caseInsCompare(const std::wstring& s1, const std::wstring& s2);

// Get folder path from a complete file path
std::string getFileFolder(const std::string& path);



/*****************************************************
	Function templates 
*****************************************************/

// string2num: used to convert from string to a numerical type (double, float, int...)
// Returns True if the conversion is succesful and False if it is not
// Found at http://www.learncpp.com/cpp-tutorial/17-2-stdstring-construction-and-destruction/
template <typename T>
inline bool string2num(const std::string& sString, T &tX)
{
	std::istringstream iStream(sString);
	return (iStream >> tX) ? true : false; // extract value into tX, return success or not
}


// readDataFromString: used to read data from the string 'inString' into the variable 'tX'
// Returns True if the conversion is succesful and False if it is not
template <typename T>
inline void readDataFromString(const std::string& inString, T &tX)
{
	std::vector<std::string> input = stringTokenize(inString, " \t");

	// 1) Verify if input contains exactly 1 element
    if ( input.size() != 1 )
	{
		throw std::runtime_error( "Unable to read data from string. More than one entry in line " + std::to_string(IO::getInLineNumber()) + ". [function readDataFromString(const std::string& inString, T &tX) in IO.h]");
	}

    // 2) Convert input from string to its corresponding numeric format (double, float, ...)
	if ( !string2num(input.at(0), tX) )
	{
        // Throw an exception if the conversion fails
		throw std::runtime_error( "Conversion from string failed. Bad data type in line " + std::to_string(IO::getInLineNumber()) + ". [function readDataFromString(const std::string& inString, T &tX) in IO.h]");
	}		
}