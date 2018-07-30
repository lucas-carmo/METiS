#pragma once

#include <string>
#include <vector>
#include <sstream>
#include <algorithm> // Defines a collection of functions especially designed to be used on ranges of elements.
#include <cctype> // This header declares a set of functions to classify and transform individual characters, like toupper
#include <cwctype> // Same thing for wide characters
#include "FOWT.h"
#include "ENVIR.h"


class IO
{
private:
	// METiS Version 
	static constexpr char m_METIS_VERSION[] = "0.0.1";

	// Members related to input
	static std::string m_inFilePath;
	static std::ifstream m_inFl;
	static unsigned int m_inLineNumber; // Stores the current line number while reading the input file

	// Members related to output of log file
	static std::string m_logFilePath;
	static std::ofstream m_logFl;
	
	// Members related to output of summary file
	static std::string m_sumFilePath;
	static std::ofstream m_sumFl;

	static std::string m_outputTimeStep; // String with the data that is output at each time step (FOWT position, hydro force components, anything that is a function of time)
	static bool m_shouldWriteOutputTimeStep;

public:
	// Set input file and output files based on the input file path
	static void setFiles(const std::string &inFlPath);

	// Functions related to Input	
	static void readLineInputFile(std::string &strInput);
	static unsigned int getInLineNumber();
	static void readInputFile(FOWT &fowt, ENVIR &envir);	

	// Functions related to output
	static std::string METiS_Header();
	static void writeLogMessage(const std::string &str);
	static void printSumFile(const FOWT &fowt, const ENVIR &envir);

	//static bool checkData(); // Depende do tipo de analise a ser feita
	
	static void print2outFile(const std::string &outFlNm);
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
inline bool readDataFromString(const std::string& inString, T &tX)
{
	std::vector<std::string> input = stringTokenize(inString, " \t");

	// 1) Verify if input contains exactly 1 element
    if ( input.size() != 1 )
	{
        std::cout << "More than one entry in line " << IO::getInLineNumber() << ". Considering only the first one. \n";		
	}

    // 2) Convert input from string to its corresponding numeric format (double, float, ...)
	if ( !string2num(input.at(0), tX) )
	{
        // Throw an exception if the conversion fails
		std::cout << "Conversion failed in readDataFromString(). Line " << IO::getInLineNumber() << ".\n";
		return false;
	}		

	return true;
}