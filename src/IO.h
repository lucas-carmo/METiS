#pragma once

#include <string>
#include <vector>
#include <sstream>
#include "FOWT.h"
#include "ENVIR.h"


class IO
{
private:
	static std::string m_header;
	static bool m_shouldWriteHeader;

	static std::string m_outputTimeStep; // String with the data that is output at each time step (FOWT position, hydro force components, anything that is a function of time)
	static bool m_shouldWriteOutputTimeStep;

public:
	static void readInputFile(const std::string &inFlNm, FOWT &fowt, ENVIR &envir);
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
}

