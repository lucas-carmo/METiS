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



/* 
Additional functions related to input/output 
*/
bool thereIsCommentInString(const std::string& str);

bool hasContent(const std::string &str);

void removeComments(std::string &str);

std::string getKeyword(const std::string &str);

std::string getData(const std::string &str);

std::vector<std::string> stringTokenize(const std::string &str, const std::string &delim = " \t");


/* 
Function templates 
*/

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

