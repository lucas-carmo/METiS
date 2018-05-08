#pragma once

#include <string>
#include <vector>
#include <sstream>


bool thereIsCommentInString(const std::string& str);

bool hasContent(const std::string &str);

void removeComments(std::string &str);

std::string getKeyword(const std::string &str);

std::string getData(const std::string &str);

std::vector<std::string> stringTokenize(const std::string &str, const std::string &delim = " \t");




/*Function templates*/

// FromString: used to convert from string to a numerical type (double, float, int...)
// Returns True if the conversion is succesful and False if it is not
// Found at http://www.learncpp.com/cpp-tutorial/17-2-stdstring-construction-and-destruction/
template <typename T>
inline bool FromString(const std::string& sString, T &tX)
{
	std::istringstream iStream(sString);
	return (iStream >> tX) ? true : false; // extract value into tX, return success or not
}