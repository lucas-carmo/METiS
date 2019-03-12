#pragma once
#include <armadillo>
#include <string>

// Define file separator for current platform
const std::string filesep =
#ifdef _WIN32
	"\\";
#else
	"/";
#endif


/*****************************************************
    Useful math/geometric operations
*****************************************************/
using namespace arma;

mat::fixed<3, 3> rotatMatrix(const vec::fixed<3> &rotation);	
mat::fixed<3, 3> rotatMatrix_deg(const vec::fixed<3> &rotation);


/*****************************************************
    String handling
*****************************************************/

// Verify whether a string contains a comment, marked by a '%'
// TODO: define comment character in the beginning of this file
bool thereIsCommentInString(const std::string& str);

// Verify whether a string has content, i.e. if it is not empty, it is not just
// white spaces or tabs, and it does not start with a comment mark ('%')
bool hasContent(const std::string &str);

// Remove comments from string, marked by a '%'
void removeComments(std::string &str);

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

// Get file name, without extension, from a complete file path
std::string getFileName(const std::string& path);

// string2num: used to convert from string to a numerical type (double, float, int...)
// Returns True if the conversion is succesful and False if it is not
// Found at http://www.learncpp.com/cpp-tutorial/17-2-stdstring-construction-and-destruction/
template <typename T>
inline bool string2num(const std::string& sString, T &tX)
{
	std::istringstream iStream(sString);
	return (iStream >> tX) ? true : false; // extract value into tX, return success or not
}

