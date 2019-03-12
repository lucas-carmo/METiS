#include "auxFunctions.h"
#include <cctype> // This header declares a set of functions to classify and transform individual characters, like toupper
#include <cwctype> // Same thing for wide characters

/*****************************************************
    Useful math/geometric operations
*****************************************************/

mat::fixed<3, 3> rotatMatrix(const vec::fixed<3> &rotation)
{
	/* Rotation matrix is RotatMat = RotatX * RotatY * RotatZ, i.e. a rotation around the Z axis,
	  followed by a rotation around the new y axis, and a rotation around the new x axis. Each rotation matrix is given by:
	
	  rotatX = { 	{ 1 ,        0        ,         0        },
			     			{ 0 , cos(rotation(3)) , -sin(rotation(3)) },
			     			{ 0 , sin(rotation(3)) ,  cos(rotation(3)) }
			   			};

	  rotatY = { 	{ cos(rotation(4))  , 0 ,  sin(rotation(4)) },
		         		{        0         , 1 ,         0        },
			     			{ -sin(rotation(4)) , 0 , cos(rotation(4)) }
			   			};

	  rotatZ = {	{ cos(rotation(5)) , -sin(rotation(5)) , 0 },			     
				 				{ sin(rotation(5)) ,  cos(rotation(5)) , 0 },
			     			{        0        ,         0        , 1 },
			   			};

	  And the resulting matrix is the one below
	*/
	mat::fixed<3, 3> rotatMatrix = { 
									{                          cos(rotation(1)) * cos(rotation(2))                               ,                          -cos(rotation(1)) * sin(rotation(2))                               ,            sin(rotation(1))          },
									{ cos(rotation(0)) * sin(rotation(2)) + sin(rotation(0)) * sin(rotation(1)) * cos(rotation(2))  ,  cos(rotation(0)) * cos(rotation(1)) - sin(rotation(0)) * sin(rotation(1)) * sin(rotation(2))  ,  -sin(rotation(0)) * cos(rotation(1)) },
									{ sin(rotation(0)) * sin(rotation(2)) - cos(rotation(0)) * sin(rotation(1)) * cos(rotation(2))  ,  sin(rotation(0)) * cos(rotation(2)) + cos(rotation(0)) * sin(rotation(1)) * sin(rotation(2))  ,   cos(rotation(0)) * cos(rotation(1)) }
								   };

	return rotatMatrix;
}

mat::fixed<3, 3> rotatMatrix_deg(const vec::fixed<3> &rotation)
{
	vec::fixed<3> rotRad;
	rotRad[0] = rotation[0] * 180 / datum::pi;
	rotRad[1] = rotation[1] * 180 / datum::pi;
	rotRad[2] = rotation[2] * 180 / datum::pi;
	return rotatMatrix(rotRad);
}


/*****************************************************
    String handling
*****************************************************/

// Verify whether a string contains a comment, marked by a '%'
bool thereIsCommentInString(const std::string& str)
{
	return (str.find("%") != std::string::npos);
}

// Verify whether a string has content, i.e. if:
// 1) it is not empty
// 2) it is not just white spaces or tabs
// 3) it does not start with a comment mark ('%')
bool hasContent(const std::string &str)
{
	// Empty strings have no content
	if (str.empty())
		return false;

	// The string has content only if it has at least one character that is neither a white-space nor a tab (\t)
	return ((str.find_first_not_of(" \t") != std::string::npos) && (str.at(0) != '%'));
}

void removeComments(std::string &str)
{
	str = str.substr(0, str.find("%", 0));
}

// Tokenize a string using a given delimiter.
// Return a std::vector with the resulting strings at the different elements
std::vector<std::string> stringTokenize(const std::string &str, const std::string &delim)
{
	std::vector<std::string> tokens;
	std::string aux = str;

	while (hasContent(aux))
	{
		if (aux.find_first_not_of(delim) == std::string::npos) // Check if there is something besides delimiters in the line
		{
			break; // Then we break the while loop and return tokens as an empty std::vector
		}

		aux = aux.substr(aux.find_first_not_of(delim)); // Get the part of aux starting at the first character that is not a delimiter
		tokens.push_back(aux.substr(0, aux.find_first_of(delim))); // Take content before next delimiter and add to tokens

		if (aux.find_first_of(delim, 0) != std::string::npos) // If there is another delimiter in the string...
		{
			aux = aux.substr(aux.find_first_of(delim, 0)); // ... keep in aux everything after the second delimiter, including the delimiter
		}
		else // Otherwise, we have already included the last element in tokens, so we can end this loop.
		{
			break;
		}
	}

	return tokens;
}


// Convert lowercase letter to uppercase and compare if they are equal
inline bool caseInsCharCompareN(char a, char b) {
	return(std::toupper(a) == std::toupper(b));
}

// Same thing for wchar
inline bool caseInsCharCompareW(wchar_t a, wchar_t b) {
	return(std::towupper(a) == std::towupper(b));
}

// Compare each character of the strings to see if they match
bool caseInsCompare(const std::string& s1, const std::string& s2) {
	return((s1.size() == s2.size()) &&
		std::equal(s1.begin(), s1.end(), s2.begin(), caseInsCharCompareN));
}

// Same thing for wstring
bool caseInsCompare(const std::wstring& s1, const std::wstring& s2) {
	return((s1.size() == s2.size()) &&
		std::equal(s1.begin(), s1.end(), s2.begin(), caseInsCharCompareW));
}


// Get folder path from a complete file path
std::string getFileFolder(const std::string& path)
{
	// Check if input string is empty
	if (path.empty())
	{
		throw std::runtime_error("Empty string passed to getFileFolder().");
	}

	std::vector<std::string> str_tokenized = stringTokenize(path, filesep);

	// If there is only one element in str_tokenized, it means that only the file name
	// was provided as an argument. Hence, the fileFolder is the current directory, "."
	if (str_tokenized.size() == 1)
	{
		return ".";
	}

	// Otherwise, return the full path until the last delimiter
	else
	{
		size_t found;
		found = path.find_last_of(filesep);
		return path.substr(0, found);
	}
}


// Get file name, without extension, from a complete file path
std::string getFileName(const std::string& path)
{
	// Check if input string is empty
	if (path.empty())
	{
		throw std::runtime_error("Empty string passed to getFileName().");
	}

	// Tokenize the passed string based on the file separator. The file name is the last part (works even if only the file name is passed as argument)
	std::vector<std::string> str_tokenized = stringTokenize(path, filesep);
	std::string flNm = str_tokenized.back();

	// We need to exclude the file extension. We get everything until the last dot and say that 
	if (flNm.find_last_not_of(".") != std::string::npos)
	{
		flNm = flNm.substr(0, flNm.find_last_of("."));
	}

	return flNm;
}