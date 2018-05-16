#include "IO.h"

#include <fstream> // Include file input/output classes
#include <iostream>
#include <cstdlib> // For exit()
#include <string>



/*****************************************************
	Defining and initializing static member variables
*****************************************************/
std::string IO::m_inFilePath = "";



/*****************************************************
    Implementation of class functions
*****************************************************/
void IO::setInputFilePath(const std::string &inFlPath)
{
	m_inFilePath = inFlPath;
}

// Function responsible for reading the input file and assigning what is read to the FOWT and ENVIR classes
void IO::readInputFile(FOWT &fowt, ENVIR &envir)
{
    std::ifstream inFl(m_inFilePath);
    if (!inFl)
    {
        std::cerr << "Unable to open " << m_inFilePath << " for reading";
		return;
        //exit(1);
    }


    while (inFl)
    {
        std::string strInput;
        std::getline(inFl, strInput);

        if ( hasContent(strInput) )
        {
            if ( thereIsCommentInString(strInput) )
            {
                removeComments(strInput);
            }
        }
        else
        {
            continue;
        }

        // Read data base on keywords
        if ( getKeyword(strInput) == "TimeStep" )
        {
            envir.readTimeStep( strInput );
            continue;
        }

        if ( getKeyword(strInput) == "TimeTotal" )
        {
            envir.readTimeTotal( strInput );
            continue;
        }

        if (getKeyword(strInput) == "TimeRamp")
        {
            envir.readTimeRamp( strInput );
            continue;
        }
    }
}


void IO::print2outFile(const std::string &outFlNm)
{}







/*****************************************************
    Additional functions related to input/output 
*****************************************************/

// Verify whether a string contains a comment, marked by a '%'
bool thereIsCommentInString(const std::string& str)
{
    return (str.find("%") != std::string::npos);
}


// Verify whether a string has content, i.e. if it is not empty, it is not just
// white spaces or tabs, and it does not start with a comment mark ('%')
bool hasContent(const std::string &str)
{
    // Empty strings have no content
    if (str.empty())
        return false;

    // The string has content only if it has at least one character that is neither a white-space nor a tab (\t)
    return ( (str.find_first_not_of(" \t") != std::string::npos) && (str.at(0) != '%') );
}



void removeComments(std::string &str)
{
    str = str.substr(0, str.find("%", 0));
}


std::string getKeyword(const std::string &str)
{
    // Get the first part of the string, the one before the first '\t' or white-space
    return str.substr(0, str.find_first_of(" \t"));
}

// Get the part of the string after the keyword, excluding the '\t' or white-space
std::string getContent(const std::string &str)
{
    std::string aux_str = str.substr(str.find_first_of(" \t", 0));
    return aux_str.substr(aux_str.find_first_not_of(" \t"));
}

// Tokenize a string using a given delimiter
std::vector<std::string> stringTokenize(const std::string &str, const std::string &delim)
{
    std::vector<std::string> tokens;
    std::string aux = str;

    while( hasContent(aux) ) {
        aux = aux.substr(aux.find_first_not_of(delim)); // Get the part of aux after the first delimiter
        tokens.push_back( aux.substr(0, aux.find_first_of(delim)) ); // Take content before next delimiter and add to tokens
        aux = aux.substr(aux.find_first_of(delim, 0)); // Keep in aux everything after the second delimiter, including the delimiter
    }

    return tokens;
}