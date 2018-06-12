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

// This is the main input function.
// It is responsible for reading the input file line by line and assigning what is read to the FOWT and ENVIR classes
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

        if ( !hasContent(strInput) )
        {
			continue;
        }
        
		if (thereIsCommentInString(strInput))
		{
			removeComments(strInput);
		}

        /* 
		Read data based on keywords
		*/

		// Read data to envir class
        if ( getKeyword(strInput) == "TimeStep" )
        {
            envir.readTimeStep( getData(strInput) );
            continue;
        }

        if ( getKeyword(strInput) == "TimeTotal" )
        {
            envir.readTimeTotal( getData(strInput) );
            continue;
        }

        if (getKeyword(strInput) == "TimeRamp")
        {
            envir.readTimeRamp( getData(strInput) );
            continue;
        }

		if (getKeyword(strInput) == "Grav")
		{
			envir.readGrav(getData(strInput));
			continue;
		}

		if (getKeyword(strInput) == "WatDens")
		{
			envir.readWatDens(getData(strInput));
			continue;
		}

		if (getKeyword(strInput) == "WatDepth")
		{
			envir.readWatDepth(getData(strInput));
			continue;
		}



		// Read data to envir class
		if (getKeyword(strInput) == "LinStiff")
		{
			fowt.readLinStiff(getData(strInput));
			continue;
		}

		if (getKeyword(strInput) == "FloaterMass")
		{
			fowt.readFloaterMass(getData(strInput));
			continue;
		}

		if (getKeyword(strInput) == "FloaterCoG")
		{
			fowt.readFloaterCoG(getData(strInput));
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
std::string getData(const std::string &str)
{
    if ( str.find_first_of(" \t", 0) != std::string::npos )
    {
        std::string aux_str = str.substr(str.find_first_of(" \t", 0));
        return aux_str.substr(aux_str.find_first_not_of(" \t"));
    }
    else
    {
        return str;
    }    
}

// Tokenize a string using a given delimiter.
// Return a std::vector with the resulting strings at the different elements
std::vector<std::string> stringTokenize(const std::string &str, const std::string &delim)
{
    std::vector<std::string> tokens;
    std::string aux = str;

    while( hasContent(aux) ) 
    {
        if ( aux.find_first_not_of(delim) == std::string::npos ) // Check if there is something besides delimiters in the line
        {
            std::cout << "The string provided has only delimiters.";
            break;
        }

        aux = aux.substr(aux.find_first_not_of(delim)); // Get the part of aux starting at the first character that is not a delimiter
        tokens.push_back( aux.substr(0, aux.find_first_of(delim)) ); // Take content before next delimiter and add to tokens

        if ( aux.find_first_of(delim, 0) != std::string::npos ) // If there is another delimiter in the string...
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