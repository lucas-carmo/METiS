#include "IO.h"

#include <fstream> // Include file input/output classes
#include <iostream>
#include <cstdlib> // For exit()
#include <string>

// Function responsible for reading the input file and assigning what is read to the FOWT and ENVIR classes
void IO::readInputFile(const std::string &inFlNm, FOWT &fowt, ENVIR &envir)
{
    std::ifstream inFl(inFlNm);
    if (!inFl)
    {
        std::cerr << "Unable to open " << inFlNm << " for reading";
        exit(1);
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







/*
    Additional functions related to input/output 
*/

// Verify whether there is a comment in a string
bool thereIsCommentInString(const std::string& str)
{
    return (str.find("%") != std::string::npos);
}


// Verify whether a string has content (in the context of data from the input file)
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
    return str.substr(0, str.find_first_of(" \t"));
}



// Pensar no caso em que a linha tem só keyword quando deveria ter dado tmbm. Que erro retornar? Como checar?
std::string getData(const std::string &str)
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