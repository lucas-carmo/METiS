#include "IO.h"
#include "auxiliarFunctions.h"

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


        if ( getKeyword(strInput) == "Time" )
        {
            envir.readTimeStepAndMaxTime( getData(strInput) );
            continue;
        }

        char cc;
        std::cin >> cc;
    }
}







void IO::print2outFile(const std::string &outFlNm)
{}





