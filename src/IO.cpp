#include "IO.h"

#include <fstream> // Include file input/output classes
#include <iostream>
#include <cstdlib> // For exit()



/*****************************************************
	Defining and initializing static member variables
*****************************************************/
std::string IO::m_inFilePath = "";
std::ifstream IO::m_inFl;
unsigned int IO::m_inLineNumber = 0;



/*****************************************************
    Implementation of class functions
*****************************************************/
void IO::setInputFile(const std::string &inFlPath)
{
	m_inFilePath = inFlPath;
	m_inFl.open(m_inFilePath);

	if (!m_inFl)
	{
		std::cerr << "Unable to open " << m_inFilePath << " for reading";
		return;
		//exit(1);
	}
}

void IO::readLineInputFile(std::string &strInput)
{
	std::getline(m_inFl, strInput);
	++m_inLineNumber; // Update counter for lineNumber
}

unsigned int IO::getInLineNumber()
{
	return m_inLineNumber;
}

// This is the main input function.
// It is responsible for reading the input file line by line and assigning what is read to the FOWT and ENVIR classes
void IO::readInputFile(FOWT &fowt, ENVIR &envir)
{
	//std::ifstream inFl(m_inFilePath);

	//if (!m_inFl)
	//{
	//	std::cerr << "Unable to open " << m_inFilePath << " for reading";
	//	return;
	//	//exit(1);
	//}


	// Classes that are members of FOWT and ENVIR
	Floater floater;		

	// Read file line by line
	while (m_inFl)
	{
		std::string strInput;		
		IO::readLineInputFile(strInput);

		if (!hasContent(strInput))
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
		if (caseInsCompare(getKeyword(strInput), "TimeStep"))
		{
			envir.readTimeStep(getData(strInput));
			continue;
		}

		if (caseInsCompare(getKeyword(strInput), "TimeTotal"))
		{
			envir.readTimeTotal(getData(strInput));
			continue;
		}

		if (caseInsCompare(getKeyword(strInput), "TimeRamp"))
		{
			envir.readTimeRamp(getData(strInput));
			continue;
		}

		if (caseInsCompare(getKeyword(strInput), "Grav"))
		{
			envir.readGrav(getData(strInput));
			continue;
		}

		if (caseInsCompare(getKeyword(strInput), "WatDens"))
		{
			envir.readWatDens(getData(strInput));
			continue;
		}

		if (caseInsCompare(getKeyword(strInput), "WatDepth"))
		{
			envir.readWatDepth(getData(strInput));
			continue;
		}



		// Read data to envir class
		if (caseInsCompare(getKeyword(strInput), "LinStiff"))
		{
			fowt.readLinStiff(getData(strInput));
			continue;
		}

		if (caseInsCompare(getKeyword(strInput), "FloaterMass"))
		{
			floater.readMass(getData(strInput));
			continue;
		}

		if (caseInsCompare(getKeyword(strInput), "FloaterCoG"))
		{
			floater.readCoG(getData(strInput));
			continue;
		}


		if (caseInsCompare(getKeyword(strInput), "Wave"))
		{
			// A list of Waves is supposed to follow the "Wave keyword"
			// 1) While the END keyword is not reached
			//Outro std::getline
			// 2) getKeyword == algum dos tipos de onda
			// 		2.1) Initializa uma "Wave wave" e passa os parâmetros dependendo da keyword
			// 3) getKeyword == END
			//		3.1) Termina
			// 4) getKeyword == Qualquer outra coisa
			//		4.1) Erro					
			//fowt.readFloaterCoG(getData(strInput));			
			IO::readLineInputFile(strInput);

			while (!caseInsCompare(getKeyword(strInput), "END") && m_inFl) 		
			{							
				if (!hasContent(strInput))
				{
					IO::readLineInputFile(strInput);
					continue;
				}

				if (thereIsCommentInString(strInput))
				{
					removeComments(strInput);
				}			

				Wave wave;	
				wave.readRegWave(strInput);	
				envir.addWave(wave); // Add this wave to the environment
				
				IO::readLineInputFile(strInput);
			}			
		}
	}

	fowt.setFloater(floater);

	IO::print2CheckVariables(fowt, envir);
}

void IO::print2outFile(const std::string &outFlNm)
{}


// Print the members of fowt and envir. Useful for debugging.
void IO::print2CheckVariables(const FOWT &fowt, const ENVIR &envir)
{
	std::cout << "FOWT:\n";
	std::cout << "Linear Stiffness:\t" << fowt.printLinStiff() << '\n';
	std::cout << "Floater:\n" << fowt.printFloater();

	std::cout << "\n\n\nENVIR:\n";
	std::cout << "Time Step:\t" << envir.printTimeStep() << '\n';
	std::cout << "Total Time:\t" << envir.printTimeTotal() << '\n';
	std::cout << "Time Ramp:\t" << envir.printTimeRamp() << '\n';
	std::cout << "Gravity:\t" << envir.printGrav() << '\n';
	std::cout << "Water Density:\t" << envir.printWatDens() << '\n';
	std::cout << "Water Depth:\t" << envir.printWatDepth() << '\n';
	std::cout << "\n" << envir.printWave() << '\n';
}



/*****************************************************
    Additional functions related to input/output 
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