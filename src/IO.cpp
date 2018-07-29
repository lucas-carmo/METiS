#include "IO.h"

#include <fstream> // Include file input/output classes
#include <iostream>
#include <stdexcept> // For std::exception
#include <sys/types.h> // Using this and the one below to check if directories and files exist
#include <sys/stat.h>
#include <cstdlib> // for exit()



/*****************************************************
	Defining and initializing static member variables
*****************************************************/
std::string IO::m_inFilePath = "";
std::ifstream IO::m_inFl;
unsigned int IO::m_inLineNumber = 0;

std::string IO::m_logFilePath = "";
std::ofstream IO::m_logFl;



/*****************************************************
    Class functions related to Input
*****************************************************/
void IO::setFiles(const std::string &inFlPath)
{
	// Check whether we are not trying to reset the input file
	if ( !m_inFilePath.empty() )
	{
		throw std::runtime_error("You can not reset the input file. Press enter to exit.\n");
	}
	m_inFilePath = inFlPath;

	// Open input file
	m_inFl.open(m_inFilePath);
	if (!m_inFl)
	{
		throw std::runtime_error("Unable to open file " + m_inFilePath + " for reading.");
	}	


	// Get path of the folder where the input file is located, as that is where the output
	// will be saved to.
	std::string folderPath = getFileFolder(m_inFilePath);

	// Path of the output folder
	std::string outputFolder = folderPath + "/output";

	// Create output folder with name "output".
	// If a folder (or file) named "output" already exists in folderPath, create a folder named output_1. If the latter exists as well, create output_2 instead.
	// Keep this process until output_n does not exist and is then created.
	struct stat info;		
	std::string outputFolder_original = outputFolder; // Original name of folderPath
	int ii = 1;	
	while ( stat(outputFolder.c_str(), &info) == 0 )
	{
		outputFolder = outputFolder_original + "_" + std::to_string(ii);			
		++ii;
	}
	system( ("mkdir " + outputFolder).c_str() );

	// Set log file
	// If we fail to open the log file, print a message to the console screen and terminate the program.
	// The other error messages will be written to the log file (as well as any other message that is fit to be output).
	m_logFilePath = outputFolder + "/log.txt";
	m_logFl.open(m_logFilePath); 
	if (!m_logFl)
	{
		throw std::runtime_error("Unable to open file " + m_logFilePath + " for writting. Press enter to exit.\n");
	}
}


// Read line from input file to string "strInput".
// The function deals with empty lines and comments using functions "hasContent" and
// "thereIsCommentInString". Besides, it updates the line number counter inLineNumber
void IO::readLineInputFile(std::string &strInput)
{
	std::getline(m_inFl, strInput); // Read next file line to string strInput
	++m_inLineNumber; // Update line number counter

	// Repeat this process until the line has some content or end of file is achieved
	while (!hasContent(strInput) && m_inFl)
	{		
		std::getline(m_inFl, strInput);
		++m_inLineNumber;
	}

	// Remove comments from line
	if (thereIsCommentInString(strInput))
	{
		removeComments(strInput);
	}				
}

unsigned int IO::getInLineNumber()
{
	return m_inLineNumber;
}


//	This is the main input function.
//	It is responsible for reading the input file line by line and assigning what is read to the FOWT and ENVIR classes
void IO::readInputFile(FOWT &fowt, ENVIR &envir)
{
	// Classes that are members of FOWT and ENVIR
	Floater floater;		

	// Read file line by line
	while (m_inFl)
	{
		std::string strInput;		
		IO::readLineInputFile(strInput);

		/**************************
		Read data based on keywords
		**************************/

		// Read data to envir
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



		// Read data to fowt
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


		if (caseInsCompare(getKeyword(strInput), "Wave")) // A list of Waves is supposed to follow the "Wave keyword"
		{			
			IO::readLineInputFile(strInput); // Read next line, since current line is just the keyword Wave

			while (!caseInsCompare(getKeyword(strInput), "END")) // The END keyword indicates the end of the list of waves
			{							
				if (!m_inFl) // Signal if the end of file is reached before the end keyword
				{
					std::cout << "\n\n\n End of file reached before END keyword in WAVE specification.\n\n";
					return;
				}

				envir.addWave( Wave(strInput) ); // Add this wave to the environment
				
				IO::readLineInputFile(strInput);
			}			
		}


		if (caseInsCompare(getKeyword(strInput), "Nodes")) // A list of Nodes is supposed to follow the "Wave keyword"
		{
			IO::readLineInputFile(strInput); // Read next line, since current line is just the keyword Nodes

			while (!caseInsCompare(getKeyword(strInput), "END"))
			{
				if (!m_inFl) // Signal if the end of file is reached before the end keyword
				{
					std::cout << "\n\n\n End of file reached before END keyword in NODES specification.\n\n";
					return;
				}

				floater.addNode(strInput); // Add this node to the floater

				IO::readLineInputFile(strInput);
			}
		}


		if (caseInsCompare(getKeyword(strInput), "Morison_circ")) // A list of circular cylinder Morison Elements is supposed to follow the "Morison_circ" keyword
		{
			IO::readLineInputFile(strInput); // Read next line, since current line is just the keyword "Morison_circ"

			while (!caseInsCompare(getKeyword(strInput), "END"))
			{
				if (!m_inFl) // Signal if the end of file is reached before the end keyword
				{
					std::cout << "\n\n\n End of file reached before END keyword in MORISON_CIRc specification.\n\n";
					return;
				}

				floater.addMorisonCirc(strInput); // Add this Morison Element to the floater

				IO::readLineInputFile(strInput);
			}
		}


		if (caseInsCompare(getKeyword(strInput), "Morison_rect")) // A list of rectangular cylinder Morison Elements is supposed to follow the "Morison_circ" keyword
		{
			IO::readLineInputFile(strInput); // Read next line, since current line is just the keyword "Morison_circ"

			while (!caseInsCompare(getKeyword(strInput), "END"))
			{
				if (!m_inFl) // Signal if the end of file is reached before the end keyword
				{
					std::cout << "\n\n\n End of file reached before END keyword in MORISON_RECT specification.\n\n";
					return;
				}

				floater.addMorisonRect(strInput); // Add this Morison Element to the floater

				IO::readLineInputFile(strInput);
			}
		}
	}

	fowt.setFloater(floater);
}


/*****************************************************
    Class functions related to Output
*****************************************************/
void IO::writeErrorMessage(const std::string &strInput)
{
	if (m_logFl) // If we are able to write to the log file, do it
	{
		m_logFl << strInput << std::endl;
	}

	// Write to the console, anyway
	std::cerr << strInput << std::endl;
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


// Get folder path from a complete file path
std::string getFileFolder(const std::string& path)
{
  size_t found;
  found = path.find_last_of("/\\");
  return path.substr(0,found);
}