#include "Wave.h"
#include "IO.h"
#include <vector>


/*****************************************************
	Setters
*****************************************************/

// READ REGULAR WAVES
// The string "wholeWaveLine" must contain the keyword that specifies how to read the wave data:
// 1) TRWave for regular waves specified by their period
// 2) FRWavefor regular waves specified by their frequency
// 3) WRWave for regular waves specified by their angular frequency
void Wave::readRegWave(const std::string &wholeWaveLine)
{	
	// Wave characteristics are divided by a space or a tab.
	// The characteristics are 
	// 1) Wave height
	// 2) Wave period OR wave frequency OR wave angular frequency
	// 3) Direction of propagation	
	std::vector<std::string> input = stringTokenize(getData(wholeWaveLine), " \t");

	// Test if the keyword is known
	if ( !caseInsCompare(getKeyword(wholeWaveLine), "TRWave") 
	  && !caseInsCompare(getKeyword(wholeWaveLine), "FRWave") 
	  && !caseInsCompare(getKeyword(wholeWaveLine), "WRWave") )
	{
		std::cout << "Error in line " << IO::getInLineNumber() 
		          << ": Unknown keyword '" << getKeyword(wholeWaveLine) << "' \n";
		return;
	}


	// Check if there are exactly three inputs (Wave height, period/frequency/angular frequency, and direction)
	if (input.size() != 3)
	{
		std::cout << "Error in line" << IO::getInLineNumber() 
		          << ": TRWave, FRWave and WRWave should be followed by three numbers"
				  << " separated by tabs ('\\t') or spaces specifying"
				  << " the wave Height, Period/Frequency/Angular Frequency, and Direction." 
				  << " However, you provided " << input.size() << " elements \n";
				  
		return;
	}


	// Finally, read wave data
	readDataFromString(input.at(0), m_height);
	readDataFromString(input.at(2), m_direction);

	if ( caseInsCompare(getKeyword(wholeWaveLine), "TRWave") )
	{
		readDataFromString(input.at(1), m_period);	
	}

	if ( caseInsCompare(getKeyword(wholeWaveLine), "FRWave") )
	{	
		double frequency{0};
		readDataFromString(input.at(1), frequency);	

		if (frequency == 0)
		{
			std::cout << "Wave frequency must be different from zero.\n";
		}
		else 
		{
			m_period = 1/frequency;
		}
	}	

	if ( caseInsCompare(getKeyword(wholeWaveLine), "WRWave") )
	{	
		double omega{0};
		readDataFromString(input.at(1), omega);	

		if (omega == 0)
		{
			std::cout << "Wave frequency must be different from zero.\n";
		}
		else
		{
			m_period = 2*M_PI/omega;
		}
	}		
}


/*****************************************************
	Getters
*****************************************************/
double Wave::height() const
{
	return m_height;
}

double Wave::period() const
{
	return m_period;
}

double Wave::direction() const
{
	return m_direction;
}



