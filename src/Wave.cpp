#include "Wave.h"
#include "IO.h"
#include <vector>


/*****************************************************
	Constructors
*****************************************************/
Wave::Wave(double height, double period, double direction)
	: m_height(height), m_period(period), m_direction(direction)
{}

// The string "wholeWaveLine" must contain the keyword that specifies how to read the wave data:
// 1) TRWave for regular waves specified by their period
// 2) FRWavefor regular waves specified by their frequency
// 3) WRWave for regular waves specified by their angular frequency
Wave::Wave(const std::string &wholeWaveLine)
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
		throw std::runtime_error( "Unknown keyword '" + getKeyword(wholeWaveLine) + "' in line " + std::to_string(IO::getInLineNumber()) + ".");
		return;
	}


	// Check if there are exactly three inputs (Wave height, period/frequency/angular frequency, and direction)
	if (input.size() != 3)
	{
		throw std::runtime_error("Unable to read the wave in line " + std::to_string(IO::getInLineNumber()) + ". Wrong number of parameters.");				  
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
			throw std::runtime_error( "Wave frequency must be different from zero. Line " + std::to_string(IO::getInLineNumber()) );
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
			throw std::runtime_error( "Wave frequency must be different from zero. Line " + std::to_string(IO::getInLineNumber()) );
		}
		else
		{
			m_period = 2*datum::pi/omega;
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



