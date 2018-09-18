#include "Wave.h"
#include "IO.h"
#include <vector>


using namespace arma;

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
		throw std::runtime_error( "Unknown keyword '" + getKeyword(wholeWaveLine) + "' in input line " + std::to_string(IO::getInLineNumber()) + ".");
		return;
	}


	// Check if there are exactly three inputs (Wave height, period/frequency/angular frequency, and direction)
	if (input.size() != 3)
	{
		throw std::runtime_error("Unable to read the wave in input line " + std::to_string(IO::getInLineNumber()) + ". Wrong number of parameters.");				  
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
			throw std::runtime_error( "Wave frequency must be different from zero. Input line " + std::to_string(IO::getInLineNumber()) );
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
			throw std::runtime_error( "Wave frequency must be different from zero. Input line " + std::to_string(IO::getInLineNumber()) );
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

double Wave::freq() const
{
	if (m_period != 0)
	{
		return 1/m_period;
	}
	else
	{
		return arma::datum::inf;
	}
}

double Wave::angFreq() const
{
	if (m_period != 0)
	{
		return 2*arma::datum::pi/m_period;
	}
	else
	{
		return arma::datum::inf;
	}
}

double Wave::waveNumber(ENVIR &envir) const
{
	double k(0);

	
	return k;
}

double Wave::waveLength(ENVIR &envir) const
{

}

/*****************************************************
	Other functions
*****************************************************/
vec::fixed<3> Wave::fluidVel(ENVIR &envir, vec &point) const
{
	vec::fixed<3> fluidVel = zeros<vec>(0);
	
	// Use a more friendly notation
	double x(point(0));
	double y(point(1));
	double z(point(2));

	// If the point is above the free surface, 
	if (z>0)
	{
		return fluidVel;
	}

	return fluidVel;
}


// vec::fixed<3> fluidAcc(vec &point) const;