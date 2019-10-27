#include "Wave.h"
#include "IO.h"
#include <vector>

using namespace arma;



/*****************************************************
	Constructors
*****************************************************/
Wave::Wave()
{
	m_height = arma::datum::nan;
	m_period = arma::datum::nan;
	m_direction = arma::datum::nan;
	m_phase = arma::datum::nan;
	m_waveNumber = arma::datum::nan;
	m_length = arma::datum::nan;
}

// The string "waveType" must contain the keyword that specifies how to read the third input:
// 1) TRWave for regular waves specified by their period
// 2) FRWavefor regular waves specified by their frequency
// 3) WRWave for regular waves specified by their angular frequency
Wave::Wave(const std::string &waveType, const double height, const double freqORperiod, const double direction, const double phase, const double watDepth, const double gravity)
{
	// Neither the wave frequency nor the period can be zero
	if (freqORperiod == 0)
	{
		throw std::runtime_error("Wave period or frequency must be different from zero. Input line " + std::to_string(IO::getInLineNumber()));
	}

	double period;
	if (caseInsCompare(waveType, "TRWave"))
	{
		period = freqORperiod;
	}
	else if (caseInsCompare(waveType, "FRWave"))
	{
		period = 1 / freqORperiod;
	}
	else if (caseInsCompare(waveType, "WRWave"))
	{
		period = 2 * datum::pi / freqORperiod;	
	}
	else
	{
		throw std::runtime_error("Unknown wave type" + waveType + " provided to Wave constructor.");
	}

	m_height = height;
	m_period = period;
	m_direction = direction;
	m_phase = phase - std::floor(phase / 360) * 360; // Make sure phase is between 0 and 2pi

	m_waveNumber = waveNumber(watDepth, gravity);
	m_length = 2 * arma::datum::pi / m_waveNumber; // The function Wave::waveNumber never outputs 0, hence this division is safe
}

/*****************************************************
	Getters
*****************************************************/
double Wave::height() const
{
	return m_height;
}

double Wave::amp() const
{
	return m_height/2;
}

double Wave::period() const
{
	return m_period;
}


double Wave::direction() const
{
	return m_direction;
}

double Wave::phase() const
{
	return m_phase;
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


double Wave::waveNumber() const
{
	if (is_finite(m_waveNumber))
	{
		return m_waveNumber;
	}	
	else
	{
		throw std::runtime_error("Tried to call Wave::waveNumber(), but wave number was not calculated yet. Try calling Wave::waveNumber(const double watDepth, const double gravity).");	
	}
}


double Wave::length() const
{
	if (is_finite(m_length))
	{
		return m_length;
	}	
	else
	{
		throw std::runtime_error("Tried to call Wave::length(), but length was not calculated yet. Try calling Wave::length(const double watDepth, const double gravity).");	
	}
}


double Wave::waveNumber(const double watDepth, const double gravity) const
{
	if (m_period == 0)
	{
		return arma::datum::inf;
	}

	// Use a more friendly notation
	double w = angFreq();
	double h = watDepth;
	double g = gravity;

	// Use the bissection method to solve the dispersion relation.
	// We want to find the root of the following equation:
	// f(x) = w^2 / g - x * tanh(h*x)
	//
	// First, we need to bracket the solution between two points a and b
	// in such a way that f(a) and f(b) have different signals.
	//
	// Since tanh(h*x) <= 1, it is easy to see that x = w^2 / g results in f(x) > 0.
	// Hence, we take a = w^2 / g and b = alpha * w^2 / g.
	// If f(b) < 0, we start the bissection method. Otherwise, we take a = b0
	// and b = alpha*b0, and so on until f(a) > 0 and f(b) < 0
	double alpha = 1.25; // This value is completely arbitrary. I should improve this whole method of calculating the wave number.
	double a = pow(w, 2) / g;
	double b = alpha*a;

	// If f(a) * f(b) < 0, it means they have different signals.
	// Otherwise, we update a and b until the solution is bracketed
	while ( (pow(w,2)/g - a*tanh(a*h)) * (pow(w,2)/g - b*tanh(b*h)) > 0 )
	{
		a = b;
		b = alpha*b;		
	}

	// x_i is the guess in the previous step and x_j in the current
	double x_i = (a+b)/2;
	double x_j(0);

	while ( std::abs(pow(w, 2) / g - x_j * tanh(x_j*h)) > m_epsWave)
	{								
		x_i = x_j;

		if ( (pow(w,2)/g - a*tanh(a*h)) * (pow(w,2)/g - x_i*tanh(x_i*h)) < 0 ) // Test with limit a
		{
			x_j = (a + x_i) / 2;
			b = x_i;
		}

		else if ( (pow(w,2)/g - b*tanh(b*h)) * (pow(w,2)/g - x_i*tanh(x_i*h)) < 0 ) // Test with limit b
		{
			x_j = (x_i + b) / 2;
			a = x_i;
		}

		else // It is very unlikely that x_i will be zero, but this possibility will be covered anyway
		{
			return x_i;
		}	
	}		

	if (x_j == 0)
	{
		throw std::runtime_error("Wave with wave number k = 0. Wave period = " + std::to_string(m_period));
	}

	return x_j;
}


double Wave::length(const double watDepth, const double gravity) const
{
	// The function Wave::waveNumber never outputs 0	
	return 2*arma::datum::pi / this->waveNumber(watDepth, gravity);
}