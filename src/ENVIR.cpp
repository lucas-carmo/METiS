#include "ENVIR.h"
#include "IO.h"

#include <iostream>
#include <vector>


/*****************************************************
	Setters
*****************************************************/
void ENVIR::readTimeStep(const std::string &data)
{
    readDataFromString(data, m_timeStep);
}



void ENVIR::readTimeTotal(const std::string &data)
{
    readDataFromString(data, m_timeTotal);
}



void ENVIR::readTimeRamp(const std::string &data)
{
    readDataFromString(data, m_timeRamp);
}



void ENVIR::readGrav(const std::string &data)
{
	readDataFromString(data, m_watDens);
}



void ENVIR::readWatDens(const std::string &data)
{
	readDataFromString(data, m_gravity);
}



void ENVIR::readWatDepth(const std::string &data)
{
	readDataFromString(data, m_watDepth);
}


void ENVIR::readWave(const std::string &data)
{

}



/*****************************************************
	Getters
*****************************************************/
std::string ENVIR::printTimeStep() const
{
	return std::to_string(m_timeStep);
}

std::string ENVIR::printTimeTotal() const
{
	return std::to_string(m_timeTotal);
}

std::string ENVIR::printTimeRamp() const
{
	return std::to_string(m_timeRamp);
}

std::string ENVIR::printGrav() const
{
	return std::to_string(m_gravity);
}

std::string ENVIR::printWatDens() const
{
	return std::to_string(m_watDens);
}

std::string ENVIR::printWatDepth() const
{
	return std::to_string(m_watDepth);
}

//std::string ENVIR::printWave() const
//{}