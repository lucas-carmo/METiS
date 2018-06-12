#include "ENVIR.h"
#include "IO.h"

#include <iostream>
#include <vector>

void ENVIR::readTimeStep(const std::string &data)
{
    readDataFromString(data, m_timeStep);
    std::cout << "Time Step: " << m_timeStep << "\n";
}



void ENVIR::readTimeTotal(const std::string &data)
{
    readDataFromString(data, m_timeTotal);
    std::cout << "Total Time: " << m_timeTotal << "\n";
}



void ENVIR::readTimeRamp(const std::string &data)
{
    readDataFromString(data, m_timeRamp);
    std::cout << "Time Ramp: " << m_timeRamp << "\n";
}



void ENVIR::readGrav(const std::string &data)
{
	readDataFromString(data, m_watDens);
	std::cout << "Gravity: " << m_watDens << "\n";
}



void ENVIR::readWatDens(const std::string &data)
{
	readDataFromString(data, m_gravity);
	std::cout << "Water Density: " << m_gravity << "\n";
}



void ENVIR::readWatDepth(const std::string &data)
{
	readDataFromString(data, m_watDepth);
	std::cout << "Water Depth: " << m_watDepth << "\n";
}




