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