#pragma once

#include <vector>
#include <string>
#include "Wave.h"


class ENVIR
{
private:
    /*
    Flags to specify the type of analysis (hydrodynamic analysis, aerodynamics, etc)
    */
    bool m_analysisFlagHydro;
    bool m_analysisFlagAero;
    bool m_analysisFlagMoor;
    bool m_analysisFlagElast;

    /*
    Data to specify the environment
    */
    double m_gravity;
    double m_watDens;
    double m_watDepth;
    std::vector<Wave> m_wave;
    double m_airDens;
    double m_windMod;
    double m_windExp;

    /*
    Data to specify the numerical analysis
    */
    double m_timeStep;
    double m_timeTotal;
    double m_timeRamp;
    double m_currentTime;
    char m_timeIntegrationMethod; // 1: 4th order Runge-Kutta

    bool m_useBEMT;
    char m_aeroSolver;
    bool m_tipLossFactor;
    bool m_hubLossFactor;
    bool m_IncTIFac;
    bool m_IncDragAIFac;
    bool m_IncDragTIFac;
    bool m_IncSkewCorr;
    bool m_TwrLoads;

public:
    // Default constructor
    // ENVIR();

    /*
    Setters used to set the data read from the input file.
    */
    void readTimeStep(const std::string &data);
    void readTimeTotal(const std::string &data);
    void readTimeRamp(const std::string &data);
    void readLinStiff(const std::string &data);
    void readGrav(const std::string &data);
    void readWatDens(const std::string &data);
    void readWatDepth(const std::string &data);

};

