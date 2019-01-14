#pragma once

#include <armadillo>
#include <vector>

using namespace arma; // For armadillo classes

class RNA
{
private:
    double rotorSpeed;
    double tilt;
    double yaw;    

    int numBlades;
    // Blade blade;
    std::vector< mat<double> > airfoils;

    double hubRadius;
    double hubHeight;
    double overhang;

public:
};

