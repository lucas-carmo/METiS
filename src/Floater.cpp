#include "Floater.h"

#include "IO.h"
#include <iostream>
#include <vector>


/*****************************************************
	Setters
*****************************************************/
void Floater::readMass(const std::string &data)
{
	readDataFromString(data, m_mass);
}

void Floater::readCoG(const std::string &data)
{
	// The coordinates of the center of gravity are separated by commas in the input string
	std::vector<std::string> input = stringTokenize(data, ",");

	if (input.size() != 3)
	{
		std::cout << "Deu ruim \n";
	}

	for (int ii = 0; ii < input.size(); ++ii)
	{
		readDataFromString(input.at(ii), m_CoG(ii));
	}
}


/*****************************************************
	Getters
*****************************************************/

/*std::string output = "";
for (int ii = 0; ii < m_linStiff.n_elem; ++ii)
{
	output = output + " " + std::to_string(m_linStiff(ii));
}
return output;*/


std::string Floater::printMass() const
{
	return std::to_string(m_mass);
}

std::string Floater::printCoG() const
{
	std::string output = "(" + std::to_string( m_CoG(0) );
	for (int ii = 1; ii < m_CoG.n_elem; ++ii)
	{
		output = output + "," + std::to_string( m_CoG(ii) );
	}
	
	return output + ")";
}

//mat Floater::rotatMat(const vec &FOWTpos)
//{}
//
//
//vec Floater::hydrodynamicForce(const ENVIR &envir, const vec &FOWTpos, const vec &FOWTvel, const vec &FOWTacc)
//{}


