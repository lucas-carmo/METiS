#include "FOWT.h"
#include "IO.h"


/*****************************************************
	Setters used to set the data read from the input file.
*****************************************************/

void FOWT::readLinStiff(const std::string &data)
{
    // The mooring line stiffness in surge, sway and yaw are separated by commas in the input string
    std::vector<std::string> input = stringTokenize( data, "," );        
    
    if ( input.size() != 3 )
    {
        std::cout << "Deu ruim \n";		
	}
    
    for ( int ii = 0; ii < 3; ++ii)
    {
        readDataFromString( input.at(ii), m_linStiff(ii) );
    }    
    
    std::cout << "Linear Stiffness: \n" << m_linStiff << "\n";
}


/*
	Floater properties
*/
void FOWT::readFloaterMass(const std::string &data)
{
	m_floater.readMass(data);
}

void FOWT::readFloaterCoG(const std::string &data)
{
	m_floater.readCoG(data);
}