#include "FOWT.h"
#include "IO.h"

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
    
    std::cout << "Linear Stiffness: " << m_linStiff << "\n";
}