#include "FOWT.h"
#include "IO.h"


/*****************************************************
	Setters
*****************************************************/

void FOWT::readLinStiff(const std::string &data)
{
    // The mooring line stiffness in surge, sway and yaw are separated by commas in the input string
    std::vector<std::string> input = stringTokenize( data, "," );        
    
    if ( input.size() != 3 )
    {
        std::cout << "Deu ruim \n";		
	}
    
    for ( int ii = 0; ii < input.size(); ++ii )
    {
        readDataFromString( input.at(ii), m_linStiff(ii) );
    }    
}


void FOWT::setFloater(const Floater &floater)
{
	m_floater = floater;
}




/*****************************************************
	Getters
*****************************************************/
std::string FOWT::printLinStiff() const
{	
	std::string output = "(" + std::to_string( m_linStiff(0) );
	for ( int ii = 1; ii < m_linStiff.n_elem; ++ii )
	{
		output = output + "," + std::to_string( m_linStiff(ii) );
	}
	return output + ")";
}


std::string FOWT::printFloater() const
{
	std::string output = "";

	output = output + "\tMass:\t" + m_floater.printMass() + "\n";
	output = output + "\tCoG:\t" + m_floater.printCoG() + "\n";
	output = output + "\tNodes:\n" + m_floater.printNodes() + "\n";

	return output;
}