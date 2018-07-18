#include "MorisonCirc.h"


using namespace arma;

/*****************************************************
	Constructors
*****************************************************/
MorisonCirc::MorisonCirc(vec cog2node1, vec cog2node2, int numIntPoints, bool botPressFlag, double diam, double CD, double CM,
						 double botDiam, double topDiam, double axialCD, double axialCa)
	: MorisonElement(cog2node1, cog2node2, numIntPoints, botPressFlag), m_diam(diam), m_CD(CD), m_CM(CM), m_botDiam(botDiam),
					 m_topDiam(topDiam), m_axialCD(axialCD), m_axialCa(axialCa)
{}




/*****************************************************
	Forces acting on the Morison Element
*****************************************************/
vec::fixed<6> MorisonCirc::hydrostaticForce(const ENVIR &envir)
{
	vec::fixed<6> force(fill::zeros); 
	return force;
}

vec::fixed<6> MorisonCirc::hydrodynamicForce(const ENVIR &envir)
{
	vec::fixed<6> force(fill::zeros);
	return force;
}



/*****************************************************
	Printing
*****************************************************/
std::string MorisonCirc::print()
{
	std::string output = "";

	output = output + "CoG_2_node1:\t(" + std::to_string(m_cog2node1(0)) + ", " + std::to_string(m_cog2node1(1)) + ", " + std::to_string(m_cog2node1(2)) + ")\n";
	output = output + "CoG_2_node1:\t(" + std::to_string(m_cog2node2(0)) + ", " + std::to_string(m_cog2node2(1)) + ", " + std::to_string(m_cog2node2(2)) + ")\n";
	output = output + "Diameter:\t" + std::to_string(m_diam) + '\n';
	output = output + "Drag Coeff.:\t" + std::to_string(m_CD) + '\n';
	output = output + "Inert. Coeff.:\t" + std::to_string(m_CM) + '\n';
	output = output + "Numb. of Int. Points:\t" + std::to_string(m_numIntPoints) + '\n';
	output = output + "Bottom diameter:\t" + std::to_string(m_botDiam) + '\n';
	output = output + "Top diameter:\t" + std::to_string(m_topDiam) + '\n';
	output = output + "Axial CD:\t" + std::to_string(m_axialCD) + '\n';
	output = output + "Axial Ca:\t" + std::to_string(m_axialCa) + '\n';
	output = output + "Bot. Press. Flag.:\t" + std::to_string(m_botPressFlag) + '\n';



	return output;
}
