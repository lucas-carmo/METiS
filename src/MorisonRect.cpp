#include "MorisonRect.h"

using namespace arma;


/*****************************************************
    Constructors
*****************************************************/
MorisonRect::MorisonRect(vec cog2node1, vec cog2node2, vec cog2node3, int numIntPoints, bool botPressFlag,
			        	double axialCD, double axialCa, double diam_X, double CD_X, double CM_X,
				        double diam_Y, double CD_Y, double CM_Y,
				        double botArea, double topArea)
            : MorisonElement(cog2node1, cog2node2, numIntPoints, botPressFlag, axialCD, axialCa),
              m_cog2node3(cog2node3), m_diam_X(diam_X), m_CD_X(CD_X), m_CM_X(CM_X),
              m_diam_Y(diam_Y), m_CD_Y(CD_Y), m_CM_Y(CM_Y), m_botArea(botArea), m_topArea(topArea)
{}


/*****************************************************
	Forces acting on the Morison Element
*****************************************************/
vec::fixed<3> MorisonRect::node1Pos(const vec::fixed<6> &floaterPos) const
{	
	// Fazer uma funcao que calcula a matriz de rotacao
	// mat::fixed<3,3> rotatMatrix(const vec::fixed<6> &FOWTpos) const
	double x = floaterPos.at(0) + m_cog2node1.at(0);
	double y = floaterPos.at(1) + m_cog2node1.at(1);
	double z = floaterPos.at(2) + m_cog2node1.at(2);

	vec::fixed<3> pos = { x, y, z };
	return pos;
}

vec::fixed<3> MorisonRect::node2Pos(const vec::fixed<6> &floaterPos) const
{	
	// Fazer uma funcao que calcula a matriz de rotacao
	// mat::fixed<3,3> rotatMatrix(const vec::fixed<6> &FOWTpos) const
	double x = floaterPos.at(0) + m_cog2node2.at(0);
	double y = floaterPos.at(1) + m_cog2node2.at(1);
	double z = floaterPos.at(2) + m_cog2node2.at(2);

	vec::fixed<3> pos = { x, y, z };
	return pos;
}

vec::fixed<3> MorisonRect::node1Vel(const vec::fixed<6> &floaterPos, const vec::fixed<6> &floaterVel) const
{
	double x = 0;
	double y = 0;
	double z = 0;

	vec::fixed<3> vel = { x, y, z };
	return vel;
}

vec::fixed<3> MorisonRect::node2Vel(const vec::fixed<6> &floaterPos, const vec::fixed<6> &floaterVel) const
{
	double x = 0;
	double y = 0;
	double z = 0;

	vec::fixed<3> vel = { x, y, z };
	return vel;
}

vec::fixed<3> MorisonRect::node1Acc(const vec::fixed<6> &floaterPos, const vec::fixed<6> &floaterVel, const vec::fixed<6> &floaterAcc) const
{
	double x = 0;
	double y = 0;
	double z = 0;

	vec::fixed<3> acc = { x, y, z };
	return acc;
}

vec::fixed<3> MorisonRect::node2Acc(const vec::fixed<6> &floaterPos, const vec::fixed<6> &floaterVel, const vec::fixed<6> &floaterAcc) const
{
	double x = 0;
	double y = 0;
	double z = 0;

	vec::fixed<3> acc = { x, y, z };
	return acc;
}





vec::fixed<6> MorisonRect::hydrostaticForce(const ENVIR &envir, const vec::fixed<6> &floaterPos) const
{
	vec::fixed<6> force(fill::zeros);
	return force;
}

vec::fixed<6> MorisonRect::hydrodynamicForce(const ENVIR &envir, const vec::fixed<6> &floaterPos, const vec::fixed<6> &floaterVel, const vec::fixed<6> &floaterAcc) const
{
	vec::fixed<6> force(fill::zeros);
	return force;
}


/*****************************************************
	Printing
*****************************************************/
std::string MorisonRect::print() const
{
	std::string output = "";

	output = output + "CoG_2_node1:\t(" + std::to_string(m_cog2node1(0)) + ", " + std::to_string(m_cog2node1(1)) + ", " + std::to_string(m_cog2node1(2)) + ")\n";
	output = output + "CoG_2_node2:\t(" + std::to_string(m_cog2node2(0)) + ", " + std::to_string(m_cog2node2(1)) + ", " + std::to_string(m_cog2node2(2)) + ")\n";
	output = output + "CoG_2_node3:\t(" + std::to_string(m_cog2node3(0)) + ", " + std::to_string(m_cog2node3(1)) + ", " + std::to_string(m_cog2node3(2)) + ")\n";
	output = output + "Diameter X:\t" + std::to_string(m_diam_X) + '\n';
	output = output + "Drag Coeff. X:\t" + std::to_string(m_CD_X) + '\n';
	output = output + "Inert. Coeff. X:\t" + std::to_string(m_CM_X) + '\n';
	output = output + "Diameter Y:\t" + std::to_string(m_diam_Y) + '\n';
	output = output + "Drag Coeff. Y:\t" + std::to_string(m_CD_Y) + '\n';
	output = output + "Inert. Coeff. Y:\t" + std::to_string(m_CM_Y) + '\n';
	output = output + "Numb. of Int. Points:\t" + std::to_string(m_numIntPoints) + '\n';
	output = output + "Bottom diameter:\t" + std::to_string(m_botArea) + '\n';
	output = output + "Top diameter:\t" + std::to_string(m_topArea) + '\n';
	output = output + "Axial CD:\t" + std::to_string(m_axialCD) + '\n';
	output = output + "Axial Ca:\t" + std::to_string(m_axialCa) + '\n';
	output = output + "Bot. Press. Flag.:\t" + std::to_string(m_botPressFlag) + '\n';

	return output;
}


