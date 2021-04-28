#include "MorisonRect.h"

using namespace arma;


/*****************************************************
    Constructors
*****************************************************/
MorisonRect::MorisonRect(const vec &node1Pos, const vec &node2Pos, const vec &node3Pos, const vec &cog, const int numIntPoints,
						 const bool botPressFlag, const double axialCD_1, const double axialCa_1, const double axialCD_2, const double axialCa_2,
						 const double diam_X, const double CD_X, const double CM_X,
						 const double diam_Y, const double CD_Y, const double CM_Y)
            : MorisonElement(node1Pos, node2Pos, cog, numIntPoints, botPressFlag, axialCD_1, axialCa_1, axialCD_2, axialCa_2),
              m_diam_X(diam_X), m_CD_X(CD_X), m_CM_X(CM_X),
              m_diam_Y(diam_Y), m_CD_Y(CD_Y), m_CM_Y(CM_Y)
{
	m_cog2node3 = node3Pos - cog;
	make_local_base(m_xvec, m_yvec_t0, m_zvec_t0, node1Pos, node2Pos);
}

void MorisonRect::evaluateQuantitiesAtBegin(const ENVIR & envir, const int hydroMode)
{
}


/*****************************************************
	Forces acting on the Morison Element
*****************************************************/
void MorisonRect::make_local_base(arma::vec::fixed<3> &xvec, arma::vec::fixed<3> &yvec, arma::vec::fixed<3> &zvec, const arma::vec::fixed<3> &n1, const arma::vec::fixed<3> &n2) const
{
}



// TODO: depois de debugar direitinho, tirar os bound checks (usar [] ao inves de () pra acessar elementos das matrizes)
mat::fixed<6, 6> MorisonRect::addedMass_perp(const double rho, const vec::fixed<3> &refP, const int hydroModet) const
{
	mat::fixed<6, 6> A(fill::zeros);

	return A;
}

double MorisonRect::A_perp(const int ii, const int jj, const vec::fixed<3> &x, const vec::fixed<3> &xG, const vec::fixed<3> &xvec, const vec::fixed<3> &yvec) const
{
	return 0;
}

mat::fixed<6, 6> MorisonRect::addedMass_paral(const double rho, const vec::fixed<3> &refPt, const int hydroMode) const
{
	mat::fixed<6, 6> A(fill::zeros);

	return A;
}

double MorisonRect::A_paral(const int ii, const int jj, const vec::fixed<3> &x, const vec::fixed<3> &xG, const vec::fixed<3> &zvec) const
{
	return 0;
}

vec::fixed<6> MorisonRect::hydrostaticForce(const double rho, const double g) const
{
	vec::fixed<6> force(fill::zeros);
	return force;
}

vec::fixed<6> MorisonRect::hydroForce_1st(const ENVIR & envir, const int hydroMode, const vec::fixed<3>& refPt) const
{
	return fill::zeros;
}

vec::fixed<6> MorisonRect::hydroForce_drag(const ENVIR & envir, const vec::fixed<3>& refPt) const
{
	return fill::zeros;
}

vec::fixed<6> MorisonRect::hydroForce_relWaveElev(const ENVIR & envir, const vec::fixed<3>& refPt) const
{
	return fill::zeros;
}

vec::fixed<6> MorisonRect::hydroForce_2ndPot(const ENVIR & envir, const vec::fixed<3>& refPt) const
{
	return fill::zeros;
}

vec::fixed<6> MorisonRect::hydroForce_convecAcc(const ENVIR & envir, const vec::fixed<3>& refPt) const
{
	return fill::zeros;
}

vec::fixed<6> MorisonRect::hydroForce_axDiverg(const ENVIR & envir, const vec::fixed<3>& refPt) const
{
	return fill::zeros;
}

vec::fixed<6> MorisonRect::hydroForce_accGradient(const ENVIR & envir, const vec::fixed<3>& refPt) const
{
	return fill::zeros;
}

vec::fixed<6> MorisonRect::hydroForce_slendBodyRot(const ENVIR & envir, const vec::fixed<3>& refPt) const
{
	return fill::zeros;
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
	output = output + "Axial CD - Node 1:\t" + std::to_string(m_axialCD_1) + '\n';
	output = output + "Axial Ca - Node 1:\t" + std::to_string(m_axialCa_1) + '\n';
	output = output + "Axial CD - Node 2:\t" + std::to_string(m_axialCD_2) + '\n';
	output = output + "Axial Ca - Node 2:\t" + std::to_string(m_axialCa_2) + '\n';
	output = output + "Bot. Press. Flag.:\t" + std::to_string(m_botPressFlag) + '\n';

	return output;
}


/*****************************************************
	Clone for creating copies of the Morison Element
*****************************************************/
MorisonRect* MorisonRect::clone() const
{
	return (new MorisonRect(*this));
}

