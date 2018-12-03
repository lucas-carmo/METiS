#include "FOWT.h"
#include "IO.h"


/*****************************************************
	Constructors
*****************************************************/
FOWT::FOWT() : m_linStiff(fill::zeros), m_mass(datum::nan), 
			   m_pos(fill::zeros), m_vel(fill::zeros), m_acc(fill::zeros)
{
	m_CoG.fill(datum::nan);	
}

/*****************************************************
	Overloaded operators
*****************************************************/
FOWT& FOWT::operator= (const FOWT &fowt)
{
	m_floater = fowt.m_floater;
	m_linStiff = fowt.m_linStiff;
	m_pos = fowt.m_pos;
	m_vel = fowt.m_vel;
	m_acc = fowt.m_acc;
    return *this;	
}

/*****************************************************
	Setters
*****************************************************/
void FOWT::readLinStiff(const std::string &data)
{
    // The mooring line stiffness in surge, sway and yaw are separated by commas in the input string
    std::vector<std::string> input = stringTokenize( data, "," );        
    
    if ( input.size() != 3 )
    {
		throw std::runtime_error("Unable to read linear stiffness in input line " + std::to_string(IO::getInLineNumber()) + ". Wrong number of parameters.");
	}
    
    for ( int ii = 0; ii < input.size(); ++ii )
    {
        readDataFromString( input.at(ii), m_linStiff(ii) );
    }    
}


void FOWT::setFloater(Floater &floater)
{
	m_floater = floater;
}




/*****************************************************
	Getters
*****************************************************/
vec::fixed<3> FOWT::CoG()
{
	// If CoG was not calculated yet, calculate it
	if (!is_finite(m_CoG))
	{
		m_CoG = m_floater.CoG();
	}
	
	return m_CoG;
}

double FOWT::mass()
{
	// If m_mass was not calculated yet, calculate it
	if (!is_finite(m_mass))
	{
		m_mass = m_floater.mass(); 
	}
	
	return m_mass;
}

vec::fixed<6> FOWT::pos() const
{
	return m_pos;
}

vec::fixed<6> FOWT::vel() const
{
	return m_vel;
}

vec::fixed<6> FOWT::acc() const
{
	return m_acc;
}


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
	output = output + "\tInertia Matrix:\t" + m_floater.printInertia() + "\n";
	output = output + "\tMorison Elements:\n" + m_floater.printMorisonElements() + "\n";

	return output;
}



/*****************************************************
	Forces, acceleration, position, etc
*****************************************************/
// Update FOWT position, velocity, acceleration and any other necessary state
void FOWT::update(const vec::fixed<6> &pos, const vec::fixed<6> &vel, const vec::fixed<6> &acc)
{
	m_pos = pos;
	m_vel = vel;
	m_acc = acc;
	m_floater.updatePosVelAcc(m_pos, m_vel, m_acc);
}

vec::fixed<6> FOWT::acceleration(const ENVIR &envir)
{
	IO::print2outLine(IO::OUTFLAG_FOWT_POS, m_pos);

	vec::fixed<6> acc;
	vec::fixed<6> inertiaMatrix = m_floater.inertia();

	vec::fixed<6> force = totalForce(envir);
	IO::print2outLine(IO::OUTFLAG_TOTAL_FORCE, force);

	acc[0] = force[0] / mass();
	acc[1] = force[1] / mass();
	acc[2] = force[2] / mass();
	acc[3] = force[3] / inertiaMatrix[0];
	acc[4] = force[4] / inertiaMatrix[1];
	acc[5] = force[5] / inertiaMatrix[2];

	IO::print2outLine(IO::OUTFLAG_TOTAL_FORCE,acc);

	return acc;
}

vec::fixed<6> FOWT::hydrodynamicForce(const ENVIR &envir)
{
	return m_floater.hydrodynamicForce(envir);
}

vec::fixed<6> FOWT::hydrostaticForce(const ENVIR &envir)
{
	return m_floater.hydrostaticForce(envir);
}

vec::fixed<6> FOWT::mooringForce()
{	
	return vec::fixed<6> {-m_linStiff[0]*m_pos[0], -m_linStiff[1]*m_pos[1], 0, 0, 0, -m_linStiff[2]*m_pos[5]};
}

vec::fixed<6> FOWT::weightForce(const ENVIR &envir)
{
	return vec::fixed<6> {0,0, -envir.gravity() * mass(), 0, 0, 0};
}

vec::fixed<6> FOWT::totalForce(const ENVIR &envir)
{
	return (hydrodynamicForce(envir) + hydrostaticForce(envir) + mooringForce() + weightForce(envir));
}