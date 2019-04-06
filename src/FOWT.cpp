#include "FOWT.h"
#include "IO.h"


/*****************************************************
	Constructors
*****************************************************/
FOWT::FOWT() : m_linStiff(fill::zeros), m_mass(datum::nan), 
			   m_disp(fill::zeros), m_vel(fill::zeros), m_acc(fill::zeros)
{
	m_CoG.fill(datum::nan);	
}

/*****************************************************
	Overloaded operators
*****************************************************/
FOWT& FOWT::operator=(const FOWT &fowt)
{
	m_floater = fowt.m_floater;
	m_linStiff = fowt.m_linStiff;
	m_disp = fowt.m_disp;
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

void FOWT::setRNA(RNA &rna)
{
	m_rna = rna;

	// Need to set the vertical distance between the hub and the CoG
	if (!arma::is_finite(m_floater.CoG()))
	{
		throw std::runtime_error( "Need to set the floater CoG before calling FOWT::setRNA(RNA &rna)" );
	}

	m_rna.setHubHeight2CoG(CoG().at(2));
}




/*****************************************************
	Getters
*****************************************************/
vec::fixed<3> FOWT::CoG()
{
	// If CoG was not calculated yet, calculate it
	if (!arma::is_finite(m_CoG))
	{
		m_CoG = m_floater.CoG();
	}
	
	return m_CoG;
}

double FOWT::mass()
{
	// If m_mass was not calculated yet, calculate it
	if (!arma::is_finite(m_mass))
	{
		m_mass = m_floater.mass(); 
	}
	
	return m_mass;
}

vec::fixed<6> FOWT::disp() const
{
	return m_disp;
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

	mat::fixed<6, 6> A = m_floater.addedMass(1);
	output = output + "\tAdded Mass for unitary density:\n";
	for (int ii = 0; ii < 6; ++ii)
	{
		output = output + "\t\t";
		for (int jj = 0; jj < 6; ++jj)
		{
			output = output + std::to_string(A(ii, jj));

			(jj == 5) ? (output = output + '\n') : (output = output + " \t\t ; \t\t");			
		}
	}	

	output = output + '\n';

	return output;
}


std::string FOWT::printRNA() const
{
	std::string output = "";

	output = output + "\tRotor Speed:\t" + std::to_string(m_rna.rotorSpeed()) + "\n";
	output = output + "\tRotor Tilt :\t" + std::to_string(m_rna.rotorTilt()) + "\n";
	output = output + "\tRotor Yaw:\t" + std::to_string(m_rna.rotorYaw()) + "\n";
	output = output + "\tNumBlades:\t" + std::to_string(m_rna.numBlades()) + "\n";
	output = output + "\tBlades Precone:\t" + std::to_string(m_rna.bladePrecone(0)) + "\n";
	output = output + "\tBlades Pitch:\t" + std::to_string(m_rna.bladePitch(0)) + "\t" + std::to_string(m_rna.bladePitch(1)) + "\t" + std::to_string(m_rna.bladePitch(2)) + "\n";
	output = output + "\tHub Radius:\t" + std::to_string(m_rna.hubRadius()) + "\n";
	output = output + "\tHub Height:\t" + std::to_string(m_rna.hubHeight()) + "\n";
	output = output + "\tOverhang:\t" + std::to_string(m_rna.overhang()) + "\n";
	output = output + "\tBlade aerodynamic properties:\n" + m_rna.printBladeAero();
	output = output + "\tAirfoils properties:\n" + m_rna.printAirfoils();

	return output;
}


/*****************************************************
	Forces, acceleration, displacement, etc
*****************************************************/
// Update FOWT displacement, velocity, acceleration and any other necessary state
void FOWT::update(const vec::fixed<6> &disp, const vec::fixed<6> &vel, const vec::fixed<6> &acc)
{
	m_disp = disp;
	m_vel = vel;
	m_acc = acc;
	m_floater.updateDispVelAcc(m_disp, m_vel, m_acc);
}

vec::fixed<6> FOWT::calcAcceleration(const ENVIR &envir)
{
	IO::print2outLine(IO::OUTFLAG_FOWT_DISP, m_disp);
	IO::print2outLine(IO::OUTFLAG_FOWT_VEL, m_vel);

	vec::fixed<6> acc(fill::zeros);	

	// Inertia matrix including added matrix
	mat::fixed<6, 6> inertiaMatrix = m_floater.addedMass(envir.watDensity()) + m_floater.inertiaMatrix();

	// Calculate the total force acting on the FOWT
	vec::fixed<6> force = totalForce(envir);

	if (!envir.isSurgeActive())
	{
		force[0] = 0;
	}

	if (!envir.isSwayActive())
	{
		force[1] = 0;
	}

	if (!envir.isHeaveActive())
	{
		force[2] = 0;
	}

	if (!envir.isRollActive())
	{
		force[3] = 0;
	}

	if (!envir.isPitchActive())
	{
		force[4] = 0;
	}

	if (!envir.isYawActive())
	{
		force[5] = 0;
	}

	IO::print2outLine(IO::OUTFLAG_TOTAL_FORCE, force);

	// Solve inertiaMatrix * acc = force
	acc = arma::solve(inertiaMatrix, force);

	// Due to the coupling effects, it is necessary to set the accelerations of the inactive DoFs to 0
	if (!envir.isSurgeActive())
	{
		acc[0] = 0;
	}

	if (!envir.isSwayActive())
	{
		acc[1] = 0;
	}

	if (!envir.isHeaveActive())
	{
		acc[2] = 0;
	}

	if (!envir.isRollActive())
	{
		acc[3] = 0;
	}

	if (!envir.isPitchActive())
	{
		acc[4] = 0;
	}

	if (!envir.isYawActive())
	{
		acc[5] = 0;
	}


	IO::print2outLine(IO::OUTFLAG_FOWT_ACC, m_acc);

	return acc;
}

vec::fixed<6> FOWT::hydrodynamicForce(const ENVIR &envir)
{
	return m_floater.hydrodynamicForce(envir);
}

vec::fixed<6> FOWT::hydrostaticForce(const double watDensity, const double gravity)
{
	return m_floater.hydrostaticForce(watDensity, gravity);
}

vec::fixed<6> FOWT::aeroForce(const ENVIR &envir)
{
	return m_rna.aeroForce(envir, m_disp + CoG(), m_vel);
}

vec::fixed<6> FOWT::mooringForce()
{	
	return vec::fixed<6> {-m_linStiff(0)*m_disp(0), -m_linStiff(1)*m_disp(1), 0, 0, 0, -m_linStiff(2)*m_disp(5)};
}

vec::fixed<6> FOWT::weightForce(const double gravity)
{
	return vec::fixed<6> {0,0, -gravity * mass(), 0, 0, 0};
}

vec::fixed<6> FOWT::totalForce(const ENVIR &envir)
{
	vec::fixed<6> aeroTest = aeroForce(envir);
	return (hydrodynamicForce(envir) + hydrostaticForce(envir.watDensity(), envir.gravity()) + mooringForce() + weightForce(envir.gravity()));
}