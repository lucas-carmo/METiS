#include "FOWT.h"
#include "IO.h"


/*****************************************************
	Constructors
*****************************************************/
FOWT::FOWT() : m_extLinStiff(fill::zeros), m_mass(datum::nan), 
			   m_disp(fill::zeros), m_vel(fill::zeros), m_acc(fill::zeros)
{
	m_CoG.fill(datum::nan);	
}

/*****************************************************
	Setters
*****************************************************/
void FOWT::setHydroMode(const int hydroMode)
{
	m_hydroMode = hydroMode;
}

void FOWT::setAeroMode(const int aeroMode)
{
	m_aeroMode = aeroMode;
}

void FOWT::setMoorMode(const int moorMode)
{
	m_moorMode = moorMode;
}

void FOWT::setDoFs(std::array<bool, 6> &dofs)
{
    m_dofs = dofs;
}


void FOWT::readExtLinStiff(const std::string &data)
{
    // The mooring line stiffness in surge, sway and yaw are separated by commas in the input string
    std::vector<std::string> input = stringTokenize( data, "," );        
    
    if ( input.size() != 3 )
    {
		throw std::runtime_error("Unable to read linear stiffness in input line " + std::to_string(IO::getInLineNumber()) + ". Wrong number of parameters.");
	}
    
    for ( int ii = 0; ii < input.size(); ++ii )
    {
        readDataFromString( input.at(ii), m_extLinStiff(ii) );
    }    
}

void FOWT::readExtConstForce(const std::string &data)
{
	// The 6 components of the external constant force are separated by commas in the input string
	std::vector<std::string> input = stringTokenize(data, ",");

	if (input.size() != 6)
	{
		throw std::runtime_error("Unable to read external constant force in input line " + std::to_string(IO::getInLineNumber()) + ". Wrong number of parameters.");
	}

	for (int ii = 0; ii < input.size(); ++ii)
	{
		readDataFromString(input.at(ii), m_extConstForce(ii));
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
int FOWT::hydroMode() const
{
	return m_hydroMode;
}

int FOWT::aeroMode() const
{
	return m_aeroMode;
}

int FOWT::moorMode() const
{
	return m_moorMode;
}

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
	std::string output = "(" + std::to_string( m_extLinStiff(0) );
	for ( int ii = 1; ii < m_extLinStiff.n_elem; ++ii )
	{
		output = output + "," + std::to_string( m_extLinStiff(ii) );
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

	mat::fixed<6, 6> A = m_floater.addedMass(1, 1);
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

	output = output + "Use Tip Loss:\t" + std::to_string(m_rna.useTipLoss()) + "\n";
	output = output + "Use Hub Loss:\t" + std::to_string(m_rna.useHubLoss()) + "\n";
	output = output + "Use Skew Correction:\t" + std::to_string(m_rna.useSkewCorr()) + "\n";
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

std::string FOWT::printHydroMode() const
{
	return std::to_string(m_hydroMode);
}

std::string FOWT::printAeroMode() const
{
	return std::to_string(m_aeroMode);
}

std::string FOWT::printMoorMode() const
{
	return std::to_string(m_moorMode);
}


std::string FOWT::printDoF() const
{
	std::string output = "";
	for (int ii = 0; ii < 6; ++ii)
	{
		output += std::to_string(m_dofs[ii]) + ' ';
	}

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
	m_floater.update(m_disp, m_vel, m_acc); // Aqui tem que passar os deslocamentos com relacao ao CoG do floater. Calcular aqui mesmo baseado na posicao do centro de referencia de movimento
}

vec::fixed<6> FOWT::calcAcceleration(const ENVIR &envir)
{
	IO::print2outLine(IO::OUTFLAG_FOWT_DISP, m_disp);
	IO::print2outLine(IO::OUTFLAG_FOWT_VEL, m_vel);

	vec::fixed<6> acc(fill::zeros);

	// Inertia matrix including added matrix
	mat::fixed<6, 6> inertiaMatrix = m_floater.addedMass(envir.watDensity(), m_hydroMode) + m_floater.inertiaMatrix();

	// Calculate the total force acting on the FOWT
	vec::fixed<6> force = totalForce(envir);
	IO::print2outLine(IO::OUTFLAG_TOTAL_FORCE, force);

	// Avoid coupling effects when a DoF is disabled and the others are not.
	// For doing so, set the calculated force to zero if the dof is deactivated.
	// Please note that the force was printed BEFORE this is done, in such a way
	// that the full 6 component force vector is printed, even if the DoF is not active.
	for (int ii = 0; ii < 6; ++ii)
	{
		if (!m_dofs[ii])
		{
			force[ii] = 0;
		}
	}

	// Solve inertiaMatrix * acc = force
	acc = arma::solve(inertiaMatrix, force);

	// Due to the coupling effects, it is necessary to set the accelerations of the inactive DoFs to 0
	for (int ii = 0; ii < 6; ++ii)
	{
		if (!m_dofs[ii])
		{
			acc[ii] = 0;
		}
	}

	IO::print2outLine(IO::OUTFLAG_FOWT_ACC, m_acc);

	return acc;
}

vec::fixed<6> FOWT::hydrodynamicForce(const ENVIR &envir)
{
	if (m_hydroMode == 0)
	{
		return vec::fixed<6> {0, 0, 0, 0, 0, 0};
	}

	return m_floater.hydrodynamicForce(envir, m_hydroMode);
}

vec::fixed<6> FOWT::hydrostaticForce(const ENVIR &envir)
{
	if (m_hydroMode == 0)
	{
		return vec::fixed<6> {0, 0, 0, 0, 0, 0};
	}

	return m_floater.hydrostaticForce(envir, m_hydroMode);
}

vec::fixed<6> FOWT::aeroForce(const ENVIR &envir)
{
	if (m_aeroMode == 1) 
	{
		return m_rna.aeroForce(envir, m_disp + join_cols(CoG(), vec::fixed<3> {0, 0, 0}), m_vel);
	}

	return vec::fixed<6> {0, 0, 0, 0, 0, 0};
}

vec::fixed<6> FOWT::mooringForce()
{	
	if (m_moorMode == 1)
	{
		return (vec::fixed<6> {-m_extLinStiff(0)*m_disp(0), -m_extLinStiff(1)*m_disp(1), 0, 0, 0, -m_extLinStiff(2)*m_disp(5)} + m_extConstForce);
	}

	return vec::fixed<6> {0, 0, 0, 0, 0, 0};
}

vec::fixed<6> FOWT::weightForce(const double gravity)
{
	return vec::fixed<6> {0,0, -gravity * mass(), 0, 0, 0};
}

vec::fixed<6> FOWT::totalForce(const ENVIR &envir)
{		
	return (hydrodynamicForce(envir) + hydrostaticForce(envir) + mooringForce() + weightForce(envir.gravity()) + aeroForce(envir));
}