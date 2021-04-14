#include "FOWT.h"
#include "IO.h"

#include <algorithm> // For std::find


/*****************************************************
	Constructors
*****************************************************/
FOWT::FOWT() : m_extLinStiff(fill::zeros), m_mass(datum::nan),
m_disp(fill::zeros), m_vel(fill::zeros),
m_disp_sd(fill::zeros), m_vel_sd(fill::zeros)
{
	m_CoG.fill(datum::nan);

	// Default values for the filter in case the user has not specified any.
	m_filterSD_omega = 0.10;
	m_filterSD_zeta = 0.20;
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


void FOWT::setExtLinStiff(const mat::fixed<6, 6> &extLinStiff)
{
	m_extLinStiff = extLinStiff;
}

void FOWT::setExtConstForce(const vec::fixed<6> &extConstForce)
{
	m_extConstForce = extConstForce;
}

void FOWT::setFilderSD(const double omega, const double zeta)
{
	m_filterSD_omega = omega;
	m_filterSD_zeta = zeta;
}


void FOWT::setFloater(Floater &floater)
{
	m_floater = floater;

	if (m_filterSD_omega < 0)
	{
		m_floater.setInstantSD(true);
	}
}

void FOWT::setRNA(RNA &rna)
{
	m_rna = rna;

	// Need to set the vertical distance between the hub and the CoG
	if (!arma::is_finite(m_floater.CoG()))
	{
		throw std::runtime_error("Need to set the floater CoG before calling FOWT::setRNA(RNA &rna)");
	}

	m_rna.setHubHeight2CoG(CoG().at(2));
}

void FOWT::setAddedMass_t0(const double density)
{
	m_floater.setAddedMass_t0(density);
}

void FOWT::setPropertiesWithIFFT(const ENVIR &envir)
{
	m_floater.setPropertiesWithIFFT(envir);
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

bool FOWT::isDoFActive(const int index)
{
	if (index >= 0 && index < 6)
		return m_dofs.at(index);
	else
		throw std::runtime_error("Invalid index in FOWT::isDoFActive(const int index).");
}

double FOWT::filterSD_omega() const
{
	return m_filterSD_omega;
}

double FOWT::filterSD_zeta() const
{
	return m_filterSD_zeta;
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

vec::fixed<6> FOWT::disp_sd() const
{
	return m_disp_sd;
}

vec::fixed<6> FOWT::constForce() const
{
	return m_extConstForce;
}

std::string FOWT::printLinStiff() const
{
	std::string output = "\n\t";
	for (int ii = 0; ii < m_extLinStiff.n_rows; ++ii)
	{
		for (int jj = 0; jj < m_extLinStiff.n_cols; ++jj)
		{
			output = output + std::to_string(m_extLinStiff.at(ii, jj));
			(jj == 5) ? (output = output + "\n\t") : (output = output + " \t ; \t");
		}
	}

	return output + "\n";
}


std::string FOWT::printFloater() const
{
	std::string output = "";

	output = output + "\tMass:\t" + m_floater.printMass() + "\n";
	output = output + "\tCoG:\t" + m_floater.printCoG() + "\n";
	output = output + "\tInertia Matrix:\t" + m_floater.printInertia() + "\n";
	output = output + "\tMorison Elements:\n" + m_floater.printMorisonElements() + "\n";

	mat::fixed<6, 6> A = m_floater.addedMass_t0();
	output = output + "\tDimensional added Mass at initial position (same units as input):\n";
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
// Update FOWT displacement, velocity, acceleration and any other necessary state.
void FOWT::update(const ENVIR &envir, const vec::fixed<6> &disp, const vec::fixed<6> &vel)
{
	m_disp = disp;
	m_vel = vel;	

	if (m_filterSD_omega < 0)
	{
		m_disp_sd = disp;
		m_vel_sd = vel;
	}

	// Aqui tem que passar os deslocamentos com relacao ao CoG do floater. Calcular aqui mesmo baseado na posicao do centro de referencia de movimento
	m_floater.update(envir, m_disp, m_vel, m_disp_sd, m_vel_sd);
}

// Evaluate (and update) the axis system that follows the slow position.
// Done so by applying a second-order low-pass filter to the displacement.
void FOWT::update_sd(const vec::fixed<6> &disp, const double dt)
{
	double wf = m_filterSD_omega;
	double zeta = m_filterSD_zeta;

	// If the filtering frequency is zero, there is nothing to be done here
	if (wf == 0)
	{
		return;
	}

	// If the filtering frequency is below zero, the slow drift position is actually 
	// equal to the instantaneous position, and it is updated in FOWT::update()
	if (wf > 0)
	{
		// Integrated using RK4. However, since the simulation time step must be small enough to 
		// capture wave frequency motions, which are way faster than slow drift motions, even an
		// explicit Euler method would work
		vec::fixed<6> acc_sd_k1 = -2 * wf*zeta*m_vel_sd - wf * wf * (m_disp_sd - m_disp);
		vec::fixed<6> vel_sd_k1 = acc_sd_k1 * dt;
		vec::fixed<6> disp_sd_k1 = m_vel_sd * dt;

		vec::fixed<6> acc_sd_k2 = -2 * wf*zeta*(m_vel_sd + vel_sd_k1 / 2) - wf * wf * (m_disp_sd + disp_sd_k1 / 2 - (m_disp + (disp - m_disp) / 2));
		vec::fixed<6> vel_sd_k2 = acc_sd_k2 * dt;
		vec::fixed<6> disp_sd_k2 = (m_vel_sd + vel_sd_k1 / 2) * dt;

		vec::fixed<6> acc_sd_k3 = -2 * wf*zeta*(m_vel_sd + vel_sd_k2 / 2) - wf * wf * (m_disp_sd + disp_sd_k2 / 2 - (m_disp + (disp - m_disp) / 2));
		vec::fixed<6> vel_sd_k3 = acc_sd_k3 * dt;
		vec::fixed<6> disp_sd_k3 = (m_vel_sd + vel_sd_k2 / 2) * dt;

		vec::fixed<6> acc_sd_k4 = -2 * wf*zeta*(m_vel_sd + vel_sd_k3) - wf * wf * (m_disp_sd + disp_sd_k3 - disp);
		vec::fixed<6> vel_sd_k4 = acc_sd_k3 * dt;
		vec::fixed<6> disp_sd_k4 = (m_vel_sd + vel_sd_k3) * dt;

		m_vel_sd += (vel_sd_k1 + 2 * vel_sd_k2 + 2 * vel_sd_k3 + vel_sd_k4) / 6;
		m_disp_sd += (disp_sd_k1 + 2 * disp_sd_k2 + 2 * disp_sd_k3 + disp_sd_k4) / 6;
	}
}


vec::fixed<6> FOWT::calcAcceleration(const ENVIR &envir)
{
	vec::fixed<6> acc(fill::zeros);	

	mat::fixed<6, 6> addedMass(fill::zeros);

	// Added mass matrix at the evaluated at slow position if the analysis is first-order. This
	// is mostly because the expressions are all expressed in terms of slow-drift variables (the ones 
	// with '_sd' in their names), since a real first-order analysis needs to consider the fixed mean position.
	//
	// If the analysis is second-order, it is reevaluated at each time step at the instantaneous position.
	//
	// However, if the slow position is fixed, just use the value that was evaluated at the beginning of the simulation
	if (m_filterSD_omega == 0 && m_hydroMode == 1)
	{
		addedMass = m_floater.addedMass_t0();
	}
	else
	{
		addedMass = m_floater.addedMass(envir.watDensity(), m_hydroMode);
	}
	
	IO::print2outLine(IO::OUTFLAG_ADDED_MASS_DIAG, addedMass.diag());

	mat::fixed<6, 6> inertiaMatrix = addedMass + m_floater.inertiaMatrix();	


	// Calculate the total force acting on the FOWT
	vec::fixed<6> force = totalForce(envir);	
	IO::print2outLine(IO::OUTFLAG_TOTAL_FORCE, force);

	// Calculate the acceleration only if at least one dof is activated
	// (i.e. if at least one element of m_dofs is equal to 'true')
	if (std::find(m_dofs.begin(), m_dofs.end(), true) != m_dofs.end())
	{
		// Avoid coupling effects when a DoF is disabled and the others are not.
		// For doing so, set the calculated force to zero if the dof is deactivated.
		// Please note that the force was printed BEFORE this is done, in such a way
		// that the full 6 component force vector is printed, even if the DoF is not active.
		for (int ii = 0; ii < 6; ++ii)
		{
			if (!m_dofs[ii])
			{
				force.at(ii) = 0;

				// Even when the dof is disabled, it is necessary to provide a non-zero inertia
				// entry at the main diagonal to avoid a singular matrix
				inertiaMatrix.row(ii).zeros();
				inertiaMatrix.col(ii).zeros();
				inertiaMatrix.at(ii, ii) = 1;
			}
		}

		// Solve inertiaMatrix * acc = force
		// Armadillo will throw its own exception if this computation fails.
		acc = arma::solve(inertiaMatrix, force);		
	}
	IO::print2outLine(IO::OUTFLAG_HD_ADD_MASS_FORCE, -(addedMass - m_floater.addedMass_t0())* acc);

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

	return m_floater.hydrostaticForce(envir);
}

vec::fixed<6> FOWT::aeroForce(const ENVIR &envir)
{
	if (m_aeroMode == 1)
	{
		vec::fixed<6> dbg = m_rna.aeroForce(envir, m_disp + join_cols(CoG(), vec::fixed<3> {0, 0, 0}), m_vel);
		IO::print2outLine(IO::OUTFLAG_DEBUG_VEC_6, dbg);
		return dbg;
	}

	return vec::fixed<6> {0, 0, 0, 0, 0, 0};
}

vec::fixed<6> FOWT::mooringForce()
{
	vec::fixed<6> force{ 0, 0, 0, 0, 0, 0 };
	if (m_moorMode == 1)
	{
		force = (-m_extLinStiff * m_disp + m_extConstForce);
	}

	IO::print2outLine(IO::OUTFLAG_MOOR_FORCE, force);

	return force;
}

vec::fixed<6> FOWT::weightForce(const double gravity)
{
	return vec::fixed<6> {0, 0, -gravity * mass(), 0, 0, 0};
}

vec::fixed<6> FOWT::totalForce(const ENVIR &envir)
{
	return (hydrodynamicForce(envir) + hydrostaticForce(envir) + mooringForce() + weightForce(envir.gravity()) + aeroForce(envir));
}
