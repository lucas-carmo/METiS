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
	m_filterSD_omega = 0;
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

void FOWT::setStiffnessMatrix(const double density, const double gravity)
{
	m_floater.setStiffnessMatrix(density, gravity);
}

void FOWT::evaluateQuantitiesAtBegin(const ENVIR &envir)
{
	if (envir.getTimeArray().size() == 0)
		throw std::runtime_error("Time array should be set before calling FOWT::setPropertiesWithIFFT.");

	m_floater.evaluateQuantitiesAtBegin(envir, m_hydroMode);
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

vec::fixed<6> FOWT::disp_1stOrd() const
{
	return m_disp_1stOrd;
}

vec::fixed<6> FOWT::vel_1stOrd() const
{
	return m_vel_1stOrd;
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

	mat::fixed<6, 6> K = m_floater.hydrostaticStiffness();
	output = output + "\tHydrostatic stiffness (same units as input):\n";
	for (int ii = 0; ii < 6; ++ii)
	{
		output = output + "\t\t";
		for (int jj = 0; jj < 6; ++jj)
		{
			output = output + std::to_string(K(ii, jj));

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
void FOWT::update(const ENVIR &envir, const vec::fixed<12> &disp, const vec::fixed<12> &vel)
{
	m_disp_1stOrd = disp.rows(0, 5);
	m_vel_1stOrd = vel.rows(0, 5);
	m_disp = disp.rows(6, 11);
	m_vel = vel.rows(6, 11);

	if (m_filterSD_omega < 0)
	{
		m_disp_sd = disp;
		m_vel_sd = vel;
	}

	// Aqui tem que passar os deslocamentos com relacao ao CoG do floater. Calcular aqui mesmo baseado na posicao do centro de referencia de movimento 
	m_floater.update(envir, m_disp, m_vel, m_disp_1stOrd, m_vel_1stOrd, m_disp_sd, m_vel_sd);
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

vec::fixed<12> FOWT::calcAcceleration(const ENVIR &envir)
{	
	mat::fixed<6, 6> addedMass(m_floater.addedMass_t0());
	IO::print2outLine(IO::OUTFLAG_ADDED_MASS_DIAG, addedMass.diag());

	mat::fixed<6, 6> inertiaMatrix = addedMass + m_floater.inertiaMatrix();	

	/*
		Calculate the force and acceleration due to first-order wave forces
	*/
	// The first order hydrodynamic forces (including hydrostatic and added mass) may be needed ahead
	// if second-order forces are required. Hence, assign the following two to avoid recalculating later.
	// The others are not necessary, as their evaluation is very fast.
	vec::fixed<6> hydrodynamicForce_1stOrd = m_floater.hydrodynamicForce_1stOrd(envir);
	vec::fixed<6> hydrostaticForce_1stOrd = m_floater.hydrostaticForce_stiffnessPart(true);
	vec::fixed<6> force_1stOrd = hydrodynamicForce_1stOrd + hydrostaticForce_1stOrd
							   + m_floater.hydrodynamicForce_drag1stOrd(envir) + mooringForce(true)
							   + m_floater.hydrostaticForce_staticBuoyancy(envir.watDensity(), envir.gravity()) + weightForce(envir.gravity());;

	// Calculate the acceleration only if at least one dof is activated 
	// (i.e. if at least one element of m_dofs is equal to 'true') 
	vec::fixed<6> acc_1stOrd(fill::zeros);
	if (std::find(m_dofs.begin(), m_dofs.end(), true) != m_dofs.end())
	{
		// Avoid coupling effects when a DoF is disabled and the others are not. 
		// Please note that the force was printed BEFORE this is done, in such a way 
		// that the full 6 component force vector is printed, even if the DoF is not active. 
		for (int ii = 0; ii < 6; ++ii)
		{
			if (!m_dofs[ii])
			{
				force_1stOrd.at(ii) = 0;

				// Even when the dof is disabled, it is necessary to provide a non-zero inertia 
				// entry at the main diagonal to avoid a singular matrix 
				inertiaMatrix.row(ii).zeros();
				inertiaMatrix.col(ii).zeros();
				inertiaMatrix.at(ii, ii) = 1;
			}
		}

		// Solve inertiaMatrix * acc = force 
		// Armadillo will throw its own exception if this computation fails. 
		acc_1stOrd = arma::solve(inertiaMatrix, force_1stOrd);
	}

	// Calculate second order forces, if necessary
	vec::fixed<6> hydroForce_2ndOrd(fill::zeros);
	if (m_hydroMode > 1)
	{
		m_floater.setNode1stAcc(acc_1stOrd);
		hydroForce_2ndOrd = m_floater.hydrodynamicForce_2ndOrd(envir, hydrodynamicForce_1stOrd + hydrostaticForce_1stOrd - m_floater.addedMass_t0()*acc_1stOrd);
	}

	// Calculate the total force acting on the FOWT
	// Remember that the weight is already included in hydroForce_1stOrd, but the drag is not!
	vec::fixed<6> force = hydrodynamicForce_1stOrd + hydroForce_2ndOrd + m_floater.hydrostaticForce_stiffnessPart(false) 
						+ m_floater.hydrodynamicForce_dragTotal(envir) + mooringForce(false) 
						+ m_floater.hydrostaticForce_staticBuoyancy(envir.watDensity(), envir.gravity()) + weightForce(envir.gravity()) + aeroForce(envir);
	IO::print2outLine(IO::OUTFLAG_TOTAL_FORCE, force);

	// Calculate the acceleration only if at least one dof is activated
	// (i.e. if at least one element of m_dofs is equal to 'true')
	vec::fixed<6> acc_total(fill::zeros);
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
		acc_total = arma::solve(inertiaMatrix, force);		
	}

	return join_cols(acc_1stOrd, acc_total);
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

vec::fixed<6> FOWT::mooringForce(bool flagUse1stOrd)
{
	vec::fixed<6> force{ 0, 0, 0, 0, 0, 0 };
	if (m_moorMode == 1)
	{
		if (flagUse1stOrd)
		{
			force = (-m_extLinStiff * m_disp_1stOrd + m_extConstForce);
		}
		else
		{
			force = (-m_extLinStiff * m_disp + m_extConstForce);
		}
	}

	IO::print2outLine(IO::OUTFLAG_MOOR_FORCE, force);

	return force;
}

vec::fixed<6> FOWT::weightForce(const double gravity)
{
	return vec::fixed<6> {0, 0, -gravity * mass(), 0, 0, 0};
}