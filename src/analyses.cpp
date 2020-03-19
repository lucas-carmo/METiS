#include "analyses.h"

#include <iostream>
#include <armadillo>
#include <iomanip> // For input/output manipulators

using namespace arma;

void timeDomainAnalysis(FOWT &fowt, ENVIR &envir)
{
	// The header of the formatted output file is written during the first time step
	IO::print2outLineHeader_turnOn();

	/*
		Variables for the RK4 method
	*/
	// FOWT state in the beginning
	vec::fixed<6> disp0(fowt.disp());
	vec::fixed<6> vel0(fowt.vel());
	vec::fixed<6> acc0(fowt.acc());

	// RK4: first estimation
	vec::fixed<6> disp_k1(arma::fill::zeros);
	vec::fixed<6> vel_k1(arma::fill::zeros);
	vec::fixed<6> acc_k1(arma::fill::zeros);

	// RK4: second estimation
	vec::fixed<6> disp_k2(arma::fill::zeros);
	vec::fixed<6> vel_k2(arma::fill::zeros);
	vec::fixed<6> acc_k2(arma::fill::zeros);

	// RK4: third estimation
	vec::fixed<6> disp_k3(arma::fill::zeros);
	vec::fixed<6> vel_k3(arma::fill::zeros);
	vec::fixed<6> acc_k3(arma::fill::zeros);

	// RK4: fourth estimation
	vec::fixed<6> disp_k4(arma::fill::zeros);
	vec::fixed<6> vel_k4(arma::fill::zeros);
	vec::fixed<6> acc_k4(arma::fill::zeros);

	// RK4: calculated values
	vec::fixed<6> disp_total(arma::fill::zeros);
	vec::fixed<6> vel_total(arma::fill::zeros);
	vec::fixed<6> acc_total(arma::fill::zeros);

	// make sure that the members of FOWT are updated
	fowt.update(disp0, vel0, acc0);

	while ( envir.time() <= envir.timeTotal() )
	{
		IO::print2outLine_turnOn();
		IO::print2outLine_decimal(envir.time()); // Time has to be printed as decimal to avoid problems with large numbers

		// FOWT state at the beginning of the time step
		disp0 = fowt.disp();
		vel0 = fowt.vel();
		acc0 = fowt.acc();

		// Print wave characteristics that may have been requested in the input file
		envir.printWaveCharact();

		// RK4: first estimation
		acc_k1 = fowt.calcAcceleration(envir);
		vel_k1 = acc_k1 * envir.timeStep();
		disp_k1 = vel0 * envir.timeStep();

		// After the first time step, we do not need to print anything else to the header of the formatted output file
		if (envir.time() == 0)
		{
			IO::print2outLineHeader_turnOff();
			IO::printOutLineHeader2outFile();
		}

		// Results are printed in the first estimation, since it is done with the state of the fowt and envir at the beginning of the time step
		IO::print2outLine_turnOff();
		IO::printOutLine2outFile();

		//
		// Other Runge-Kutta estimations for evaluating the next time step:
		//

		// RK4: second estimation
		// Update fowt and environment
		fowt.update( disp0 + disp_k1/2 , vel0 + vel_k1/2 , acc_k1);
		envir.stepTime(envir.timeStep()/2);

		acc_k2 = fowt.calcAcceleration(envir);
		vel_k2 = acc_k2 * envir.timeStep();
		disp_k2 = (vel0 + vel_k1/2) * envir.timeStep();

		// RK4: third estimation
		// Update only fowt. Environment is already at t+dt/2
		fowt.update( disp0 + disp_k2/2 , vel0 + vel_k2/2 , acc_k2);

		acc_k3 = fowt.calcAcceleration(envir);
		vel_k3 = acc_k3 * envir.timeStep();
		disp_k3 = (vel0 + vel_k2/2) * envir.timeStep();

		// RK4: fourth estimation
		// Update fowt and environment, which needs to be in t+dt (hence, just need to add dt/2)
		fowt.update( disp0 + disp_k3 , vel0 + vel_k3 , acc_k3);
		envir.stepTime(envir.timeStep()/2);

		acc_k4 = fowt.calcAcceleration(envir);
		vel_k4 = acc_k4 * envir.timeStep();
		disp_k4 = (vel0 + vel_k3) * envir.timeStep();

		// Calculate new state of the FOWT
		acc_total = (acc_k1 + 2*acc_k2 + 2*acc_k3 + acc_k4) / 6;
		vel_total = vel0 + (vel_k1 + 2*vel_k2 + 2*vel_k3 + vel_k4) / 6;
		disp_total = disp0 + (disp_k1 + 2*disp_k2 + 2*disp_k3 + disp_k4) / 6;

		// We only need to update fowt, as envir was already updated during the RK4 steps
		fowt.update(disp_total, vel_total, acc_total);

		// Print progress to the screen only after integer seconds.
		// The criterion is to verify if the current time is different
		// from its floor value by at most one time step.
		if (almostEqual(envir.time(), std::floor(envir.time()), envir.timeStep()))
		{
		std::cout << "   Progress: " << std::setprecision(0) << envir.time() << " of " << envir.timeTotal() << " seconds -- " << std::fixed << std::setprecision(1) << 100 * envir.time() / envir.timeTotal() << "%" << '\r';
		std::fflush(stdout);
		}
	}
}
