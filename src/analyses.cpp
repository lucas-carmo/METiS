#include "analyses.h"

#include <iostream>
#include <armadillo>
#include <iomanip> // For input/output manipulators

using namespace arma;

void timeDomainAnalysis(FOWT &fowt, ENVIR &envir)
{
	// Avoid printing to output files before the simulation begins
	IO::print2outLine_turnOff();
	IO::print2outLineHeader_turnOff();

	/*
		Variables for the RK4 method
	*/
	// Cash-Karp Parameters for Embedded Runga-Kutta Method
	// From 'Numerical Recipes in C, Second edition, Press et al', Section 16.2
	double a2{ 0.2 }, a3{ 0.3 }, a4{ 0.6 }, a5{ 1.0 }, a6{ 7 / 8. };
	double b21{ 0.2 }, b31{ 0.075 }, b41{ 0.3 }, b51{ -11 / 54. }, b61{ 1631 / 55296. };
	double b32{ 0.225 }, b42{ -0.9 }, b52{ 2.5 }, b62{ 175 / 512. };
	double b43{ 1.2 }, b53{ -70 / 27. }, b63{ 575 / 13824. };
	double b54{ 35 / 27. }, b64{ 44275 / 110592. }, b65{ 253 / 4096. };
	double c1{ 37 / 378. }, c2{ 0 }, c3{ 250 / 621. }, c4{ 125 / 594. }, c5{ 0 }, c6{ 512 / 1771. };
	double c1s{ 2825 / 27648. }, c2s{ 0 }, c3s{ 18575 / 48384. }, c4s{ 13525 / 55296. }, c5s{ 277 / 14336. }, c6s{ 0.25 };

	// FOWT state in the beginning
	vec::fixed<6> disp0(fowt.disp());
	vec::fixed<6> vel0(fowt.vel());

	// RK4 estimations
	vec::fixed<6> disp_k1(arma::fill::zeros);
	vec::fixed<6> vel_k1(arma::fill::zeros);
	vec::fixed<6> acc_k1(arma::fill::zeros);

	vec::fixed<6> disp_k2(arma::fill::zeros);
	vec::fixed<6> vel_k2(arma::fill::zeros);
	vec::fixed<6> acc_k2(arma::fill::zeros);

	vec::fixed<6> disp_k3(arma::fill::zeros);
	vec::fixed<6> vel_k3(arma::fill::zeros);
	vec::fixed<6> acc_k3(arma::fill::zeros);

	vec::fixed<6> disp_k4(arma::fill::zeros);
	vec::fixed<6> vel_k4(arma::fill::zeros);
	vec::fixed<6> acc_k4(arma::fill::zeros);

	vec::fixed<6> disp_k5(arma::fill::zeros);
	vec::fixed<6> vel_k5(arma::fill::zeros);
	vec::fixed<6> acc_k5(arma::fill::zeros);

	vec::fixed<6> disp_k6(arma::fill::zeros);
	vec::fixed<6> vel_k6(arma::fill::zeros);
	vec::fixed<6> acc_k6(arma::fill::zeros);

	// RK4: calculated values
	vec::fixed<6> disp_total(arma::fill::zeros);
	vec::fixed<6> vel_total(arma::fill::zeros);
	vec::fixed<6> acc_total(arma::fill::zeros);

	// make sure that the members of FOWT are updated
	fowt.update(envir, disp0, vel0);

	// The header of the formatted output file is written during the first time step and is then turned off
	IO::print2outLineHeader_turnOn();

	double h = envir.timeStep();
	while ( envir.time() <= envir.timeTotal() )
	{

		IO::print2outLine_turnOn();
		IO::print2outLine_decimal(envir.time()); // Time has to be printed as decimal to avoid problems with large numbers

		// FOWT state at the beginning of the time step
		disp0 = fowt.disp();
		vel0 = fowt.vel();

		// Output wave and FOWT characteristics that may have been requested in the input file
		envir.printWaveCharact();
		IO::print2outLine(IO::OUTFLAG_FOWT_DISP, disp0);
		IO::print2outLine(IO::OUTFLAG_FOWT_VEL, vel0);
		IO::print2outLine(IO::OUTFLAG_FOWT_ACC, acc_total); // Acc calculated from previous time step
		IO::print2outLine(IO::OUTFLAG_FOWT_DISP_SD, fowt.disp_sd());


		// RK4: first estimation
		acc_k1 = fowt.calcAcceleration(envir);
		vel_k1 = acc_k1 * h;
		disp_k1 = vel0 * h;

		// After the first time step, we do not need to print anything else to the header of the formatted output file
		if (envir.time() == 0)
		{
			IO::print2outLineHeader_turnOff();
			IO::printOutLineHeader2outFile();
		}		
		// Results are printed in the first estimation, since it is done with the state of the fowt and envir at the beginning of the time step
		IO::print2outLine_turnOff();
		IO::printOutLine2outFile();

		/*
		Other Runge-Kutta estimations for evaluating the next time step:
		*/

		// RK4: second estimation
		envir.stepTime(a2*h);
		fowt.update(envir, disp0 + b21*disp_k1, vel0 + b21*vel_k1);				

		acc_k2 = fowt.calcAcceleration(envir);
		vel_k2 = acc_k2 * h;
		disp_k2 = (vel0 + b21*vel_k1) * h;

		// RK4: third estimation
		envir.stepTime(-a2*h + a3*h); // To step only a3*h from initial time
		fowt.update(envir, disp0 + b31 * disp_k1 + b32 * disp_k2, vel0 + b31 * vel_k1 + b32 * vel_k2);

		acc_k3 = fowt.calcAcceleration(envir);
		vel_k3 = acc_k3 * h;
		disp_k3 = (vel0 + b31 * vel_k1 + b32 * vel_k2) * h;		

		// RK4: fourth estimation
		envir.stepTime(-a3 * h + a4 * h); // To step only a4*h from initial time
		fowt.update(envir, disp0 + b41 * disp_k1 + b42 * disp_k2 + b43 * disp_k3, vel0 + b41 * vel_k1 + b42 * vel_k2 + b43 * vel_k3);

		acc_k4 = fowt.calcAcceleration(envir);
		vel_k4 = acc_k4 * h;
		disp_k4 = (vel0 + b41 * vel_k1 + b42 * vel_k2 + b43 * vel_k3) * h;

		// RK4: fifth estimation
		envir.stepTime(-a4 * h + a5 * h); // To step only a5*h from initial time
		fowt.update(envir, disp0 + b51 * disp_k1 + b52 * disp_k2 + b53 * disp_k3 + b54 * disp_k4, vel0 + b51 * vel_k1 + b52 * vel_k2 + b53 * vel_k3 + b54 * vel_k4);

		acc_k5 = fowt.calcAcceleration(envir);
		vel_k5 = acc_k5 * h;
		disp_k5 = (vel0 + b51 * vel_k1 + b52 * vel_k2 + b53 * vel_k3 + b54 * vel_k4) * h;

		// RK4: sixth estimation
		envir.stepTime(-a5 * h + a6 * h); // To step only a6*h from initial time
		fowt.update(envir, disp0 + b61 * disp_k1 + b62 * disp_k2 + b63 * disp_k3 + b64 * disp_k4 + b65 * disp_k5, vel0 + b61 * vel_k1 + b62 * vel_k2 + b63 * vel_k3 + b64 * vel_k4 + b65 * vel_k5);

		acc_k6 = fowt.calcAcceleration(envir);
		vel_k6 = acc_k6 * h;
		disp_k6 = (vel0 + b61 * vel_k1 + b62 * vel_k2 + b63 * vel_k3 + b64 * vel_k4 + b65 * vel_k5) * h;

		// Calculate new state of the FOWT
		acc_total = c1 * acc_k1 + c2 * acc_k2 + c3 * acc_k3 + c4 * acc_k4 + c5 * acc_k5 + c6*acc_k6;
		vel_total = vel0 + c1*vel_k1 + c2*vel_k2 + c3*vel_k3 + c4*vel_k4 + c5*vel_k5 + c6*vel_k6;
		disp_total = disp0 + c1*disp_k1 + c2*disp_k2 + c3*disp_k3 + c4*disp_k4 + c5*disp_k5 + c6*disp_k6;

		// Update for next time step
		envir.stepTime(-a6 * h + h);
		fowt.update(envir, disp_total, vel_total);
		fowt.update_sd(disp_total, h);

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
