#include "analyses.h"


#include <cmath>
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
	// - Components from 0 to 5: quantities obtained from the solution of the first-order hydrodynamic problem
	// - Components from 6 to 11: total quantities obtained considering all the forces
	vec::fixed<12> disp0(join_cols(fowt.disp_1stOrd(), fowt.disp()));
	vec::fixed<12> vel0(join_cols(fowt.vel_1stOrd(), fowt.vel()));

	// RK4 estimations
	vec::fixed<12> disp_k1(arma::fill::zeros);
	vec::fixed<12> vel_k1(arma::fill::zeros);
	vec::fixed<12> acc_k1(arma::fill::zeros);

	vec::fixed<12> disp_k2(arma::fill::zeros);
	vec::fixed<12> vel_k2(arma::fill::zeros);
	vec::fixed<12> acc_k2(arma::fill::zeros);

	vec::fixed<12> disp_k3(arma::fill::zeros);
	vec::fixed<12> vel_k3(arma::fill::zeros);
	vec::fixed<12> acc_k3(arma::fill::zeros);

	vec::fixed<12> disp_k4(arma::fill::zeros);
	vec::fixed<12> vel_k4(arma::fill::zeros);
	vec::fixed<12> acc_k4(arma::fill::zeros);

	vec::fixed<12> disp_k5(arma::fill::zeros);
	vec::fixed<12> vel_k5(arma::fill::zeros);
	vec::fixed<12> acc_k5(arma::fill::zeros);

	vec::fixed<12> disp_k6(arma::fill::zeros);
	vec::fixed<12> vel_k6(arma::fill::zeros);
	vec::fixed<12> acc_k6(arma::fill::zeros);

	// RK4: calculated values
	vec::fixed<12> disp_total(arma::fill::zeros);
	vec::fixed<12> vel_total(arma::fill::zeros);
	vec::fixed<12> acc_total(arma::fill::zeros);

	// Check if this is a moving or fixed body
	bool movBody = false;
	for (int ii = 0; ii < 6; ++ii)
	{
		if (fowt.isDoFActive(ii))
		{
			movBody = true;
			break;
		}
	}

	// Fifth-order Runge-Kutta method with adaptive stepsize
	//
	// If the error 'err' between the RK5 and the embedded 4th order method is 
	// greater than what is required by delta0, the time step is reduced. Otherwise, it is increased.
	double h = envir.timeStep();
	double h_aux = h; // Extra time step used to print every print step without losing the time step calculated by the adaptive stepsize algorithm
	double epsRel{ 0.005 };
	double epsAbs{ 1e-6 };
	vec::fixed<12> delta0{ 0 };
	vec::fixed<12> delta1{ 0 };
	vec::fixed<12> dispErr(fill::zeros);
	double safFact = 0.8;
	double factor{ 1 };
	bool condition = true;

	// make sure that the members of FOWT are updated
	fowt.update(envir, disp0, vel0);

	// The header of the formatted output file is written during the first time step and is then turned off
	IO::print2outLineHeader_turnOn();
	double currentPrintStep{ 0 };
	double nextPrintStep{ 0 };
	bool h_from_aux{ false };
	while (envir.time() <= envir.timeTotal())
	{
		double remainder = std::fmod(envir.time(), envir.printStep());
		bool flagPrintStep(almostEqual(remainder, 0, epsAbs) || almostEqual(remainder, envir.printStep(), epsAbs));
		if (flagPrintStep)
		{
			currentPrintStep = envir.time();
			IO::print2outLine_turnOn();
			IO::print2outLine_decimal(envir.time()); // Time has to be printed as decimal to avoid problems with large numbers
		}

		// FOWT state at the beginning of the time step
		disp0 = join_cols(fowt.disp_1stOrd(), fowt.disp());
		vel0 = join_cols(fowt.vel_1stOrd(), fowt.vel());

		// Output wave and FOWT characteristics that may have been requested in the input file
		envir.printWaveCharact();
		IO::print2outLine(IO::OUTFLAG_FOWT_DISP_1ST, disp0.rows(0, 5));
		IO::print2outLine(IO::OUTFLAG_FOWT_VEL_1ST, vel0.rows(0, 5));
		IO::print2outLine(IO::OUTFLAG_FOWT_ACC_1ST, acc_total.rows(0, 5));
		IO::print2outLine(IO::OUTFLAG_FOWT_DISP, disp0.rows(6,11));
		IO::print2outLine(IO::OUTFLAG_FOWT_VEL, vel0.rows(6, 11));
		IO::print2outLine(IO::OUTFLAG_FOWT_ACC, acc_total.rows(6, 11));		
		IO::print2outLine(IO::OUTFLAG_FOWT_DISP_SD, fowt.disp_sd());

		// Print progress to the screen
		if (flagPrintStep)
		{
			std::cout << "   Progress: " << std::setprecision(2) << envir.time() << " of " << envir.timeTotal() << " seconds -- " << std::fixed << std::setprecision(1) << 100 * envir.time() / envir.timeTotal() << "%" << '\r';
			std::fflush(stdout);
		}

		// Acceleration evaluated considering the previous time step does not change
		// in the loop for the adaptive stepsize 
		acc_k1 = fowt.calcAcceleration(envir);

		// After the first time step, we do not need to print anything else to the header of the formatted output file
		if (envir.time() == 0)
		{
			IO::printOutLineHeader2outFile();
			IO::print2outLineHeader_turnOff();
		}
		// Results are printed in the first estimation, since it is done with the state of the fowt and envir at the beginning of the time step
		IO::printOutLine2outFile();
		IO::print2outLine_turnOff();

		// Loop for the adaptive stepsize control
		condition = true;		
		while (condition && movBody)
		{			
			// RK5: first estimation
			vel_k1 = acc_k1 * h;
			disp_k1 = vel0 * h;

			// RK5: second estimation
			envir.stepTime(a2*h);
			fowt.update(envir, disp0 + b21 * disp_k1, vel0 + b21 * vel_k1);

			acc_k2 = fowt.calcAcceleration(envir);
			vel_k2 = acc_k2 * h;
			disp_k2 = (vel0 + b21 * vel_k1) * h;

			// RK5: third estimation
			envir.stepTime(-a2 * h + a3 * h); // To step only a3*h from initial time
			fowt.update(envir, disp0 + b31 * disp_k1 + b32 * disp_k2, vel0 + b31 * vel_k1 + b32 * vel_k2);

			acc_k3 = fowt.calcAcceleration(envir);
			vel_k3 = acc_k3 * h;
			disp_k3 = (vel0 + b31 * vel_k1 + b32 * vel_k2) * h;

			// RK5: fourth estimation
			envir.stepTime(-a3 * h + a4 * h); // To step only a4*h from initial time
			fowt.update(envir, disp0 + b41 * disp_k1 + b42 * disp_k2 + b43 * disp_k3, vel0 + b41 * vel_k1 + b42 * vel_k2 + b43 * vel_k3);

			acc_k4 = fowt.calcAcceleration(envir);
			vel_k4 = acc_k4 * h;
			disp_k4 = (vel0 + b41 * vel_k1 + b42 * vel_k2 + b43 * vel_k3) * h;

			// RK5: fifth estimation
			envir.stepTime(-a4 * h + a5 * h); // To step only a5*h from initial time
			fowt.update(envir, disp0 + b51 * disp_k1 + b52 * disp_k2 + b53 * disp_k3 + b54 * disp_k4, vel0 + b51 * vel_k1 + b52 * vel_k2 + b53 * vel_k3 + b54 * vel_k4);

			acc_k5 = fowt.calcAcceleration(envir);
			vel_k5 = acc_k5 * h;
			disp_k5 = (vel0 + b51 * vel_k1 + b52 * vel_k2 + b53 * vel_k3 + b54 * vel_k4) * h;

			// RK5: sixth estimation
			envir.stepTime(-a5 * h + a6 * h); // To step only a6*h from initial time
			fowt.update(envir, disp0 + b61 * disp_k1 + b62 * disp_k2 + b63 * disp_k3 + b64 * disp_k4 + b65 * disp_k5, vel0 + b61 * vel_k1 + b62 * vel_k2 + b63 * vel_k3 + b64 * vel_k4 + b65 * vel_k5);

			acc_k6 = fowt.calcAcceleration(envir);
			vel_k6 = acc_k6 * h;
			disp_k6 = (vel0 + b61 * vel_k1 + b62 * vel_k2 + b63 * vel_k3 + b64 * vel_k4 + b65 * vel_k5) * h;

			// Results obtained with the RK5 and the associated error			
			disp_total = c1 * disp_k1 + c2 * disp_k2 + c3 * disp_k3 + c4 * disp_k4 + c5 * disp_k5 + c6 * disp_k6;
			dispErr = (c1 - c1s)*disp_k1 + (c2 - c2s)*disp_k2 + (c3 - c3s)*disp_k3 + (c4 - c4s)*disp_k4 + (c5 - c5s)*disp_k5 + (c6 - c6s)*disp_k6;
			/*disp_total = c1 * vel_k1 + c2 * vel_k2 + c3 * vel_k3 + c4 * vel_k4 + c5 * vel_k5 + c6 * vel_k6;
			dispErr = (c1 - c1s)*vel_k1 + (c2 - c2s)*vel_k2 + (c3 - c3s)*vel_k3 + (c4 - c4s)*vel_k4 + (c5 - c5s)*vel_k5 + (c6 - c6s)*vel_k6;*/

			// Evaluate the error
			delta1 = arma::abs(dispErr);
			for (int ii = 0; ii < delta0.size(); ++ii)
			{				
				delta0.at(ii) = epsRel * ((std::abs(disp_total.at(ii)) > epsAbs) ? std::abs(disp_total.at(ii)) : epsAbs);
			}
			factor = arma::min(delta0 / delta1);

			if (!arma::is_finite(factor) && !h_from_aux)
			{
				throw std::runtime_error("The adaptive stepsize RK5 diverged.");
			}

			

			if (factor >= 1)
			{
				// Update for next time step.
				acc_total = c1 * acc_k1 + c2 * acc_k2 + c3 * acc_k3 + c4 * acc_k4 + c5 * acc_k5 + c6 * acc_k6;
				vel_total = vel0 + c1 * vel_k1 + c2 * vel_k2 + c3 * vel_k3 + c4 * vel_k4 + c5 * vel_k5 + c6 * vel_k6;
				disp_total += disp0;

				// When the time step is calculated in order to match the print step,
				// there are cases where the resulting 'h' is so small that it is not 
				// possible to step due to roundoff errors. Hence, set it to the desired value directly.
				if (h < epsAbs*epsRel && h_from_aux)
				{
					envir.setCurrentTime(nextPrintStep);
				}
				else
				{
					envir.stepTime(-a6 * h + h);
				}
		
				fowt.update(envir, disp_total, vel_total);
				fowt.update_sd(disp_total.rows(6,11), h);

				condition = false;

				double aux = safFact * std::pow(factor, 0.2);
				if (aux > 1 && !h_from_aux) // Avoid to reduce the time step when it was supposed to increase
				{
					h *= aux;
				}
				else if (h_from_aux)
				{
					h = h_aux;
				}
			}
			else
			{
				envir.stepTime(-a6 * h);
				h *= safFact * std::pow(factor, 0.25);
			}

			if (h > envir.printStep())
			{
				h = envir.printStep();
			}
		}				

		// If the next print step is going to be crossed in this time step, 
		// update time to the print step instead
		nextPrintStep = currentPrintStep + envir.printStep();
		h_from_aux = false;
		if ((envir.time()+h > nextPrintStep + epsAbs) && envir.time() < nextPrintStep - epsAbs)
		{
			h_aux = h;
			h_from_aux = true;
			h = nextPrintStep - envir.time();
		}

		if (!movBody)
		{
			envir.stepTime(h);
			h = envir.timeStep();
		}
	}
}
