#include "analyses.h"

#include <iostream>
#include <armadillo>
#include <iomanip> // For input/output manipulators

using namespace arma;

void timeDomainAnalysis_FOWT(FOWT &fowt, ENVIR &envir)
{    
    IO::print2outLineHeader_turnOn(); // The header of the formatted output file is written during the first time step

    /*
        Variables for the RK4 method
    */
    // FOWT state in the beginning
    vec::fixed<6> pos0(fowt.pos());
    vec::fixed<6> vel0(fowt.vel());
    vec::fixed<6> acc0(fowt.acc());

    // RK4: first estimation
    vec::fixed<6> pos_k1(arma::fill::zeros);
    vec::fixed<6> vel_k1(arma::fill::zeros);
    vec::fixed<6> acc_k1(arma::fill::zeros);

    // RK4: second estimation
    vec::fixed<6> pos_k2(arma::fill::zeros);
    vec::fixed<6> vel_k2(arma::fill::zeros);
    vec::fixed<6> acc_k2(arma::fill::zeros);

    // RK4: third estimation
    vec::fixed<6> pos_k3(arma::fill::zeros);
    vec::fixed<6> vel_k3(arma::fill::zeros);
    vec::fixed<6> acc_k3(arma::fill::zeros);

    // RK4: fourth estimation
    vec::fixed<6> pos_k4(arma::fill::zeros);
    vec::fixed<6> vel_k4(arma::fill::zeros);
    vec::fixed<6> acc_k4(arma::fill::zeros);

    // RK4: calculated values
    vec::fixed<6> pos_total(arma::fill::zeros);
    vec::fixed<6> vel_total(arma::fill::zeros);
    vec::fixed<6> acc_total(arma::fill::zeros);

	// make sure that the members of FOWT are updated
	fowt.update(pos0, vel0, acc0);

    while ( envir.time() <= envir.timeTotal() )    
    {          
        IO::print2outLine_turnOn();	
        IO::print2outLine(envir.time());

        // FOWT state at the beginning of the time step
        pos0 = fowt.pos();	
        vel0 = fowt.vel();
        acc0 = fowt.acc();   

        // RK4: first estimation
        acc_k1 = fowt.calcAcceleration(envir);
        vel_k1 = acc_k1 * envir.timeStep();
        pos_k1 = vel0 * envir.timeStep();        

        // Printing is done only in the first estimation, since it is done with the previous position
        IO::print2outLine_turnOff();
        
        // After the first time step, we do not need to print anything else to the header of the formatted output file
        if (envir.time() == 0)
        {
			IO::print2outLineHeader_turnOff();
			IO::printOutLineHeader2outFile();
        }


        // RK4: second estimation
        // Update fowt and environment 
        fowt.update( pos0 + pos_k1/2 , vel0 + vel_k1/2 , acc_k1); 
        envir.stepTime(envir.timeStep()/2);

        acc_k2 = fowt.calcAcceleration(envir);
        vel_k2 = acc_k2 * envir.timeStep();
        pos_k2 = (vel0 + vel_k1/2) * envir.timeStep();     


        // RK4: third estimation
        // Update only fowt. Environment is already at t+dt/2
        fowt.update( pos0 + pos_k2/2 , vel0 + vel_k2/2 , acc_k2); 

        acc_k3 = fowt.calcAcceleration(envir);
        vel_k3 = acc_k3 * envir.timeStep();
        pos_k3 = (vel0 + vel_k2/2) * envir.timeStep(); 
        

        // RK4: fourth estimation
        // Update fowt and environment, which needs to be in t+dt (hence, just need to add dt/2)
        fowt.update( pos0 + pos_k3 , vel0 + vel_k3 , acc_k3);
        envir.stepTime(envir.timeStep()/2);

        acc_k4 = fowt.calcAcceleration(envir);
        vel_k4 = acc_k4 * envir.timeStep();
        pos_k4 = (vel0 + vel_k3) * envir.timeStep(); 


        // Calculate new state of the FOWT
        acc_total = (acc_k1 + 2*acc_k2 + 2*acc_k3 + acc_k4) / 6;
        vel_total = vel0 + (vel_k1 + 2*vel_k2 + 2*vel_k3 + vel_k4) / 6;
        pos_total = pos0 + (pos_k1 + 2*pos_k2 + 2*pos_k3 + pos_k4) / 6;      

        // We only need to update fowt, as envir was already updated during the RK4 steps
        fowt.update(pos_total, vel_total, acc_total);

		IO::printOutLine2outFile();
        
		// Print progress to the screen
        std::cout << round(100 * envir.time() / envir.timeTotal()) << "%" << '\r';
        std::fflush(stdout);        
    }
}