//
//  main.cpp
//  dust
//
//  Created by Giancarlo Venturato on 28/11/24.
//
#include "simulation.hpp"
#include "root.h"

int main() {
    
    //First run the simulation for only the advection equation, particle size 1e-6, time domain 5Myrs
    alpha=0; a_d = 1e-6; R=1; Lt = 5e6;
    study_dust_evolution(R, alpha, "R1");
    
    
    //plot the formation time for all radii
    study_form_time();
    
    
    //study the height of the disk after 1Mrys for every alpha
    Lt = 1.5e6; study_diffusion();
    
    //Run the simulation with the diffusion implemented, time domain 2Myrs (twice the formation time at R=100)
    alpha=1e-4; a_d = 1e-6; R=100; Lt = 2e6;
    study_dust_evolution(R, alpha, "R100_diffusion");
    
    
    //study dust grow for rays of the gas giant
    cout << "\n\nI'm studying the particle size growth...\n\n";
    //to make it efficient you should keep Lt * a_d = 10 (for Nz = 1000)
    alpha = 1e-4; R = 5.2; Lt = 10; a_d = 1; //you need to init a_d as the smallest power of 10 that forms planetsimals
    double a_est = check_dust_grow(R, a_d);
    cout << "I estimate at R = " << R << ", that we need a_d = " << a_est << " for planetsimals formation" << endl;
    
    R = 9.6; Lt = 10; a_d = 1;
    a_est = check_dust_grow(R, a_d);
    cout << "I estimate at R = " << R << ", that we need a_d = " << a_est << " for planetsimals formation" << endl;

    R = 19; Lt = 100; a_d = 0.1;
    a_est = check_dust_grow(R, a_d);
    cout << "I estimate at R = " << R << ", that we need a_d = " << a_est << " for planetsimals formation" << endl;
    
    R = 30; Lt = 100;
    a_d = 0.1;
    a_est = check_dust_grow(R, a_d);
    cout << "I estimate at R = " << R << ", that we need a_d = " << a_est << " for planetsimals formation" << endl;
     
  
    return 0;
}
