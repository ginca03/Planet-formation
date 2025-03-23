

#ifndef simulation_h
#define simulation_h

//C++ libraries
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>


using namespace std;

//Physics constants
const double G = 6.67430e-11;   // Gravitational constant
const double M_sun = 2e30;      // Sun Mass (kg)
const double AU = 1.5e11;       // UM (m)
const double kb = 1.38e-23;     // Boltzmann constant (J/K)
const double mH = 1.67e-27;     // Hidrogen mass (kg)

//Physics parameters
double a_d_0 = 1e-6;          // particle width  (m)
const double rho_in = 1000.0;   // intern particle mass (kg/m^3)
const double rho_g_0 = 5e-7;     //density of gas at z=0 and R=0 (kg m^-3)
const double T0 = 264.0;        // temperature at 1 AU (K)
const double year = 60*60*24*365; //a year in seconds ≈ 1e-7

//Function prototypes for parameters
double T(double); double Omega(double); double c_s(double); double rho_g(double, double); double t_s(double, double); double u(double, double);

//Root function prototypes
void plot_rho();
void save_rho(const vector<vector<double>>& rho_d, const char* filename);
void scatter_plot(const vector<vector<double>> &data, string s);
void scatter_plot_spline(const vector<vector<double>> &data, string s);


//Simulation parameters
int Nz = 2000;         // Grid points z
const double Lz = 4.0; // total height [H(R)]

double Lt = 2e6; // default time domain = 2Myears
double dt = 100;     // default time step (years)

double dz = Lz / (Nz-1);    //z-step (fractions of H(R)
int Nt = static_cast<int>(ceil(Lt/dt)); // number of time steps

const double NR = 20;           //number of R steps
const double R_max = 100; const double R_min = 1;

//Scaling: R is in AU, t is in years, z in fractions of H(R) rho is in rho_g_0

//Declaring global variables
vector<vector<double> > rho_d;
double R = 1;
double a_d = a_d_0;          // particle width  (m)
double alpha = 0.0;


double H(double R){return c_s(R)/Omega(R); } //Scale Height of the disk at R (m)

//Temperature at R (K)
double T(double R){return T0 / sqrt(R); }

//Keplerian angular velocity (s^-1)
double Omega(double R) {
    return sqrt(G * M_sun / pow(R * AU,3)); }

//sound velocity (m/s)
double c_s(double R) {
    return sqrt(kb * T(R) / (2 * mH));}

//gas density at R and z (rho_g_0) // z is in H(R)
double rho_g(double R, double z){
    return pow(R,-3) *  exp(-0.5 * pow(z, 2)); }

//Settling timescale (s)
double t_s(double R, double z, double a_d) {
    return sqrt(M_PI / 8.0) * (rho_in * a_d) / (rho_g_0 * rho_g(R, z) * c_s(R)); }

//Dust velocity (H(R)/year)
double u(double R, double z, double a_d){return -pow(Omega(R)*year,2) * (z) * (t_s(R, z, a_d)/year) ;}


double choose_dt(double R, double alpha, double a_d){
    
    dt = 0.8 * dz / abs(u(R,Lz/2,a_d)); // choose the timesteps so that we don't have problems with advection
    
    //calculate diffusion coefficient
    double D = alpha * pow(c_s(R),2) / Omega(R) * (year/(H(R)*H(R))); // D is in H(R)^2/years
    //cout << "D = " << D << " H(R)^2/years" << endl; //debug
    
    //if time-step of advection is too big we choose one smaller fitted for the diffusion
    if((D * dt / (dz*dz)) > 0.5){
        dt = (dz*dz) / (3*D);}
    
    return dt;
    
}


//solve using dynamic values of dt
void solve_dyn(double R, double alpha){
    //lambda exressions
    auto u_z = [R](int i){return u(R, (i-Nz/2) * dz, a_d);};
    auto rho_gas = [R](int i) { return rho_g(R, (i-Nz/2) * dz); };
    
    dt = choose_dt(R, alpha, a_d);
    Nt = static_cast<int>(ceil(Lt/dt)); //choose number of steps so that we cover the domain
    cout << "Solving for R = " << R << ", alpha = " << alpha << ", a_d = " << a_d << ",\tNt = " << Nt << ",\tLt = " << Lt << endl;
    
    //calculate diffusion coefficient
    double D = alpha * pow(c_s(R),2) / Omega(R) * (year/(H(R)*H(R))); // D is in H(R)^2/years
    //cout << "D = " << D << " H(R)^2/years" << endl; //debug
    
    //resize rho_d;
    rho_d.clear(); //bring rho_d to default
    rho_d.resize(Nt, vector<double> (Nz, 0.0));
    
    //Init with uniform distrubution
    for(int i=0; i < Nz; i++){
        rho_d[0][i] = 0.01 * rho_g(R, (i-Nz/2) * dz);}
    
    if((D * dt / (dz*dz)) > 0.5) {cerr << "Diffusion CFL > 1/2 for R = " << R << "AU,\talpha = " << alpha << "\t d = " << D * dt / (dz*dz) << endl; return; }
    
    
    //SOLVE
    for(int n=0; n < Nt-1; n++){ //iterate along time
        
        //Boundary conditions: no dust in the borders
        rho_d[n][0] = 0.0; rho_d[n][Nz-1] = 0.0;
        
        for(int i=1; i < Nz-1; i++){ //iterate along z
            
            
            //ADVECTION PART
            double adv_term = 0.0;
            
            if(abs(u_z(i) * dt / dz) > 1) {cerr << "Advection CFL > 1 for R = " << R << ", a_d = " << a_d << ", at z-step: " << i << "\ta = " << abs(u_z(i) * dt / dz) << endl; return; }
            
            //upwind method
            if (u_z(i) > 0) {
                adv_term = -(rho_d[n][i] * u_z(i) - rho_d[n][i-1] * u_z(i-1)) / dz; // backwards difference
            } else {
                adv_term = -(rho_d[n][i+1] *  u_z(i+1) - rho_d[n][i] * u_z(i)) / dz; // forward difference
            }

            
            //DIFFUSIVE PART
            double diffusive_term = 0.0;
            
            //Euler method for diffusive term
            //using dy/dx = (y_(i+1)-y_(i-1))/2dx
            
            //first approximate d(rho_d/rho_g)/dz = grad
            auto grad = [&, n, rho_gas](int i) {
                return (rho_d[n][i + 1] / rho_gas(i + 1) - rho_d[n][i - 1] / rho_gas(i - 1)) / (2 * dz);
            };
            //cout << "Grad(i=" << i << ") = " << grad(i) << endl; //debug
            
            //the dust condensate in the center after enough time so it becomes really narrow and near the center drho/dz goes to infinity so we lower the grad to grad_max for numerical stability
            double grad_max = 1e300;
            
            //clip grad between [-grad_max, + grad_max]
            auto reg_grad = [&, grad](int i) {
                return (abs(grad(i)) < grad_max) ? grad(i) : copysign(grad_max, grad(i));
            };

            //then approximate diffusive_term = D * d(rho_g*grad) using finite central difference
            if(i > 2 && i < Nz -2){//avoid out of bound
                diffusive_term = D * (rho_gas(i + 1) * reg_grad(i + 1) - rho_gas(i - 1) * reg_grad(i - 1)) / (2 * dz);}
            
            //DEBUG
            //if( (n % 1000 == 0 || (n % 100 == 0 && n < 1000) || (n % 10 == 0 && n < 100)) && (i<50 || i>950 || i % 100 == 0 || i == 498 || i == 499 || i == 501 || i == 502)){//cout << "reg_grad(n=" << n << ", i=" << i << ") = " << reg_grad(i);
                
               // cout << "(n=" << n << ", i=" << i << ")" << "\tadv_term = " << adv_term << "\tdiffusive_term = " << diffusive_term << "\trho_d = " << rho_d[n][i] << endl;}
            
            //|| (i>85 && i<105) || (i>885 && i<915)
             
            // Add new rho term
            rho_d[n+1][i] = rho_d[n][i] + dt * (adv_term + diffusive_term);
        }
    }
    //reset values of dt and Nt
    dt = 100;
    Nt = static_cast<int>(ceil(Lt/dt)); // number of time steps
}


void debug_evolution() {
    dt = choose_dt(R, alpha, a_d);
    Nt = static_cast<int>(ceil(Lt/dt));
    
    cout << "Omega = " << Omega(R) << " [1/s]" << endl;
    cout << "c_s = " << c_s(R) << " [m/s]" << endl;
    cout << "H = " << H(R) << " [m]" << endl;
    cout << "z_i [H(R)]= ";
    for (int i = 0; i < Nz; i += Nz / 10) { // Print out every tenth z_i
        cout << (i-Nz/2) * dz << " ";} cout << Nz/2 * dz << endl;
    
    cout << "rho_g(R, z_i) [rho_g_0] = ";
    for (int i = 0; i < Nz; i += Nz / 10) { // Print out rho_g for every tenth z_i
        cout << rho_g(R, (i-Nz/2) * dz) << " ";} cout << rho_g(R, Nz/2 * dz) << endl;
    
    cout << "t_s(R, z_i, a_d) [years] = ";
    for (int i = 0; i < Nz; i += Nz / 10) { // Print out rho_g for every tenth z_i
        cout << t_s(R, (i-Nz/2) * dz, a_d) / year << " ";} cout << t_s(R, Nz/2 * dz, a_d) / year << endl;
    
    cout << "u(z_i) [H(R)/year] = ";
    for (int i = 0; i < Nz; i += Nz / 10) { // Print out velocity u for every tenth z_i
        cout << u(R, (i-Nz/2) * dz, a_d) << " ";} cout << u(R, Nz/2 * dz, a_d) << endl;
    
    cout << "a(z_i) = ";
    for (int i = 0; i < Nz; i += Nz / 10) { // Print out velocity u for every tenth z_i
        cout << u(R, (i-Nz/2) * dz, a_d) * dt / dz << " ";} cout << u(R, Nz/2 * dz, a_d) * dt / dz << endl;
    
    
    for (int n = 0; n < Nt; n += Nt / 10) { // Print every tenth time step
        cout << "Time=  " << n*dt << ":\n";
        cout << "Rho_d = ";
        for (int i = 0; i < Nz; i += Nz / 10) { // Print every tenth spacial step
            cout << rho_d[n][i] << " ";}
        cout << rho_d[n][Nz-1] << endl;
    }
    //reset values of dt and Nt
    dt = 100;
    Nt = static_cast<int>(ceil(Lt/dt)); // number of time steps
}

//Studies evolution of dust density
void study_dust_evolution(double R_chosen, double alpha_chosen, string name){
    R = R_chosen; alpha = alpha_chosen;
    solve_dyn(R, alpha);
    
    cout << "\n\nValues for R = " << R << " AU" << ",\t alpha = " << alpha << ",\t a_d = " << a_d << "\n";
    debug_evolution();
    string filename = name + ".root";
    const char* fname = filename.c_str();
    save_rho(rho_d, fname);
    //plot_rho();
}


//Finds index at wich local dust is more than local gas
int check_form(vector<vector<double>> rho_d, double R_chosen) { //HIGHLY UNEFFICIENT :(
    R = R_chosen;
    dt = choose_dt(R, alpha, a_d);
    Nt = static_cast<int>(ceil(Lt/dt));
    for (int n = 0; n < Nt-1; n++) {
        for (int i = 0; i < Nz; i++) { 
            if (rho_d[n][i] > 100 * rho_d[0][Nz/2]) {return n;
                cout << "R = " << R << "\tform_time = " << n * dt << "years\n";}
        }}
    cerr << "I could not detect formation of planetesimals at R = " << R << ", alpha = " << alpha << ", a_d = " << a_d << endl;
    return -1;
}


void study_form_time(){
    vector<vector<double>> plan_sim;
    double dR = (R_max-R_min)/NR; alpha=0.0;
    
    cout << "\n\nHow long does it take for a planetsimal to form? ...I am slow please wait :)" << endl;
    for(int j=0; j < NR; j++){ //Check time formation of platesimals for every R
        R = R_min + j * dR;
        solve_dyn(R, alpha);
        dt = choose_dt(R, alpha, a_d);
        double form_time = dt * check_form(rho_d, R);
        plan_sim.push_back({R, form_time});
    }
    
   //Print rays (UA) and time to form planetsimals
    /*
    for (int i = 0; i < NR; i += NR / 10) { // Print out R for every tenth R_i
        cout << "R = " << plan_sim[i][0] << " AU\t time = " << plan_sim[i][1]  << " years"<< endl;} */

    scatter_plot_spline(plan_sim, "Time of plantesimals formation in relation to radius with spline fit;R[AU]; Time [years]");
    
}

//return the z-index (assuming symmetry in z) in which the density is non zero (numerically not super small respect to the center)
int check_height(vector<vector<double>> rho_d){
    dt = choose_dt(R, alpha, a_d);
    Nt = static_cast<int>(ceil(Lt/dt));
    for (int n = (int)1e6/dt; n < Nt; n++) { //iterate time starting from 1Myr
        for (int i = 0; i < Nz; i++) { //iterate in z
            if (rho_d[n][i] > 1e-3 * pow(R,-3)) {return i;}

        }}
    cerr << "Height of the disk could not be checked for R = " << R << ", alpha = " << alpha << endl; return -1;
}

//extimates alpha assuming that at R=100 planesimals take more than 1Myr to form
void study_diffusion(){
    R = 100;    //STUDY DIFFUSION AT R = 100 AU
    int Nalpha = 10; //steps of the iteration
    double alpha_0 = 1e-13; double height = 0.0;
    vector<vector<double>> alpha_h;
    cout << "\n\nI'm studying the diffusion at R = " << R << " AU\t from alpha= " << alpha_0 << " to alpha= " << alpha_0*pow(10,Nalpha) << "\t...\n";
    
    //first check alpha = 0;
    alpha = 0; solve_dyn(R, alpha);
    height = H(R) * dz * (Nz/2 - check_height(rho_d)) / AU; //Height is in AU for confrontation with the picture
    cout << "For aplha = " << alpha << "\tdisk height = " << height << "AU" << endl;
    
    for(int i=0 ; i < Nalpha; i++){ //start from alpha_0
        alpha = pow(10, i) * alpha_0; //multiply by 10 each iteration
        solve_dyn(R, alpha);
        height = H(R) * dz * (Nz/2 - check_height(rho_d)) / AU; //Height is in AU for confrontation with the picture
        cout << "For aplha = " << alpha << "\tdisk height = " << height << "AU" << endl;
        alpha_h.push_back({alpha, height});
    }
    
    scatter_plot(alpha_h, "Growth of disk height increasing the coefficient of diffusion; alpha; height [AU]");
}


//gives an extimation of dimension of the dust particle to form a planetsimal at ray R
double check_dust_grow(double R, double a_0){
    alpha = 1e-4;
    
    int N_ad = 10; double dad = 0.1 * a_0; //Number of steps and value of a_d-step
    double a_est = 0; //stores estimated value
    
    for(int i=0 ; i < N_ad; i++){
        a_d = a_0 - i*dad; //multiply by 10 each iteration
        
        dt = choose_dt(R, alpha, a_d);
        Nt = static_cast<int>(ceil(Lt/dt));
        
        //reinitialise rho_d to a vector of the right dimension
        rho_d.clear();
        rho_d.resize(Nt, vector<double> (Nz, 0.0));
        
        solve_dyn(R, alpha);
        
        double form_time = dt * check_form(rho_d, R);
        if(form_time < 0) {a_est = a_0 - (i-1) * dad; break;} //the a_d before was the last one valid
        
    }
    
    dt = 100;
    Nt = static_cast<int>(ceil(Lt/dt));
    
    rho_d.clear(); //bring rho_d to default
    rho_d.resize(Nt, vector<double> (Nz, 0.0));
    return a_est;
}


#endif /* simulation_h */
