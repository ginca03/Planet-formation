
#ifndef root_h
#define root_h

//Root Libraries
#include <TCanvas.h>
#include <TGraph.h>
#include <TSpline.h>
#include <TH2F.h>
#include <TApplication.h>
#include <TFile.h>
#include <TTree.h>

//USE WITH SOLVE_DYN
void plot_rho() {
    TCanvas* c_rho = new TCanvas("canvas", "Dust Density Evolution", 800, 600);
    TH2F* hist = new TH2F("hist", "Dust Density [kg/m^3] Evolution;z [H(R)];t [years]", Nz, -H(R) / 2, H(R) / 2, Nt, 0, Nt * dt );
    
    dt = choose_dt(R, alpha, a_d); // choose the timesteps dynamically
    Nt = static_cast<int>(ceil(Lt/dt));

    
    //Fill hist
    for (int n = 0; n < Nt; n++) {
        for (int i = 0; i < Nz; i++) {
            double z = (i - Nz / 2) * dz;  // z in H(R)
            double t = n * dt;   // time in years

            // Trova l'indice del bin
            int bin_x = hist->GetXaxis()->FindBin(z);
            int bin_y = hist->GetYaxis()->FindBin(t);

            // Imposta il contenuto direttamente
            double density_value = rho_d[n][i]; // density in rho_g_0;
            hist->SetBinContent(bin_x, bin_y, density_value);
        }
    }
    
    //Cosmetics
    hist->SetStats(0);  hist->SetContour(100);
    
    hist->Draw("COLZ");
    c_rho->Update();
    
    //reset values of dt and Nt
    dt = 100;
    Nt = static_cast<int>(ceil(Lt/dt)); // number of time steps
}


void save_rho(const vector<vector<double>>& rho_d, const char* filename) {
    TFile* file = new TFile(filename, "RECREATE");
    TTree* tree = new TTree(filename, "Dust density evolution");
    
    dt = 0.8 * dz / abs(u(R, 2, a_d)); // choose the timesteps dynamically
    Nt = static_cast<int>(ceil(Lt/dt));
    
    double D = D = alpha * pow(c_s(R),2) / Omega(R) * (year/(H(R)*H(R))); // D is in H(R)^2/years
    if((D * dt / (dz*dz)) > 0.5){
        dt = (dz*dz) / (3*D);
        Nt = static_cast<int>(ceil(Lt/dt));}
    

    double z, t, rho;  // Temporary variables
    tree->Branch("z", &z, "z/D");
    tree->Branch("t", &t, "t/D");
    tree->Branch("rho", &rho, "rho/D");

    for (int n = 0; n < Nt; n++) {
        t = n * dt;  // time in years
        for (int i = 0; i < Nz; i++) {
            z = (i - Nz / 2) * dz;  // z in H(R)
            rho = rho_d[n][i];  // density in rho_g_0
            tree->Fill();
        }
    }

    tree->Write();
    file->Close();
    cout << "I saved rho_d in file: " << filename << endl;
    
    //reset values of dt and Nt
    dt = 100;
    Nt = static_cast<int>(ceil(Lt/dt)); // number of time steps
}




//USE WITH SOLVE_STAT

void scatter_plot(const vector<vector<double>> &data, string s) {
    int n = data.size();
    vector<double> x(n), y(n);
    for (int i = 0; i < n; ++i) {
        x[i] = data[i][0]; y[i] = data[i][1];}

    TGraph *graph = new TGraph(n, x.data(), y.data());

    // Cosmetics
    graph->SetTitle(s.c_str());
    graph->SetMarkerStyle(20); graph->SetMarkerSize(1.0); graph->SetMarkerColor(kBlue);

    TCanvas *canvas_scatter = new TCanvas("canvas_scatter", "Scatter Plot", 800, 600);
    graph->Draw("AP");
    canvas_scatter->Update();
}

void scatter_plot_spline(const vector<vector<double>> &data, string s) {
    int n = data.size();
    vector<double> x(n), y(n);
    for (int i = 0; i < n; ++i) {
        x[i] = data[i][0];
        y[i] = data[i][1];
    }

    TGraph *graph = new TGraph(n, x.data(), y.data());

    // Cosmetics
    graph->SetTitle(s.c_str());
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1.0);
    graph->SetMarkerColor(kBlue);

    // Spline interpolation
    TSpline3 *spline = new TSpline3("spline", graph);
    spline->SetLineColor(kRed);
    spline->SetLineWidth(2);

    TCanvas *canvas_scatter_spline = new TCanvas("canvas_scatter_spline", "Scatter Plot with spline", 800, 600);
    graph->Draw("AP");
    spline->Draw("same");
    canvas_scatter_spline->Update();
}


#endif /* root_h */
