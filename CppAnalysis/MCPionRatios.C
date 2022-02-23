
#include <iostream>
#include <fstream>
#include <string>

#include <vector>
#include <cstdlib>
#include <chrono>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include <time.h>

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TVector3.h>

#include "../Auxiliary/bank.h"
#include "../Auxiliary/BBand.h"
#include "../Auxiliary/BEvent.h"
#include "../Auxiliary/constants.h"
#include "../Auxiliary/bandhit.cpp"
#include "../Auxiliary/clashit.cpp"
#include "../Auxiliary/genpart.cpp"

TH1D* mmiss_hist;
TH1D* Q2_hist;
TH1D* xB_hist;
TH1D* Y_hist;
TH1D* W_hist;
TH1D* Vez_hist;
TH1D* eP_hist;
TH1D* theta_e_hist;
TH1D* phi_e_hist;

TH1D* p_pi_hist;
TH1D* theta_pi_hist;
TH1D* phi_pi_hist;
TH1D* phi_trento_hist;
TH1D* pt_hist;
TH1D* z_hist;

TH2D* Q2_xb_hist;
TH2D* Q2_W_hist;
TH2D* xb_W_hist;
TH2D* E_theta_e_hist;
TH2D* E_phi_e_hist;
TH2D* theta_phi_e_hist;

TH2D* theta_p_e_hist;
TH2D* theta_p_pi_hist;
TH2D* z_pt_hist;
TH2D* pt_x_hist;

TH3D* xb_z_pt_hist;

TH1D* e_dc_s_hist;
TH1D* pi_dc_s_hist;
TH1D* pi_no_hist;

TH2D* e_xy_hist_S[6];
TH2D* e_theta_phi_hist_S[6];
TH2D* pi_theta_phi_hist_S[6];
TH2D* pi_p_theta_hist_S[6];
TH2D* z_pt_hist_S[6];

TF1* theta_momentum_cut[6][2][2];

// sector, pion, top/bottom, parameter
const double theta_momentum_cut_params[6][2][2][2] = 
    {{{{4.48444148, 9.17214918}, {6.13646294, 33.37018884}}, {{7.20468119, 14.62704809}, {6.93859495, 32.13113786}}},
    {{{4.3350149, 9.7613775}, {6.25178768, 32.71433338}}, {{7.04899446, 14.82707565}, {6.58857115, 33.00689677}}},
    {{{4.52162786, 8.48187808}, {6.17111499, 32.53055689}}, {{6.99496389, 15.32228833}, {6.64405866, 33.6622875 }}},
    {{{4.27388498, 9.01677273}, {6.17205936, 32.62168359}}, {{7.01765029, 14.82602243}, {6.35632293, 34.36439859}}},
    {{{4.02349973, 9.73539597}, {6.05465832, 33.12697781}}, {{6.79515088, 15.30515795}, {6.58713963, 33.28105369}}},
    {{{4.2116445,  9.68950529}, {5.41861825, 34.35713488}}, {{7.21262506, 14.53342703}, {6.96618552, 31.76762161}}}};

using namespace std;

void   InitializeHists();
void   InitializeCuts();
double sq(double x) {return x*x;}

void MCPionRatios(Int_t NeventsMax=-1, Int_t pionCode=211, Int_t neutrons=0, TString infilename="", TString outfileName="") {
    InitializeHists();
    InitializeCuts();

    TTree *tree;

    if (neutrons == 0) {
        TFile *infile      = new TFile(infilename);
        tree               = (TTree*)infile->Get("tree");
    }
    else if (neutrons == 1) {
        TFile *infile      = new TFile(infilename);
        tree               = (TTree*)infile->Get("T");
    }

    TFile *outfile = new TFile(outfileName, "RECREATE");

    TLorentzVector* e = 0;
    TVector3* Ve = 0;
    TClonesArray* pi = new TClonesArray("TLorentzVector");
    TClonesArray* Vpi = new TClonesArray("TVector3");
    Double_t Ebeam, xB, Q2, omega, W;
    Double_t e_DC_x[3], e_DC_y[3], e_DC_sector;
    Double_t pi_DC_sector[20];
    Int_t Npis;

    tree->SetBranchAddress("Ebeam", &Ebeam);
    tree->SetBranchAddress("xB", &xB);
    tree->SetBranchAddress("Q2", &Q2);
    tree->SetBranchAddress("omega", &omega);
    tree->SetBranchAddress("W", &W);

    tree->SetBranchAddress("e", &e);
    tree->SetBranchAddress("Ve", &Ve);
    tree->SetBranchAddress("pi", &pi);
    tree->SetBranchAddress("Vpi", &Vpi);

    tree->SetBranchAddress("e_DC_x", &e_DC_x);
    tree->SetBranchAddress("e_DC_y", &e_DC_y);
    tree->SetBranchAddress("e_DC_sector", &e_DC_sector);
    tree->SetBranchAddress("pi_DC_sector", &pi_DC_sector);

    if (pionCode == 211) {
        tree->SetBranchAddress("Npips", &Npis);
    }
    else {
        tree->SetBranchAddress("Npims", &Npis);
    }

    const TVector3 l(0, 0, 10.2);

    double y;
    TLorentzVector* ppi;
    TVector3* Vxpi;

    for (int sEvent = 0; sEvent < tree->GetEntries(); sEvent++) {
        if (sEvent == NeventsMax) break;

        tree->GetEntry(sEvent);
        Ebeam = 0;
        xB = 0;
        Q2 = 0;
        omega = 0;
        W = 0;
        y = 0;
        e_DC_sector = 0;

        tree->GetEntry(sEvent);

        TVector3 pbeam(0, 0, Ebeam);
        TVector3 lp = e->Vect();
        y = omega/Ebeam;
        
        if (Q2 <= 1 || W <= 2 || y >= 0.75 || e->Theta() <= 0.0872665 || e->Theta() >= 0.610865) continue;

        pi_no_hist->Fill(Npis);

        for (int piNo = 0; piNo < Npis; piNo++) {
            ppi = (TLorentzVector*)pi->ConstructedAt(piNo);
            Vxpi = (TVector3*)Vpi->ConstructedAt(piNo);

            TVector3 q = l - lp;
            TVector3 Ph = ppi->Vect();
            double mmiss = sqrt(sq(omega + mN - ppi->E()) - (q - Ph).Mag2());
            double z = ppi->E() / omega;
            int sector = (int)pi_DC_sector[piNo]-1;
            double theta = ppi->Theta()*180/M_PI;

            // if (Vxpi->z() == 0) continue;
            if (ppi->P() < 1.25 || ppi->P() > 5 || ppi->Theta() <= 0.0872665 || ppi->Theta() >= 0.610865) continue;
            if (z <= 0.3) continue;
            if (sector < 0) continue;
            // if (ppi->E() / omega <= 0.3) continue;
            
            double pip_tp_min = theta_momentum_cut[sector][0][0]->Eval(ppi->P());
            double pip_tp_max = theta_momentum_cut[sector][0][1]->Eval(ppi->P());
            double pim_tp_min = theta_momentum_cut[sector][1][0]->Eval(ppi->P());
            double pim_tp_max = theta_momentum_cut[sector][1][1]->Eval(ppi->P());

            if (theta < pip_tp_min || theta < pim_tp_min || theta > pip_tp_max || theta > pim_tp_max) continue;

            mmiss_hist->Fill(mmiss);

            if (mmiss <= 1.5) continue;

            TVector3 qu = q.Unit();

            TVector3 quxl = qu.Cross(l);
            TVector3 quxPh = qu.Cross(Ph);
            TVector3 lxPh = l.Cross(Ph);

            double phi_tr_y = lxPh.Dot(qu) / (quxl.Mag() * quxPh.Mag());
            double phi_tr_x = quxl.Unit().Dot(quxPh.Unit());
            double phi_tr_tan = atan2(phi_tr_y, phi_tr_x);
            phi_tr_tan *= 180 / M_PI;
            if (phi_tr_tan < 0) phi_tr_tan += 360;
            double phi_tr = ((q.Cross(l).Dot(Ph)) / abs(q.Cross(l).Dot(Ph))) * (q.Cross(l).Dot(q.Cross(Ph)) / q.Cross(l).Mag() * q.Cross(Ph).Mag());

            double pt = ppi->Perp(q);

            theta_p_e_hist->Fill(e->P(), e->Theta() * 180 / M_PI);

            Q2_hist->Fill(Q2);
            xB_hist->Fill(xB);
            Y_hist->Fill(y);
            W_hist->Fill(W);
            Vez_hist->Fill(Ve->z());
            eP_hist->Fill(e->P());
            theta_e_hist->Fill(e->Theta()*180/3.14);
            phi_e_hist->Fill(e->Phi()*180/3.14);

            Q2_xb_hist->Fill(xB, Q2);
            Q2_W_hist->Fill(W, Q2);
            xb_W_hist->Fill(W, xB);
            E_theta_e_hist->Fill(e->Theta()*180/3.14, e->E());
            E_phi_e_hist->Fill(e->Phi()*180/3.14, e->E());
            theta_phi_e_hist->Fill(e->Phi()*180/3.14, e->Theta()*180/3.14);

            theta_p_pi_hist->Fill(ppi->P(), ppi->Theta() * 180 / M_PI);
            p_pi_hist->Fill(ppi->P());
            theta_pi_hist->Fill(ppi->Theta() * 180 / M_PI);
            phi_pi_hist->Fill(ppi->Phi() * 180 / M_PI);
            pt_hist->Fill(pt);
            z_hist->Fill(ppi->E() / omega);
            z_pt_hist->Fill(pt, ppi->E() / omega);
            pt_x_hist->Fill(xB, pt);

            phi_trento_hist->Fill(phi_tr_tan);
     
            xb_z_pt_hist->Fill(xB, z, pt);

            e_dc_s_hist->Fill(e_DC_sector);
            pi_dc_s_hist->Fill(pi_DC_sector[piNo]);

            e_xy_hist_S[(int)e_DC_sector-1]->Fill(e_DC_x[0], e_DC_y[0]);
            e_theta_phi_hist_S[(int)e_DC_sector-1]->Fill(e->Theta()*180/M_PI, e->Phi()*180/M_PI);
            pi_theta_phi_hist_S[(int)pi_DC_sector[piNo]-1]->Fill(ppi->Theta()*180/M_PI, ppi->Phi()*180/M_PI);
            pi_p_theta_hist_S[(int)pi_DC_sector[piNo]-1]->Fill(ppi->P(), ppi->Theta()*180/M_PI);
            z_pt_hist_S[(int)pi_DC_sector[piNo]-1]->Fill(z, pt);
        }
    }

    outfile->cd();

    Q2_hist->Write();
    xB_hist->Write();
    Y_hist->Write();
    W_hist->Write();
    Vez_hist->Write();
    eP_hist->Write();
    theta_e_hist->Write();
    phi_e_hist->Write();
    p_pi_hist->Write();
    theta_pi_hist->Write();
    phi_pi_hist->Write();
    phi_trento_hist->Write();
    pt_hist->Write();
    mmiss_hist->Write();
    z_hist->Write();
    Q2_xb_hist->Write();
    Q2_W_hist->Write();
    xb_W_hist->Write();
    E_theta_e_hist->Write();
    E_phi_e_hist->Write();
    theta_phi_e_hist->Write();
    theta_p_e_hist->Write();
    theta_p_pi_hist->Write();
    z_pt_hist->Write();
    pt_x_hist->Write();
    xb_z_pt_hist->Write();
    e_dc_s_hist->Write();
    pi_dc_s_hist->Write();
    pi_no_hist->Write();

    for (int s = 1; s <= 6; s++) {
        e_xy_hist_S[s-1]->Write();
        e_theta_phi_hist_S[s-1]->Write();
        pi_theta_phi_hist_S[s-1]->Write();
        pi_p_theta_hist_S[s-1]->Write();
        z_pt_hist_S[s-1]->Write();
    }

    outfile->Close();
    
    cout << "Complete" << endl;
}

void InitializeHists() {
    Q2_hist = new TH1D("Q2", "Q2; Q2 [GeV^2]", 100, 0, 12);
    xB_hist = new TH1D("xB", "xB; xB", 100, 0, 0.8);
    Y_hist = new TH1D("Y", "Y; Y", 100, 0, 1);
    W_hist = new TH1D("W", "W; W[GeV]", 100, 1, 5);
    Vez_hist = new TH1D("VeZ", "Electron Longitudinal Vertex; Ve_z", 100, -15, 15);
    eP_hist = new TH1D("eP", "Electron Momentum; Momentum [GeV/c]", 100, 0, 10);
    theta_e_hist = new TH1D("theta_e", "Electron Theta; Theta [degrees]", 180, 0, 40);
    phi_e_hist = new TH1D("phi_e", "Electron Phi; Phi [degrees]", 360, -180, 180);

    p_pi_hist = new TH1D("p", "Pion Momentum; p [GeV]", 100, 0, 6);
    theta_pi_hist = new TH1D("theta", "Pion angle; theta [deg]", 160, 0, 40);
    phi_pi_hist = new TH1D("phi", "Pion azimuthal angle; phi [deg]", 360, -180, -180);
    phi_trento_hist = new TH1D("phi_t", "Pion Phi Trento; phi trento [deg]", 360, 0, 360);
    pt_hist = new TH1D("pt", "Pion Transverse Momentum; Momentum [GeV/c]", 100, 0, 1.4);
    mmiss_hist = new TH1D("mmiss", "Missing mass; M_miss [GeV]", 200, 0, 4);
    z_hist = new TH1D("z", "z; z", 100, 0, 1);

    Q2_xb_hist = new TH2D("Q2_xb", "Q2 vs xB; xB; Q2 [GeV^2]", 100, 0, 1, 100, 0, 12);
    Q2_W_hist = new TH2D("Q2_W", "Q2 vs W; W [GeV]; Q2 [GeV^2]", 100, 1, 5, 100, 0, 12);
    xb_W_hist = new TH2D("xb_W", "xB vs W; W [GeV]; xB", 100, 1, 5, 100, 0, 1);
    E_theta_e_hist = new TH2D("E_theta", "E vs theta; theta [deg]; E [GeV]", 50, 0, 50, 100, 0, 12);
    E_phi_e_hist = new TH2D("E_phi", "E vs phi; phi [deg]; E [GeV]", 360, -180, 180, 100, 0, 12);
    theta_phi_e_hist = new TH2D("theta_phi", "theta vs phi; phi [deg]; theta [deg]", 360, -180, 180, 40, 0, 40);

    theta_p_e_hist = new TH2D("theta_p_e", "Theta vs momentum - electron; p_e [GeV]; theta_e [deg]", 100, 2, 9, 35, 5, 40);
    theta_p_pi_hist = new TH2D("theta_p_pi", "Theta vs momentum - pion; p_pi [GeV]; theta_pi [deg]", 100, 1, 5, 35, 5, 40);
    z_pt_hist = new TH2D("z_pi", "Z vs Pt; Pt [GeV]; z", 100, 0, 1.6, 100, 0, 1);
    pt_x_hist = new TH2D("pt_x", "Pt vs x; x; Pt[GeV]", 100, 0, 1, 100, 0, 1.6);

    xb_z_pt_hist = new TH3D("xb_z_pt", "xB vs z vs Pt", 100, 0, 1, 100, 0.3, 1.6, 100, 0, 1.6);

    e_dc_s_hist = new TH1D("e_dc_s", "e_dc_s", 7, 0, 6);
    pi_dc_s_hist = new TH1D("pi_dc_s", "pi_dc_s", 7, 0, 6);
    pi_no_hist = new TH1D("pi_no", "pi_no", 21, 0, 20);

    for (int s = 1; s <= 6; s++) {
        e_xy_hist_S[s-1] = new TH2D(Form("e_xy_S_%d",s), Form("sector %d",s), 100, -300, 300, 100, -300, 300);
        e_theta_phi_hist_S[s-1] = new TH2D(Form("e_theta_phi_S_%d",s), Form("sector %d",s), 60, 5, 35, 360, -180, 180);
        pi_theta_phi_hist_S[s-1] = new TH2D(Form("pi_theta_phi_S_%d",s), Form("sector %d",s), 30, 5, 35, 360, -180, 180);
        pi_p_theta_hist_S[s-1] = new TH2D(Form("pi_p_theta_S_%d",s), Form("sector %d",s), 100, 1, 5, 160, 5, 35);
        z_pt_hist_S[s-1] = new TH2D(Form("z_pt_S_%d",s), Form("sector %d",s), 100, 0.3, 1.6, 100, 0, 1);
    }
}

void InitializeCuts() {
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                theta_momentum_cut[i][j][k] = new TF1(Form("pion_cut_%d_%d_%d", i, j, k),Form("%f+%f/TMath::Power(x,1.)", theta_momentum_cut_params[i][j][k][0], theta_momentum_cut_params[i][j][k][1]),0,5.);
            }
        }
    }
}
