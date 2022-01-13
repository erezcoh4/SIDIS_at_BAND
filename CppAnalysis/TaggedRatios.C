
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

TString DataPath             = "/volatile/clas12/users/akiral/BAND/";
std::string GoodRunsFilename = "/work/clas12/users/akiral/SIDIS_at_BAND/macros/runlists/good_runs_10-2.txt";

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

TH1D* theta_pi_root_hist;
TH1D* phi_pi_root_hist;

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

using namespace std;

void   InitializeHists();
double sq(double x) {return x*x;}

void TaggedRatios(Int_t RunNumber=6420, Int_t NeventsMax=-1, Int_t pionCode=211, Int_t neutrons=1) {
    InitializeHists();

    TString pionType;
    if (pionCode == 211) {
        pionType = "plus";
    }
    else if (pionCode == -211) {
        pionType = "minus";
    }
    else {
        cout << "Unknown particle code for pion: " << pionCode << endl;
        return;
    }

    char RunNumberStr[20];
    sprintf( RunNumberStr, "00%d", RunNumber );

    TString outfileName;
    TTree *tree;

    if (neutrons == 0) {
        TString infilename = DataPath + "SIDIS_skimming/skimmed_SIDIS_inc_" + RunNumberStr + "_e_pi" + pionType + ".root";
        outfileName        = DataPath + "histograms/pion_ratio_" + RunNumberStr + "_" + pionType + ".root";

        TFile *infile      = new TFile(infilename);
        tree               = (TTree*)infile->Get("tree");
    }
    else if (neutrons == 1) {
        TString infilename = DataPath + "merged_SIDIS_and_BAND_skimming/skimmed_SIDIS_and_BAND_inc_" + RunNumberStr + "_e_pi" + pionType + "_n.root";
        outfileName        = DataPath + "histograms/tagged/tagged_ratio_" + RunNumberStr + "_" + pionType + ".root";

        TFile *infile      = new TFile(infilename);
        tree               = (TTree*)infile->Get("T");
    }

    TFile *outfile = new TFile(outfileName, "RECREATE");

    clashit* eHit = 0;
    TLorentzVector* e = 0;
    TVector3* Ve = 0;
    Double_t Ebeam, xB, Q2, omega, W, W2;
    Double_t pipx[20], pipy[20], pipz[20], pie[20];
    Int_t Npis;

    tree->SetBranchAddress("eHit", &eHit);
    tree->SetBranchAddress("Ebeam", &Ebeam);
    tree->SetBranchAddress("xB", &xB);
//    tree->SetBranchAddress("Q2", &Q2);
    tree->SetBranchAddress("omega", &omega);
//    tree->SetBranchAddress("W", &W);

    tree->SetBranchAddress("e", &e);
    tree->SetBranchAddress("Ve", &Ve);

    if (pionCode == 211) {
        tree->SetBranchAddress("Npips", &Npis);
        tree->SetBranchAddress("piplus_Px", &pipx);
        tree->SetBranchAddress("piplus_Py", &pipy);
        tree->SetBranchAddress("piplus_Pz", &pipz);
        tree->SetBranchAddress("piplus_E", &pie);
    }
    else {
        tree->SetBranchAddress("Npims", &Npis);
        tree->SetBranchAddress("piminus_Px", &pipx);
        tree->SetBranchAddress("piminus_Py", &pipy);
        tree->SetBranchAddress("piminus_Pz", &pipz);
        tree->SetBranchAddress("piminus_E", &pie);
    }

    const TVector3 l(0, 0, 10.2);

    double y;
    TLorentzVector* ppi;

    for (int sEvent = 0; sEvent < tree->GetEntries(); sEvent++) {
        if (sEvent == NeventsMax) break;

        tree->GetEntry(sEvent);
        Ebeam = 0;
        xB = 0;
        Q2 = 0;
        omega = 0;
        W2 = 0;
        y = 0;

        tree->GetEntry(sEvent);

        Q2 = eHit->getQ2();
        W2 = eHit->getW2();
        W = sqrt(W2);

        TVector3 pbeam(0, 0, Ebeam);
        TVector3 lp = e->Vect();
        y = omega/Ebeam;

        if (Q2 <= 1 || W <= 2 || y >= 0.75 || e->Theta() <= 0.0872665 || e->Theta() >= 0.610865) continue;

        for (int piNo = 0; piNo < Npis; piNo++) {
            ppi = new TLorentzVector(pipx[piNo], pipy[piNo], pipz[piNo], pie[piNo]);

            TVector3 q = l - lp;
            TVector3 Ph = ppi->Vect();
            double mmiss = sqrt(sq(omega + mN - ppi->E()) - (q - Ph).Mag2());

            if (ppi->P() < 1.25 || ppi->P() > 5 || ppi->Theta() <= 0.0872665 || ppi->Theta() >= 0.610865) continue;
            // if (ppi->E() / omega <= 0.3) continue;

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
    theta_pi_root_hist->Write();
    phi_pi_root_hist->Write();
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

    outfile->Close();
    
    cout << "Complete" << endl;
}

void InitializeHists() {
    Q2_hist = new TH1D("Q2", "Q2; Q2 [GeV^2]", 100, 0, 12);
    xB_hist = new TH1D("xB", "xB; xB", 100, 0, 0.8);
    Y_hist = new TH1D("Y", "Y; Y", 100, 0, 1);
    W_hist = new TH1D("W", "W; W[GeV]", 100, 1, 5);
    Vez_hist = new TH1D("VeZ", "Electron Longitudinal Vertex; Ve_z", 100, -20, 20);
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

    theta_pi_root_hist = new TH1D("theta_root", "Pion angle; theta [deg]", 101, -10000, 100);
    phi_pi_root_hist = new TH1D("phi_root", "Pion azimuthal angle; phi [deg]", 101, -10000, 100);

    Q2_xb_hist = new TH2D("Q2_xb", "Q2 vs xB; xB; Q2 [GeV^2]", 100, 0, 1, 100, 0, 10);
    Q2_W_hist = new TH2D("Q2_W", "Q2 vs W; W [GeV]; Q2 [GeV^2]", 100, 0, 5, 100, 0, 15);
    xb_W_hist = new TH2D("xb_W", "xB vs W; W [GeV]; xB", 100, 0, 5, 100, 0, 1);
    E_theta_e_hist = new TH2D("E_theta", "E vs theta; theta [deg]; E [GeV]", 50, 0, 50, 100, 0, 10);
    E_phi_e_hist = new TH2D("E_phi", "E vs phi; phi [deg]; E [GeV]", 360, -180, 180, 100, 0, 10);
    theta_phi_e_hist = new TH2D("theta_phi", "theta vs phi; phi [deg]; theta [deg]", 360, -180, 180, 50, 0, 50);

    theta_p_e_hist = new TH2D("theta_p_e", "Theta vs momentum - electron; p_e [GeV]; theta_e [deg]", 100, 2, 9, 30, 5, 35);
    theta_p_pi_hist = new TH2D("theta_p_pi", "Theta vs momentum - pion; p_pi [GeV]; theta_pi [deg]", 100, 1, 5, 30, 5, 35);
    z_pt_hist = new TH2D("z_pi", "Z vs Pt; Pt [GeV]; z", 100, 0, 1.6, 100, 0, 1);
    pt_x_hist = new TH2D("pt_x", "Pt vs x; x; Pt[GeV]", 100, 0, 1, 100, 0, 1.6);
}
