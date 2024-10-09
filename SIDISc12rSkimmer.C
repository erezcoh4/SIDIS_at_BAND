// last edit June-22, 2023

#ifdef __CINT__
#pragma link C++ class std::vector<TVector3>+;
#pragma link C++ class vector<TVector3>+;
#pragma link C++ class std::vector<TLorentzVector>+;
#pragma link C++ class vector<TLorentzVector>+;
#endif

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TH1.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include "clas12reader.h"
#include "Auxiliary/DCfid_SIDIS.cpp"
#include "Auxiliary/csv_reader.h"
#include "Auxiliary/SIDISatBAND_auxiliary.cpp"

using namespace clas12;
SIDISatBAND_auxiliary aux;

// Results in CSV file
TString csvheader = ( (TString)"status,runnum,evnum,beam_helicity,"
                     +(TString)"e_P,e_Theta,e_Phi,e_Vz,"
                     +(TString)"pi_P,pi_Theta,pi_Phi,pi_Vz,"
                     +(TString)"Q2,xB,omega,y,"
                     +(TString)"e_DC_sector,pi_DC_sector,"
                     +(TString)"e_Theta_qFrame,e_Phi_qFrame,"
                     +(TString)"pi_qFrame_Theta,pi_qFrame_Phi,"
                     +(TString)"pi_qFrame_pT,pi_qFrame_pL,"
                     +(TString)"Zpi,Zpi_LC,"
                     +(TString)"W,M_x,"
                     +(TString)"xF,eta_pi,"
                     +(TString)"W_d,M_x_d,"
                     +(TString)"q,qStar,"
                     );

// addition to csv for GEMC simulations
TString csvheader_GEMCaddition = ( (TString)"e_P_g,e_Theta_g,e_Phi_g,e_Vz_g,"
                                  +(TString)"pi_P_g,pi_Theta_g,pi_Phi_g,pi_Vz_g,"
                                  +(TString)"Q2_g,xB_g,omega_g,y_g,"
                                  );

std::vector<int> csvprecisions = {0,0,0,0,9,9,9,9,9,9,9,9,9,9,9,9,0,0,};



// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
// start clock
auto start = std::chrono::high_resolution_clock::now();
// declare methods
TVector3                GetParticleVertex (clas12::region_part_ptr rp);
void                     SetLorentzVector (TLorentzVector &p4, clas12::region_part_ptr rp);
void                      OpenOutputFiles (TString csvfilename);
void                     CloseOutputFiles (TString OutDataPath, TString outfilename);
void                      SetOutputTTrees ();
bool   CheckIfElectronPassedSelectionCuts (Double_t e_PCAL_x, Double_t e_PCAL_y,
                                           Double_t e_PCAL_W,Double_t e_PCAL_V,
                                           Double_t e_E_PCAL,
                                           Double_t e_E_ECIN, Double_t e_E_ECOUT,
                                           TLorentzVector e,
                                           TVector3 Ve,
                                           Double_t e_DC_sector,
                                           Double_t e_DC_x[3],
                                           Double_t e_DC_y[3],
                                           Double_t e_DC_z[3],
                                           int torusBending);
bool      CheckIfPionPassedSelectionCuts (TString pionCharge, // "pi+" or "pi-"
                                          Double_t DC_sector,
                                          Double_t DC_x[3], Double_t DC_y[3], Double_t DC_z[3],
                                          Double_t chi2PID, Double_t p,
                                          TVector3 Ve,      TVector3 Vpi,
                                          int fdebug);
int                       GetBeamHelicity (event_ptr p_event, int runnum, int fdebug);
double                      GetBeamEnergy (int fdebug);
void                InitializeFileReading (int NeventsMax,int c12Nentries, int fdebug);
void                  InitializeVariables ();
void                      OpenResultFiles (TString outfilepath, TString outfilename );
void           ExtractElectronInformation (int fdebug);
void              ExtractPionsInformation (int fdebug);
void               ExtractPipsInformation (int pipsIdx, int fdebug );
void               ExtractPimsInformation (int pimsIdx, int fdebug );
void                ExtractKpsInformation (int kIdx, int fdebug );
void                ExtractKmsInformation (int kIdx, int fdebug );
void            ComputeElectronKinematics ();
void                ComputePionKinematics (TLorentzVector pi, TLorentzVector pi_qFrame);
void                   WriteEventToOutput (int fdebug);
void                        FinishProgram (TString outfilepath, TString outfilename);
void                   GetParticlesByType (int evnum, int fdebug );
void              Stream_e_pi_line_to_CSV (TString pionCharge, int piIdx,
                                           bool passed_cuts_e_pi,
                                           bool passed_cuts_e_pi_kinematics,
                                           int fdebug );
void               Stream_e_K_line_to_CSV (TString KCharge, int KIdx,
                                           bool passed_cuts_e_K,
                                           bool passed_cuts_e_K_kinematics,
                                           int fdebug );
TVector3            RotateVectorTo_qFrame (TVector3 V);
void                        MoveTo_qFrame (int fdebug);
void                          SetDataPath (TString fDataPath, Double_t fEbeam) ;
void                             SetSimPi (TString fSimPi);
void                          SetSkimming (TString fSkimming) ;
void                         SetInclusive ( int fInclusive );
void                             SetEbeam ( double fEbeam );
void                              SetIsMC ( bool fIsMC = false );
void                         SetVerbosity ( int _fdebug_ = 0 );
// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.

// globals
TString Skimming = "", DataPath = "", prefix = "", SimPi = "", SimK = "";
TString SpecificFilePath = "", SpecificFilename = "";
auto db = TDatabasePDG::Instance();
double        Pe_phi, q_phi, q_theta; // "q-frame" parameters
double                   Zpi, Zpi_LC;
double                   ZK,  ZK_LC;

bool    ePastCutsInEvent = false,       EventPassedCuts = false;
bool pipsPastCutsInEvent = false, eepipsPastCutsInEvent = false;
bool pimsPastCutsInEvent = false, eepimsPastCutsInEvent = false;
bool  KpsPastCutsInEvent = false,  eeKpsPastCutsInEvent = false;
bool  KmsPastCutsInEvent = false,  eeKmsPastCutsInEvent = false;

bool IsMC = false; // GEMC simulations
// meta-data
int fdebug = 0;
int torusBending = -1; // -1 for In-bending, +1 for Out-bending
int DC_layers[3] = {6,18,36};// Region 1 is denoted at DC detector 6, Region 2 is denoted 18, Region 3 - as 36
int DC_layer, runnum=0, evnum=0, beam_helicity=0;
// helicity of the electron +1 along the beam and -1 opposite to it
int status, NeventsMaxToProcess, Nevents_processed, Nevents_passed_e_cuts;
int Nevents_passed_pips_cuts,               Nevents_passed_e_pips_cuts;
int Nevents_passed_e_pips_kinematics_cuts,  Nevents_passed_pims_cuts;
int Nevents_passed_e_pims_cuts,             Nevents_passed_e_pims_kinematics_cuts;
int Nevents_passed_Kps_cuts,                Nevents_passed_e_Kps_cuts;
int Nevents_passed_e_Kps_kinematics_cuts;
int Nevents_passed_Kms_cuts,                Nevents_passed_e_Kms_cuts;
int Nevents_passed_e_Kms_kinematics_cuts;
int inclusive; // tag to look at inclusive run - all the events with no selection

// number of particles per event
int Ne, Nn, Np, Npips, Npims, NKps, NKms, Ngammas, Nd;

// leading electron
// electron energy deposit in PCAL [GeV], in ECAL_in [GeV], in ECAL_out [GeV]...
double e_E_PCAL, e_E_ECIN, e_E_ECOUT, e_PCAL_W, e_PCAL_V, e_PCAL_x, e_PCAL_y, e_PCAL_z;
double e_PCAL_sector, e_DC_sector, e_DC_Chi2N, e_DC_x[3], e_DC_y[3], e_DC_z[3];

// positive pions
bool pipsPastSelectionCuts[NMAXPIONS], eepipsPastKinematicalCuts[NMAXPIONS];
int            pips_region[NMAXPIONS];
double        pips_chi2PID[NMAXPIONS];
double         pips_PCAL_W[NMAXPIONS];
double         pips_PCAL_V[NMAXPIONS];
double         pips_PCAL_x[NMAXPIONS];
double         pips_PCAL_y[NMAXPIONS];
double         pips_PCAL_z[NMAXPIONS];
double    pips_PCAL_sector[NMAXPIONS];
double      pips_DC_sector[NMAXPIONS];
double          pips_Chi2N[NMAXPIONS];
double        pips_DC_x[NMAXPIONS][3];
double        pips_DC_y[NMAXPIONS][3];
double        pips_DC_z[NMAXPIONS][3];
double         pips_E_PCAL[NMAXPIONS];
double         pips_E_ECIN[NMAXPIONS];
double        pips_E_ECOUT[NMAXPIONS];
double               Zpips[NMAXPIONS]; // hadron rest-frame energy
double             ZpipsLC[NMAXPIONS]; // hadron z on the light-cone
double           piplus_Px[NMAXPIONS];
double           piplus_Py[NMAXPIONS];
double           piplus_Pz[NMAXPIONS];
double            piplus_E[NMAXPIONS];
double           Vpiplus_X[NMAXPIONS];
double           Vpiplus_Y[NMAXPIONS];
double           Vpiplus_Z[NMAXPIONS];
double       piplus_qFrame_pT[NMAXPIONS]; // transverse momentum relative to q
double       piplus_qFrame_pL[NMAXPIONS]; // longitudinal momentum relative to q
double    piplus_qFrame_Theta[NMAXPIONS];
double      piplus_qFrame_Phi[NMAXPIONS];


// negative pions
bool     pimsPastSelectionCuts[NMAXPIONS];
bool eepimsPastKinematicalCuts[NMAXPIONS];
int            pims_region[NMAXPIONS];
double        pims_chi2PID[NMAXPIONS];
double         pims_PCAL_W[NMAXPIONS];
double         pims_PCAL_V[NMAXPIONS];
double         pims_PCAL_x[NMAXPIONS];
double         pims_PCAL_y[NMAXPIONS];
double         pims_PCAL_z[NMAXPIONS];
double    pims_PCAL_sector[NMAXPIONS];
double      pims_DC_sector[NMAXPIONS];
double          pims_Chi2N[NMAXPIONS];
double        pims_DC_x[NMAXPIONS][3];
double        pims_DC_y[NMAXPIONS][3];
double        pims_DC_z[NMAXPIONS][3];
double         pims_E_PCAL[NMAXPIONS];
double         pims_E_ECIN[NMAXPIONS];
double        pims_E_ECOUT[NMAXPIONS];
double               Zpims[NMAXPIONS]; // hadron rest-frame energy
double             ZpimsLC[NMAXPIONS]; // hadron z on the light-cone
double              piminus_Px[NMAXPIONS];
double              piminus_Py[NMAXPIONS];
double              piminus_Pz[NMAXPIONS];
double               piminus_E[NMAXPIONS];
double              Vpiminus_X[NMAXPIONS];
double              Vpiminus_Y[NMAXPIONS];
double              Vpiminus_Z[NMAXPIONS];
double       piminus_qFrame_pT[NMAXPIONS]; // transverse momentum relative to q
double       piminus_qFrame_pL[NMAXPIONS]; // longitudinal momentum relative to q
double    piminus_qFrame_Theta[NMAXPIONS];
double      piminus_qFrame_Phi[NMAXPIONS];

// K+
bool KpsPastSelectionCuts[NMAXKAONS], eeKpsPastKinematicalCuts[NMAXKAONS];
int            Kps_region[NMAXKAONS];
double        Kps_chi2PID[NMAXKAONS];
double         Kps_PCAL_W[NMAXKAONS];
double         Kps_PCAL_V[NMAXKAONS];
double         Kps_PCAL_x[NMAXKAONS];
double         Kps_PCAL_y[NMAXKAONS];
double         Kps_PCAL_z[NMAXKAONS];
double    Kps_PCAL_sector[NMAXKAONS];
double      Kps_DC_sector[NMAXKAONS];
double          Kps_Chi2N[NMAXKAONS];
double        Kps_DC_x[NMAXKAONS][3];
double        Kps_DC_y[NMAXKAONS][3];
double        Kps_DC_z[NMAXKAONS][3];
double         Kps_E_PCAL[NMAXKAONS];
double         Kps_E_ECIN[NMAXKAONS];
double        Kps_E_ECOUT[NMAXKAONS];
double               ZKps[NMAXKAONS]; // hadron rest-frame energy
double             ZKpsLC[NMAXKAONS]; // hadron z on the light-cone
double           Kplus_Px[NMAXKAONS];
double           Kplus_Py[NMAXKAONS];
double           Kplus_Pz[NMAXKAONS];
double            Kplus_E[NMAXKAONS];
double           VKplus_X[NMAXKAONS];
double           VKplus_Y[NMAXKAONS];
double           VKplus_Z[NMAXKAONS];
double       Kplus_qFrame_pT[NMAXKAONS]; // transverse momentum relative to q
double       Kplus_qFrame_pL[NMAXKAONS]; // longitudinal momentum relative to q
double    Kplus_qFrame_Theta[NMAXKAONS];
double      Kplus_qFrame_Phi[NMAXKAONS];


// K-
bool     KmsPastSelectionCuts[NMAXKAONS];
bool eeKmsPastKinematicalCuts[NMAXKAONS];
int                Kms_region[NMAXKAONS];
double            Kms_chi2PID[NMAXKAONS];
double             Kms_PCAL_W[NMAXKAONS];
double             Kms_PCAL_V[NMAXKAONS];
double             Kms_PCAL_x[NMAXKAONS];
double         Kms_PCAL_y[NMAXKAONS];
double         Kms_PCAL_z[NMAXKAONS];
double    Kms_PCAL_sector[NMAXKAONS];
double      Kms_DC_sector[NMAXKAONS];
double          Kms_Chi2N[NMAXKAONS];
double        Kms_DC_x[NMAXKAONS][3];
double        Kms_DC_y[NMAXKAONS][3];
double        Kms_DC_z[NMAXKAONS][3];
double         Kms_E_PCAL[NMAXKAONS];
double         Kms_E_ECIN[NMAXKAONS];
double        Kms_E_ECOUT[NMAXKAONS];
double               ZKms[NMAXKAONS]; // hadron rest-frame energy
double             ZKmsLC[NMAXKAONS]; // hadron z on the light-cone
double              Kminus_Px[NMAXKAONS];
double              Kminus_Py[NMAXKAONS];
double              Kminus_Pz[NMAXKAONS];
double               Kminus_E[NMAXKAONS];
double              VKminus_X[NMAXKAONS];
double              VKminus_Y[NMAXKAONS];
double              VKminus_Z[NMAXKAONS];
double       Kminus_qFrame_pT[NMAXKAONS]; // transverse momentum relative to q
double       Kminus_qFrame_pL[NMAXKAONS]; // longitudinal momentum relative to q
double    Kminus_qFrame_Theta[NMAXKAONS];
double      Kminus_qFrame_Phi[NMAXKAONS];



// Output root file and tree
TFile * outFile_e_piplus, * outFile_e_piminus;
TTree * outTree_e_piplus, * outTree_e_piminus;
TFile * outFile_e_Kplus,  * outFile_e_Kminus;
TTree * outTree_e_Kplus,  * outTree_e_Kminus;
TFile * outFile_e_piplus_no_cuts, * outFile_e_piminus_no_cuts;
TTree * outTree_e_piplus_no_cuts, * outTree_e_piminus_no_cuts;

// Output CSV file
std::ofstream   CSVfile_e_piplus,  SelectedEventsCSVfile_e_piplus,  SelectedEventsCSVfile_e_piplus_kinematics;
std::ofstream   CSVfile_e_piminus, SelectedEventsCSVfile_e_piminus, SelectedEventsCSVfile_e_piminus_kinematics;
std::ofstream   CSVfile_e_Kplus,   SelectedEventsCSVfile_e_Kplus,   SelectedEventsCSVfile_e_Kplus_kinematics;
std::ofstream   CSVfile_e_Kminus,  SelectedEventsCSVfile_e_Kminus,  SelectedEventsCSVfile_e_Kminus_kinematics;
// vectors in lab-frame
TLorentzVector   Beam, target, e, q, pi, K;
TLorentzVector              d_rest, p_rest;
std::vector<TLorentzVector>         piplus; // positive pions
std::vector<TLorentzVector>        piminus; // negative pions
std::vector<TLorentzVector>          Kplus; // positive Kaons
std::vector<TLorentzVector>         Kminus; // negative Kaons
// reconstructed vertex position
TVector3                                Ve;
std::vector<TVector3>              Vpiplus;
std::vector<TVector3>             Vpiminus;
std::vector<TVector3>               VKplus;
std::vector<TVector3>              VKminus;

// kinematics
Double_t Ebeam, omega, y, xB, Q2, xF, eta_pi;
Double_t W, W2, W_d, W2_d, M_x, M_x_d, qStar;


// MC information,
// generated quantities
TLorentzVector  P_mc_particle;
TLorentzVector  e_g,    pi_g, q_g;
TVector3 V_mc_particle, Ve_g, Vpi_g;
Double_t e_P_g, e_Theta_g, e_Phi_g, e_Vz_g;
Double_t pi_P_g, pi_Theta_g, pi_Phi_g, pi_Vz_g, Q2_g, xB_g, omega_g, y_g;

// vectors in q-frame
TLorentzVector               e_qFrame, q_qFrame;
std::vector<TLorentzVector>       piplus_qFrame;
std::vector<TLorentzVector>      piminus_qFrame;
std::vector<TLorentzVector>       Kplus_qFrame;
std::vector<TLorentzVector>      Kminus_qFrame;

// auxiliary
DCfid_SIDIS dcfid;
std::vector<region_part_ptr>  electrons, neutrons, protons, gammas, deuterons;
std::vector<region_part_ptr>  pipluses, piminuses;
std::vector<region_part_ptr>  Kpluses,  Kminuses;

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void SetVerbosity( int _fdebug_ ){
    fdebug = _fdebug_;
    aux.SetVerbosity( fdebug );
}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void SetIsMC( bool fIsMC ){
    IsMC = fIsMC;
    if (Skimming=="p_uniform_distribution"){
        IsMC = true;
    }
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void SetSimPi (TString fSimPi) {
    SimPi   = fSimPi; // Simulated pion species
    SimK    = "";
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void SetSimK (TString fSimK) {
    SimPi  = "";
    SimK   = fSimK; // Simulated Kaon species
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void SetDataPath (TString fDataPath, Double_t fEbeam) {
    prefix   = "sidisdvcs_"; // default
    
    if (fDataPath=="" || fDataPath=="sidisdvcs" || fDataPath=="sidis dvcs"){
        // sidis-dvcs train files, used since July 2022
        // (the 'usual' train files)
        if (fEbeam==10.2){
            DataPath = "/cache/clas12/rg-b/production/recon/spring2019/torus-1/pass1/v0/dst/train/sidisdvcs/";
        } else if (fEbeam==10.4){
            DataPath = "/cache/clas12/rg-b/production/recon/spring2020/torus-1/pass1/v1/dst/train/sidisdvcs/";
        } else if (fEbeam==10.6){
            DataPath = "/cache/clas12/rg-b/production/recon/spring2019/torus-1/pass1/v0/dst/train/sidisdvcs/";
        }
        prefix   = "sidisdvcs_";
    }
    else if (fDataPath=="inclusive" || fDataPath=="inc"){
        // inclusive train files, used until July 2022
        // (inclusive train files were only generated in the beginning of RGB without any backup)
        DataPath = "/volatile/clas12/rg-b/production/recon/spring2019/torus-1/pass1/v0/dst/train_20200610/inc/";
        prefix   = "inc_";
    }
    else if (fDataPath=="nSidis" || fDataPath=="nsidis"){
        // free-p data from RGA data
        // For RGA we use nSidis, they key difference is sidisdvcs has e_p > 1 GeV and nSidis has e_p > 2 GeV.
        DataPath = "/cache/clas12/rg-a/production/recon/spring2019/torus-1/pass1/v1/dst/train/nSidis/";
        prefix   = "nSidis_";
    }
    else if (fDataPath=="AcceptanceCorrection"){
        // GEMC simulations of "white" spectra
        // i.e. (e,e'π) events with no physics generator
        DataPath = "/volatile/clas12/users/ecohen/GEMC/hipo/10.2/AcceptanceCorrection/";
        prefix = "p_uniform_distribution";
    }
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void SetSkimming (TString fSkimming) {
    if (fSkimming=="" || fSkimming=="SIDIS_skimming" ){
        // d(e,e'π) files from RGB data
        Skimming = "SIDIS_skimming";
    }
    else if (fSkimming=="RGA_Free_proton"){
        // p(e,e'π) files from RGA data
        Skimming = "RGA_Free_proton";
    }
    else if (fSkimming=="p_uniform_distribution"){
        // (e,e'π) generated uniformly from "white" spectra with no physics
        Skimming = "p_uniform_distribution";
    }
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void SetEbeam (double fEbeam) {
    
    // [GeV]
    // RGA the enrgy was 10.6
    // RGB Spring-2019 the enrgy was 10.2
    // RGB Fall-2019 the enrgy was 10.4096
    Ebeam = fEbeam;
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
bool CheckIfElectronPassedSelectionCuts(Double_t e_PCAL_x, Double_t e_PCAL_y,
                                        Double_t e_PCAL_W, Double_t e_PCAL_V,
                                        Double_t e_E_PCAL,
                                        Double_t e_E_ECIN, Double_t e_E_ECOUT,
                                        TLorentzVector e,
                                        TVector3 Ve,
                                        Double_t e_DC_sector,
                                        Double_t e_DC_x[3],
                                        Double_t e_DC_y[3],
                                        Double_t e_DC_z[3],
                                        int torusBending){
    
    // decide if electron in event passes event selection cuts
    
    // DC - fiducial cuts on DC
    // from bandsoft_tools/skimmers/electrons.cpp,
    // where eHit.getDC_x1() - x position in first region of the drift chamber
    // same for y1,x2,y2,...
    // eHit.getDC_sector() - sector
    // checking DC Fiducials
    // torusBending         torus magnet bending:   ( 1 = inbeding, -1 = outbending    )
    
    // sometimes the readout-sector is 0. This is funny
    // Justin B. Estee (June-21): I also had this issue. I am throwing away sector 0. The way you check is plot the (x,y) coordinates of the sector and you will not see any thing. Double check me but I think it is 0.
    if (e_DC_sector == 0) return false;
    
    for (int regionIdx=0; regionIdx<3; regionIdx++) {
        // DC_e_fid:
        // sector:  1-6
        // layer:   1-3
        // bending: 0(out)/1(in)
        // std::cout << "e_DC_sector: " << e_DC_sector << ", regionIdx: " << regionIdx << std::endl;
        int bending  = 1 ? (torusBending==-1) : 0;
        bool DC_fid  = dcfid.DC_fid_xy_sidis(11,                 // particle PID,
                                             e_DC_x[regionIdx],  // x
                                             e_DC_y[regionIdx],  // y
                                             e_DC_sector,        // sector
                                             regionIdx+1,        // layer
                                             bending);           // torus bending
        if (DC_fid == false) {
            return false;
        }
    }
    
    
    if(!(true
         // fiducial cuts on PCAL
         //fabs(e_PCAL_x)>0
         //&&  fabs(e_PCAL_y)>0
         &&  e_PCAL_W > aux.cutValue_e_PCAL_W
         &&  e_PCAL_V > aux.cutValue_e_PCAL_V
         
         // Electron Identification Refinement  - PCAL Minimum Energy Deposition Cut
         &&  e_E_PCAL > aux.cutValue_e_E_PCAL
         
         // Sampling fraction cut
         && ((e_E_PCAL + e_E_ECIN + e_E_ECOUT)/e.P()) > aux.cutValue_SamplingFraction_min
         && (e_E_ECIN/e.P() > aux.cutValue_PCAL_ECIN_SF_min - e_E_PCAL/e.P()) // RGA AN puts "<" here mistakenly
         
         // Cut on z-vertex position: in-bending torus field -13.0 cm < Vz < +12.0 cm
         // Spring 19 and Spring 2020 in-bending.
         // Fall 2019 (without low-energy-run) was out-bending.
         &&  ((aux.cutValue_Vz_min < Ve.Z()) && (Ve.Z() < aux.cutValue_Vz_max))
         )) return false;
    
    return true;
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
bool CheckIfPionPassedSelectionCuts(TString pionCharge, // "pi+" or "pi-"
                                    Double_t DC_sector,
                                    Double_t DC_x[3], Double_t DC_y[3], Double_t DC_z[3],
                                    Double_t chi2PID, Double_t p,
                                    TVector3 Ve,
                                    TVector3 Vpi,
                                    int fdebug){
    
    // decide if pion (pi+ or pi-) passed event selection cuts
    //
    // input:
    // --------
    // DC_x, DC_y   pi drift-chamber coordinates
    // chi2PID      pi chi2PID     (pips_chi2PID)
    // p            pi momentum    (pi.P())
    //
    // comments
    // ---------------
    // DC - fiducial cuts on DC
    if (fdebug>3) {
        std::cout << "CheckIfPionPassedSelectionCuts()" << std::endl;
    }
    if (DC_sector == 0) { if (fdebug>2){std::cout << "DC_sector=0 (funny...)" << std::endl;} return false;}
    
    int PDGcode;
    double    C;
    if (pionCharge=="pi+"){
        PDGcode = 211;
        C       = 0.88;
    } else if (pionCharge=="pi-") {
        PDGcode = -211;
        C       = 0.93;
    } else {
        std::cout << "π charge ill-defined, returning false" << std::endl;
        return false;
    }
    
    for (int regionIdx=0; regionIdx<3; regionIdx++) {
        // DC_e_fid:
        // sector:  1-6
        // layer:   1-3
        // bending: 0(out)/1(in)
        int bending  = 1 ? (torusBending==-1) : 0;
        // new version Aug-11,2021
        if (fdebug>3) {
            std::cout << "dcfid.DC_fid_th_ph_sidis(): "
            << DC_x[regionIdx] <<     ","
            << DC_y[regionIdx] <<     ","
            << DC_z[regionIdx] <<     ","
            << DC_sector       <<     ","
            << regionIdx+1     <<     ","
            << bending         <<     ","
            << std::endl;
            
        }
        bool DC_fid  = dcfid.DC_fid_th_ph_sidis(PDGcode,            // particle PID
                                                DC_x[regionIdx],    // x
                                                DC_y[regionIdx],    // y
                                                DC_z[regionIdx],    // z
                                                DC_sector,          // sector
                                                regionIdx+1,        // layer
                                                bending);           // torus bending
        if (DC_fid == false) {
            return false;
        }
    }
    
    if (fdebug>3) {
        std::cout << "in CheckIfPionPassedSelectionCuts()"<< std::endl
        << "pion charge: "          << pionCharge               << ","
        << "DC_x[0]: "              << DC_x[0]                  << ","
        << "chi2PID:"               << chi2PID                  << ","
        << "Chi2PID_pion_lowerBound( p="<<p<<", C="<<C<<" ): "
        << aux.Chi2PID_pion_lowerBound( p, C ) << ","
        << "Chi2PID_pion_upperBound( p="<<p<<", C="<<C<<" ): "
        << aux.Chi2PID_pion_upperBound( p, C ) << ","
        << "fabs((Ve-Vpi).Z()): "   << fabs((Ve-Vpi).Z())       << ","
        << std::endl;
    }
    if(!
       // pi+ Identification Refinement - chi2PID vs. momentum
       (( aux.Chi2PID_pion_lowerBound( p, C ) < chi2PID
         &&
         chi2PID < aux.Chi2PID_pion_upperBound( p , C ) )
        
        // Cut on the z-Vertex Difference Between Electrons and Hadrons.
        &&  ( fabs((Ve-Vpi).Z()) < aux.cutValue_Ve_Vpi_dz_max )
        )) {
        return false;
    }
    if (fdebug>3) { std::cout << "succesfully passed CheckIfPionPassedSelectionCuts(), return true" << std::endl; }
    return true;
}



// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
TVector3 GetParticleVertex(clas12::region_part_ptr rp){
    TVector3 V(rp->par()->getVx(),
               rp->par()->getVy(),
               rp->par()->getVz());
    return V;
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void SetLorentzVector (TLorentzVector &p4,clas12::region_part_ptr rp){
    p4.SetXYZM(rp->par()->getPx(),
               rp->par()->getPy(),
               rp->par()->getPz(),
               p4.M());
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void OpenOutputFiles (TString outfilename){
    
    
    //    // Create output tree
    //    outFile_e_piplus_no_cuts  = new TFile( outfilename + "_e_piplus_no_cuts.root"  ,"RECREATE");
    //    outTree_e_piplus_no_cuts  = new TTree( "tree" , "(e,e'pi+) event information - no cuts");
    //    outFile_e_piminus_no_cuts = new TFile( outfilename + "_e_piminus_no_cuts.root"  ,"RECREATE");
    //    outTree_e_piminus_no_cuts = new TTree( "tree" , "(e,e'pi+) event information - no cuts");
    TString header = csvheader;
    
    if (!IsMC){
        outFile_e_piplus  = new TFile( outfilename + "_e_piplus.root"  ,"RECREATE");
        outTree_e_piplus  = new TTree( "tree" , "(e,e'pi+) event information");
        outFile_e_piminus = new TFile( outfilename + "_e_piminus.root" ,"RECREATE");
        outTree_e_piminus = new TTree( "tree" , "(e,e'pi-) event  information");
        
        SelectedEventsCSVfile_e_piplus_kinematics.open( outfilename + "_e_piplus_selected_eepi_kinematics.csv" );
        SelectedEventsCSVfile_e_piplus_kinematics << header << std::endl;
        SelectedEventsCSVfile_e_piminus_kinematics.open( outfilename + "_e_piminus_selected_eepi_kinematics.csv" );
        SelectedEventsCSVfile_e_piminus_kinematics << header << std::endl;
    } else if (IsMC) {
        // GEMC simulation
        // Create output csv files
        header += csvheader_GEMCaddition;
        
        if (SimPi=="piplus"){
            outFile_e_piplus  = new TFile( outfilename + "_e_piplus.root"  ,"RECREATE");
            outTree_e_piplus  = new TTree( "tree" , "(e,e'pi+) event information");
            SelectedEventsCSVfile_e_piplus_kinematics.open( outfilename + "_e_piplus_selected_eepi_kinematics.csv" );
            SelectedEventsCSVfile_e_piplus_kinematics << header << std::endl;
        } else if (SimPi=="piminus"){
            outFile_e_piminus = new TFile( outfilename + "_e_piminus.root" ,"RECREATE");
            outTree_e_piminus = new TTree( "tree" , "(e,e'pi-) event  information");
            SelectedEventsCSVfile_e_piminus_kinematics.open( outfilename + "_e_piminus_selected_eepi_kinematics.csv" );
            SelectedEventsCSVfile_e_piminus_kinematics << header << std::endl;
        }
    }
    if (fdebug>1) std::cout << "Done OpenOutputFiles( " << outfilename << ")" << std::endl;
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void CloseOutputFiles (TString OutDataPath, TString outfilename){
    // close output ROOT and CSV files
    std::cout
    << "Done processesing "  <<  Nevents_processed          << " events,"
    << std::setprecision(3)
    << (float)Nevents_passed_e_cuts/Nevents_processed       << " events passed e cuts,"
    << "output files ready in root/csv formats in "
    << std::endl;
    
    
    //    CSVfile_e_piplus                .close();
    //    SelectedEventsCSVfile_e_piplus  .close();
    //    CSVfile_e_piminus               .close();
    //    SelectedEventsCSVfile_e_piminus .close();
    
    
    // close output ROOT
    //    outFile_e_piplus_no_cuts->cd();
    //    outTree_e_piplus_no_cuts->Write();
    //    outFile_e_piplus_no_cuts->Close();
    //
    //    outFile_e_piminus_no_cuts->cd();
    //    outTree_e_piminus_no_cuts->Write();
    //    outFile_e_piminus_no_cuts->Close();
    
    if ((!IsMC)
        ||
        ((IsMC) && (SimPi=="piplus"))
        ){
        SelectedEventsCSVfile_e_piplus_kinematics  .close();
        int Nentires_e_piplus  = outTree_e_piplus  -> GetEntries();
        
        outFile_e_piplus->cd();
        outTree_e_piplus->Write();
        outFile_e_piplus->Close();
        
        
        std::cout
        << (float)Nevents_passed_pips_cuts/Nevents_processed    << " events passed π+ cuts,"
        << std::endl
        << "\t" << (float)Nevents_passed_e_pips_cuts/Nevents_processed  << " passed (e,e'π+) cuts,"
        << std::endl
        << "\t\t" << (float)Nevents_passed_e_pips_kinematics_cuts/Nevents_processed  << " also passed kinematical cuts,"
        << std::endl;
        
        std::cout << "wrote "  << Nentires_e_piplus  << " to (e,e'π+) root file, "
        << std::endl << outFile_e_piplus -> GetName()
        << "and "  << Nevents_passed_e_pips_kinematics_cuts << " to (e,e'π-) csv file (after kinematical cuts) "
        << std::endl << OutDataPath + outfilename + "_e_piplus_selected_eepi_kinematics.csv"
        << std::endl;
    }
    
    if ((!IsMC)
        ||
        ((IsMC) && (SimPi=="piminus"))
        ){
        SelectedEventsCSVfile_e_piminus_kinematics .close();
        int Nentires_e_piminus = outTree_e_piminus -> GetEntries();
        
        outFile_e_piminus->cd();
        outTree_e_piminus->Write();
        outFile_e_piminus->Close();
        
        std::cout
        << (float)Nevents_passed_pims_cuts/Nevents_processed    << " events passed π- cuts,"
        << std::endl
        << "\t" << (float)Nevents_passed_e_pims_cuts/Nevents_processed  << " passed (e,e'π-) cuts,"
        << std::endl
        <<  "\t\t" << (float)Nevents_passed_e_pims_kinematics_cuts/Nevents_processed  << " also passed kinematical cuts,"
        << std::endl;
        
        std::cout << "output files ready in root/csv formats in " << std::endl
        << std::endl
        << "wrote "  << Nentires_e_piminus << " to (e,e'π-) root file. "
        << std::endl << outFile_e_piminus -> GetName()
        << "and "  << Nevents_passed_e_pims_kinematics_cuts << " to (e,e'π-) csv file (after kinematical cuts) "
        << std::endl << OutDataPath + outfilename + "_e_piminus_selected_eepi_kinematics.csv"
        << std::endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
int GetBeamHelicity( event_ptr p_event, int runnum, int fdebug ){
    // deprecated as of Aug-11, 2021,
    //    since
    // however we keep it here for more data
    // get beam helicity (+1 along the beam and -1 opposite to it)
    // [Christopher Dilks <dilks@jlab.org>, email from Aug-5, 2021]
    // for more items
    // [https://github.com/JeffersonLab/clas12root/blob/master/AccesssingBankDataInCpp.txt]
    
    //// helFlip: if true, REC::Event.helicity has opposite sign from reality
    //def helFlip
    //if(RG=="RGA") helFlip = true
    //else if(RG=="RGB") {
    //  helFlip = true
    //  if(runnum>=11093 && runnum<=11283) helFlip = false // fall, 10.4 GeV period only
    //  else if(runnum>=11323 && runnum<=11571) helFlip = false // winter
    //};
    //else if(RG=="RGK") helFlip = false
    //else if(RG=="RGF") helFlip = true
    if (fdebug>5) std::cout << "beam_helicity = c12.event()->getHelicity()" << std::endl;
    
    beam_helicity = p_event->getHelicity();
    
    if (fdebug>5) std::cout << "check spin flip" << std::endl;
    // we are working here on RGB data
    bool helFlip = true;
    if      (runnum>=11093 && runnum<=11283)    helFlip = false; // falls, 10.4 GeV period only
    else if (runnum>=11323 && runnum<=11571)    helFlip = false; // winter
    
    if (helFlip) {
        beam_helicity = -1 * beam_helicity;
    }
    if (fdebug>5) std::cout << "done GetBeamHelicity() " << std::endl;
    return beam_helicity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SetOutputTTrees(){
    // pi+
    
    if (fdebug>1) std::cout << "SetOutputTTrees()" << std::endl;
    
    
    // pi+
    if ((!IsMC)
        ||
        (IsMC && SimPi=="piplus")){
        
        outTree_e_piplus->Branch("eventnumber"          ,&evnum                 );
        outTree_e_piplus->Branch("runnum"               ,&runnum                );
        outTree_e_piplus->Branch("inclusive"            ,&inclusive             );
        outTree_e_piplus->Branch("e_E_PCAL"             ,&e_E_PCAL              );
        outTree_e_piplus->Branch("e_E_ECIN"             ,&e_E_ECIN              );
        outTree_e_piplus->Branch("e_E_ECOUT"            ,&e_E_ECOUT             );
        outTree_e_piplus->Branch("e_PCAL_W"             ,&e_PCAL_W              );
        outTree_e_piplus->Branch("e_PCAL_V"             ,&e_PCAL_V              );
        outTree_e_piplus->Branch("e_PCAL_x"             ,&e_PCAL_x              );
        outTree_e_piplus->Branch("e_PCAL_y"             ,&e_PCAL_y              );
        outTree_e_piplus->Branch("e_PCAL_z"             ,&e_PCAL_z              );
        outTree_e_piplus->Branch("e_PCAL_sector"        ,&e_PCAL_sector         );
        outTree_e_piplus->Branch("e_DC_sector"          ,&e_DC_sector           );
        outTree_e_piplus->Branch("e_DC_Chi2N"           ,&e_DC_Chi2N            );
        outTree_e_piplus->Branch("e_DC_x"               ,&e_DC_x                , "e_DC_x[3]/D"         );
        outTree_e_piplus->Branch("e_DC_y"               ,&e_DC_y                , "e_DC_y[3]/D"         );
        outTree_e_piplus->Branch("e_DC_z"               ,&e_DC_z                , "e_DC_z[3]/D"         );
        outTree_e_piplus->Branch("Npi"                  ,&Npips                 );
        
        outTree_e_piplus->Branch("pi_region"            ,&pips_region           , "pi_region[Npi]/I"   );
        outTree_e_piplus->Branch("pi_chi2PID"           ,&pips_chi2PID          , "pi_chi2PID[Npi]/D"    );
        outTree_e_piplus->Branch("pi_PCAL_x"            ,&pips_PCAL_x           , "pi_PCAL_x[Npi]/D"     );
        outTree_e_piplus->Branch("pi_PCAL_y"            ,&pips_PCAL_y           , "pi_PCAL_y[Npi]/D"     );
        outTree_e_piplus->Branch("pi_PCAL_z"            ,&pips_PCAL_z           , "pi_PCAL_z[Npi]/D"     );
        outTree_e_piplus->Branch("pi_PCAL_sector"       ,&pips_PCAL_sector      , "pi_PCAL_sector[Npi]/D");
        outTree_e_piplus->Branch("pi_DC_sector"         ,&pips_DC_sector        , "pi_DC_sector[Npi]/D"  );
        outTree_e_piplus->Branch("pi_Chi2N"             ,&pips_Chi2N            , "pi_Chi2N[Npi]/D"      );
        outTree_e_piplus->Branch("pi_DC_x"              ,&pips_DC_x             , "pi_DC_x[Npi][3]/D"    );
        outTree_e_piplus->Branch("pi_DC_y"              ,&pips_DC_y             , "pi_DC_y[Npi][3]/D"    );
        outTree_e_piplus->Branch("pi_DC_z"              ,&pips_DC_z             , "pi_DC_z[Npi][3]/D"    );
        outTree_e_piplus->Branch("pi_E_PCAL"            ,&pips_E_PCAL           , "pi_E_PCAL[Npi]/D"     );
        outTree_e_piplus->Branch("pi_E_ECIN"            ,&pips_E_ECIN           , "pi_E_ECIN[Npi]/D"     );
        outTree_e_piplus->Branch("pi_E_ECIN"            ,&pips_E_ECIN           , "pi_E_ECIN[Npi]/D"     );
        outTree_e_piplus->Branch("pi_E_ECOUT"           ,&pips_E_ECOUT          , "pi_E_ECOUT[Npi]/D"    );
        outTree_e_piplus->Branch("DC_layers"            ,&DC_layers             , "DC_layers[3]/I"      );
        outTree_e_piplus->Branch("e"                    ,&e                     );
        outTree_e_piplus->Branch("Ve"                   ,&Ve                    );
        outTree_e_piplus->Branch("pi"                   ,&piplus                );
        outTree_e_piplus->Branch("Vpi"                  ,&Vpiplus               );
        outTree_e_piplus->Branch("Beam"                 ,&Beam                  );
        outTree_e_piplus->Branch("beam_helicity"        ,&beam_helicity         );
        outTree_e_piplus->Branch("q"                    ,&q                     );
        outTree_e_piplus->Branch("Ebeam"                ,&Ebeam                 );
        outTree_e_piplus->Branch("xB"                   ,&xB                    );
        outTree_e_piplus->Branch("Q2"                   ,&Q2                    );
        outTree_e_piplus->Branch("omega"                ,&omega                 );
        outTree_e_piplus->Branch("W_d"                  ,&W_d                   );
        outTree_e_piplus->Branch("W"                    ,&W                     );
        outTree_e_piplus->Branch("Z"                    ,Zpips                  );
        outTree_e_piplus->Branch("Z_LC"                 ,ZpipsLC                );
        outTree_e_piplus->Branch("y"                    ,&y                     );
        outTree_e_piplus->Branch("EventPassedCuts"      ,&EventPassedCuts       );
        outTree_e_piplus->Branch("ePastCutsInEvent"     ,&ePastCutsInEvent      );
        outTree_e_piplus->Branch("eepipsPastKinematicalCuts",&eepipsPastKinematicalCuts ,"eepipsPastKinematicalCuts[Npi]/O"  );
        outTree_e_piplus->Branch("piPastCutsInEvent"    ,&pipsPastCutsInEvent   ,"piPastCutsInEvent/O"  );
        outTree_e_piplus->Branch("eepipsPastCutsInEvent",&eepipsPastCutsInEvent ,"eepipsPastCutsInEvent/O"  );
        outTree_e_piplus->Branch("Npips"                ,&Npips                 );
        outTree_e_piplus->Branch("Npims"                ,&Npims                 );
        outTree_e_piplus->Branch("Nelectrons"           ,&Ne                    );
        outTree_e_piplus->Branch("Ngammas"              ,&Ngammas               );
        outTree_e_piplus->Branch("Nprotons"             ,&Np                    );
        outTree_e_piplus->Branch("Nneutrons"            ,&Nn                    );
        outTree_e_piplus->Branch("piplus_Px"                ,&piplus_Px              , "piplus_Px[Npi]/D"    );
        outTree_e_piplus->Branch("piplus_Py"                ,&piplus_Py              , "piplus_Py[Npi]/D"    );
        outTree_e_piplus->Branch("piplus_Pz"                ,&piplus_Pz              , "piplus_Pz[Npi]/D"    );
        outTree_e_piplus->Branch("piplus_E"                 ,&piplus_E               , "piplus_E[Npi]/D"    );
        outTree_e_piplus->Branch("Vpiplus_X"                ,&Vpiplus_X              , "Vpiplus_X[Npi]/D"    );
        outTree_e_piplus->Branch("Vpiplus_Y"                ,&Vpiplus_Y              , "Vpiplus_Y[Npi]/D"    );
        outTree_e_piplus->Branch("Vpiplus_Z"                ,&Vpiplus_Z              , "Vpiplus_Z[Npi]/D"    );
        outTree_e_piplus->Branch("piplus_qFrame_pT"         ,&piplus_qFrame_pT       , "piplus_qFrame_pT[Npi]/D");
        outTree_e_piplus->Branch("piplus_qFrame_pL"         ,&piplus_qFrame_pL       , "piplus_qFrame_pL[Npi]/D");
        outTree_e_piplus->Branch("piplus_qFrame_Theta"      ,&piplus_qFrame_Theta    , "piplus_qFrame_Theta[Npi]/D");
        outTree_e_piplus->Branch("piplus_qFrame_Phi"        ,&piplus_qFrame_Phi      , "piplus_qFrame_Phi[Npi]/D");
    }
    
    if ((!IsMC)
        ||
        (IsMC && SimPi=="piminus")){
        
        outTree_e_piminus->Branch("eventnumber"          ,&evnum                 );
        outTree_e_piminus->Branch("runnum"               ,&runnum                );
        outTree_e_piminus->Branch("inclusive"            ,&inclusive             );
        outTree_e_piminus->Branch("e_E_PCAL"             ,&e_E_PCAL              );
        outTree_e_piminus->Branch("e_E_ECIN"             ,&e_E_ECIN              );
        outTree_e_piminus->Branch("e_E_ECOUT"            ,&e_E_ECOUT             );
        outTree_e_piminus->Branch("e_PCAL_W"             ,&e_PCAL_W              );
        outTree_e_piminus->Branch("e_PCAL_V"             ,&e_PCAL_V              );
        outTree_e_piminus->Branch("e_PCAL_x"             ,&e_PCAL_x              );
        outTree_e_piminus->Branch("e_PCAL_y"             ,&e_PCAL_y              );
        outTree_e_piminus->Branch("e_PCAL_z"             ,&e_PCAL_z              );
        outTree_e_piminus->Branch("e_PCAL_sector"        ,&e_PCAL_sector         );
        outTree_e_piminus->Branch("e_DC_sector"          ,&e_DC_sector           );
        outTree_e_piminus->Branch("e_DC_Chi2N"           ,&e_DC_Chi2N            );
        outTree_e_piminus->Branch("e_DC_x"               ,&e_DC_x                , "e_DC_x[3]/D"            );
        outTree_e_piminus->Branch("e_DC_y"               ,&e_DC_y                , "e_DC_y[3]/D"            );
        outTree_e_piminus->Branch("e_DC_z"               ,&e_DC_z                , "e_DC_z[3]/D"            );
        
        outTree_e_piminus->Branch("Npi"                  ,&Npims                 );
        outTree_e_piminus->Branch("pi_region"            ,&pims_region           , "pi_region[Npi]/I"   );
        outTree_e_piminus->Branch("pi_chi2PID"           ,&pims_chi2PID          , "pi_chi2PID[Npi]/D"       );
        outTree_e_piminus->Branch("pi_PCAL_x"            ,&pims_PCAL_x           , "pi_PCAL_x[Npi]/D"        );
        outTree_e_piminus->Branch("pi_PCAL_y"            ,&pims_PCAL_y           , "pi_PCAL_y[Npi]/D"        );
        outTree_e_piminus->Branch("pi_PCAL_z"            ,&pims_PCAL_z           , "pi_PCAL_z[Npi]/D"        );
        outTree_e_piminus->Branch("pi_PCAL_sector"       ,&pims_PCAL_sector      , "pi_PCAL_sector[Npi]/D"   );
        outTree_e_piminus->Branch("pi_DC_sector"         ,&pims_DC_sector        , "pi_DC_sector[Npi]/D"     );
        outTree_e_piminus->Branch("pi_Chi2N"             ,&pims_Chi2N            , "pi_Chi2N[Npi]/D"         );
        outTree_e_piminus->Branch("pi_DC_x"              ,&pims_DC_x             , "pi_DC_x[Npi][3]/D"       );
        outTree_e_piminus->Branch("pi_DC_y"              ,&pims_DC_y             , "pi_DC_y[Npi][3]/D"       );
        outTree_e_piminus->Branch("pi_DC_z"              ,&pims_DC_z             , "pi_DC_z[Npi][3]/D"       );
        outTree_e_piminus->Branch("pi_E_PCAL"            ,&pims_E_PCAL           , "pi_E_PCAL[Npi]/D"        );
        outTree_e_piminus->Branch("pi_E_ECIN"            ,&pims_E_ECIN           , "pi_E_ECIN[Npi]/D"        );
        outTree_e_piminus->Branch("pi_E_ECIN"            ,&pims_E_ECIN           , "pi_E_ECIN[Npi]/D"        );
        outTree_e_piminus->Branch("pi_E_ECOUT"           ,&pims_E_ECOUT          , "pi_E_ECOUT[Npi]/D"       );
        outTree_e_piminus->Branch("DC_layers"            ,&DC_layers             , "DC_layers[3]"           );
        outTree_e_piminus->Branch("e"                   ,&e                     );
        outTree_e_piminus->Branch("pi"                  ,&piminus               );
        outTree_e_piminus->Branch("Ve"                  ,&Ve                    );
        outTree_e_piminus->Branch("Vpi"                 ,&Vpiminus              );
        outTree_e_piminus->Branch("Beam"                ,&Beam                  );
        outTree_e_piminus->Branch("beam_helicity"       ,&beam_helicity         );
        outTree_e_piminus->Branch("q"                   ,&q                     );
        outTree_e_piminus->Branch("Ebeam"               ,&Ebeam                 );
        outTree_e_piminus->Branch("xB"                  ,&xB                    );
        outTree_e_piminus->Branch("Q2"                  ,&Q2                    );
        outTree_e_piminus->Branch("omega"               ,&omega                 );
        outTree_e_piminus->Branch("W_d"                 ,&W_d                   );
        outTree_e_piminus->Branch("W"                   ,&W                     );
        outTree_e_piminus->Branch("Z"                   ,Zpims                  );
        outTree_e_piminus->Branch("Z_LC"                ,ZpimsLC                );
        outTree_e_piminus->Branch("y"                   ,&y                     );
        outTree_e_piminus->Branch("EventPassedCuts"      ,&EventPassedCuts       );
        outTree_e_piminus->Branch("ePastCutsInEvent"     ,&ePastCutsInEvent      );
        outTree_e_piminus->Branch("eepimsPastKinematicalCuts",&eepimsPastKinematicalCuts ,"eepimsPastKinematicalCuts[Npi]/O"  );
        outTree_e_piminus->Branch("piPastCutsInEvent"    ,&pimsPastCutsInEvent   ,"piPastCutsInEvent/O" );
        outTree_e_piminus->Branch("eepimsPastCutsInEvent",&eepimsPastCutsInEvent ,"eepimsPastCutsInEvent/O"  );
        outTree_e_piminus->Branch("Npips"                ,&Npips                 );
        outTree_e_piminus->Branch("Npims"                ,&Npims                 );
        outTree_e_piminus->Branch("Nelectrons"           ,&Ne                    );
        outTree_e_piminus->Branch("Ngammas"              ,&Ngammas               );
        outTree_e_piminus->Branch("Nprotons"             ,&Np                    );
        outTree_e_piminus->Branch("Nneutrons"            ,&Nn                    );
        
        outTree_e_piminus->Branch("piminus_Px"           ,&piminus_Px              , "piminus_Px[Npi]/D"    );
        outTree_e_piminus->Branch("piminus_Py"                ,&piminus_Py              , "piminus_Py[Npi]/D"    );
        outTree_e_piminus->Branch("piminus_Pz"                ,&piminus_Pz              , "piminus_Pz[Npi]/D"    );
        outTree_e_piminus->Branch("piminus_E"                 ,&piminus_E               , "piminus_E[Npi]/D"    );
        outTree_e_piminus->Branch("Vpiminus_X"                ,&Vpiminus_X              , "Vpiminus_X[Npi]/D"    );
        outTree_e_piminus->Branch("Vpiminus_Y"                ,&Vpiminus_Y              , "Vpiminus_Y[Npi]/D"    );
        outTree_e_piminus->Branch("Vpiminus_Z"                ,&Vpiminus_Z              , "Vpiminus_Z[Npi]/D"    );
        outTree_e_piminus->Branch("piminus_qFrame_pT"         ,&piminus_qFrame_pT       , "piminus_qFrame_pT[Npi]/D");
        outTree_e_piminus->Branch("piminus_qFrame_pL"       ,&piminus_qFrame_pL       , "piminus_qFrame_pL[Npi]/D");
        outTree_e_piminus->Branch("piminus_qFrame_Theta"    ,&piminus_qFrame_Theta    , "piminus_qFrame_Theta[Npi]/D");
        outTree_e_piminus->Branch("piminus_qFrame_Phi"      ,&piminus_qFrame_Phi      , "piminus_qFrame_Phi[Npi]/D");
    }
    if (fdebug>1) std::cout << "Done SetOutputTTrees()" << std::endl;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void InitializeFileReading(int NeventsMax, int c12Nentries, int fdebug){
    if (fdebug>1) {
        std::cout << "InitializeFileReading( " << NeventsMax << " , " << c12Nentries << " , " << fdebug << ")" << std::endl;
    }
    // Ebeam   = GetBeamEnergy ( fdebug );
    Beam    .SetPxPyPzE (0, 0, Ebeam, Ebeam );
    target  .SetXYZM    (0, 0, 0,     aux.Md    );
    d_rest  .SetXYZM    (0, 0, 0,     aux.Md    );
    p_rest  .SetXYZM    (0, 0, 0,     aux.Mp    );
    
    NeventsMaxToProcess = NeventsMax;
    if (NeventsMax<0) NeventsMaxToProcess = c12Nentries;
    Nevents_processed           = 0;
    Nevents_passed_e_cuts       = 0;
    Nevents_passed_pips_cuts    = 0;
    Nevents_passed_pims_cuts    = 0;
    Nevents_passed_e_pips_cuts  = 0;
    Nevents_passed_e_pims_cuts  = 0;
    Nevents_passed_e_pips_kinematics_cuts = 0;
    Nevents_passed_e_pims_kinematics_cuts = 0;
    if (fdebug>1) {
        std::cout << "NeventsMaxToProcess =  " << NeventsMaxToProcess << "" << std::endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void InitializeVariables(){
    e = TLorentzVector(0,0,0,db->GetParticle( 11   )->Mass());
    qStar                               = -9999;
    xB          = Q2        = omega     = -9999;
    xF          = y                     = -9999;
    xB_g        = Q2_g      = omega_g   = -9999;
    y_g                                 = -9999;
    W_d         = W                     = -9999;
    M_x         = M_x_d                 = -9999;
    e_E_ECIN    = e_E_ECOUT = e_E_PCAL  = -9999;
    e_PCAL_W    = e_PCAL_V              = -9999;
    e_PCAL_x    = e_PCAL_y  = e_PCAL_z  = -9999;
    e_PCAL_sector                       = -9999;
    e_DC_sector = e_DC_Chi2N            = -9999;
    for (int regionIdx=0; regionIdx<3; regionIdx++) {
        e_DC_x[regionIdx]               = -9999;
        e_DC_y[regionIdx]               = -9999;
        e_DC_z[regionIdx]               = -9999;
    }
    Pe_phi = q_phi = q_theta            = 0;
    
    Ve                                  = TVector3();
    ePastCutsInEvent                    = false;
    
    electrons   .clear();
    neutrons    .clear();
    protons     .clear();
    gammas      .clear();
    
    piplus          .clear();
    piminus         .clear();
    piplus_qFrame   .clear();
    piminus_qFrame  .clear();
    Vpiplus     .clear();
    Vpiminus    .clear();
    pipluses    .clear();
    piminuses   .clear();

    Kplus          .clear();
    Kminus         .clear();
    Kplus_qFrame   .clear();
    Kminus_qFrame  .clear();
    VKplus     .clear();
    VKminus    .clear();
    Kpluses    .clear();
    Kminuses   .clear();
    
    for (int piIdx=0; piIdx<NMAXPIONS; piIdx++) {
        pips_region[piIdx]                          = -9999;
        pips_chi2PID[piIdx]                         = -9999;
        pips_DC_sector[piIdx]                       = -9999;
        pips_PCAL_sector[piIdx]                     = -9999;
        pips_PCAL_W[piIdx] = pips_PCAL_V[piIdx]     = -9999;
        pips_PCAL_x[piIdx] = pips_PCAL_y[piIdx]     = -9999;
        pips_PCAL_z[piIdx]                          = -9999;
        pips_E_PCAL[piIdx]                          = -9999;
        pips_E_ECIN[piIdx] = pips_E_ECOUT[piIdx]    = -9999;
        
        pims_region[piIdx]                          = -9999;
        pims_chi2PID[piIdx]                         = -9999;
        pims_DC_sector[piIdx]                       = -9999;
        pims_PCAL_sector[piIdx]                     = -9999;
        pims_PCAL_W[piIdx] = pims_PCAL_V[piIdx]     = -9999;
        pims_PCAL_x[piIdx] = pims_PCAL_y[piIdx]     = -9999;
        pims_PCAL_z[piIdx]                          = -9999;
        pims_E_PCAL[piIdx]                          = -9999;
        pims_E_ECIN[piIdx] = pims_E_ECOUT[piIdx]    = -9999;
        for (int regionIdx=0; regionIdx<3; regionIdx++) {
            pips_DC_x[piIdx][regionIdx]= pips_DC_y[piIdx][regionIdx]    = -9999;
            pips_DC_z[piIdx][regionIdx]                                 = -9999;
            pims_DC_x[piIdx][regionIdx]= pims_DC_y[piIdx][regionIdx]    = -9999;
            pims_DC_z[piIdx][regionIdx]                                 = -9999;
        }
        piplus       .push_back( TLorentzVector(0,0,0,db->GetParticle( 211 )->Mass()) );
        piplus_qFrame.push_back( TLorentzVector(0,0,0,db->GetParticle( 211 )->Mass()) );
        Vpiplus .push_back( TVector3() );
        pipsPastSelectionCuts[piIdx]                = false;
        eepipsPastKinematicalCuts[piIdx]            = false;
        
        piminus       .push_back( TLorentzVector(0,0,0,db->GetParticle( -211 )->Mass()) );
        piminus_qFrame.push_back( TLorentzVector(0,0,0,db->GetParticle( -211 )->Mass()) );
        Vpiminus.push_back( TVector3() );
        pimsPastSelectionCuts[piIdx]                = false;
        eepimsPastKinematicalCuts[piIdx]            = false;
        
        piplus_Px[piIdx]    = piplus_Py[piIdx]  = piplus_Pz[piIdx]  = piplus_E[piIdx]   = -9999;
        piminus_Px[piIdx]   = piminus_Py[piIdx] = piminus_Pz[piIdx] = piminus_E[piIdx]  = -9999;
        Vpiplus_X[piIdx]    = Vpiplus_Y[piIdx]  = Vpiplus_Z[piIdx]                      = -9999;
        Vpiminus_X[piIdx]   = Vpiminus_Y[piIdx] = Vpiminus_Z[piIdx]                     = -9999;
        piplus_qFrame_pT[piIdx]    = piplus_qFrame_pL[piIdx]                            = -9999;
        piminus_qFrame_pT[piIdx]   = piminus_qFrame_pL[piIdx]                           = -9999;
        piplus_qFrame_Theta[piIdx] = piplus_qFrame_Phi[piIdx]                           = -9999;
        piminus_qFrame_Theta[piIdx]= piminus_qFrame_Phi[piIdx]                          = -9999;
        
    }
    
    for (int KIdx=0; KIdx<NMAXKAONS; KIdx++) {
        Kps_region[KIdx]                          = -9999;
        Kps_chi2PID[KIdx]                         = -9999;
        Kps_DC_sector[KIdx]                       = -9999;
        Kps_PCAL_sector[KIdx]                     = -9999;
        Kps_PCAL_W[KIdx] = Kps_PCAL_V[KIdx]     = -9999;
        Kps_PCAL_x[KIdx] = Kps_PCAL_y[KIdx]     = -9999;
        Kps_PCAL_z[KIdx]                          = -9999;
        Kps_E_PCAL[KIdx]                          = -9999;
        Kps_E_ECIN[KIdx] = Kps_E_ECOUT[KIdx]    = -9999;
        
        Kms_region[KIdx]                          = -9999;
        Kms_chi2PID[KIdx]                         = -9999;
        Kms_DC_sector[KIdx]                       = -9999;
        Kms_PCAL_sector[KIdx]                     = -9999;
        Kms_PCAL_W[KIdx] = Kms_PCAL_V[KIdx]     = -9999;
        Kms_PCAL_x[KIdx] = Kms_PCAL_y[KIdx]     = -9999;
        Kms_PCAL_z[KIdx]                          = -9999;
        Kms_E_PCAL[KIdx]                          = -9999;
        Kms_E_ECIN[KIdx] = Kms_E_ECOUT[KIdx]    = -9999;
        for (int regionIdx=0; regionIdx<3; regionIdx++) {
            Kps_DC_x[KIdx][regionIdx]= Kps_DC_y[KIdx][regionIdx]    = -9999;
            Kps_DC_z[KIdx][regionIdx]                                 = -9999;
            Kms_DC_x[KIdx][regionIdx]= Kms_DC_y[KIdx][regionIdx]    = -9999;
            Kms_DC_z[KIdx][regionIdx]                                 = -9999;
        }
        Kplus       .push_back( TLorentzVector(0,0,0,db->GetParticle( 211 )->Mass()) );
        Kplus_qFrame.push_back( TLorentzVector(0,0,0,db->GetParticle( 211 )->Mass()) );
        VKplus .push_back( TVector3() );
        KpsPastSelectionCuts[KIdx]                = false;
        eeKpsPastKinematicalCuts[KIdx]            = false;
        
        Kminus       .push_back( TLorentzVector(0,0,0,db->GetParticle( -211 )->Mass()) );
        Kminus_qFrame.push_back( TLorentzVector(0,0,0,db->GetParticle( -211 )->Mass()) );
        VKminus.push_back( TVector3() );
        KmsPastSelectionCuts[KIdx]                = false;
        eeKmsPastKinematicalCuts[KIdx]            = false;
        
        Kplus_Px[KIdx]    = Kplus_Py[KIdx]  = Kplus_Pz[KIdx]  = Kplus_E[KIdx]   = -9999;
        Kminus_Px[KIdx]   = Kminus_Py[KIdx] = Kminus_Pz[KIdx] = Kminus_E[KIdx]  = -9999;
        VKplus_X[KIdx]    = VKplus_Y[KIdx]  = VKplus_Z[KIdx]                      = -9999;
        VKminus_X[KIdx]   = VKminus_Y[KIdx] = VKminus_Z[KIdx]                     = -9999;
        Kplus_qFrame_pT[KIdx]    = Kplus_qFrame_pL[KIdx]                            = -9999;
        Kminus_qFrame_pT[KIdx]   = Kminus_qFrame_pL[KIdx]                           = -9999;
        Kplus_qFrame_Theta[KIdx] = Kplus_qFrame_Phi[KIdx]                           = -9999;
        Kminus_qFrame_Theta[KIdx]= Kminus_qFrame_Phi[KIdx]                          = -9999;
        
    }
    DC_layer                                        = -9999;
    status                                          = 1; // 0 is good...
    
    pipsPastCutsInEvent                             = false;
    eepipsPastCutsInEvent                           = false;
    pimsPastCutsInEvent                             = false;
    eepimsPastCutsInEvent                           = false;
    
    KpsPastCutsInEvent                             = false;
    eeKpsPastCutsInEvent                           = false;
    KmsPastCutsInEvent                             = false;
    eeKmsPastCutsInEvent                           = false;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void OpenResultFiles( TString outfilepath, TString outfilename ){
    OpenOutputFiles( outfilepath + outfilename );
    // output tree branches
    SetOutputTTrees();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExtractElectronInformation(int fdebug){
    // ------------------------------------------------------------------------------------------------
    // extract electron information
    // ------------------------------------------------------------------------------------------------
    // find leading electron as the one with highest energy
    double  leading_e_E;
    int     leading_e_index = 0;
    SetLorentzVector(e,electrons[0]);
    TLorentzVector e_tmp(0,0,0,db->GetParticle(11)->Mass());
    for (int eIdx=0; eIdx < Ne; eIdx++) {
        SetLorentzVector(e_tmp  ,electrons[eIdx]);
        double Ee = e_tmp.E();
        if (Ee > leading_e_E) {
            leading_e_index = eIdx;
            leading_e_E     = Ee;
        }
    }
    // set leading electron 4-momentum
    SetLorentzVector(e , electrons[leading_e_index]);
    // set leading electron vertex
    Ve              = GetParticleVertex( electrons[leading_e_index] );
    
    // detector information on electron
    auto e_PCAL_info= electrons[leading_e_index]->cal(PCAL);
    e_E_PCAL        = e_PCAL_info->getEnergy();
    e_PCAL_sector   = e_PCAL_info->getSector();
    e_PCAL_V        = e_PCAL_info->getLv();
    e_PCAL_W        = e_PCAL_info->getLw();
    e_E_ECIN        = electrons[leading_e_index]->cal(ECIN)->getEnergy();
    e_E_ECOUT       = electrons[leading_e_index]->cal(ECOUT)->getEnergy();
    
    // hit position in PCAL
    e_PCAL_x        = e_PCAL_info->getX();
    e_PCAL_y        = e_PCAL_info->getY();
    e_PCAL_z        = e_PCAL_info->getZ();
    
    // Drift Chamber tracking system
    auto e_DC_info  = electrons[leading_e_index]->trk(DC);
    e_DC_sector     = e_DC_info->getSector(); // tracking sector
    e_DC_Chi2N      = e_DC_info->getChi2N();  // tracking chi^2/NDF
    
    for (int regionIdx=0; regionIdx<3; regionIdx++) {
        int DC_layer = DC_layers[regionIdx];
        e_DC_x[regionIdx] = electrons[leading_e_index]->traj(DC,DC_layer)->getX();
        e_DC_y[regionIdx] = electrons[leading_e_index]->traj(DC,DC_layer)->getY();
        e_DC_z[regionIdx] = electrons[leading_e_index]->traj(DC,DC_layer)->getZ();
    }
    if (fdebug > 2) std::cout << "extracted electron information and computed kinematics" << std::endl;
    
    // ------------------------------------------------------------------------------------------------
    // now, check if electron passed event selection requirements
    // ------------------------------------------------------------------------------------------------
    ePastCutsInEvent = CheckIfElectronPassedSelectionCuts(e_PCAL_x, e_PCAL_y,
                                                          e_PCAL_W, e_PCAL_V,
                                                          e_E_PCAL, e_E_ECIN,
                                                          e_E_ECOUT,
                                                          e, Ve,
                                                          e_PCAL_sector, // e_PCAL_sector should be consistent with e_DC_sector
                                                          e_DC_x, e_DC_y, e_DC_z,
                                                          torusBending );
    if (inclusive)         ePastCutsInEvent = true;
    if (ePastCutsInEvent)  Nevents_passed_e_cuts++ ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExtractPionsInformation(int fdebug){
    
    // positive pions)
    for (int pipsIdx=0; pipsIdx < Npips; pipsIdx++) {
        ExtractPipsInformation( pipsIdx, fdebug );
    }
    // negative pions
    for (int pimsIdx=0; pimsIdx < Npims; pimsIdx++) {
        ExtractPimsInformation( pimsIdx, fdebug );
    }
    
    // move to q-frame and define the pion momentum with respect to q
    MoveTo_qFrame( fdebug );
    
    // done
    if (fdebug > 2) std::cout << "done extracting pion information" << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExtractKaonsInformation(int fdebug){
    
    // positive pions)
    for (int Kidx=0; Kidx < NKps; Kidx++) {
        ExtractKpsInformation( Kidx, fdebug );
    }
    // negative pions
    for (int Kidx=0; Kidx < NKms; Kidx++) {
        ExtractKmsInformation( Kidx, fdebug );
    }
    
    // move to q-frame and define the pion momentum with respect to q
    MoveTo_qFrame( fdebug );
    
    // done
    if (fdebug > 2) std::cout << "done extracting Kaon information" << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void WriteEventToOutput(int fdebug){
    // (Maybe) write this event to "selected events csv-file"
    bool            IsSelected_eepi = false;
    bool             IsSelected_eeK = false;
    
    //ePastCutsInEvent = true;
    if (inclusive == 1) {
        pipsPastCutsInEvent     = true;
        pimsPastCutsInEvent     = true;
        eepipsPastCutsInEvent   = true;
        eepimsPastCutsInEvent   = true;
        KpsPastCutsInEvent      = true;
        KmsPastCutsInEvent      = true;
        eeKpsPastCutsInEvent    = true;
        eeKmsPastCutsInEvent    = true;
    }
    
    // Fill "no-cuts" TTrees even if nothing passed criteria
    // outTree_e_piplus_no_cuts -> Fill();
    // outTree_e_piminus_no_cuts -> Fill();
    
    if (ePastCutsInEvent && pipsPastCutsInEvent) {
        if ((!IsMC) || (IsMC && SimPi=="piplus")){
            IsSelected_eepi = true;
            outTree_e_piplus -> Fill();
            if (fdebug>3) std::cout << "Filling (e,e'pi+) TTree with this event!" << std::endl;
            
            Nevents_passed_e_pips_cuts ++ ;
            if (eepipsPastCutsInEvent) Nevents_passed_e_pips_kinematics_cuts ++;
            
            for (int pipsIdx=0; pipsIdx<Npips; pipsIdx++) {
                Stream_e_pi_line_to_CSV( "pi+", pipsIdx,
                                        pipsPastSelectionCuts[pipsIdx], eepipsPastKinematicalCuts[pipsIdx],
                                        fdebug );
            }
        }
    }
    
    if (ePastCutsInEvent && pimsPastCutsInEvent) {
        if ((!IsMC) || (IsMC && SimPi=="piminus")){
            IsSelected_eepi = true;
            outTree_e_piminus -> Fill();
            if (fdebug>3) std::cout << "Filling (e,e'pi-) TTree with this event!" << std::endl;
            Nevents_passed_e_pims_cuts ++ ;
            if (eepimsPastCutsInEvent) Nevents_passed_e_pims_kinematics_cuts ++;
            
            for (int pimsIdx=0; pimsIdx<Npims; pimsIdx++) {
                Stream_e_pi_line_to_CSV( "pi-", pimsIdx,
                                        pimsPastSelectionCuts[pimsIdx], eepimsPastKinematicalCuts[pimsIdx],
                                        fdebug );
            }
        }
    }
    
    if (ePastCutsInEvent && KmsPastCutsInEvent) {
        if ((!IsMC) || (IsMC && SimPi=="Kminus")){
            IsSelected_eeK = true;
            outTree_e_Kminus -> Fill();
            if (fdebug>3) std::cout << "Filling (e,e'K-) TTree with this event!" << std::endl;
            Nevents_passed_e_Kms_cuts ++ ;
            if (eeKmsPastCutsInEvent) Nevents_passed_e_Kms_kinematics_cuts ++;
            
            for (int Kidx=0; Kidx<NKms; Kidx++) {
                Stream_e_K_line_to_CSV( "K-", Kidx,
                                       KmsPastSelectionCuts[Kidx], eeKmsPastKinematicalCuts[Kidx],
                                       fdebug );
            }
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void FinishProgram(TString outfilepath, TString outfilename){
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    
    std::cout << "Done. Elapsed time: " << elapsed.count() << std::endl;
    CloseOutputFiles( outfilepath, outfilename );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ComputeElectronKinematics(){
    // compute event kinematics (from e-only information)
    q       = Beam - e;
    Q2      = -q.Mag2();
    omega   = q.E();
    xB      = Q2/(2. * aux.Mp * q.E());
    y       = omega / Ebeam;
    W       = sqrt((p_rest + q).Mag2());
    W_d     = sqrt((d_rest + q).Mag2());
    
    if (IsMC){
        q_g       = Beam - e_g;
        Q2_g      = -q_g.Mag2();
        omega_g   = q_g.E();
        xB_g      = Q2_g/(2. * aux.Mp * q_g.E());
        y_g       = omega_g / Ebeam;
        if (fdebug>1)
            std::cout
            << "ComputeElectronKinematics( for MC )"
            << std::endl
            << "Ebeam: "   << Ebeam     << ","
            << "omega_g: " << omega_g   << ","
            << std::endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ComputePionKinematics(TLorentzVector pi, TLorentzVector pi_qFrame){
    
    // pion energy fraction
    Zpi          = pi.E()/omega;
    Zpi_LC       = (pi_qFrame.E() + pi_qFrame.Pz()) / (q.E() + q.P());
    
    // additional kinematical variables
    // assuming scattering off a proton at rest
    M_x     = ( q + p_rest - pi ).Mag();
    M_x_d   = ( q + d_rest - pi ).Mag();
    eta_pi  = 0.5 * log((pi_qFrame.E()+pi_qFrame.Pz()) /
                        (pi_qFrame.E()-pi_qFrame.Pz()));
    
    qStar   = calcQStar( e_qFrame.Vect(), pi_qFrame.Vect(), Ebeam );

    //    xF      = 2. * (pi.Vect().Dot(q.Vect())) / (q.P() * W);
    
    // Compute xF following Chris Dilks
    // [https://github.com/c-dilks/dispin/blob/a48674568ac4c92af4c10d633397332ba727d20c/src/Dihadron.cxx#L100-L121]
    
    // lab frame 4-vectors
    TLorenzVector  vecElectron = e;
    TLorenzVector      vecBeam = TLorentzVector(0.0,0.0,TMath::Sqrt(TMath::Power(Ebeam,2)-TMath::Power(aux.Me,2)),Ebeam );
    TLorenzVector    vecTarget = TLorentzVector(0.0,0.0,0.0, aux.Mp );
    TLorenzVector         vecQ = vecBeam - vecElectron;
    TLorenzVector         vecW = vecBeam + vecTarget - vecElectron;
    TLorenzVector  boostvecCom = vecQ + vecTarget;
    TVector3          ComBoost = -1 * boostvecCom.BoostVector();

    // -- boost to CoM frame
    TLorenzVector    vecPh_com = pi_qFrame; // P+q COM frame Ph
    TLorenzVector  disVecQ_com = q_qFrame;  // P+q COM frame Q
    vecPh_com                  .Boost( ComBoost );
    disVecQ_com                .Boost( ComBoost );
    TVector3           pPh_com = vecPh_com.Vect();
    TVector3            pQ_com = disVecQ_com.Vect();
    // compute xF
    xF                         = 2 * pPh_com.Dot(pQ_com) / (vecW.M() * pQ_com.Mag());
    
    if (fdebug>3){
        std::cout
        << "ComputePionKinematics()"
        << std::endl
        << "Zpi: "                  << Zpi      << ","
        << "M_x: "                  << M_x      << " GeV/c2,"
        << "W: "                    << W        << " GeV/c2, "
        << std::endl
        << "pi_qFrame.Vect().Mag()" << pi_qFrame.Vect().Mag() << ", "
        << "xF: "                   << xF       << ", "
        << "2. * (pi.Vect().Dot(q.Vect())) / (q.P() * W): "                   << 2. * (pi.Vect().Dot(q.Vect())) / (q.P() * W)       << ", "
        << std::endl
        << "eta_pi: "               << eta_pi   << ","
        << "qStar: "                << qStar    << ","
        << std::endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExtractPipsInformation( int pipsIdx, int fdebug ){
    if (fdebug>2)
        std::cout << "ExtractPipsInformation( pipsIdx=" << pipsIdx << ", fdebug=" << fdebug << " )" << std::endl;
    
    
    // Extract positive pion information
    pips_region[pipsIdx] = pipluses[pipsIdx]->getRegion();
    // First - we restrict ourselves to pions only from the central detector
    if( pipluses[pipsIdx]->getRegion() != FD ){
        if (fdebug>2){
            SetLorentzVector(pi  ,pipluses[pipsIdx]);
            std::cout
            << "piplus [" << pipsIdx << "] not from FD (from " << pipluses[pipsIdx]->getRegion()
            << ") "
            << "p = "       << pi.P() << " GeV/c, "
            << "theta = "   << pi.Theta()*180./3.1415 << "˚"
            << ", not extracting information..."
            << std::endl;
        }
        return;
    }else{
        if (fdebug>2) std::cout << "piplus ["<<pipsIdx<<"] is from FD ("<<pipluses[pipsIdx]->getRegion()<< "), extracting information..." << std::endl;
    }
    
    // Now we extract the information on this pion
    SetLorentzVector(piplus[pipsIdx]  ,pipluses[pipsIdx]);
    Zpips[pipsIdx]              = piplus[pipsIdx].E() / omega;
    Vpiplus[pipsIdx]            = GetParticleVertex( pipluses[pipsIdx] );
    pips_chi2PID[pipsIdx]       = pipluses[pipsIdx]->par()->getChi2Pid();
    
    // EC in and out
    pips_E_ECIN[pipsIdx]        = pipluses[pipsIdx]->cal(ECIN)->getEnergy();
    pips_E_ECOUT[pipsIdx]       = pipluses[pipsIdx]->cal(ECOUT)->getEnergy();
    // PCAL
    auto pips_PCAL_info         = pipluses[pipsIdx]->cal(PCAL);
    pips_E_PCAL[pipsIdx]        = pips_PCAL_info->getEnergy();
    pips_PCAL_sector[pipsIdx]   = pips_PCAL_info->getSector();
    pips_PCAL_V[pipsIdx]        = pips_PCAL_info->getLv();
    pips_PCAL_W[pipsIdx]        = pips_PCAL_info->getLw();
    pips_PCAL_x[pipsIdx]        = pips_PCAL_info->getX();
    pips_PCAL_y[pipsIdx]        = pips_PCAL_info->getY();
    pips_PCAL_z[pipsIdx]        = pips_PCAL_info->getZ();
    // DC
    auto pips_DC_info           = pipluses[pipsIdx]->trk(DC);
    pips_DC_sector[pipsIdx]     = pips_DC_info->getSector(); // tracking sector
    pips_Chi2N[pipsIdx]         = pips_DC_info->getChi2N();  // tracking chi^2/NDF
    for (int regionIdx=0; regionIdx<3; regionIdx++) {
        DC_layer = DC_layers[regionIdx];
        pips_DC_x[pipsIdx][regionIdx] = pipluses[pipsIdx]->traj(DC,DC_layer)->getX();
        pips_DC_y[pipsIdx][regionIdx] = pipluses[pipsIdx]->traj(DC,DC_layer)->getY();
        pips_DC_z[pipsIdx][regionIdx] = pipluses[pipsIdx]->traj(DC,DC_layer)->getZ();
        if (fdebug>3) {
            std::cout
            << "pips_DC_sector[pipsIdx="<<pipsIdx<<"]="
            << pips_DC_sector[pipsIdx]
            << ", DC_layer = " << DC_layer
            << std::endl
            << "pips_DC_x[pipsIdx="<<pipsIdx<<"][regionIdx="<<regionIdx<<"]="
            << pips_DC_x[pipsIdx][regionIdx]
            << std::endl
            << "pips_DC_y[pipsIdx="<<pipsIdx<<"][regionIdx="<<regionIdx<<"]="
            << pips_DC_y[pipsIdx][regionIdx]
            << std::endl
            << "pips_DC_z[pipsIdx="<<pipsIdx<<"][regionIdx="<<regionIdx<<"]="
            << pips_DC_z[pipsIdx][regionIdx]
            << std::endl;
        }
    }
    // ------------------------------------------------------------------------------------------------
    // now, check if pion passed event selection requirements
    // ------------------------------------------------------------------------------------------------
    pipsPastSelectionCuts[pipsIdx] = CheckIfPionPassedSelectionCuts("pi+",
                                                                    pips_DC_sector[pipsIdx],
                                                                    pips_DC_x[pipsIdx],
                                                                    pips_DC_y[pipsIdx],
                                                                    pips_DC_z[pipsIdx],
                                                                    pips_chi2PID[pipsIdx],
                                                                    piplus[pipsIdx].P(),
                                                                    Ve,
                                                                    Vpiplus[pipsIdx],
                                                                    fdebug);
    eepipsPastKinematicalCuts[pipsIdx] = aux.eepiPassedKinematicalCriteria(Ebeam,
                                                                           omega,
                                                                           Q2,
                                                                           y,
                                                                           W,
                                                                           piplus[pipsIdx],
                                                                           e);
    
    if (pipsPastSelectionCuts[pipsIdx]) {
        pipsPastCutsInEvent = true;
        Nevents_passed_pips_cuts ++;
        if (eepipsPastKinematicalCuts[pipsIdx]) {
            eepipsPastCutsInEvent = true;
        }
    }
    
    piplus_Px[pipsIdx] = piplus[pipsIdx].Px();
    piplus_Py[pipsIdx] = piplus[pipsIdx].Py();
    piplus_Pz[pipsIdx] = piplus[pipsIdx].Pz();
    piplus_E[pipsIdx]  = piplus[pipsIdx].E();
    Vpiplus_X[pipsIdx] = Vpiplus[pipsIdx].X();
    Vpiplus_Y[pipsIdx] = Vpiplus[pipsIdx].Y();
    Vpiplus_Z[pipsIdx] = Vpiplus[pipsIdx].Z();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExtractPimsInformation( int pimsIdx, int fdebug ){
    // Extract negative pion information
    pims_region[pimsIdx] = piminuses[pimsIdx]->getRegion();
    // First - we restrict ourselves to pions only from the central detector
    if( piminuses[pimsIdx]->getRegion() != FD ){
        if (fdebug>2){
            SetLorentzVector(pi  ,piminuses[pimsIdx] );
            
            std::cout << "piminus ["<<pimsIdx<<"] not from FD (from "
            << piminuses[pimsIdx]->getRegion()
            << ") "
            << "p = "       << pi.P() << " GeV/c, "
            << "theta = "   << pi.Theta()*180./3.1415 << "˚"
            << ", not extracting information..." << std::endl;
        }
        return;
    }
    else{
        if (fdebug>2) std::cout << "piminus ["<<pimsIdx<<"] is from FD ("<<piminuses[pimsIdx]->getRegion()<< "), extracting information..." << std::endl;
    }
    
    // Now we extract pion information
    SetLorentzVector(piminus[pimsIdx]  ,piminuses[pimsIdx]);
    Zpims[pimsIdx]              = piminus[pimsIdx].E() / omega;
    Vpiminus[pimsIdx]           = GetParticleVertex( piminuses[pimsIdx] );
    pims_chi2PID[pimsIdx]       = piminuses[pimsIdx]->par()->getChi2Pid();
    
    // EC in and out
    pims_E_ECIN[pimsIdx]        = piminuses[pimsIdx]->cal(ECIN)->getEnergy();
    pims_E_ECOUT[pimsIdx]       = piminuses[pimsIdx]->cal(ECOUT)->getEnergy();
    // PCAL
    auto pims_PCAL_info         = piminuses[pimsIdx]->cal(PCAL);
    pims_E_PCAL[pimsIdx]        = pims_PCAL_info->getEnergy();
    pims_PCAL_sector[pimsIdx]   = pims_PCAL_info->getSector();
    pims_PCAL_V[pimsIdx]        = pims_PCAL_info->getLv();
    pims_PCAL_W[pimsIdx]        = pims_PCAL_info->getLw();
    pims_PCAL_x[pimsIdx]        = pims_PCAL_info->getX();
    pims_PCAL_y[pimsIdx]        = pims_PCAL_info->getY();
    pims_PCAL_z[pimsIdx]        = pims_PCAL_info->getZ();
    // DC
    auto pims_DC_info           = piminuses[pimsIdx]->trk(DC);
    pims_DC_sector[pimsIdx]     = pims_DC_info->getSector(); // tracking sector
    pims_Chi2N[pimsIdx]         = pims_DC_info->getChi2N();  // tracking chi^2/NDF
    for (int regionIdx=0; regionIdx<3; regionIdx++) {
        DC_layer = DC_layers[regionIdx];
        pims_DC_x[pimsIdx][regionIdx] = piminuses[pimsIdx]->traj(DC,DC_layer)->getX();
        pims_DC_y[pimsIdx][regionIdx] = piminuses[pimsIdx]->traj(DC,DC_layer)->getY();
        pims_DC_z[pimsIdx][regionIdx] = piminuses[pimsIdx]->traj(DC,DC_layer)->getZ();
    }
    // ------------------------------------------------------------------------------------------------
    // now, check if pion passed event selection requirements
    // ------------------------------------------------------------------------------------------------
    pimsPastSelectionCuts[pimsIdx] = CheckIfPionPassedSelectionCuts("pi-",
                                                                    pims_DC_sector[pimsIdx],
                                                                    pims_DC_x[pimsIdx],
                                                                    pims_DC_y[pimsIdx],
                                                                    pims_DC_z[pimsIdx],
                                                                    pims_chi2PID[pimsIdx],  piminus[pimsIdx].P(),
                                                                    Ve,
                                                                    Vpiminus[pimsIdx],
                                                                    fdebug);
    eepimsPastKinematicalCuts[pimsIdx] = aux.eepiPassedKinematicalCriteria(Ebeam,
                                                                           omega,
                                                                           Q2,
                                                                           y,
                                                                           W,
                                                                           piminus[pimsIdx],
                                                                           e);
    if (pimsPastSelectionCuts[pimsIdx]) {
        pimsPastCutsInEvent = true;
        Nevents_passed_pims_cuts ++;
        if (eepimsPastKinematicalCuts[pimsIdx]) {
            eepimsPastCutsInEvent = true;
        }
    }
    
    piminus_Px[pimsIdx]          = piminus[pimsIdx].Px();
    piminus_Py[pimsIdx]          = piminus[pimsIdx].Py();
    piminus_Pz[pimsIdx]          = piminus[pimsIdx].Pz();
    piminus_E[pimsIdx]           = piminus[pimsIdx].E();
    Vpiminus_X[pimsIdx]          = Vpiminus[pimsIdx].X();
    Vpiminus_Y[pimsIdx]          = Vpiminus[pimsIdx].Y();
    Vpiminus_Z[pimsIdx]          = Vpiminus[pimsIdx].Z();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExtractKpsInformation( int KpsIdx, int fdebug ){
    if (fdebug>2)
        std::cout << "ExtractKpsInformation( KpsIdx=" << KpsIdx << ", fdebug=" << fdebug << " )" << std::endl;
    // CONTINUE HERE! COMPLETE THIS FUNCITON AND THEN ALSO FOR K-....
    
    // Extract positive pion information
    Kps_region[KpsIdx] = Kpluses[KpsIdx]->getRegion();
    // First - we restrict ourselves to pions only from the central detector
    if( pipluses[KpsIdx]->getRegion() != FD ){
        if (fdebug>2){
            SetLorentzVector(K  ,Kpluses[KpsIdx]);
            std::cout
            << "piplus [" << KpsIdx << "] not from FD (from " << Kpluses[KpsIdx]->getRegion()
            << ") "
            << "p = "       << K.P() << " GeV/c, "
            << "theta = "   << K.Theta()*180./3.1415 << "˚"
            << ", not extracting information..."
            << std::endl;
        }
        return;
    }else{
        if (fdebug>2) std::cout << "piplus ["<<KpsIdx<<"] is from FD ("<<Kpluses[KpsIdx]->getRegion()<< "), extracting information..." << std::endl;
    }
    
    // Now we extract the information on this Kaon
    SetLorentzVector(Kplus[KpsIdx]  ,Kpluses[KpsIdx]);
    ZKps[KpsIdx]              = Kplus[KpsIdx].E() / omega;
    VKplus[KpsIdx]            = GetParticleVertex( Kpluses[KpsIdx] );
    Kps_chi2PID[KpsIdx]       = Kpluses[KpsIdx]->par()->getChi2Pid();
    
    // EC in and out
    Kps_E_ECIN[KpsIdx]        = Kpluses[KpsIdx]->cal(ECIN)->getEnergy();
    Kps_E_ECOUT[KpsIdx]       = Kpluses[KpsIdx]->cal(ECOUT)->getEnergy();
    // PCAL
    auto Kps_PCAL_info        = Kpluses[KpsIdx]->cal(PCAL);
    Kps_E_PCAL[KpsIdx]        = Kps_PCAL_info->getEnergy();
    Kps_PCAL_sector[KpsIdx]   = Kps_PCAL_info->getSector();
    Kps_PCAL_V[KpsIdx]        = Kps_PCAL_info->getLv();
    Kps_PCAL_W[KpsIdx]        = Kps_PCAL_info->getLw();
    Kps_PCAL_x[KpsIdx]        = Kps_PCAL_info->getX();
    Kps_PCAL_y[KpsIdx]        = Kps_PCAL_info->getY();
    Kps_PCAL_z[KpsIdx]        = Kps_PCAL_info->getZ();
    // DC
    auto Kps_DC_info          = Kpluses[KpsIdx]->trk(DC);
    Kps_DC_sector[KpsIdx]     = Kps_DC_info->getSector(); // tracking sector
    Kps_Chi2N[KpsIdx]         = Kps_DC_info->getChi2N();  // tracking chi^2/NDF
    for (int regionIdx=0; regionIdx<3; regionIdx++) {
        DC_layer = DC_layers[regionIdx];
        Kps_DC_x[KpsIdx][regionIdx] = Kpluses[KpsIdx]->traj(DC,DC_layer)->getX();
        Kps_DC_y[KpsIdx][regionIdx] = Kpluses[KpsIdx]->traj(DC,DC_layer)->getY();
        Kps_DC_z[KpsIdx][regionIdx] = Kpluses[KpsIdx]->traj(DC,DC_layer)->getZ();
    }
    //    // ------------------------------------------------------------------------------------------------
    //    // now, check if pion passed event selection requirements
    //    // ------------------------------------------------------------------------------------------------
    //    pipsPastSelectionCuts[pipsIdx] = CheckIfPionPassedSelectionCuts("pi+",
    //                                                                    pips_DC_sector[pipsIdx],
    //                                                                    pips_DC_x[pipsIdx],
    //                                                                    pips_DC_y[pipsIdx],
    //                                                                    pips_DC_z[pipsIdx],
    //                                                                    pips_chi2PID[pipsIdx],
    //                                                                    piplus[pipsIdx].P(),
    //                                                                    Ve,
    //                                                                    Vpiplus[pipsIdx],
    //                                                                    fdebug);
    //    eepipsPastKinematicalCuts[pipsIdx] = aux.eepiPassedKinematicalCriteria(Ebeam,
    //                                                                           omega,
    //                                                                           Q2,
    //                                                                           y,
    //                                                                           W,
    //                                                                           piplus[pipsIdx],
    //                                                                           e);
    //
    //    if (pipsPastSelectionCuts[pipsIdx]) {
    //        pipsPastCutsInEvent = true;
    //        Nevents_passed_pips_cuts ++;
    //        if (eepipsPastKinematicalCuts[pipsIdx]) {
    //            eepipsPastCutsInEvent = true;
    //        }
    //    }
    //
    //    piplus_Px[pipsIdx] = piplus[pipsIdx].Px();
    //    piplus_Py[pipsIdx] = piplus[pipsIdx].Py();
    //    piplus_Pz[pipsIdx] = piplus[pipsIdx].Pz();
    //    piplus_E[pipsIdx]  = piplus[pipsIdx].E();
    //    Vpiplus_X[pipsIdx] = Vpiplus[pipsIdx].X();
    //    Vpiplus_Y[pipsIdx] = Vpiplus[pipsIdx].Y();
    //    Vpiplus_Z[pipsIdx] = Vpiplus[pipsIdx].Z();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GetParticlesByType (int evnum, int fdebug){
    // get particles by type
    Ne      = electrons .size();
    Nn      = neutrons  .size();
    Np      = protons   .size();
    Npips   = pipluses  .size();
    Npims   = piminuses .size();
    Ngammas = gammas    .size();
    Nd      = deuterons.size();
    if (fdebug>2){
        std::cout
        << "Particles in event "     << evnum        << ": "
        << "N(electrons): "                 << Ne           <<  ","
        << "N(protons): "                   << Np           <<  ","
        << "N(neutrons): "                  << Nn           <<  ","
        << "N(pi+): "                       << Npips        <<  ","
        << "N(pi-): "                       << Npims        <<  ","
        << "N(gammas): "                    << Ngammas      <<  ","
        << "N(deuterons): "                 << Nd           <<  ","
        << std::endl;
        
        if (fdebug>6){
            std::cout
            << "size(piplus): "         << piplus.size()        << ","
            << "size(piminus): "        << piminus.size()       << ","
            << std::endl;
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Stream_e_pi_line_to_CSV( TString pionCharge, int piIdx,
                             bool passed_cuts_e_pi,
                             bool passed_cuts_e_pi_kinematics,
                             int fdebug ){
    TLorentzVector         pi;
    TLorentzVector  pi_qFrame;
    TVector3              Vpi;
    int          pi_DC_sector;
    
    if (pionCharge=="pi+") {
        pi              = piplus       .at(piIdx);
        pi_qFrame       = piplus_qFrame.at(piIdx);
        Vpi             = Vpiplus         [piIdx];
        pi_DC_sector    = pips_DC_sector  [piIdx];
    }
    else if (pionCharge=="pi-") {
        pi           = piminus       .at(piIdx);
        pi_qFrame    = piminus_qFrame.at(piIdx);
        Vpi          = Vpiminus         [piIdx];
        pi_DC_sector = pims_DC_sector   [piIdx];
    }
    else {
        std::cout << "Stream_e_pi_line_to_CSV(): bad pi charge! " << std::endl;
        return;
    }
    
    ComputePionKinematics( pi, pi_qFrame );
    
    // now stream data to CSV file
    std::vector<double> variables =
    {   (double)status, (double)runnum,     (double)evnum,      (double)beam_helicity,
        e.P(),          e.Theta(),          e.Phi(),            Ve.Z(),
        pi.P(),         pi.Theta(),         pi.Phi(),           Vpi.Z(),
        Q2,             xB,                 omega,              y,
        (double)e_DC_sector,                (double)pi_DC_sector,
        e_qFrame.Theta(),                   e_qFrame.Phi(),
        pi_qFrame.Theta(),                  pi_qFrame.Phi(),
        pi_qFrame.Pt(),                     pi_qFrame.Pz(),
        Zpi,            Zpi_LC,
        W,              M_x,
        xF,             eta_pi,
        W_d,            M_x_d,
        q.P(),          qStar,
    };
    
    if (IsMC){
        // GEMC simulations - add truth information
        
        std::vector<double> MCvariables ={
            e_g.P(),  e_g.Theta(),          e_g.Phi(),            Ve_g.Z(),
            pi_g.P(), pi_g.Theta(),         pi_g.Phi(),           Vpi_g.Z(),
            Q2_g,     xB_g,                 omega_g,              y_g     };
        
        variables.insert(variables.end(), std::begin(MCvariables), std::end(MCvariables));
        
    }
    
    // decide which file to write...
    if (pionCharge=="pi+") {
        
        // do not write empty pi+ lines to (e,e'π+) file if simulation is (e,e'π-)
        if ((!IsMC)
            ||
            ((IsMC) && (SimPi=="piplus"))){
            
            if (passed_cuts_e_pi && passed_cuts_e_pi_kinematics) {
                aux.StreamToCSVfile(SelectedEventsCSVfile_e_piplus_kinematics,
                                    variables,
                                    csvprecisions );
            }
        }
    } else if (pionCharge=="pi-") {
        
        // do not write empty pi- lines to (e,e'π-) file if simulation is (e,e'π+)
        if ((!IsMC)
            ||
            ((IsMC) && (SimPi=="piminus"))){
            
            
            if (passed_cuts_e_pi && passed_cuts_e_pi_kinematics) {
                aux.StreamToCSVfile(SelectedEventsCSVfile_e_piminus_kinematics,
                                    variables,
                                    csvprecisions );
            }
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MoveTo_qFrame(int fdebug){
    if (fdebug>1){
        std::cout << "Moving to q-Frame" <<std::endl;
    }
    //    Move to the "q-frame" and define the pion momentum in this frame
    //    q-frame is defined as follows:
    //    z axis is defined by the q - parallel to q
    //    x axis is defined by the e' - such that p(e') resides in the x-z plane
    
    // (1) define q-angles
    q_phi   = q.Phi();
    q_theta = q.Theta();
    
    // (2) rotate Pe and q according to q angles
    TVector3 Pe = e.Vect();
    TVector3 Pq = q.Vect();
    
    Pe       .RotateZ(-q_phi);
    Pe       .RotateY(-q_theta);
    Pe_phi = Pe.Phi();
    Pe       .RotateZ(-Pe_phi);
    Pq       .RotateZ(-q_phi);
    Pq       .RotateY(-q_theta);
    
    
    // (3) verify on q and Pe that the frame-change is done correctly
    //    RotateVectorTo_qFrame( &Pe );
    e_qFrame.SetVectM( Pe, aux.Me );
    //    RotateVectorTo_qFrame( &Pq );
    q_qFrame.SetVectM( Pq, q.M() );
    
    if (fdebug>2){
        aux.Print4Vector( e, "e" );
        aux.Print4Vector( e_qFrame, "e in q-Frame" );
        aux.Print4Vector( q , "q");
        aux.Print4Vector( q_qFrame, "q in q-Frame" );
    }
    
    
    // (4) rotate pions to this q-frame
    // pi+ int the q-Frame
    for (int piIdx=0; piIdx<Npips; piIdx++) {
        TVector3 Ppiplus = RotateVectorTo_qFrame( piplus.at(piIdx).Vect() );
        piplus_qFrame.at(piIdx).SetVectM( Ppiplus, aux.Mpi  );
        TLorentzVector Ppi_q = piplus_qFrame.at(piIdx);
        // fill variables that later go to TTree
        piplus_qFrame_pT[piIdx]   = Ppi_q.Pt();
        piplus_qFrame_pL[piIdx]   = Ppi_q.Pz();
        piplus_qFrame_Theta[piIdx]= Ppi_q.Theta();
        piplus_qFrame_Phi[piIdx]  = Ppi_q.Phi();
        
        // z on the light-cone
        ZpipsLC[piIdx]            = (Ppi_q.E() - Ppi_q.Pz()) / (q.E() - q.P());
        
        if (fdebug>1) aux.Print4Vector( Ppi_q, "pi+(" + std::to_string(piIdx) + ")" );
        
    }
    // pi- int the q-Frame
    for (int piIdx=0; piIdx<Npims; piIdx++) {
        TVector3 Ppiminus = RotateVectorTo_qFrame( piminus.at(piIdx).Vect() );
        piminus_qFrame.at(piIdx).SetVectM( Ppiminus, aux.Mpi );
        TLorentzVector Ppi_q = piminus_qFrame.at(piIdx);
        // fill variables that later go to TTree
        piminus_qFrame_pT[piIdx]   = Ppi_q.Pt();
        piminus_qFrame_pL[piIdx]   = Ppi_q.Pz();
        piminus_qFrame_Theta[piIdx]= Ppi_q.Theta();
        piminus_qFrame_Phi[piIdx]  = Ppi_q.Phi();
        // z on the light-cone
        ZpimsLC[piIdx]             = (Ppi_q.E() - Ppi_q.Pz()) / (q.E() - q.P());
        
        if (fdebug>1) aux.Print4Vector( Ppi_q, "pi-(" + std::to_string(piIdx) + ")" );
    }
    // (5) rotate Kaons to this q-frame
    // K+ int the q-Frame
    for (int KIdx=0; KIdx<NKps; KIdx++) {
        TVector3 PKplus = RotateVectorTo_qFrame( Kplus.at(KIdx).Vect() );
        Kplus_qFrame.at(KIdx).SetVectM( PKplus, aux.MK  );
        TLorentzVector PK_q = Kplus_qFrame.at(KIdx);
        // fill variables that later go to TTree
        Kplus_qFrame_pT[KIdx]   = PK_q.Pt();
        Kplus_qFrame_pL[KIdx]   = PK_q.Pz();
        Kplus_qFrame_Theta[KIdx]= PK_q.Theta();
        Kplus_qFrame_Phi[KIdx]  = PK_q.Phi();
        // z on the light-cone
        ZKpsLC[KIdx]            = (PK_q.E() - PK_q.Pz()) / (q.E() - q.P());
        if (fdebug>1) aux.Print4Vector( PK_q, "K+(" + std::to_string(KIdx) + ")" );
    }
    // K- int the q-Frame
    for (int KIdx=0; KIdx<NKms; KIdx++) {
        TVector3 PKminus = RotateVectorTo_qFrame( Kminus.at(KIdx).Vect() );
        Kminus_qFrame.at(KIdx).SetVectM( PKminus, aux.MK  );
        TLorentzVector PK_q = Kminus_qFrame.at(KIdx);
        // fill variables that later go to TTree
        Kminus_qFrame_pT[KIdx]   = PK_q.Pt();
        Kminus_qFrame_pL[KIdx]   = PK_q.Pz();
        Kminus_qFrame_Theta[KIdx]= PK_q.Theta();
        Kminus_qFrame_Phi[KIdx]  = PK_q.Phi();
        // z on the light-cone
        ZKmsLC[KIdx]            = (PK_q.E() - PK_q.Pz()) / (q.E() - q.P());
        if (fdebug>1) aux.Print4Vector( PK_q, "K-(" + std::to_string(KIdx) + ")" );
    }
    
    
    if (fdebug>2){
        std::cout
        << "size(piplus_qFrame): "  << piplus_qFrame.size() << ","
        << "size(piminus_qFrame): " << piminus_qFrame.size()<< ","
        << std::endl;
    }
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TVector3 RotateVectorTo_qFrame( TVector3 v ){
    // move to q-Pe system: q is the z axis, Pe is in x-z plane: Pe=(Pe[x],0,Pe[q])
    v.RotateZ( -q_phi  );
    v.RotateY( -q_theta);
    v.RotateZ( -Pe_phi );
    return v;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SetInclusive( int fInclusive ){
    inclusive = fInclusive;
    if (inclusive == 1) std::cout << "Running as inclusive (all data registered event regardless of the application of event selection cuts)" << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SetSpecificFilename (TString fSpecificFilename=""){
    SpecificFilename = fSpecificFilename;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SetSpecificFilePath (TString fSpecificFilePath=""){
    SpecificFilePath = fSpecificFilePath;
}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
// Main functionality
// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void SIDISc12rSkimmer(int RunNumber    = 6420   ,
                      int NeventsMax   = -1     ,
                      int fdebug       = 1      ,
                      int PrintProgress= 50000  ,
                      // first event to analyze
                      int FirstEvent   = 0      ,
                      // skimming:
                      // "SIDIS_skimming"  , "RGA_Free_proton" , "p_uniform_distribution"
                      TString fSkimming= "SIDIS_skimming",
                      // data-path symbolic name
                      // "sidisdvcs", "inc", "nSidis", "AcceptanceCorrection"
                      TString fDataPath= "sidisdvcs",
                      double fEbeam    = 10.2, // [GeV]
                      // simulation - π+ / π- 
                      TString fSimPi   = "piplus", //""/"piplus" / "piminus"
                      // inclusive means do not apply cuts
                      int   fInclusive = 0,
                      //simulation - K+ / K-
                      TString fSimK    = "", // ""/"Kplus" / "Kminus"
                      // specific file name
                      TString fSpecificFilePath = "",
                      TString fSpecificFilename = "",
                      bool       fIsMC = false
                      ){
    
    
    SetVerbosity        ( fdebug     );
    SetDataPath         ( fDataPath, fEbeam );
    SetSkimming         ( fSkimming  );
    SetEbeam            ( fEbeam     );
    SetSimPi            ( fSimPi     );
    SetSimK             ( fSimK      );
    SetInclusive        ( fInclusive );
    SetIsMC             ( fIsMC      );
    SetSpecificFilename ( fSpecificFilename );
    SetSpecificFilePath ( fSpecificFilePath );
    
    TString RunNumberStr = aux.GetRunNumberSTR(RunNumber, Skimming);
    
    // read cut values
    aux.loadCutValues("macros/cuts/BANDcutValues.csv",torusBending);
    
    // define input filename
    TString infilename, outfilepath, outfilename;
    
    if (SpecificFilePath=="" || SpecificFilename=="") {
        infilename  = DataPath + prefix + RunNumberStr + ".hipo";
        outfilepath = "/volatile/clas12/users/ecohen/BAND/" + Skimming + "/";
        outfilename = "skimmed_SIDIS_" + prefix + RunNumberStr;
    }
    else {
        infilename  = SpecificFilePath + "/" + SpecificFilename + ".hipo";
        outfilepath = SpecificFilePath + "/";
        outfilename = "skimmed_SIDIS_" + SpecificFilename;
    }
    
    
    
    // Simulation - define input and output file names for GEMC simulations slightly different
    // if we read off a simulation file - the simulation was done
    // for an (e,e'π+) or an (e,e'π-) data set
    // and this is labeled in the filename
    if (IsMC){
        if (SimPi != ""){
            infilename = (DataPath
                          + "/"   + SimPi
                          + "/ee" + SimPi + "_" + prefix
                          + "_" + RunNumberStr + "_reco.hipo");
            outfilename = "skimmed_SIDIS_" + SimPi + "_" + prefix + "_" + RunNumberStr;
        }
        else if (SimK != ""){
            infilename = (DataPath
                          + "/"   + SimK
                          + "/ee" + SimK + "_" + prefix
                          + "_" + RunNumberStr + "_reco.hipo");
            outfilename = "skimmed_SIDIS_" + SimK + "_" + prefix + "_" + RunNumberStr;
        }
    }
    
    if (fdebug>1){
        std::cout
        << "Input file name: " << std::endl
        << infilename          << std::endl
        << "Output file name: "<< std::endl
        << outfilepath + outfilename
        << std::endl;
    }
    
    // open input file
    TChain fake("hipo");
    fake.Add(infilename.Data());
    //get the hipo data
    auto files = fake.GetListOfFiles();
    
    // open output files
    OpenResultFiles( outfilepath, outfilename );
    
    // start analysis
    // step over events and extract information....
    for(Int_t i=0;i<files->GetEntries();i++){
        
        //create the event reader
        if (fdebug) std::cout << "reading file " << i << std::endl;
        clas12reader c12(files->At(i)->GetTitle(),{0});
        InitializeFileReading( NeventsMax, c12.getReader().getEntries(), fdebug );
        
        int event = 0;
        
        // process the events...
        while((c12.next()==true) && (event < NeventsMaxToProcess)){
            
            event++;
            if (event > FirstEvent) {
                
                runnum = c12.runconfig()->getRun();
                evnum  = c12.runconfig()->getEvent();
                if (fdebug>2) {
                    std::cout
                    << ".............................................."
                    << std::endl
                    << "begin analysis of event " << evnum
                    << " (run " << runnum << ")"
                    << std::endl;
                }
                
                GetBeamHelicity    ( c12.event() , runnum, fdebug );
                InitializeVariables();
                // Get Particles By Type
                electrons   = c12.getByID( 11   );
                neutrons    = c12.getByID( 2112 );
                protons     = c12.getByID( 2212 );
                pipluses    = c12.getByID( 211  );
                piminuses   = c12.getByID(-211  );
                gammas      = c12.getByID( 22   );
                deuterons   = c12.getByID( 1000010020 );
                Kpluses     = c12.getByID( 321  );
                Kminuses    = c12.getByID(-321  );
                GetParticlesByType ( evnum, fdebug );
                
                if (fdebug>8) {
                    std::cout
                    << "Ne: "           << Ne        << ","
                    << "inclusive: "    << inclusive << ","
                    << "Npips: "        << Npips     << ","
                    << "Npims: "        << Npims     << ","
                    << "NKps: "         << NKps      << ","
                    << "NKms: "         << NKms      << ","
                    << std::endl;
                }
                // filter events, extract information, and compute event kinematics:
                // we keep only events of;
                // d(e,e’π+)X, d(e,e’π-)X, d(e,e’K+)X, d(e,e’K-)X
                if(  0 < Ne // after studying some MC and data, we need to kill events with more than 1 electron
                   &&
                   ((inclusive == 1) || (0 < Npips) || (0 < Npims)|| (0 < NKps) || (0 < NKms)) // "inclusive" means (e,e') events
                   &&
                   (Npips < NMAXPIONS) && (Npims < NMAXPIONS) // we don't want to crash the memeory
                   &&
                   (NKps < NMAXKAONS)  && (NKms < NMAXKAONS) // we don't want to crash the memeory
                   ){
                    
                    ExtractElectronInformation  (fdebug);
                    ComputeElectronKinematics   ();
                    ExtractPionsInformation     (fdebug);
                    ExtractKaonsInformation     (fdebug);
                    WriteEventToOutput          (fdebug);
                    
                } else {
                    if (fdebug>1) {
                        std::cout << "Skipped computations in this event as there are not enough particles: "
                        << "Ne = " << Ne
                        << ",N(π+) = " << Npips << ", N(π-) = " << Npims
                        << ",N(K+) = " << NKps  << ", N(K-) = " << NKms
                        << std::endl ;
                    }
                }
                if (fdebug>1) {
                    std::cout << "done processing event " << evnum
                    << " (" << event << "/" << NeventsMaxToProcess<< ") "
                    << std::endl << "------------------------------------------------------------" << std::endl << std::endl ;
                }
                
                if (IsMC){
                    // add truth-information,
                    // i.e. generated electron and generated pion information
                    auto mcpbank = c12.mcparts();
                    const Int_t Ngen=mcpbank->getRows();
                    if (fdebug>1) std::cout << "Grabbing truth-information of " << Ngen << " particles" << std::endl;
                    
                    for( Int_t i_mc =0; i_mc< Ngen ; i_mc++){
                        mcpbank -> setEntry(i_mc);
                        
                        P_mc_particle.SetXYZM( mcpbank->getPx() , mcpbank->getPy() , mcpbank->getPz() , mcpbank->getMass() );
                        V_mc_particle.SetXYZ( mcpbank->getVx() , mcpbank->getVy() , mcpbank->getVz() );
                        auto pid = mcpbank->getPid();
                        
                        if (fdebug>2){
                            std::cout << "MC particle PDG code " << pid
                            << std::setprecision(4)
                            << ", p: "    << P_mc_particle.P()          << " GeV/c, "
                            << ", theta: "<< P_mc_particle.Theta()      << ", " << P_mc_particle.Theta()*r2d  << " deg, "
                            << ", phi: "  << P_mc_particle.Phi()        << ", "<< P_mc_particle.Phi()*r2d    << " deg, "
                            << ", V(z): " << V_mc_particle.Z()          << " cm"
                            << std::endl;
                        }
                        
                        if ( pid==11 ) {
                            e_g = P_mc_particle;
                            Ve_g = V_mc_particle;
                        }
                        else if ( (pid==211) && (SimPi=="piplus") ) {
                            pi_g = P_mc_particle;
                            Vpi_g = V_mc_particle;
                        }
                        else if ( (pid==-211) && (SimPi=="piminus") ) {
                            pi_g = P_mc_particle;
                            Vpi_g = V_mc_particle;
                        }
                        else if ( (pid==321) && (SimPi=="Kplus") ) {
                            pi_g = P_mc_particle;
                            Vpi_g = V_mc_particle;
                        }
                        else if ( (pid==-321) && (SimPi=="Kminus") ) {
                            pi_g = P_mc_particle;
                            Vpi_g = V_mc_particle;
                        }
                        else  {
                            if (fdebug>2){
                                std::cout << "MC particle PDG code " << pid << " do not match generated particles: e (" << 11 << ") + ";
                                if ( SimPi=="piplus")
                                    std::cout << " π+ (" << 211;
                                else if ( SimPi=="piminus")
                                    std::cout << " π- (" << -211;
                                else if ( SimK=="Kplus")
                                    std::cout << " K+ (" << 321;
                                else if ( SimK=="Kminus")
                                    std::cout << " K- (" << -321;
                                std::cout << ")" << std::endl;
                            }
                        }
                    }
                }
                Nevents_processed++;
            }
            if (fdebug && event%PrintProgress==0 && (event > FirstEvent)){
                std::cout
                << std::setprecision(1)
                << event << "/" << NeventsMaxToProcess
                << std::endl;
            }
        } // end event loop
        
    } // end file loop
    
    
    FinishProgram( outfilepath, outfilename);
}
