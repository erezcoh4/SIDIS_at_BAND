#include <vector>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TH1.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include <time.h>

#include "Auxiliary/bank.h"
#include "Auxiliary/BBand.h"
#include "Auxiliary/BEvent.h"
#include "Auxiliary/constants.h"
#include "Auxiliary/bandhit.cpp"
#include "Auxiliary/clashit.cpp"
#include "Auxiliary/genpart.cpp"

#include "Auxiliary/SIDISatBAND_auxiliary.cpp"

#define NMAXEVENTS 5000000
#define NMAXNML 20 // maximal number of XXX

SIDISatBAND_auxiliary aux;

// time
clock_t tStart = clock();

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
                     +(TString)"W_d,");

std::vector<int> csvprecisions = {
    0,0,0,0,
    9,9,9,9,
    9,9,9,9,
    9,9,9,9,
    0,0,
};

TString                   DataPath;
TString                   Filename;
TString       CSVFilename_e_e_pips;
TString       CSVFilename_e_e_pims;


int               Nevents_e_e_pips;
int               Nevents_e_e_pims;
int                         fdebug;
int               NeventsInRGATree;
Long64_t                    runnum;
Long64_t                     evnum;
int                     NMaxEvents;
int                            nml;
//int                            nPi;
int                           npiP;
int                           npiM;
int                      inclusive; // tag to look at inclusive run - all the events with no selection
int       Ne, Nn, Np, Npips, Npims;
int                          e_idx; // leading electron index

TFile                    * RGAFile;
TTree                    * RGATree;

// kinematics
Double_t           Ebeam, omega, y, xB, Q2;
Double_t                                xF;
Double_t                            eta_pi;
Double_t                             W, W2;
Double_t                         W_d, W2_d;
Double_t                               M_x;
Double_t                                w2; // omega^2
double                                 Zpi;
double                              Zpi_LC;


// Input TTree variables
Float_t                Ep[NMAXNML]; // electron momentum
Float_t               Epx[NMAXNML]; // electron momentum
Float_t               Epy[NMAXNML]; // electron momentum
Float_t               Epz[NMAXNML]; // electron momentum
Float_t               Evz[NMAXNML]; // electron vertex in z-direction
Int_t                Epid[NMAXNML]; // electron PID
Int_t             Esector[NMAXNML]; // electron sector
Float_t            Etheta[NMAXNML];
Float_t              Ephi[NMAXNML];
Int_t          Efiducial1[NMAXNML];
Int_t          Efiducial2[NMAXNML];
Int_t          Efiducial3[NMAXNML];
Float_t      ePCAL_energy[NMAXNML];
Float_t           ePCAL_x[NMAXNML];
Float_t           ePCAL_y[NMAXNML];
Float_t           ePCAL_z[NMAXNML];
Float_t          ePCAL_lu[NMAXNML];
Float_t          ePCAL_lv[NMAXNML];
Float_t          ePCAL_lw[NMAXNML];
Float_t      eECIN_energy[NMAXNML];
Float_t     eECOUT_energy[NMAXNML];


Float_t                piPp[NMAXNML]; // electron momentum
Float_t               piPpx[NMAXNML]; // electron momentum
Float_t               piPpy[NMAXNML]; // electron momentum
Float_t               piPpz[NMAXNML]; // electron momentum
Float_t               piPvz[NMAXNML]; // electron vertex in z-direction
Int_t                piPpid[NMAXNML]; // electron PID
Int_t             piPsector[NMAXNML]; // electron sector
Float_t            piPtheta[NMAXNML];
Float_t              piPphi[NMAXNML];
Int_t          piPfiducial1[NMAXNML];
Int_t          piPfiducial2[NMAXNML];
Int_t          piPfiducial3[NMAXNML];
Float_t      piPPCAL_energy[NMAXNML];
Float_t           piPPCAL_x[NMAXNML];
Float_t           piPPCAL_y[NMAXNML];
Float_t           piPPCAL_z[NMAXNML];
Float_t          piPPCAL_lu[NMAXNML];
Float_t          piPPCAL_lv[NMAXNML];
Float_t          piPPCAL_lw[NMAXNML];
Float_t      piPECIN_energy[NMAXNML];
Float_t     piPECOUT_energy[NMAXNML];



Float_t                piMp[NMAXNML]; // electron momentum
Float_t               piMpx[NMAXNML]; // electron momentum
Float_t               piMpy[NMAXNML]; // electron momentum
Float_t               piMpz[NMAXNML]; // electron momentum
Float_t               piMvz[NMAXNML]; // electron vertex in z-direction
Int_t                piMpid[NMAXNML]; // electron PID
Int_t             piMsector[NMAXNML]; // electron sector
Float_t            piMtheta[NMAXNML];
Float_t              piMphi[NMAXNML];
Int_t          piMfiducial1[NMAXNML];
Int_t          piMfiducial2[NMAXNML];
Int_t          piMfiducial3[NMAXNML];
Float_t      piMPCAL_energy[NMAXNML];
Float_t           piMPCAL_x[NMAXNML];
Float_t           piMPCAL_y[NMAXNML];
Float_t           piMPCAL_z[NMAXNML];
Float_t          piMPCAL_lu[NMAXNML];
Float_t          piMPCAL_lv[NMAXNML];
Float_t          piMPCAL_lw[NMAXNML];
Float_t      piMECIN_energy[NMAXNML];
Float_t     piMECOUT_energy[NMAXNML];
Float_t                Nu[NMAXNML];
Float_t                qx[NMAXNML];
Float_t                qy[NMAXNML];
Float_t                qz[NMAXNML];



// Our SIDIS analysis variables
// cut values
std::vector<std::pair<std::string, double>> cutValues;
double               cutValue_Vz_min;
double               cutValue_Vz_max;
double             cutValue_e_PCAL_W;
double             cutValue_e_PCAL_V;
double             cutValue_e_E_PCAL;
double cutValue_SamplingFraction_min;
double     cutValue_PCAL_ECIN_SF_min;
double        cutValue_Ve_Vpi_dz_max;
double               cutValue_Q2_min;
double               cutValue_Q2_max;
double                cutValue_W_min;
double                cutValue_y_max;
double          cutValue_e_theta_min;
double          cutValue_e_theta_max;
double         cutValue_pi_theta_min;
double         cutValue_pi_theta_max;
double              cutValue_Ppi_min;
double              cutValue_Ppi_max;
double               cutValue_Pe_min;
double               cutValue_Pe_max;
double              cutValue_Zpi_min;
double              cutValue_Zpi_max;
double        Pe_phi, q_phi, q_theta; // "q-frame" parameters

bool        ePastCutsInEvent = false;
bool     pipsPastCutsInEvent = false;
bool   eepipsPastCutsInEvent = false;
bool     pimsPastCutsInEvent = false;
bool         EventPassedCuts = false;
bool   eepimsPastCutsInEvent = false;

// meta-data
int           torusBending = -1; // -1 for In-bending, +1 for Out-bending
int    DC_layers[3] = {6,18,36};// Region 1 is denoted at DC detector 6, Region 2 is denoted 18, Region 3 - as 36
int                    DC_layer;
int               beam_helicity; // helicity of the electron +1 along the beam and -1 opposite to it
int                      status;
int         NeventsMaxToProcess;
int           Nevents_processed;
int       Nevents_passed_e_cuts;
int    Nevents_passed_pips_cuts;
int  Nevents_passed_e_pips_cuts;
int Nevents_passed_e_pips_kinematics_cuts;
int    Nevents_passed_pims_cuts;
int  Nevents_passed_e_pims_cuts;
int Nevents_passed_e_pims_kinematics_cuts;
int               PrintProgress;

// leading electron
double            e_E_PCAL; // electron energy deposit in PCAL [GeV]
double            e_E_ECIN; // electron energy deposit in ECAL_in [GeV]
double           e_E_ECOUT; // electron energy deposit in ECAL_out [GeV]
double            e_PCAL_W;
double            e_PCAL_V;
double            e_PCAL_x;
double            e_PCAL_y;
double            e_PCAL_z;
double       e_PCAL_sector;
double         e_DC_sector;
double          e_DC_Chi2N;
double           e_DC_x[3];
double           e_DC_y[3];
double           e_DC_z[3];

// positive pions
bool     pipsPastSelectionCuts[NMAXPIONS];
bool eepipsPastKinematicalCuts[NMAXPIONS];
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
double             ZpipsLC[NMAXPIONS]; // hadron rest-frame energy on the light cone

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
double             ZpimsLC[NMAXPIONS]; // hadron rest-frame energy on the light cone

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


// Output root file and tree
//TFile * outFile_e_piplus, * outFile_e_piminus;
//TTree * outTree_e_piplus, * outTree_e_piminus;
//TFile * outFile_e_piplus_no_cuts, * outFile_e_piminus_no_cuts;
//TTree * outTree_e_piplus_no_cuts, * outTree_e_piminus_no_cuts;

// Output CSV file
std::ofstream   CSVfile_e_e_pips, CSVfile_e_e_pims;


// vectors in lab-frame
TLorentzVector          Beam, target, e, q;
TLorentzVector              d_rest, p_rest;
std::vector<TLorentzVector>         piplus; // positive pions
std::vector<TLorentzVector>        piminus; // negative pions
// reconstructed vertex position
TVector3                                Ve;
std::vector<TVector3>              Vpiplus;
std::vector<TVector3>             Vpiminus;

// vectors in q-frame
TLorentzVector               e_qFrame, q_qFrame;
std::vector<TLorentzVector>       piplus_qFrame;
std::vector<TLorentzVector>      piminus_qFrame;



// methods
void                            PrintDone ();
void                        OpenInputFile ();
void                       CloseInputFile ();
void                      OpenOutputFiles ();
void                     CloseOutputFiles ();
void                InitializeFileReading ();
void                         SetVerbosity (int ffdebug )          {
    fdebug = ffdebug;
    aux.SetVerbosity(fdebug);
} ;
void                          SetDataPath (TString fDataPath )    {DataPath = fDataPath;} ;
void                          SetFileName (TString fFileName )    {Filename = fFileName;} ;
void                        SetNMaxEVents (int fNMAXevents )      {NMaxEvents = fNMAXevents;};
void                         SetInclusive (int fInclusive )       {inclusive = fInclusive;};
void                           ReadEvents ();
void                        SetInputTTree ();
void                  InitializeVariables ();
void                         ReadRGAEvent (int event);
void                    ComputeKinematics ();
void           ExtractElectronInformation ();
void              ExtractPionsInformation ();
void               ExtractPipsInformation (int pipsIdx);
void               ExtractPimsInformation (int pimsIdx);
void                   WriteEventToOutput ();
bool   CheckIfElectronPassedSelectionCuts (Int_t fid1,
                                           Int_t fid2,
                                           Int_t fid3,
                                           TLorentzVector e,
                                           TVector3 Ve);

bool      CheckIfPionPassedSelectionCuts (TString pionCharge, // "pi+" or "pi-"
                                          Int_t fid1,
                                          Int_t fid2,
                                          Int_t fid3,
                                          Double_t chi2PID,
                                          Double_t p,
                                          TVector3 Ve,
                                          TVector3 Vpi);
void             Stream_e_pi_line_to_CSV (TString pionCharge,
                                          int piIdx,
                                          bool passed_cuts_e_pi,
                                          bool passed_cuts_e_pi_kinematics);
void               ComputePionKinematics (TLorentzVector pi,
                                          TLorentzVector pi_qFrame);

void                       MoveTo_qFrame ();
TVector3           RotateVectorTo_qFrame (TVector3 v);
void                       GetBeamEnergy ();
void                 PrintRGAInformation ();
void                    SetPrintProgress ( int fPrintProgress ){
    PrintProgress = fPrintProgress;
};

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
// Main functionality
// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void ReadIgorRGAFile(TString fFileName="ntupleNew",
                     int fNMAXevents=-1,
                     int ffdebug=1,
                     int fPrintProgress=50000,
                     TString fDataPath="/Users/erezcohen/Desktop/data/BAND/RGA_Free_proton/",
                     int fInclusive=0,
                     int torusBending=-1
                     ){
    
    SetVerbosity          ( ffdebug );
    SetDataPath           ( fDataPath );
    SetFileName           ( fFileName );
    SetNMaxEVents         ( fNMAXevents );
    SetInclusive          ( fInclusive );
    SetPrintProgress      ( fPrintProgress );
    // read cut values
    aux.loadCutValues     ("macros/cuts/BANDcutValues.csv",torusBending);
    
    InitializeFileReading ();
    OpenInputFile         ();
    OpenOutputFiles       ();
    ReadEvents            ();
    CloseOutputFiles      ();
    CloseInputFile        ();
    PrintDone             ();
}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void InitializeVariables (){
    xB          = Q2        = omega     = -9999;
    W           = W_d                   = -9999;
    //    xF          = y         = M_X       = -9999;
    //    e_E_ECIN    = e_E_ECOUT = e_E_PCAL  = -9999;
    //    e_PCAL_W    = e_PCAL_V              = -9999;
    //    e_PCAL_x    = e_PCAL_y  = e_PCAL_z  = -9999;
    //    e_PCAL_sector                       = -9999;
    //    e_DC_sector = e_DC_Chi2N            = -9999;
    //    for (int regionIdx=0; regionIdx<3; regionIdx++) {
    //        e_DC_x[regionIdx]               = -9999;
    //        e_DC_y[regionIdx]               = -9999;
    //        e_DC_z[regionIdx]               = -9999;
    //    }
    //    Pe_phi = q_phi = q_theta            = 0;
    //
    //    Ve                                  = TVector3();
    //    ePastCutsInEvent                    = false;
    //
    piplus          .clear();
    piminus         .clear();
    piplus_qFrame   .clear();
    piminus_qFrame  .clear();
    //    Vpiplus     .clear();
    //    Vpiminus    .clear();
    //    pipluses    .clear();
    //    piminuses   .clear();
    //    electrons   .clear();
    //    neutrons    .clear();
    //    protons     .clear();
    //    gammas      .clear();
    for (int piIdx=0; piIdx<NMAXPIONS; piIdx++) {
        //        pips_chi2PID[piIdx]                         = -9999;
        //        pips_DC_sector[piIdx]                       = -9999;
        //        pips_PCAL_sector[piIdx]                     = -9999;
        //        pips_PCAL_W[piIdx] = pips_PCAL_V[piIdx]     = -9999;
        //        pips_PCAL_x[piIdx] = pips_PCAL_y[piIdx]     = -9999;
        //        pips_PCAL_z[piIdx]                          = -9999;
        //        pips_E_PCAL[piIdx]                          = -9999;
        //        pips_E_ECIN[piIdx] = pips_E_ECOUT[piIdx]    = -9999;
        //
        //        pims_chi2PID[piIdx]                         = -9999;
        //        pims_DC_sector[piIdx]                       = -9999;
        //        pims_PCAL_sector[piIdx]                     = -9999;
        //        pims_PCAL_W[piIdx] = pims_PCAL_V[piIdx]     = -9999;
        //        pims_PCAL_x[piIdx] = pims_PCAL_y[piIdx]     = -9999;
        //        pims_PCAL_z[piIdx]                          = -9999;
        //        pims_E_PCAL[piIdx]                          = -9999;
        //        pims_E_ECIN[piIdx] = pims_E_ECOUT[piIdx]    = -9999;
        //        for (int regionIdx=0; regionIdx<3; regionIdx++) {
        //            pips_DC_x[piIdx][regionIdx]= pips_DC_y[piIdx][regionIdx]    = -9999;
        //            pips_DC_z[piIdx][regionIdx]                                 = -9999;
        //            pims_DC_x[piIdx][regionIdx]= pims_DC_y[piIdx][regionIdx]    = -9999;
        //            pims_DC_z[piIdx][regionIdx]                                 = -9999;
        //        }
        piplus       .push_back( TLorentzVector(0,0,0,aux.Mpips) );
        piplus_qFrame.push_back( TLorentzVector(0,0,0,aux.Mpips) );
        Vpiplus .push_back( TVector3() );
        pipsPastSelectionCuts[piIdx]                = false;
        eepipsPastKinematicalCuts[piIdx]            = false;
        
        piminus       .push_back( TLorentzVector(0,0,0,aux.Mpims) );
        piminus_qFrame.push_back( TLorentzVector(0,0,0,aux.Mpims) );
        Vpiminus.push_back( TVector3() );
        pimsPastSelectionCuts[piIdx]                = false;
        eepimsPastKinematicalCuts[piIdx]            = false;
        
        piplus_Px[piIdx]    = piplus_Py[piIdx]  = piplus_Pz[piIdx]  = piplus_E[piIdx]   = -9999;
        piminus_Px[piIdx]   = piminus_Py[piIdx] = piminus_Pz[piIdx] = piminus_E[piIdx]  = -9999;
        Vpiplus_X[piIdx]    = Vpiplus_Y[piIdx]  = Vpiplus_Z[piIdx]                      = -9999;
        Vpiminus_X[piIdx]   = Vpiminus_Y[piIdx] = Vpiminus_Z[piIdx]                     = -9999;
        piplus_qFrame_pT[piIdx]    = piplus_qFrame_pL[piIdx]                            = -9999;
        piminus_qFrame_pT[piIdx]   = piminus_qFrame_pL[piIdx]                           = -9999;
        
    }
    //    DC_layer                                        = -9999;
    status                                          = 1; // 0 is good...
    
    pipsPastCutsInEvent                             = false;
    eepipsPastCutsInEvent                           = false;
    pimsPastCutsInEvent                             = false;
    eepimsPastCutsInEvent                           = false;
}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void ReadRGAEvent(int event){
    RGATree->GetEntry(event);
    Ne    = nml;
    Npips = npiP;
    Npims = npiM;
    
    // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    // Aug-2022
    // problem in Igor' file: no run number and no event number
    //
    // to solve this temporarily I assign the
    // event number to be the TTree entry number here
    //
    if (evnum==0) evnum=event;
    //
    //
    // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ComputeKinematics(){
    // compute event kinematics (from e-only information)
    // assign kinematics to leading electron
    q       .SetPxPyPzE( qx[e_idx], qy[e_idx], qz[e_idx], Nu[e_idx] );
    GetBeamEnergy();
    Beam    .SetPxPyPzE (0, 0, Ebeam, Ebeam );
    
    Q2      = -q.Mag2();
    omega   = q.E();
    xB      = Q2/(2. * aux.Mp * q.E());
    
    //    W2      = (target + q).Mag2(); //    W2      = Mp2 - Q2 + 2. * omega * Mp;
    //    W       = sqrt(W2);
    y       = omega / Ebeam;
    
    W       = sqrt((p_rest + q).Mag2());
    W_d     = sqrt((d_rest + q).Mag2());
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void GetBeamEnergy (){
    Ebeam = Ep[e_idx] + Nu[e_idx];
}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void ReadEvents (){
    if (fdebug>1) {
        std::cout << "------------------------------------------------" << std::endl;
        std::cout << "ReadEvents ("<<NeventsInRGATree<<" events)" << std::endl;
        std::cout << "------------------------------------------------" << std::endl;
    }
    
    for (int event=0; event < NeventsInRGATree ; event++) {
        if ( (NMaxEvents>0) && (event>=NMaxEvents) ){
            if (fdebug>1) {
                std::cout
                << std::endl
                << "Stop reading events at event " << event << std::endl;
            }
            break;
        }
        
        
        if (fdebug>1) {
            // start event
            std::cout
            << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
            << std::endl;
        }
        
        InitializeVariables ();
        ReadRGAEvent        (event);
        
        if (fdebug>1){
            std::cout
            << "run "    << runnum      << ","
            << "event "  << event       << "/" << NeventsInRGATree << ","
            << "evnum "  << evnum       << ","
            << std::endl
            << "N(π+): " << Npips       << ","
            << "N(π-): " << Npims       << ","
            << std::setprecision(3)
            << "E(beam): " << Ebeam     << " GeV,"
            << std::endl;
        }
        // filter events, extract information, and compute event kinematics:
        // we keep only d(e,e’pi+)X and d(e,e’pi-)X events
        if(  0 < Ne // after studying some MC and data, we need to kill events with more than 1 electron
           &&
           (inclusive == 1 || (0 < Npips) || (0 < Npims)) // "inclusive" means (e,e') events
           &&
           (Npips < NMAXPIONS) && (Npims < NMAXPIONS) // we don't want to crash the memeory
           ){
            
            ExtractElectronInformation  ();
            ComputeKinematics           ();
            ExtractPionsInformation     ();
            MoveTo_qFrame               ();
            WriteEventToOutput          ();
            
            if (fdebug>1){
                if (fdebug>3){
                    PrintRGAInformation();
                }
                
                if (fdebug>4){
                    std::cout
                    << "q: ("    << q.Px() << "," << q.Py()<< "," << q.Pz() << ";" << q.E()<<") "
                    << std::endl
                    << "Q2: " << Q2 << " (GeV/c)2, "
                    << "xB: " << xB << ", "
                    << std::endl;
                }
            }
        } else {
            if (fdebug>1) {
                std::cout << "Skipped computations in this event as there are not enough particles: "
                << "Ne = " << Ne << ",Npips = " << Npips << ",Npims = " << Npims << std::endl ;
            }
        }
        if (fdebug>1) {
            // finish event
            std::cout
            << "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv"
            << std::endl<< std::endl;
        }
        if (fdebug && event%PrintProgress==0){
            std::cout << std::setprecision(1) << " event " << event << std::endl;
        }
    }
    if (fdebug>1) std::cout << "Done ReadEvents ()" << std::endl;
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void SetInputTTree (){
    NeventsInRGATree = RGATree->GetEntries();
    RGATree  -> SetBranchAddress("runNum"           ,&runnum            );
    RGATree  -> SetBranchAddress("evNum"            ,&evnum             );
    RGATree  -> SetBranchAddress("nml"              ,&nml               );
    
    RGATree  -> SetBranchAddress("Epid"             ,Epid              );
    RGATree  -> SetBranchAddress("Ep"               ,Ep                );
    RGATree  -> SetBranchAddress("Epx"              ,Epx               );
    RGATree  -> SetBranchAddress("Epy"              ,Epy               );
    RGATree  -> SetBranchAddress("Epz"              ,Epz               );
    RGATree  -> SetBranchAddress("Evz"              ,Evz               );
    RGATree  -> SetBranchAddress("Etheta"           ,Etheta            );
    RGATree  -> SetBranchAddress("Ephi"             ,Ephi              );
    RGATree  -> SetBranchAddress("Efiducial1"       ,Efiducial1        );
    RGATree  -> SetBranchAddress("Efiducial2"       ,Efiducial2        );
    RGATree  -> SetBranchAddress("Efiducial3"       ,Efiducial3        );
    RGATree  -> SetBranchAddress("Efiducial3"       ,Efiducial3        );
    RGATree  -> SetBranchAddress("ePCAL_energy"     ,ePCAL_energy      );
    RGATree  -> SetBranchAddress("ePCAL_x"          ,ePCAL_x           );
    RGATree  -> SetBranchAddress("ePCAL_y"          ,ePCAL_y           );
    RGATree  -> SetBranchAddress("ePCAL_z"          ,ePCAL_z           );
    RGATree  -> SetBranchAddress("ePCAL_lu"         ,ePCAL_lu          );
    RGATree  -> SetBranchAddress("ePCAL_lv"         ,ePCAL_lv          );
    RGATree  -> SetBranchAddress("ePCAL_lw"         ,ePCAL_lw          );
    RGATree  -> SetBranchAddress("eECIN_energy"     ,eECIN_energy      );
    RGATree  -> SetBranchAddress("eECOUT_energy"    ,eECOUT_energy     );
    RGATree  -> SetBranchAddress("Esector"          ,Esector           );

    // RGATree  -> SetBranchAddress("nPi"              ,&nPi              );
    
    RGATree  -> SetBranchAddress("npiP"               ,&npiP               );
    RGATree  -> SetBranchAddress("piPpid"             ,piPpid              );
    RGATree  -> SetBranchAddress("piPp"               ,piPp                );
    RGATree  -> SetBranchAddress("piPpx"              ,piPpx               );
    RGATree  -> SetBranchAddress("piPpy"              ,piPpy               );
    RGATree  -> SetBranchAddress("piPpz"              ,piPpz               );
    RGATree  -> SetBranchAddress("piPvz"              ,piPvz               );
    RGATree  -> SetBranchAddress("piPtheta"           ,piPtheta            );
    RGATree  -> SetBranchAddress("piPphi"             ,piPphi              );
    RGATree  -> SetBranchAddress("piPfiducial1"       ,piPfiducial1        );
    RGATree  -> SetBranchAddress("piPfiducial2"       ,piPfiducial2        );
    RGATree  -> SetBranchAddress("piPfiducial3"       ,piPfiducial3        );
    RGATree  -> SetBranchAddress("piPfiducial3"       ,piPfiducial3        );
    RGATree  -> SetBranchAddress("piPPCAL_energy"     ,piPPCAL_energy      );
    RGATree  -> SetBranchAddress("piPPCAL_x"          ,piPPCAL_x           );
    RGATree  -> SetBranchAddress("piPPCAL_y"          ,piPPCAL_y           );
    RGATree  -> SetBranchAddress("piPPCAL_z"          ,piPPCAL_z           );
    RGATree  -> SetBranchAddress("piPPCAL_lu"         ,piPPCAL_lu          );
    RGATree  -> SetBranchAddress("piPPCAL_lv"         ,piPPCAL_lv          );
    RGATree  -> SetBranchAddress("piPPCAL_lw"         ,piPPCAL_lw          );
    RGATree  -> SetBranchAddress("piPECIN_energy"     ,piPECIN_energy      );
    RGATree  -> SetBranchAddress("piPECOUT_energy"    ,piPECOUT_energy     );
    RGATree  -> SetBranchAddress("piPsector"          ,piPsector           );
    
    RGATree  -> SetBranchAddress("npiM"               ,&npiM               );
    RGATree  -> SetBranchAddress("piMpid"             ,piMpid              );
    RGATree  -> SetBranchAddress("piMp"               ,piMp                );
    RGATree  -> SetBranchAddress("piMpx"              ,piMpx               );
    RGATree  -> SetBranchAddress("piMpy"              ,piMpy               );
    RGATree  -> SetBranchAddress("piMpz"              ,piMpz               );
    RGATree  -> SetBranchAddress("piMvz"              ,piMvz               );
    RGATree  -> SetBranchAddress("piMtheta"           ,piMtheta            );
    RGATree  -> SetBranchAddress("piMphi"             ,piMphi              );
    RGATree  -> SetBranchAddress("piMfiducial1"       ,piMfiducial1        );
    RGATree  -> SetBranchAddress("piMfiducial2"       ,piMfiducial2        );
    RGATree  -> SetBranchAddress("piMfiducial3"       ,piMfiducial3        );
    RGATree  -> SetBranchAddress("piMfiducial3"       ,piMfiducial3        );
    RGATree  -> SetBranchAddress("piMPCAL_energy"     ,piMPCAL_energy      );
    RGATree  -> SetBranchAddress("piMPCAL_x"          ,piMPCAL_x           );
    RGATree  -> SetBranchAddress("piMPCAL_y"          ,piMPCAL_y           );
    RGATree  -> SetBranchAddress("piMPCAL_z"          ,piMPCAL_z           );
    RGATree  -> SetBranchAddress("piMPCAL_lu"         ,piMPCAL_lu          );
    RGATree  -> SetBranchAddress("piMPCAL_lv"         ,piMPCAL_lv          );
    RGATree  -> SetBranchAddress("piMPCAL_lw"         ,piMPCAL_lw          );
    RGATree  -> SetBranchAddress("piMECIN_energy"     ,piMECIN_energy      );
    RGATree  -> SetBranchAddress("piMECOUT_energy"    ,piMECOUT_energy     );
    RGATree  -> SetBranchAddress("piMsector"          ,piMsector           );
    
    
    RGATree  -> SetBranchAddress("Nu"               ,Nu                );
    RGATree  -> SetBranchAddress("qx"               ,qx                );
    RGATree  -> SetBranchAddress("qy"               ,qy                );
    RGATree  -> SetBranchAddress("qz"               ,qz                );
    
    //    ******************************************************************************
    //    *Tree    :T         : e'pion                                                 *
    //    *Br    1 :runNum    : runNum/L                                               *
    //    *Br    2 :evNum     : evNum/L                                                *
    //    *Br    4 :Epid      : Epid[nml]/I                                            *
    //    *Br    5 :Echarge   : Echarge[nml]/I                                         *
    //    *Br    6 :Ep        : Ep[nml]/F                                              *
    //    *Br    7 :Epx       : Epx[nml]/F                                             *
    //    *Br    8 :Epy       : Epy[nml]/F                                             *
    //    *Br    9 :Epz       : Epz[nml]/F                                             *
    //    *Br   10 :Etheta    : Etheta[nml]/F                                          *
    //    *Br   11 :Ephi      : Ephi[nml]/F                                            *
    //    *Br   12 :Etime     : Etime[nml]/F                                           *
    //    *Br   13 :Epath     : Epath[nml]/F                                           *
    //    *Br   14 :Evz       : Evz[nml]/F                                             *
    //    *Br   15 :Estat     : Estat[nml]/I                                           *
    //    *Br   16 :EdE       : EdE[nml]/F                                             *
    //    *Br   17 :Ebeta     : Ebeta[nml]/F                                           *
    //    *Br   18 :Esector   : Esector[nml]/I                                         *
    //    *Br   19 :Eregion   : Eregion[nml]/I                                         *
    //    *Br   20 :Efiducial1 : Efiducial1[nml]/I                                     *
    //    *Br   21 :Efiducial2 : Efiducial2[nml]/I                                     *
    //    *Br   22 :Efiducial3 : Efiducial3[nml]/I                                     *
    //    *Br   23 :ePCAL_time : ePCAL_time[nml]/F                                     *
    //    *Br   24 :ePCAL_energy : ePCAL_energy[nml]/F                                 *
    //    *Br   25 :ePCAL_path : ePCAL_path[nml]/F                                     *
    //    *Br   26 :ePCAL_sector : ePCAL_sector[nml]/I                                 *
    //    *Br   27 :ePCAL_x   : ePCAL_x[nml]/F                                         *
    //    *Br   28 :ePCAL_y   : ePCAL_y[nml]/F                                         *
    //    *Br   29 :ePCAL_z   : ePCAL_z[nml]/F                                         *
    //    *Br   30 :ePCAL_lu  : ePCAL_lu[nml]/F                                        *
    //    *Br   31 :ePCAL_lv  : ePCAL_lv[nml]/F                                        *
    //    *Br   32 :ePCAL_lw  : ePCAL_lw[nml]/F                                        *
    //    *Br   33 :ePCAL_status : ePCAL_status[nml]/I                                 *
    //    *Br   34 :eECIN_time : eECIN_time[nml]/F                                     *
    //    *Br   35 :eECIN_energy : eECIN_energy[nml]/F                                 *
    //    *Br   36 :eECIN_path : eECIN_path[nml]/F                                     *
    //    *Br   37 :eECIN_sector : eECIN_sector[nml]/I                                 *
    //    *Br   38 :eECIN_x   : eECIN_x[nml]/F                                         *
    //    *Br   39 :eECIN_y   : eECIN_y[nml]/F                                         *
    //    *Br   40 :eECIN_z   : eECIN_z[nml]/F                                         *
    //    *Br   41 :eECIN_lu  : eECIN_lu[nml]/F                                        *
    //    *Br   42 :eECIN_lv  : eECIN_lv[nml]/F                                        *
    //    *Br   43 :eECIN_lw  : eECIN_lw[nml]/F                                        *
    //    *Br   44 :eECIN_status : eECIN_status[nml]/I                                 *
    //    *Br   45 :eECOUT_time : eECOUT_time[nml]/F                                   *
    //    *Br   46 :eECOUT_energy : eECOUT_energy[nml]/F                               *
    //    *Br   47 :eECOUT_path : eECOUT_path[nml]/F                                   *
    //    *Br   48 :eECOUT_sector : eECOUT_sector[nml]/I                               *
    //    *Br   49 :eECOUT_x  : eECOUT_x[nml]/F                                        *
    //    *Br   50 :eECOUT_y  : eECOUT_y[nml]/F                                        *
    //    *Br   51 :eECOUT_z  : eECOUT_z[nml]/F                                        *
    //    *Br   52 :eECOUT_lu : eECOUT_lu[nml]/F                                       *
    //    *Br   53 :eECOUT_lv : eECOUT_lv[nml]/F                                       *
    //    *Br   54 :eECOUT_lw : eECOUT_lw[nml]/F                                       *
    //    *Br   55 :eECOUT_status : eECOUT_status[nml]/I                               *
    //    *Br   56 :EdcX      : EdcX[nml][3]/F                                         *
    //    *Br   57 :EdcY      : EdcY[nml][3]/F                                         *
    //    *Br   58 :EdcZ      : EdcZ[nml][3]/F                                         *
    //    *Br   59 :npiP      : npiP/I                                                 *
    //    *Br   60 :piPpid    : piPpid[npiP]/I                                         *
    //    *Br   61 :piPcharge : piPcharge[npiP]/I                                      *
    //    *Br   62 :piPp      : piPp[npiP]/F                                           *
    //    *Br   63 :piPpx     : piPpx[npiP]/F                                          *
    //    *Br   64 :piPpy     : piPpy[npiP]/F                                          *
    //    *Br   65 :piPpz     : piPpz[npiP]/F                                          *
    //    *Br   66 :piPtheta  : piPtheta[npiP]/F                                       *
    //    *Br   67 :piPphi    : piPphi[npiP]/F                                         *
    //    *Br   68 :piPtime   : piPtime[npiP]/F                                        *
    //    *Br   69 :piPpath   : piPpath[npiP]/F                                        *
    //    *Br   70 :piPvz     : piPvz[npiP]/F                                          *
    //    *Br   71 :piPstat   : piPstat[npiP]/I                                        *
    //    *Br   72 :piPdE     : piPdE[npiP]/F                                          *
    //    *Br   73 :piPbeta   : piPbeta[npiP]/F                                        *
    //    *Br   74 :piPsector : piPsector[npiP]/I                                      *
    //    *Br   75 :piPregion : piPregion[npiP]/I                                      *
    //    *Br   76 :piPfiducial1 : piPfiducial1[npiP]/I                                *
    //    *Br   77 :piPfiducial2 : piPfiducial2[npiP]/I                                *
    //    *Br   78 :piPfiducial3 : piPfiducial3[npiP]/I                                *
    //    *Br   79 :piPPCAL_time : piPPCAL_time[npiP]/F                                *
    //    *Br   80 :piPPCAL_energy : piPPCAL_energy[npiP]/F                            *
    //    *Br   81 :piPPCAL_path : piPPCAL_path[npiP]/F                                *
    //    *Br   82 :piPPCAL_sector : piPPCAL_sector[npiP]/I                            *
    //    *Br   83 :piPPCAL_x : piPPCAL_x[npiP]/F                                      *
    //    *Br   84 :piPPCAL_y : piPPCAL_y[npiP]/F                                      *
    //    *Br   85 :piPPCAL_z : piPPCAL_z[npiP]/F                                      *
    //    *Br   86 :piPPCAL_lu : piPPCAL_lu[npiP]/F                                    *
    //    *Br   87 :piPPCAL_lv : piPPCAL_lv[npiP]/F                                    *
    //    *Br   88 :piPPCAL_lw : piPPCAL_lw[npiP]/F                                    *
    //    *Br   89 :piPPCAL_status : piPPCAL_status[npiP]/I                            *
    //    *Br   90 :piPECIN_time : piPECIN_time[npiP]/F                                *
    //    *Br   91 :piPECIN_energy : piPECIN_energy[npiP]/F                            *
    //    *Br   92 :piPECIN_path : piPECIN_path[npiP]/F                                *
    //    *Br   93 :piPECIN_sector : piPECIN_sector[npiP]/I                            *
    //    *Br   94 :piPECIN_x : piPECIN_x[npiP]/F                                      *
    //    *Br   95 :piPECIN_y : piPECIN_y[npiP]/F                                      *
    //    *Br   96 :piPECIN_z : piPECIN_z[npiP]/F                                      *
    //    *Br   97 :piPECIN_lu : piPECIN_lu[npiP]/F                                    *
    //    *Br   98 :piPECIN_lv : piPECIN_lv[npiP]/F                                    *
    //    *Br   99 :piPECIN_lw : piPECIN_lw[npiP]/F                                    *
    //    *Br  100 :piPECIN_status : piPECIN_status[npiP]/I                            *
    //    *Br  101 :piPECOUT_time : piPECOUT_time[npiP]/F                              *
    //    *Br  102 :piPECOUT_energy : piPECOUT_energy[npiP]/F                          *
    //    *Br  103 :piPECOUT_path : piPECOUT_path[npiP]/F                              *
    //    *Br  104 :piPECOUT_sector : piPECOUT_sector[npiP]/I                          *
    //    *Br  105 :piPECOUT_x : piPECOUT_x[npiP]/F                                    *
    //    *Br  106 :piPECOUT_y : piPECOUT_y[npiP]/F                                    *
    //    *Br  107 :piPECOUT_z : piPECOUT_z[npiP]/F                                    *
    //    *Br  108 :piPECOUT_lu : piPECOUT_lu[npiP]/F                                  *
    //    *Br  109 :piPECOUT_lv : piPECOUT_lv[npiP]/F                                  *
    //    *Br  110 :piPECOUT_lw : piPECOUT_lw[npiP]/F                                  *
    //    *Br  111 :piPECOUT_status : piPECOUT_status[npiP]/I                          *
    //    *Br  112 :piPdcX    : piPdcX[npiP][3]/F                                      *
    //    *Br  113 :piPdcY    : piPdcY[npiP][3]/F                                      *
    //    *Br  114 :piPdcZ    : piPdcZ[npiP][3]/F                                      *
    //    *Br  115 :npiM      : npiM/I                                                 *
    //    *Br  116 :piMpid    : piPpid[npiM]/I                                         *
    //    *Br  117 :piMcharge : piPcharge[npiM]/I                                      *
    //    *Br  118 :piMp      : piMp[npiM]/F                                           *
    //    *Br  119 :piMpx     : piMpx[npiM]/F                                          *
    //    *Br  120 :piMpy     : piMpy[npiM]/F                                          *
    //    *Br  121 :piMpz     : piMpz[npiM]/F                                          *
    //    *Br  122 :piMtheta  : piMtheta[npiM]/F                                       *
    //    *Br  123 :piMphi    : piMphi[npiM]/F                                         *
    //    *Br  124 :piMtime   : piMtime[npiM]/F                                        *
    //    *Br  125 :piMpath   : piMpath[npiM]/F                                        *
    //    *Br  126 :piMvz     : piMvz[npiM]/F                                          *
    //    *Br  127 :piMstat   : piMstat[npiM]/I                                        *
    //    *Br  128 :piMdE     : piMdE[npiM]/F                                          *
    //    *Br  129 :piMbeta   : piMbeta[npiM]/F                                        *
    //    *Br  130 :piMsector : piMsector[npiM]/I                                      *
    //    *Br  131 :piMregion : piMregion[npiM]/I                                      *
    //    *Br  132 :piMfiducial1 : piMfiducial1[npiM]/I                                *
    //    *Br  133 :piMfiducial2 : piMfiducial2[npiM]/I                                *
    //    *Br  134 :piMfiducial3 : piMfiducial3[npiM]/I                                *
    //    *Br  135 :piMPCAL_time : piMPCAL_time[npiM]/F                                *
    //    *Br  136 :piMPCAL_energy : piMPCAL_energy[npiM]/F                            *
    //    *Br  137 :piMPCAL_path : piMPCAL_path[npiM]/F                                *
    //    *Br  138 :piMPCAL_sector : piMPCAL_sector[npiM]/I                            *
    //    *Br  139 :piMPCAL_x : piMPCAL_x[npiM]/F                                      *
    //    *Br  140 :piMPCAL_y : piMPCAL_y[npiM]/F                                      *
    //    *Br  141 :piMPCAL_z : piMPCAL_z[npiM]/F                                      *
    //    *Br  142 :piMPCAL_lu : piMPCAL_lu[npiM]/F                                    *
    //    *Br  143 :piMPCAL_lv : piMPCAL_lv[npiM]/F                                    *
    //    *Br  144 :piMPCAL_lw : piMPCAL_lw[npiM]/F                                    *
    //    *Br  145 :piMPCAL_status : piMPCAL_status[npiM]/I                            *
    //    *Br  146 :piMECIN_time : piMECIN_time[npiM]/F                                *
    //    *Br  147 :piMECIN_energy : piMECIN_energy[npiM]/F                            *
    //    *Br  148 :piMECIN_path : piMECIN_path[npiM]/F                                *
    //    *Br  149 :piMECIN_sector : piMECIN_sector[npiM]/I                            *
    //    *Br  150 :piMECIN_x : piMECIN_x[npiM]/F                                      *
    //    *Br  151 :piMECIN_y : piMECIN_y[npiM]/F                                      *
    //    *Br  152 :piMECIN_z : piMECIN_z[npiM]/F                                      *
    //    *Br  153 :piMECIN_lu : piMECIN_lu[npiM]/F                                    *
    //    *Br  154 :piMECIN_lv : piMECIN_lv[npiM]/F                                    *
    //    *Br  155 :piMECIN_lw : piMECIN_lw[npiM]/F                                    *
    //    *Br  156 :piMECIN_status : piMECIN_status[npiM]/I                            *
    //    *Br  157 :piMECOUT_time : piMECOUT_time[npiM]/F                              *
    //    *Br  158 :piMECOUT_energy : piMECOUT_energy[npiM]/F                          *
    //    *Br  159 :piMECOUT_path : piMECOUT_path[npiM]/F                              *
    //    *Br  160 :piMECOUT_sector : piMECOUT_sector[npiM]/I                          *
    //    *Br  161 :piMECOUT_x : piMECOUT_x[npiM]/F                                    *
    //    *Br  162 :piMECOUT_y : piMECOUT_y[npiM]/F                                    *
    //    *Br  163 :piMECOUT_z : piMECOUT_z[npiM]/F                                    *
    //    *Br  164 :piMECOUT_lu : piMECOUT_lu[npiM]/F                                  *
    //    *Br  165 :piMECOUT_lv : piMECOUT_lv[npiM]/F                                  *
    //    *Br  166 :piMECOUT_lw : piMECOUT_lw[npiM]/F                                  *
    //    *Br  167 :piMECOUT_status : piMECOUT_status[npiM]/I                          *
    //    *Br  168 :piMdcX    : piMdcX[npiM][3]/F                                      *
    //    *Br  169 :piMdcY    : piMdcY[npiM][3]/F                                      *
    //    *Br  170 :piMdcZ    : piMdcZ[npiM][3]/F                                      *
    //    *Br  171 :Q2        : Q2[nml]/F                                              *
    //    *Br  172 :Nu        : Nu[nml]/F                                              *
    //    *Br  173 :q         : q[nml]/F                                               *
    //    *Br  174 :qx        : qx[nml]/F                                              *
    //    *Br  175 :qy        : qy[nml]/F                                              *
    //    *Br  176 :qz        : qz[nml]/F                                              *
    //    *Br  177 :xB        : xB[nml]/F                                              *
    //    *Br  178 :W2        : W2[nml]/F                                              *
    
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void PrintDone(){
    std::cout << "Done reading RGA file " << Filename << ", execution time: "
    << double(clock() - tStart) / (double)CLOCKS_PER_SEC
    << " sec "<< std::endl;
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void OpenInputFile (){
    
    TString InputROOTFilename = DataPath + "/" + Filename + ".root";
    if (fdebug>2) std::cout << "Opening " << InputROOTFilename << std::endl;
    RGAFile = new TFile( InputROOTFilename );
    
    RGATree = (TTree*)RGAFile->Get("T");
    
    if (fdebug>4) RGATree->Print();
    
    NeventsInRGATree = RGATree->GetEntries();
    SetInputTTree();
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void OpenOutputFiles (){
    
    CSVFilename_e_e_pips = DataPath + "/" + Filename + "_e_e_piplus";
    CSVFilename_e_e_pims = DataPath + "/" + Filename + "_e_e_piminus";
    
    if (fdebug>2) {
        std::cout << "Opening output files: "
        << CSVFilename_e_e_pips
        << " and "
        << CSVFilename_e_e_pims
        << std::endl
        << csvheader
        << std::endl;
    }
    
    // Write csv header output csv files
    CSVfile_e_e_pips.open( CSVFilename_e_e_pips + "_selected_eepi_kinematics.csv" );
    CSVfile_e_e_pips << csvheader << std::endl;
    
    // Write csv header output csv files
    CSVfile_e_e_pims.open( CSVFilename_e_e_pims + "_selected_eepi_kinematics.csv" );
    CSVfile_e_e_pims << csvheader << std::endl;
    
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void InitializeFileReading (){
    Nevents_e_e_pips = Nevents_e_e_pims = 0;
    NeventsInRGATree = 0;
    target  .SetXYZM    (0, 0, 0,     aux.Md    );
    d_rest  .SetXYZM    (0, 0, 0,     aux.Md    );
    p_rest  .SetXYZM    (0, 0, 0,     aux.Mp    );
    
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void CloseInputFile (){
    RGAFile->Close();
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void CloseOutputFiles (){
    
    // close output CSV
    CSVfile_e_e_pips.close();
    CSVfile_e_e_pims.close();
    
    std::cout << "Done producing output file. They are ready in " << std::endl << DataPath << std::endl;
    std::cout << "Wrote "
    << Nevents_e_e_pips << " (e,e'π+) events"
    << " and "
    << Nevents_e_e_pims << " (e,e'π-) events." << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExtractElectronInformation(){
    // ------------------------------------------------------------------------------------------------
    // extract electron information
    // ------------------------------------------------------------------------------------------------
    
    // use only leading (fastest) electron
    // find the leading electron
    e_idx = 0;
    if (nml>1){
        Double_t Pe_lead = 0;
        for (int e_idx_tmp=0; e_idx_tmp<nml; e_idx_tmp++){
            if (Ep[e_idx_tmp] > Pe_lead){
                e_idx   = e_idx_tmp;
                Pe_lead = Ep[e_idx_tmp];
            }
        }
    }
    e.SetXYZM( Epx[e_idx], Epy[e_idx], Epz[e_idx], aux.Me );
    Ve = TVector3( 0, 0, Evz[e_idx] );
    e_DC_sector = Esector[e_idx];
    
    if (fdebug>2) {
        std::cout
        << "p(e): "         << e.P()        << ","
        << "Vz(e): "        << Ve.Z()       << ","
        << "DC-sector(e): " << e_DC_sector  << ","
        << std::endl;
    }
    
    //    // detector information on electron
    //    auto e_PCAL_info= electrons[leading_e_index]->cal(PCAL);
    //    e_E_PCAL        = e_PCAL_info->getEnergy();
    //    e_PCAL_sector   = e_PCAL_info->getSector();
    //    e_PCAL_V        = e_PCAL_info->getLv();
    //    e_PCAL_W        = e_PCAL_info->getLw();
    //    e_E_ECIN        = electrons[leading_e_index]->cal(ECIN)->getEnergy();
    //    e_E_ECOUT       = electrons[leading_e_index]->cal(ECOUT)->getEnergy();
    //
    //    // hit position in PCAL
    //    e_PCAL_x        = e_PCAL_info->getX();
    //    e_PCAL_y        = e_PCAL_info->getY();
    //    e_PCAL_z        = e_PCAL_info->getZ();
    //
    //    // Drift Chamber tracking system
    //    auto e_DC_info  = electrons[leading_e_index]->trk(DC);
    //    e_DC_sector     = e_DC_info->getSector(); // tracking sector
    //    e_DC_Chi2N      = e_DC_info->getChi2N();  // tracking chi^2/NDF
    //
    //    for (int regionIdx=0; regionIdx<3; regionIdx++) {
    //        int DC_layer = DC_layers[regionIdx];
    //        e_DC_x[regionIdx] = electrons[leading_e_index]->traj(DC,DC_layer)->getX();
    //        e_DC_y[regionIdx] = electrons[leading_e_index]->traj(DC,DC_layer)->getY();
    //        e_DC_z[regionIdx] = electrons[leading_e_index]->traj(DC,DC_layer)->getZ();
    //    }
    
    //    if (fdebug > 2) std::cout << "extracted electron information and computed kinematics" << std::endl;
    
    // ------------------------------------------------------------------------------------------------
    // now, check if electron passed event selection requirements
    // ------------------------------------------------------------------------------------------------
    //    ePastCutsInEvent = CheckIfElectronPassedSelectionCuts(e_PCAL_x, e_PCAL_y,
    //                                                              e_PCAL_W, e_PCAL_V,
    //                                                              e_E_PCAL, e_E_ECIN,
    //                                                              e_E_ECOUT,
    //                                                              e, Ve,
    //                                                              e_PCAL_sector, // e_PCAL_sector should be consistent with e_DC_sector
    //                                                              e_DC_x, e_DC_y, e_DC_z,
    //                                                              torusBending );
    ePastCutsInEvent = CheckIfElectronPassedSelectionCuts(Efiducial1[e_idx],
                                                          Efiducial2[e_idx],
                                                          Efiducial3[e_idx],
                                                          e,
                                                          Ve);
    
    if (ePastCutsInEvent)  Nevents_passed_e_cuts++ ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExtractPionsInformation(){
    
    // positive pions
    for (int pipsIdx=0; pipsIdx < Npips; pipsIdx++) {
        ExtractPipsInformation( pipsIdx );
    }
    // negative pions
    for (int pimsIdx=0; pimsIdx < Npims; pimsIdx++) {
        ExtractPimsInformation( pimsIdx );
    }
    
    
    // done
    if (fdebug > 2) std::cout << "done extracting pion information" << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MoveTo_qFrame(){
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
    
    if (fdebug>3){
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
        
        if (fdebug>1) aux.Print4Vector( Ppi_q, "π+(" + std::to_string(piIdx) + ")" );
        
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
        
        if (fdebug>1) aux.Print4Vector( Ppi_q, "π-(" + std::to_string(piIdx) + ")" );
    }
    if (fdebug>4){
        std::cout
        << "size(piplus_qFrame): "  << piplus_qFrame.size() << ","
        << "size(piminus_qFrame): " << piminus_qFrame.size()<< ","
        << std::endl;
    }
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void WriteEventToOutput(){
    if (fdebug>1) std::cout << "WriteEventToOutput()" << std::endl;
    // (Maybe) write this event to "selected events csv-file"
    bool IsSelected_eepi = false;
    
    //ePastCutsInEvent = true;
    if (inclusive == 1) {
        pipsPastCutsInEvent = true;
        pimsPastCutsInEvent = true;
        eepipsPastCutsInEvent = true;
        eepimsPastCutsInEvent = true;
    }
    
    //    if (fdebug>1) std::cout << "Fill 'no-cuts' TTrees even if nothing passed criteria" << std::endl;
    //    // Fill "no-cuts" TTrees even if nothing passed criteria
    //    outTree_e_piplus_no_cuts -> Fill();
    //    outTree_e_piminus_no_cuts -> Fill();
    
    if (fdebug>1) {
        std::cout
        << "ePastCutsInEvent: "    << ePastCutsInEvent << ","
        << "pipsPastCutsInEvent: " << pipsPastCutsInEvent
        << std::endl;
    }
    
    
    if (ePastCutsInEvent && pipsPastCutsInEvent) {
        IsSelected_eepi = true;
        //        outTree_e_piplus -> Fill();
        if (fdebug>3) std::cout << "Filling (e,e'π+) TTree with this event!" << std::endl;
        
        Nevents_passed_e_pips_cuts ++ ;
        if (eepipsPastCutsInEvent) Nevents_passed_e_pips_kinematics_cuts ++;
        
        for (int pipsIdx=0; pipsIdx<Npips; pipsIdx++) {
            Stream_e_pi_line_to_CSV( "pi+", pipsIdx,
                                    pipsPastSelectionCuts[pipsIdx],
                                    eepipsPastKinematicalCuts[pipsIdx]
                                    );
        }
    }
    
    if (ePastCutsInEvent && pimsPastCutsInEvent) {
        IsSelected_eepi = true;
        //        outTree_e_piminus -> Fill();
        if (fdebug>3) std::cout << "Filling (e,e'π-) TTree with this event!" << std::endl;
        Nevents_passed_e_pims_cuts ++ ;
        if (eepimsPastCutsInEvent) Nevents_passed_e_pims_kinematics_cuts ++;
        
        for (int pimsIdx=0; pimsIdx<Npims; pimsIdx++) {
            Stream_e_pi_line_to_CSV( "pi-", pimsIdx,
                                    pimsPastSelectionCuts[pimsIdx],
                                    eepimsPastKinematicalCuts[pimsIdx]
                                    );
        }
    }
    if (fdebug>1) std::cout << "Done WriteEventToOutput()" << std::endl;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExtractPipsInformation( int pipsIdx ){
    if (fdebug>2)
        std::cout << "ExtractPipsInformation(pipsIdx=" << pipsIdx << ")" << std::endl;
    
    
    // extract positive pion information
    piplus.at(pipsIdx).SetXYZM(piPpx[pipsIdx],
                               piPpy[pipsIdx],
                               piPpz[pipsIdx],
                               aux.Mpips);
    
    Zpips[pipsIdx]          = piplus[pipsIdx].E() / omega;
    Vpiplus[pipsIdx]        = TVector3( 0, 0 , piPvz[pipsIdx] );
    pips_DC_sector[pipsIdx] = piPsector[pipsIdx];
    
    if (fdebug>2) {
        std::cout
        << "p(π+ "          <<pipsIdx<<"): " << piplus[pipsIdx].P()     << ","
        << "z(π+ "          <<pipsIdx<<"): " << Zpips[pipsIdx]          << ","
        << "Vz(π+ "         <<pipsIdx<<"): " << Vpiplus[pipsIdx].Z()    << ","
        << "DC-sector(π+ "  <<pipsIdx<<"): " << pips_DC_sector[pipsIdx] << ","
        << std::endl;
    }
    
    //    Vpiplus[pipsIdx]            = GetParticleVertex( pipluses[pipsIdx] );
    //    pips_chi2PID[pipsIdx]       = pipluses[pipsIdx]->par()->getChi2Pid();
    //
    //    // EC in and out
    //    pips_E_ECIN[pipsIdx]        = pipluses[pipsIdx]->cal(ECIN)->getEnergy();
    //    pips_E_ECOUT[pipsIdx]       = pipluses[pipsIdx]->cal(ECOUT)->getEnergy();
    //    // PCAL
    //    auto pips_PCAL_info         = pipluses[pipsIdx]->cal(PCAL);
    //    pips_E_PCAL[pipsIdx]        = pips_PCAL_info->getEnergy();
    //    pips_PCAL_sector[pipsIdx]   = pips_PCAL_info->getSector();
    //    pips_PCAL_V[pipsIdx]        = pips_PCAL_info->getLv();
    //    pips_PCAL_W[pipsIdx]        = pips_PCAL_info->getLw();
    //    pips_PCAL_x[pipsIdx]        = pips_PCAL_info->getX();
    //    pips_PCAL_y[pipsIdx]        = pips_PCAL_info->getY();
    //    pips_PCAL_z[pipsIdx]        = pips_PCAL_info->getZ();
    //    // DC
    //    auto pips_DC_info           = pipluses[pipsIdx]->trk(DC);
    //    pips_DC_sector[pipsIdx]     = pips_DC_info->getSector(); // tracking sector
    //    pips_Chi2N[pipsIdx]         = pips_DC_info->getChi2N();  // tracking chi^2/NDF
    //    for (int regionIdx=0; regionIdx<3; regionIdx++) {
    //        DC_layer = DC_layers[regionIdx];
    //        pips_DC_x[pipsIdx][regionIdx] = pipluses[pipsIdx]->traj(DC,DC_layer)->getX();
    //        pips_DC_y[pipsIdx][regionIdx] = pipluses[pipsIdx]->traj(DC,DC_layer)->getY();
    //        pips_DC_z[pipsIdx][regionIdx] = pipluses[pipsIdx]->traj(DC,DC_layer)->getZ();
    //        if (fdebug>3) {
    //            std::cout
    //            << "pips_DC_sector[pipsIdx="<<pipsIdx<<"]="
    //            << pips_DC_sector[pipsIdx]
    //            << ", DC_layer = " << DC_layer
    //            << std::endl
    //            << "pips_DC_x[pipsIdx="<<pipsIdx<<"][regionIdx="<<regionIdx<<"]="
    //            << pips_DC_x[pipsIdx][regionIdx]
    //            << std::endl
    //            << "pips_DC_y[pipsIdx="<<pipsIdx<<"][regionIdx="<<regionIdx<<"]="
    //            << pips_DC_y[pipsIdx][regionIdx]
    //            << std::endl
    //            << "pips_DC_z[pipsIdx="<<pipsIdx<<"][regionIdx="<<regionIdx<<"]="
    //            << pips_DC_z[pipsIdx][regionIdx]
    //            << std::endl;
    //        }
    //    }
    // ------------------------------------------------------------------------------------------------
    // now, check if pion passed event selection requirements
    // ------------------------------------------------------------------------------------------------
    //    pipsPastSelectionCuts[pipsIdx] = CheckIfPionPassedSelectionCuts("pi+",
    //                                                                     pips_DC_sector[pipsIdx],
    //                                                                     pips_DC_x[pipsIdx],
    //                                                                     pips_DC_y[pipsIdx],
    //                                                                     pips_DC_z[pipsIdx],
    //                                                                     pips_chi2PID[pipsIdx],  piplus[pipsIdx].P(),
    //                                                                     Ve,
    //                                                                     Vpiplus[pipsIdx],
    //                                                                     fdebug);
    
    pipsPastSelectionCuts[pipsIdx] = CheckIfPionPassedSelectionCuts("pi+",
                                                                    piPfiducial1[pipsIdx],
                                                                    piPfiducial2[pipsIdx],
                                                                    piPfiducial3[pipsIdx],
                                                                    piPpid[pipsIdx],
                                                                    piplus[pipsIdx].P(),
                                                                    Ve,
                                                                    Vpiplus[pipsIdx]);
    
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
    
    //    piplus_Px[pipsIdx] = piplus[pipsIdx].Px();
    //    piplus_Py[pipsIdx] = piplus[pipsIdx].Py();
    //    piplus_Pz[pipsIdx] = piplus[pipsIdx].Pz();
    //    piplus_E[pipsIdx]  = piplus[pipsIdx].E();
    //    Vpiplus_X[pipsIdx] = Vpiplus[pipsIdx].X();
    //    Vpiplus_Y[pipsIdx] = Vpiplus[pipsIdx].Y();
    //    Vpiplus_Z[pipsIdx] = Vpiplus[pipsIdx].Z();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExtractPimsInformation( int pimsIdx ){
    // extract negative pion information
    piminus.at(pimsIdx).SetXYZM(piMpx[pimsIdx],
                                piMpy[pimsIdx],
                                piMpz[pimsIdx],
                                aux.Mpims);
    
    
    Zpims[pimsIdx]          = piminus[pimsIdx].E() / omega;
    Vpiminus[pimsIdx]       = TVector3( 0, 0 , piMvz[pimsIdx] );
    pims_DC_sector[pimsIdx] = piMsector[pimsIdx];
    
    if (fdebug>2) {
        std::cout
        << "p(π+ "          <<pimsIdx<<"): " << piminus[pimsIdx].P()     << ","
        << "z(π+ "          <<pimsIdx<<"): " << Zpims[pimsIdx]           << ","
        << "Vz(π+ "         <<pimsIdx<<"): " << Vpiminus[pimsIdx].Z()    << ","
        << "DC-sector(π+ "  <<pimsIdx<<"): " << pims_DC_sector[pimsIdx]  << ","
        << std::endl;
    }
    
    //    Vpiminus[pimsIdx]           = GetParticleVertex( piminuses[pimsIdx] );
    //    pims_chi2PID[pimsIdx]       = piminuses[pimsIdx]->par()->getChi2Pid();
    //
    //    // EC in and out
    //    pims_E_ECIN[pimsIdx]        = piminuses[pimsIdx]->cal(ECIN)->getEnergy();
    //    pims_E_ECOUT[pimsIdx]       = piminuses[pimsIdx]->cal(ECOUT)->getEnergy();
    //    // PCAL
    //    auto pims_PCAL_info         = piminuses[pimsIdx]->cal(PCAL);
    //    pims_E_PCAL[pimsIdx]        = pims_PCAL_info->getEnergy();
    //    pims_PCAL_sector[pimsIdx]   = pims_PCAL_info->getSector();
    //    pims_PCAL_V[pimsIdx]        = pims_PCAL_info->getLv();
    //    pims_PCAL_W[pimsIdx]        = pims_PCAL_info->getLw();
    //    pims_PCAL_x[pimsIdx]        = pims_PCAL_info->getX();
    //    pims_PCAL_y[pimsIdx]        = pims_PCAL_info->getY();
    //    pims_PCAL_z[pimsIdx]        = pims_PCAL_info->getZ();
    //    // DC
    //    auto pims_DC_info           = piminuses[pimsIdx]->trk(DC);
    //    pims_DC_sector[pimsIdx]     = pims_DC_info->getSector(); // tracking sector
    //    pims_Chi2N[pimsIdx]         = pims_DC_info->getChi2N();  // tracking chi^2/NDF
    //    for (int regionIdx=0; regionIdx<3; regionIdx++) {
    //        DC_layer = DC_layers[regionIdx];
    //        pims_DC_x[pimsIdx][regionIdx] = piminuses[pimsIdx]->traj(DC,DC_layer)->getX();
    //        pims_DC_y[pimsIdx][regionIdx] = piminuses[pimsIdx]->traj(DC,DC_layer)->getY();
    //        pims_DC_z[pimsIdx][regionIdx] = piminuses[pimsIdx]->traj(DC,DC_layer)->getZ();
    //    }
    // ------------------------------------------------------------------------------------------------
    // now, check if pion passed event selection requirements
    // ------------------------------------------------------------------------------------------------
    //    pimsPastSelectionCuts[pimsIdx] = CheckIfPionPassedSelectionCuts("pi-",
    //                                                                     pims_DC_sector[pimsIdx],
    //                                                                     pims_DC_x[pimsIdx],
    //                                                                     pims_DC_y[pimsIdx],
    //                                                                     pims_DC_z[pimsIdx],
    //                                                                     pims_chi2PID[pimsIdx],  piminus[pimsIdx].P(),
    //                                                                     Ve,
    //                                                                     Vpiminus[pimsIdx],
    //                                                                     fdebug);
    
    
    pimsPastSelectionCuts[pimsIdx] = CheckIfPionPassedSelectionCuts("pi-",
                                                                    piMfiducial1[pimsIdx],
                                                                    piMfiducial2[pimsIdx],
                                                                    piMfiducial3[pimsIdx],
                                                                    piMpid[pimsIdx],
                                                                    piminus[pimsIdx].P(),
                                                                    Ve,
                                                                    Vpiminus[pimsIdx]);
    
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
    
    //    piminus_Px[pimsIdx]          = piminus[pimsIdx].Px();
    //    piminus_Py[pimsIdx]          = piminus[pimsIdx].Py();
    //    piminus_Pz[pimsIdx]          = piminus[pimsIdx].Pz();
    //    piminus_E[pimsIdx]           = piminus[pimsIdx].E();
    //    Vpiminus_X[pimsIdx]          = Vpiminus[pimsIdx].X();
    //    Vpiminus_Y[pimsIdx]          = Vpiminus[pimsIdx].Y();
    //    Vpiminus_Z[pimsIdx]          = Vpiminus[pimsIdx].Z();
}




// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
bool CheckIfElectronPassedSelectionCuts(Int_t fid1,
                                        Int_t fid2,
                                        Int_t fid3,
                                        TLorentzVector e,
                                        TVector3 Ve){
    
    if (fdebug>2) {
        std::cout
        << "CheckIfElectronPassedSelectionCuts("
        << fid1 << "," << fid2 << "," << fid3 << ","
        << Ve.Z()
        << ")"
        << std::endl;
        
        if ((aux.cutValue_Vz_min < Ve.Z()) && (Ve.Z() < aux.cutValue_Vz_max))
            std::cout << "passed z-cut" << std::endl;
        
        if (fid1 && fid2 && fid3 == false)
            std::cout << "Did not pass fiducial cuts" << std::endl;
        
        if(!(true
             // Cut on z-vertex position: in-bending torus field -13.0 cm < Vz < +12.0 cm
             // Spring 19 and Spring 2020 in-bending.
             // Fall 2019 (without low-energy-run) was out-bending.
             &&  ((aux.cutValue_Vz_min < Ve.Z()) && (Ve.Z() < aux.cutValue_Vz_max))
             ))
            std::cout << "Did not pass electron cuts" << std::endl;
        else
            std::cout << "Succesfully passed electron cuts" << std::endl;
    }
    bool DC_fid = fid1 && fid2 && fid3;
    if (DC_fid == false) return false;
    
    if(!(true
         // Cut on z-vertex position: in-bending torus field -13.0 cm < Vz < +12.0 cm
         // Spring 19 and Spring 2020 in-bending.
         // Fall 2019 (without low-energy-run) was out-bending.
         &&  ((aux.cutValue_Vz_min < Ve.Z()) && (Ve.Z() < aux.cutValue_Vz_max))
         )) return false;
    
    return true;
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
bool CheckIfPionPassedSelectionCuts(TString pionCharge, // "pi+" or "pi-"
                                    Int_t fid1,
                                    Int_t fid2,
                                    Int_t fid3,
                                    Double_t chi2PID,
                                    Double_t p,
                                    TVector3 Ve,
                                    TVector3 Vpi){
    
    // decide if pion (pi+ or pi-) passed event selection cuts
    //
    // input:
    // --------
    // p            pi momentum    (pi.P())
    //
    // comments
    // ---------------
    // DC - fiducial cuts on DC
    
    if (fdebug>2) {
        std::cout
        << "CheckIfPionPassedSelectionCuts("
        << pionCharge   << ","
        << fid1         << ","
        << fid2         << ","
        << fid3         << ","
        << std::setprecision(4)
        << chi2PID      << ","
        << std::setprecision(3)
        << p            << ","
        << "|z(e)-z(π)| = "<< (Ve-Vpi).Z() << ")"
        << std::endl;
    }
    
    bool DC_fid = fid1 && fid2 && fid3;
    if (DC_fid == false) return false;
    
    
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
    
    
    
    
    // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    // Aug-25, 2022
    // There is a problem here with chi2PID
    // and the following block bypasses it temporarily
    // so that at this point there is effectively no cut
    // on the chi2pid of the pion.
    //
    // We need to solve it in the future.
    //
    // In RGB SIDIS studies we filter pions based on their chi2pid information:
    //
    // region_part_ptr::par()->getChi2Pid()
    //
    // This variable represents the goodness of fit
    // for the hypothesis that it is a pion.
    // In the RGA datafile that Igor prepared,
    // there is a variable called *piPpid* for π+,
    // and similarly *piMpid* for π-,
    // but it only stores the PDG code of the particles, 211 and -211 respectively.
    //
    // On Aug-25,2022
    // I asked Igor to add the chi2pid information to the RGA data-files
    // and when he does this the following block should be corrected
    //
    // The reason why it is written as a block with the awd if-statement
    // is to add a precaution in case I forget to make the correction:
    // the bypass is done only in case where
    // chi2pid == PDGcode for the particle
    //
    if (chi2PID == 211 || chi2PID == -211){
        chi2PID = (aux.Chi2PID_pion_lowerBound( p, C )
                   +
                   aux.Chi2PID_pion_upperBound( p , C ) )/2;
    }
    //
    //
    //
    // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    
    if (fdebug>3) {
        std::cout << "in CheckIfPionPassedSelectionCuts()"<< std::endl
        << "pion charge: "          << pionCharge               << ","
        << "chi2PID:"               << chi2PID                  << ","
        << std::endl
        << "Chi2PID_pion_lowerBound( p="<<p<<", C="<<C<<" ): "
        << aux.Chi2PID_pion_lowerBound( p, C ) << ","
        << std::endl
        << "Chi2PID_pion_upperBound( p="<<p<<", C="<<C<<" ): "
        << aux.Chi2PID_pion_upperBound( p, C ) << ","
        << std::endl
        << "|z(e)-z(π)| = "         << fabs((Ve-Vpi).Z())       << ","
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
        if (fdebug>3) {
            std::cout
            << "Did not pass CheckIfPionPassedSelectionCuts(), return false"
            << std::endl;
        }
        return false;
    }
    if (fdebug>3) {
        std::cout
        << "succesfully passed CheckIfPionPassedSelectionCuts(), return true"
        << std::endl;
    }
    return true;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Stream_e_pi_line_to_CSV( TString pionCharge, int piIdx,
                             bool passed_cuts_e_pi,
                             bool passed_cuts_e_pi_kinematics){
    TLorentzVector  pi;
    TLorentzVector  pi_qFrame;
    TVector3        Vpi;
    int    pi_DC_sector;
    
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
        W_d,
    };
    
    // decide which file to write...
    if (pionCharge=="pi+") {
        if (passed_cuts_e_pi && passed_cuts_e_pi_kinematics) {
            Nevents_e_e_pips ++ ;
            aux.StreamToCSVfile(CSVfile_e_e_pips,
                                variables,
                                csvprecisions );
        }
    }
    else if (pionCharge=="pi-") {
        if (passed_cuts_e_pi && passed_cuts_e_pi_kinematics) {
            Nevents_e_e_pims ++ ;
            aux.StreamToCSVfile(CSVfile_e_e_pims,
                                variables,
                                csvprecisions );
        }
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
    xF      = 2. * (pi.Vect().Dot(q.Vect())) / (q.P() * W);
    eta_pi  = 0.5 * log((pi_qFrame.E()+pi_qFrame.Pz()) /
                        (pi_qFrame.E()-pi_qFrame.Pz()));
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
void PrintRGAInformation(){
for (int e_idx=0; e_idx<nml; e_idx++) {
    std::cout
    << "Ep["<<e_idx<<"] "               << Ep[e_idx]            << ","
    << std::endl
    << "Epx["<<e_idx<<"] "              << Epx[e_idx]           << ","
    << "Epy["<<e_idx<<"] "              << Epy[e_idx]           << ","
    << "Epz["<<e_idx<<"] "              << Epz[e_idx]           << ","
    << std::endl
    << "Epid["<<e_idx<<"] "             << Epid[e_idx]          << ","
    << "Esector["<<e_idx<<"] "          << Esector[e_idx]       << ","
    << std::endl
    << "Etheta["<<e_idx<<"] "           << Etheta[e_idx]        << ","
    << "Ephi["<<e_idx<<"] "             << Ephi[e_idx]          << ","
    << std::endl
    << "Efiducial1["<<e_idx<<"] "       << Efiducial1[e_idx]    << ","
    << "Efiducial2["<<e_idx<<"] "       << Efiducial2[e_idx]    << ","
    << "Efiducial3["<<e_idx<<"] "       << Efiducial3[e_idx]    << ","
    << std::endl
    << "ePCAL_energy["<<e_idx<<"] "     << ePCAL_energy[e_idx]  << ","
    << "eECIN_energy["<<e_idx<<"] "     << eECIN_energy[e_idx]  << ","
    << "eECOUT_energy["<<e_idx<<"] "    << eECOUT_energy[e_idx] << ","
    << std::endl
    << "Nu["<<e_idx<<"] "               << Nu[e_idx]            << ","
    << std::endl
    << "qx["<<e_idx<<"] "              << qx[e_idx]           << ","
    << "qy["<<e_idx<<"] "              << qy[e_idx]           << ","
    << "qz["<<e_idx<<"] "              << qz[e_idx]           << ","
    << std::endl;
}
for (int pips_idx=0; pips_idx < npiP; pips_idx++) {
    std::cout
    << "piPp["<<pips_idx<<"] "               << piPp[pips_idx]            << ","
    << std::endl
    << "piPpx["<<pips_idx<<"] "              << piPpx[pips_idx]           << ","
    << "piPpy["<<pips_idx<<"] "              << piPpy[pips_idx]           << ","
    << "piPpz["<<pips_idx<<"] "              << piPpz[pips_idx]           << ","
    << std::endl
    << "piPpid["<<pips_idx<<"] "             << piPpid[pips_idx]          << ","
    << "piPsector["<<pips_idx<<"] "          << piPsector[pips_idx]       << ","
    << std::endl
    << "piPtheta["<<pips_idx<<"] "           << piPtheta[pips_idx]        << ","
    << "piPphi["<<pips_idx<<"] "             << piPphi[pips_idx]          << ","
    << std::endl
    << "piPfiducial1["<<pips_idx<<"] "       << piPfiducial1[pips_idx]    << ","
    << "piPfiducial2["<<pips_idx<<"] "       << piPfiducial2[pips_idx]    << ","
    << "piPfiducial3["<<pips_idx<<"] "       << piPfiducial3[pips_idx]    << ","
    << std::endl
    << "piPPCAL_energy["<<pips_idx<<"] "     << piPPCAL_energy[pips_idx]  << ","
    << "piPECIN_energy["<<pips_idx<<"] "     << piPECIN_energy[pips_idx]  << ","
    << "piPECOUT_energy["<<pips_idx<<"] "    << piPECOUT_energy[pips_idx] << ","
    << std::endl;
}
for (int pims_idx=0; pims_idx < npiM; pims_idx++) {
    std::cout
    << "piMp["<<pims_idx<<"] "               << piMp[pims_idx]            << ","
    << std::endl
    << "piMpx["<<pims_idx<<"] "              << piMpx[pims_idx]           << ","
    << "piMpy["<<pims_idx<<"] "              << piMpy[pims_idx]           << ","
    << "piMpz["<<pims_idx<<"] "              << piMpz[pims_idx]           << ","
    << std::endl
    << "piMpid["<<pims_idx<<"] "             << piMpid[pims_idx]          << ","
    << "piMsector["<<pims_idx<<"] "          << piMsector[pims_idx]       << ","
    << std::endl
    << "piMtheta["<<pims_idx<<"] "           << piMtheta[pims_idx]        << ","
    << "piMphi["<<pims_idx<<"] "             << piMphi[pims_idx]          << ","
    << std::endl
    << "piMfiducial1["<<pims_idx<<"] "       << piMfiducial1[pims_idx]    << ","
    << "piMfiducial2["<<pims_idx<<"] "       << piMfiducial2[pims_idx]    << ","
    << "piMfiducial3["<<pims_idx<<"] "       << piMfiducial3[pims_idx]    << ","
    << std::endl
    << "piMPCAL_energy["<<pims_idx<<"] "     << piMPCAL_energy[pims_idx]  << ","
    << "piMECIN_energy["<<pims_idx<<"] "     << piMECIN_energy[pims_idx]  << ","
    << "piMECOUT_energy["<<pims_idx<<"] "    << piMECOUT_energy[pims_idx] << ","
    << std::endl;
}
}
