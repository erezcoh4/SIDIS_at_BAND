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
//#define NMAXPIONS 5 // maximal allowed number of pions
//#define r2d 180./3.1415 // radians to degrees
SIDISatBAND_auxiliary aux;


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
// Globals
// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
TString csvheader = ((TString)"status,runnum,evnum,beam_helicity,"
                     +(TString)"e_P,e_Theta,e_Phi,e_Vz,"
                     +(TString)"pi_P,pi_Theta,pi_Phi,pi_Vz,"
                     +(TString)"n_P,n_Theta,n_Phi,n_Vz,"
                     +(TString)"n_HitPos_X,n_HitPos_Y,n_HitPos_Z,"
                     +(TString)"n_E,alpha_n,"
                     +(TString)"Q2,xB,omega,y,"
                     +(TString)"e_DC_sector,pi_DC_sector,"
                     +(TString)"pi_qFrame_Theta,pi_qFrame_Phi,"
                     +(TString)"pi_qFrame_pT,pi_qFrame_pL,"
                     +(TString)"n_Theta_qFrame,n_Phi_qFrame,"
                     +(TString)"n_pT_qFrame,n_pL_qFrame,"
                     +(TString)"E_init,P_init,alpha_init,v_init,"
                     +(TString)"Zpi,Zpi_LC,"
                     +(TString)"Zpi_LC_Prime,"
                     +(TString)"W,M_x,"
                     +(TString)"xF,eta_pi,"
                     +(TString)"W_Prime,M_x_Prime,"
                     +(TString)"xF_Prime,"
                     +(TString)"xB_Prime"
                     );
std::vector<int> csvprecisions = {0,0,0,0,
    9,9,9,9,
    9,9,9,9,
    9,9,9,9,
    9,9,9,
    9,9,
    9,9,9,9,
    0,0,
};



TString DataPath = "/volatile/clas12/users/ecohen/BAND/";
TString                            prefix = "sidisdvcs_";
TString                              skimmedBANDFilename;
TString                             skimmedSIDISFilename;
TString                                       pionCharge; // "pi+" or "pi-"
TString                                          pionStr;

// Input root files and trees
TFile * SIDISFile, * BANDFile;
TTree * SIDISTree, * BANDTree;

// Output root file and tree
TFile * MergedFile;
TTree * MergedTree;

// Output CSV file
std::ofstream   CSVfile_e_pi_n;

// time
clock_t tStart = clock();

// globals
//Double_t       Me  = 0.00051099895; // GeV/c2
//Double_t            Mpi = 0.139570; // GeV/c2
//Double_t       Mpims  = 0.13957039; // GeV/c2
//Double_t       Mpips  = 0.13957039; // GeV/c2
//Double_t            Mp  = 0.938272; // GeV/c2
//Double_t             Mn = 0.939565; // GeV/c2
//Double_t                Md = 1.875; // GeV/c2
//Double_t             Mp2 = Mp * Mp;
Double_t                        xB; // Bjorken x
Double_t                         y; // QE y-variable
Double_t                        Q2;
Double_t                     omega;
Double_t                        w2; // omega^2
Double_t          Zpips[NMAXPIONS]; // energy fraction rest frame
Double_t       Zpips_LC[NMAXPIONS]; // energy fraction rest frame on the light cone
Double_t          Zpims[NMAXPIONS]; // energy fraction rest frame
Double_t       Zpims_LC[NMAXPIONS]; // energy fraction rest frame on the light cone
Double_t                     W, W2; // energy of the hadronic system
//Double_t              W_p_rest;
//Double_t              W_d_rest;
//Double_t                     Emiss;
Double_t                    E_init; // proton initial state energy 
Double_t                   alpha_n; // light cone fraction of momentum of the recoil neutron
Double_t                alpha_init; // light cone fraction of momentum of the initial proton
Double_t                    v_init; // virtuality of the initial proton
Double_t          W_Prime, W2prime; // moving proton
Double_t                        xF; // Feynman x for a standing proton
Double_t                  xF_Prime; // Feynman x for a moving proton
Double_t                    eta_pi; // rapidity


Double_t                   xB_Prime; // moving proton: x' = Q2 / (W'^2 - mN^2 + Q2)
//Double_t                   xB_Prime2; // moving proton defined as
// x' = Q2 / (2. * ((Md - Es) * omega + Pn_Vect*q->Vect() ));

Double_t                        En; // (spectator) neutron energy
Double_t                        Es; // spectator energy
Double_t                        Ps; // spectator momentum
Double_t                  theta_sq; // spectator angle with respect to momentum transfer
Double_t            M_x, M_x_Prime;
Double_t                     n_ToF;
Double_t Zpi, Zpi_LC, Zpi_LC_Prime;



Int_t                               NeventsToMerge;
Int_t                       BANDrunID, BANDeventID;
Int_t                     SIDISrunID, SIDISeventID;
Int_t                  EventIDsToMerge[NMAXEVENTS];
Int_t                         Npions, Npips, Npims;
Int_t                          Ne, Np, Nn, Ngammas;
Int_t                                       status;
std::vector<Int_t>                    BANDeventIDs;
std::vector<Int_t>                   SIDISeventIDs;
std::vector<Int_t>                BANDeventIndices;
std::vector<Int_t>               SIDISeventIndices;
std::vector<Int_t>             EventNumbersToMerge;
std::vector<std::size_t>   BANDEventIndicesToMerge;
std::vector<std::size_t>  SIDISEventIndicesToMerge;

// SIDIS Tree
//TLorentzVector     *target = new TLorentzVector(0, 0, 0, aux.Md );
TLorentzVector     *Beam=0;
TLorentzVector        *e=0;
TLorentzVector        *q=0;
TLorentzVector          Pn; // neutron momentum
TLorentzVector Band_data_e; // electron information in BAND TTree
TLorentzVector      p_init;
TLorentzVector d_rest(0, 0, 0, aux.Md );
TLorentzVector p_rest(0, 0, 0, aux.Mp );

// "q-frame" parameters
double              Pe_phi;
double      q_phi, q_theta;
// vectors in q-frame
TLorentzVector       e_qFrame;
TLorentzVector       q_qFrame;
TLorentzVector      pi_qFrame;
TLorentzVector      Pn_qFrame;
TLorentzVector  p_init_qFrame;

// reconstructed vertex position
TVector3             *Ve=0;
TVector3           Pn_Vect; // neutron 3-momentum
TVector3       Band_e_Vect; // electron 3-momentum as written in BAND file
TVector3          n_HitPos;

std::vector<TVector3>               Vpiplus;
std::vector<TVector3>              Vpiminus;
std::vector<TLorentzVector>          piplus;
std::vector<TLorentzVector>         piminus;
std::vector<TLorentzVector>   piplus_qFrame;
std::vector<TLorentzVector>  piminus_qFrame;

bool                     eepiPastCutsInEvent;
bool                   eepipsPastCutsInEvent;
bool                   eepimsPastCutsInEvent;
bool                     goodneutron = false;
bool    eepipsPastKinematicalCuts[NMAXPIONS];
bool        pimsPastSelectionCuts[NMAXPIONS];
bool    eepimsPastKinematicalCuts[NMAXPIONS];

int             fdebug = 1;
int        eventnumber = 0;
int        nleadindex = -1;
int              nMult = 0;
int            genMult = 0;
int          beam_helicity;

double           Ebeam = 0;
double    gated_charge = 0;
double        livetime = 0;
double       starttime = 0;
double         current = 0;
double          weight = 0;
                                
double                    e_DC_sector;

double           piplus_Px[NMAXPIONS];
double           piplus_Py[NMAXPIONS];
double           piplus_Pz[NMAXPIONS];
double            piplus_E[NMAXPIONS];
double           Vpiplus_X[NMAXPIONS];
double           Vpiplus_Y[NMAXPIONS];
double           Vpiplus_Z[NMAXPIONS];
double      pips_DC_sector[NMAXPIONS];
double    piplus_qFrame_pT[NMAXPIONS]; // transverse momentum with respect to q
double    piplus_qFrame_pL[NMAXPIONS]; // longitudinal momentum with respect to q
double piplus_qFrame_Theta[NMAXPIONS];
double   piplus_qFrame_Phi[NMAXPIONS];

double           piminus_Px[NMAXPIONS];
double           piminus_Py[NMAXPIONS];
double           piminus_Pz[NMAXPIONS];
double            piminus_E[NMAXPIONS];
double           Vpiminus_X[NMAXPIONS];
double           Vpiminus_Y[NMAXPIONS];
double           Vpiminus_Z[NMAXPIONS];
double       pims_DC_sector[NMAXPIONS];
double    piminus_qFrame_pT[NMAXPIONS]; // transverse momentum with respect to q
double    piminus_qFrame_pL[NMAXPIONS]; // longitudinal momentum with respect to q
double piminus_qFrame_Theta[NMAXPIONS];
double   piminus_qFrame_Phi[NMAXPIONS];

clashit *eHit             = new clashit;
TClonesArray      * nHits = new TClonesArray("bandhit"); // BAND neutrons in BAND analysis
TClonesArray    * mcParts = new TClonesArray("genpart");


void                      OpenInputFiles (TString RunStr);
void                     OpenOutputFiles (TString RunStr);
//void                     StreamToCSVfile (std::vector<Double_t> observables);
void                     CloseInputFiles ();
void                    CloseOutputFiles (TString OutDataPath);
void            MergeSIDISandBANDevents  (int NMAXeventsToMerge=10,
                                          int PrintProgress=5000);
Int_t          CreateListOfEventsToMerge (TTree * BANDTree,
                                          TTree * SIDISTree,
                                          int NMAXeventsToMerge=-1);
void           Stream_e_pi_n_line_to_CSV (int piIdx,
                                          bool passed_cuts_e_pi_kinematics,
                                          bool passed_cuts_n);

void                 InitializeVariables ();
void                        GetSIDISData ( int MergedEvtId );
void                         GetBANDData ( int MergedEvtId );
void                      MergeEventData ();
void             SetInputAndOutputTTrees ();
//void           ComputeElectronKinematics ();
//void               ComputePionKinematics (TLorentzVector pi,
//                                          TLorentzVector pi_qFrame);
void                           PrintDone ();
//void                    PrintMonitorHello ();
void                       SetPionCharge ( TString fpionCharge ) {
    pionCharge = fpionCharge;
    if (pionCharge=="pi+") {
        pionStr = "_e_piplus";
    } else if (pionCharge=="pi-") {
        pionStr = "_e_piminus";
    } else {
        pionStr = "_no_pion_charge_info";
    }
};
void                         SetDataPath ( TString fDataPath )   {DataPath = fDataPath + "/";};
void                           SetPrefix ( TString fPrefix )     {prefix = fPrefix;};
void                        SetVerbosity ( int ffdebug )         {fdebug = ffdebug;};
void                           PrintTime ( TString prefix ){
    std::cout << prefix << ", after "
    << double(clock() - tStart) / (double)CLOCKS_PER_SEC
    << " sec "<< std::endl;
}
void                       MoveTo_qFrame ();
TVector3           RotateVectorTo_qFrame ( TVector3 v );
void                   ComputeKinematics (TLorentzVector pi);
//void               Print4Vector ( TLorentzVector v, std::string label );


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
// Main functionality
// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void MergeSIDISandBANDSkimmers(int RunNumber=6420,
                               TString fpionCharge="pi+", // "pi+" or "pi-"
                               int NMAXeventsToMerge=-1,
                               int ffdebug=1,
                               int PrintProgress=1000,
                               TString fDataPath="/volatile/clas12/users/ecohen/BAND/",
                               TString fPrefix="sidisdvcs_"){
    
    SetPionCharge    ( fpionCharge );
    SetVerbosity     ( ffdebug );
    SetDataPath      ( fDataPath );
    SetPrefix        ( fPrefix );
    char RunNumberStr[20];
    sprintf( RunNumberStr, "00%d", RunNumber );
    OpenInputFiles   ( (TString)RunNumberStr );
    OpenOutputFiles  ( (TString)RunNumberStr );
    MergeSIDISandBANDevents( NMAXeventsToMerge, PrintProgress );
    CloseOutputFiles (DataPath + "merged_SIDIS_and_BAND_skimming/");
    CloseInputFiles  ();
    PrintDone        ();
}



// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void MergeSIDISandBANDevents (int NMAXeventsToMerge=10,
                              int PrintProgress=1000){
    
    Int_t   NeventsBAND  = BANDTree->GetEntries();
    Int_t   NeventsSIDIS = SIDISTree->GetEntries();
    // Create a list of events to merge
    NeventsToMerge = CreateListOfEventsToMerge(BANDTree, SIDISTree, NMAXeventsToMerge);
    if ((NMAXeventsToMerge>0) && (NeventsToMerge>NMAXeventsToMerge)) NeventsToMerge = NMAXeventsToMerge;
        
    if (fdebug>2){
        PrintTime((TString)"Created list of (" + (TString)std::to_string(NeventsToMerge) + (TString)") events to merge");
        std::cout << "................................................................" << std::endl;
    }
    
    // assign TTree branches to variables
    SetInputAndOutputTTrees ();
        
    // step over list of events-to-merge and merge them...
    for (int MergedEvtId=0; MergedEvtId < NeventsToMerge; MergedEvtId++) {

        // initialize
        InitializeVariables ();

        // grab electron and pion information from SIDIS TTree
        GetSIDISData( MergedEvtId );
        
        // grab neturon information from BAND
        GetBANDData ( MergedEvtId );
                
        // move to q-frame and define the pion and neutron momentum with respect to q
        MoveTo_qFrame       ();

        // merge information about the event from both sources
        MergeEventData      ();
        
        if (fdebug && MergedEvtId%PrintProgress==0) std::cout << MergedEvtId << "/" << NeventsToMerge << std::endl;
    } // end merged event loop
        
    if (fdebug>2) std::cout << "Merged " << NeventsToMerge << " SIDIS and BAND events." << std::endl;
    
}



// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void OpenInputFiles (TString RunStr){
    
    // Skimmed neutron-BAND files and "Final skim files" from
    // /volatile/clas12/users/segarrae/BAND/v3.1/10.2/final/tagged
    skimmedBANDFilename = (DataPath + "neutron_skimming/"
                           + "final_tagged_5MeV_250mevc_fiducial_thetaCut_allBars_"
                           + RunStr + ".root");
    
    // previous:
    //  + "ncalibration_shiftedskim_"  + RunStr + ".root"); // Efrain' file
    //    + "ncalibrationtest_0option_"  + RunStr + ".root"); // Florian' file
    //    + "skimmed_neutrons_inc_"  + RunStr + ".root"); // My file (global time shifts not working...)
    // Sep-21, "ncalibration_newclass" skimmer Tree name is "calib"
    //    BANDTree                      = (TTree*)BANDFile->Get("calib");
    if (fdebug>2) std::cout << "Opening " << skimmedBANDFilename << std::endl;
    BANDFile                      = new TFile( skimmedBANDFilename );
    // Nov-23, "final_tagged" skimmer Tree name is "tagged"
    BANDTree                      = (TTree*)BANDFile->Get("tagged");
    
    skimmedSIDISFilename = (DataPath + "SIDIS_skimming/"
                            + "skimmed_SIDIS_" + prefix + RunStr + pionStr + ".root");
    if (fdebug>2) std::cout << "Opening " << skimmedSIDISFilename << std::endl;
    SIDISFile                     = new TFile( skimmedSIDISFilename );
    SIDISTree                     = (TTree*)SIDISFile->Get("tree");
}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void OpenOutputFiles (TString RunStr){
    
    TString skimmedMergedFilename = (DataPath + "merged_SIDIS_and_BAND_skimming/"
                                     + "skimmed_SIDIS_and_BAND_" + prefix + RunStr + pionStr + "_n" );
    
    if (fdebug>2) {
        std::cout
        << "Opening output file: " << skimmedMergedFilename
        << ".root/csv " << std::endl;
    }
    
    // Create output tree
    MergedFile = new TFile( skimmedMergedFilename + ".root" ,"RECREATE");
    MergedTree = new TTree( "T" , "Event information from merged SIDIS and BAND skimmers");
    
    // Create output csv files
    CSVfile_e_pi_n.open( skimmedMergedFilename + ".csv" );
    CSVfile_e_pi_n << csvheader << std::endl;
    
}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void CloseInputFiles (){
    SIDISFile->Close();
    BANDFile->Close();
}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void CloseOutputFiles (TString OutDataPath){
    
    // close output CSV
    CSVfile_e_pi_n.close();
    int Nevents = MergedTree->GetEntries();
    
    // close output ROOT
    MergedFile->cd();
    MergedTree->Write();
    MergedFile->Close();
    
    std::cout << "Done producing output files. They are ready in " << std::endl << OutDataPath << std::endl
    << "merged " << Nevents << " SIDIS (e,e'" << pionCharge << ") and BAND (neutron) events" << std::endl;
}


//// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
//void StreamToCSVfile (std::vector<Double_t> observables){
//    
//    
////    for (auto v:observables) {
////        CSVfile_e_pi_n << std::setprecision(9) << std::fixed << v << ",";
////    }
////    CSVfile_e_pi_n << std::endl;
//    
//    for (int j=0; j < observables.size(); j++){
//        int precision = 9;
//        if (j < csvprecisions.size()) precision = csvprecisions.at(j);
//        auto v = observables.at(j);
//        CSVfile_e_pi_n << std::setprecision(precision) << std::fixed << v << ",";
//    }
//    CSVfile_e_pi_n << std::endl;
//    
//    
//    if (fdebug>1) {
//        std::cout << std::fixed;
//        std::cout << "StreamToCSVfile()" << std::endl;
//        std::cout << csvheader << std::endl;
//        for (auto v:observables) std::cout << std::fixed << v << ",";
//        std::cout << std::endl;
//    }
//}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
Int_t CreateListOfEventsToMerge(TTree * BANDTree,
                                TTree * SIDISTree,
                                int NMAXeventsToMerge){
    
    // Decide which event-indices to merge from the two TTrees
    // Since this functionallity takes the most time,
    // we need to be very creative to make it less time-consuming
    if (fdebug>2) {
        std::cout << "CreateListOfEventsToMerge()" << std::endl;
    }
    
    Int_t                    NmergedEvents = 0;
    Int_t                      NeventsBAND = 0;
    Int_t                     NeventsSIDIS = 0;
    bool              DefineEventIDvectors = true;
    bool                 ListEventsToMerge = true;
    BANDeventIDs                .clear();
    BANDeventIndices            .clear();
    BANDEventIndicesToMerge     .clear();
    SIDISeventIDs               .clear();
    SIDISeventIndices           .clear();
    SIDISEventIndicesToMerge    .clear();
    EventNumbersToMerge         .clear();
    
    
    // First, define two vectors that contain the event IDs in each TTree
    if (DefineEventIDvectors){
        
        //        TTreeReader BANDReader("calib", BANDFile); // for ncalibration skimming
        TTreeReader BANDReader("tagged", BANDFile); // for tagged skimming
        TTreeReaderValue<Int_t>     fBANDeventID(BANDReader, "eventnumber");
        TTreeReaderValue<bool>  fBANDgoodneutron(BANDReader, "goodneutron");
        int BANDeventIndex = 0;
        while (BANDReader.Next()) {
            if (*fBANDgoodneutron==true){
                BANDeventIDs.push_back(*fBANDeventID);
                BANDeventIndices.push_back( NeventsBAND );
                NeventsBAND++;
            }
            BANDeventIndex++;
        }
        
        TTreeReader SIDISReader("tree", SIDISFile);
        TTreeReaderValue<Int_t>   fSIDISeventID(SIDISReader, "eventnumber");
        int SIDISeventIndex = 0;
        if (pionCharge=="pi+"){
            TTreeReaderValue<bool>  fSIDISeepipsPastCutsInEvent(SIDISReader, "eepipsPastCutsInEvent");
            while (SIDISReader.Next()) {
                if (*fSIDISeepipsPastCutsInEvent==true){
                    SIDISeventIDs       .push_back(*fSIDISeventID);
                    SIDISeventIndices   .push_back( SIDISeventIndex );
                    NeventsSIDIS++;
                }
                SIDISeventIndex++;
            }
        }
        else if (pionCharge=="pi-"){ // same for pi-
            TTreeReaderValue<bool>  fSIDISeepimsPastCutsInEvent(SIDISReader, "eepimsPastCutsInEvent");
            while (SIDISReader.Next()) {
                if (*fSIDISeepimsPastCutsInEvent==true){
                    SIDISeventIDs       .push_back(*fSIDISeventID);
                    SIDISeventIndices   .push_back( SIDISeventIndex );
                    NeventsSIDIS++;
                }
                SIDISeventIndex++;
            }
        }

        // Print-outs
        if (fdebug>3){
            PrintTime ( "done filling BANDeventIDs and SIDISeventIDs" );
            if (fdebug>4){
                std::cout << "BANDeventIDs: " << std::endl << std::endl;
                for (int BANDevent=0; BANDevent<50; BANDevent++){
                    std::cout << BANDeventIDs[BANDevent] << ",";
                }
                std::cout << "..." << std::endl << std::endl;
                std::cout << "SIDISeventIDs: "  << std::endl << std::endl;
                for (int SIDISevent=0; SIDISevent < std::min(50,NeventsSIDIS); SIDISevent++) {
                    std::cout << SIDISeventIDs[SIDISevent] << ",";
                }
                if (NeventsSIDIS > 50) {
                    std::cout << "..." << std::endl << std::endl;
                } else {
                    std::cout << std::endl;
                }
            }
        }
    }
    
    // second, list the events that we want to merge
    if (ListEventsToMerge){
        
        std::vector<int> tmpBANDeventIDs(BANDeventIDs);
        std::sort(std::begin(tmpBANDeventIDs), std::end(tmpBANDeventIDs));
        BANDEventIndicesToMerge.reserve(BANDeventIDs.size());
        
        std::vector<int> tmpSIDISeventIDs(SIDISeventIDs);
        std::sort(std::begin(tmpSIDISeventIDs), std::end(tmpSIDISeventIDs));
        SIDISEventIndicesToMerge.reserve(SIDISeventIDs.size());

        for (size_t i = 0; i < SIDISeventIDs.size(); ++i) {
            if ( (NMAXeventsToMerge>0) && (SIDISEventIndicesToMerge.size()>NMAXeventsToMerge-1) ) {
                break;
            }
            if (std::binary_search(std::begin(tmpBANDeventIDs),
                                   std::end(tmpBANDeventIDs),
                                   SIDISeventIDs[i])) {
                SIDISEventIndicesToMerge.push_back( SIDISeventIndices.at(i) );
            }
        }
        for (size_t i = 0; i < BANDeventIDs.size(); ++i) {
            if ( (NMAXeventsToMerge>0) && (BANDEventIndicesToMerge.size()>NMAXeventsToMerge-1) ) {
                break;
            }
            if (std::binary_search(std::begin(tmpSIDISeventIDs),
                                   std::end(tmpSIDISeventIDs),
                                   BANDeventIDs[i])) {
                BANDEventIndicesToMerge.push_back( BANDeventIndices.at(i));
                EventNumbersToMerge.push_back( BANDeventIDs[i] );
            }
        }
        NmergedEvents = EventNumbersToMerge.size();
        
        // print-outs
        if (fdebug>4){
            std::cout << NmergedEvents << " events to merge: ";
            for(int i=0; i < NmergedEvents ; i++ ) std::cout << BANDeventIDs.at(BANDEventIndicesToMerge.at(i)) << ' ';
            std::cout << std::endl<< std::endl;

            std::cout << "indices in BANDEventIndicesToMerge: ";
            for(std::vector<int>::size_type i : BANDEventIndicesToMerge)
                std::cout << i << ' ';

            std::cout << std::endl << std::endl;

            std::cout << "indices in SIDISEventIndicesToMerge: ";
            for(std::vector<int>::size_type i : SIDISEventIndicesToMerge)
                std::cout << i << ' ';
            std::cout << std::endl << std::endl;

            PrintTime("Done creating a list of " + std::to_string(NmergedEvents) + " events to merge");
        }
    }
    return NmergedEvents;
}

//
//// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
//Double_t ComputeLightConeFraction( TLorentzVector p ){
//    // compute light-cone momentum fraction
//    Double_t m = p.Mag();
//    Double_t alpha = (p.E() - p.Z())/m;
//    return alpha;
//}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void SetInputAndOutputTTrees (){
    if (fdebug>3) {
        std::cout << "SetInputAndOutputTTrees()" << std::endl;
    }
    
    // SIDIS input tree branches
    SIDISTree  -> SetBranchAddress("eventnumber"            ,&SIDISeventID          );
    SIDISTree  -> SetBranchAddress("runnum"                 ,&SIDISrunID            );
    SIDISTree  -> SetBranchAddress("e"                      ,&e                     );
    SIDISTree  -> SetBranchAddress("Ve"                     ,&Ve                    );
    SIDISTree  -> SetBranchAddress("Beam"                   ,&Beam                  );
    SIDISTree  -> SetBranchAddress("beam_helicity"          ,&beam_helicity         );
    SIDISTree  -> SetBranchAddress("q"                      ,&q                     );
    SIDISTree  -> SetBranchAddress("Nelectrons"             ,&Ne                     );
    SIDISTree  -> SetBranchAddress("Ngammas"                ,&Ngammas                );
    SIDISTree  -> SetBranchAddress("Nprotons"               ,&Np                     );
    SIDISTree  -> SetBranchAddress("Nneutrons"              ,&Nn                     );
    SIDISTree  -> SetBranchAddress("Npips"                  ,&Npips                  );
    SIDISTree  -> SetBranchAddress("Npims"                  ,&Npims                  );
    SIDISTree  -> SetBranchAddress("e_DC_sector"            ,&e_DC_sector            );
    SIDISTree  -> SetBranchAddress("omega"                  ,&omega                  );
    SIDISTree  -> SetBranchAddress("xB"                     ,&xB                     );
    SIDISTree  -> SetBranchAddress("Q2"                     ,&Q2                     );
    SIDISTree  -> SetBranchAddress("y"                      ,&y                      );
    SIDISTree  -> SetBranchAddress("W"                         ,&W                      );
    
    // branches that depend on pion charge
    if (pionCharge=="pi+") {
        SIDISTree  -> SetBranchAddress("eepipsPastCutsInEvent"
                                                            ,&eepipsPastCutsInEvent    );
        SIDISTree  -> SetBranchAddress("eepipsPastKinematicalCuts"
                                                            ,&eepipsPastKinematicalCuts);
        SIDISTree -> SetBranchAddress("piplus_Px"           ,&piplus_Px                );
        SIDISTree -> SetBranchAddress("piplus_Py"           ,&piplus_Py                );
        SIDISTree -> SetBranchAddress("piplus_Pz"           ,&piplus_Pz                );
        SIDISTree -> SetBranchAddress("piplus_E"            ,&piplus_E                 );
        SIDISTree -> SetBranchAddress("Vpiplus_X"           ,&Vpiplus_X                );
        SIDISTree -> SetBranchAddress("Vpiplus_Y"           ,&Vpiplus_Y                );
        SIDISTree -> SetBranchAddress("Vpiplus_Z"           ,&Vpiplus_Z                );
        SIDISTree -> SetBranchAddress("pi_DC_sector"        ,&pips_DC_sector           );
        SIDISTree -> SetBranchAddress("piplus_qFrame_pT"    ,&piplus_qFrame_pT         );
        SIDISTree -> SetBranchAddress("piplus_qFrame_pL"    ,&piplus_qFrame_pL         );
        SIDISTree -> SetBranchAddress("piplus_qFrame_Theta" ,&piplus_qFrame_Theta      );
        SIDISTree -> SetBranchAddress("piplus_qFrame_Phi"   ,&piplus_qFrame_Phi        );
        SIDISTree -> SetBranchAddress("Z"                   ,&Zpips        );
        SIDISTree -> SetBranchAddress("Z_LC"                ,&Zpips_LC        );
        
    } else if (pionCharge=="pi-") {
        SIDISTree  -> SetBranchAddress("eepimsPastCutsInEvent"
                                                            ,&eepimsPastCutsInEvent     );
        SIDISTree  -> SetBranchAddress("eepimsPastKinematicalCuts"
                                                            ,&eepimsPastKinematicalCuts );
        SIDISTree -> SetBranchAddress("piminus_Px"          ,&piminus_Px                );
        SIDISTree -> SetBranchAddress("piminus_Py"          ,&piminus_Py                );
        SIDISTree -> SetBranchAddress("piminus_Pz"          ,&piminus_Pz                );
        SIDISTree -> SetBranchAddress("piminus_E"           ,&piminus_E                 );
        SIDISTree -> SetBranchAddress("Vpiminus_X"          ,&Vpiminus_X                );
        SIDISTree -> SetBranchAddress("Vpiminus_Y"          ,&Vpiminus_Y                );
        SIDISTree -> SetBranchAddress("Vpiminus_Z"          ,&Vpiminus_Z                );
        SIDISTree -> SetBranchAddress("pi_DC_sector"        ,&pims_DC_sector            );
        SIDISTree -> SetBranchAddress("piminus_qFrame_pT"   ,&piminus_qFrame_pT         );
        SIDISTree -> SetBranchAddress("piminus_qFrame_pL"   ,&piminus_qFrame_pL         );
        SIDISTree -> SetBranchAddress("piminus_qFrame_Theta",&piminus_qFrame_Theta      );
        SIDISTree -> SetBranchAddress("piminus_qFrame_Phi"  ,&piminus_qFrame_Phi        );
        SIDISTree -> SetBranchAddress("Z"                   ,&Zpims        );
        SIDISTree -> SetBranchAddress("Z_LC"                ,&Zpims_LC        );
    }
    
    
    // BAND input tree branches
    BANDTree   -> SetBranchAddress("eventnumber"  ,&BANDeventID);
    BANDTree   -> SetBranchAddress("Runno"        ,&BANDrunID);
    BANDTree   -> SetBranchAddress("Ebeam"        ,&Ebeam);
    BANDTree   -> SetBranchAddress("gated_charge" ,&gated_charge);
    BANDTree   -> SetBranchAddress("livetime"     ,&livetime);
    BANDTree   -> SetBranchAddress("starttime"    ,&starttime);
    BANDTree   -> SetBranchAddress("current"      ,&current);
    BANDTree   -> SetBranchAddress("nMult"        ,&nMult);
    BANDTree   -> SetBranchAddress("nHits"        ,&nHits); // BAND neutrons in BAND skimming
    BANDTree   -> SetBranchAddress("eHit"         ,&eHit); // CLAS12 electrons in BAND skimming
    BANDTree   -> SetBranchAddress("goodneutron"  ,&goodneutron);
    BANDTree   -> SetBranchAddress("nleadindex"   ,&nleadindex);
    
    
    // Merged Tree - containing all variables...
    // run and event number (ID) have to be consistent in the merged tree,
    // so it does not matter from where we take them...
    MergedTree->Branch("Runno"                  ,&SIDISrunID            );
    MergedTree->Branch("eventnumber"            ,&SIDISeventID          );
    MergedTree->Branch("Ebeam"                  ,&Ebeam                 );
    MergedTree->Branch("gated_charge"           ,&gated_charge          );
    MergedTree->Branch("livetime"               ,&livetime              );
    MergedTree->Branch("starttime"              ,&starttime             );
    MergedTree->Branch("current"                ,&current               );
    MergedTree->Branch("weight"                 ,&weight                );
    MergedTree->Branch("nMult"                  ,&nMult                 );
    MergedTree->Branch("nHits"                  ,&nHits                 ); // BAND neutrons in BAND skimming
    MergedTree->Branch("eHit"                   ,&eHit                  ); // CLAS12 electrons in BAND skimming
    MergedTree->Branch("goodneutron"            ,&goodneutron           );
    MergedTree->Branch("nleadindex"             ,&nleadindex            );
    MergedTree->Branch("e"                      ,&e                     );
    MergedTree->Branch("Ve"                     ,&Ve                    );
    MergedTree->Branch("Beam"                   ,&Beam                  );
    MergedTree->Branch("q"                      ,&q                     );
    MergedTree->Branch("xB"                     ,&xB                    );
    MergedTree->Branch("Q2"                     ,&Q2                    );
    MergedTree->Branch("omega"                  ,&omega                 );
    MergedTree->Branch("eepiPastCutsInEvent"    ,&eepiPastCutsInEvent   );
    MergedTree->Branch("Npips"                  ,&Npips                 );
    MergedTree->Branch("Npims"                  ,&Npims                 );
    MergedTree->Branch("Nelectrons"             ,&Ne                    );
    MergedTree->Branch("Ngammas"                ,&Ngammas               );
    MergedTree->Branch("Nprotons"               ,&Np                    );
    MergedTree->Branch("Nneutrons"              ,&Nn                    );
    MergedTree->Branch("y"                      ,&y                     );

    MergedTree->Branch("Npips"                    ,&Npips                  , "Npips/I"    );
    MergedTree->Branch("piplus_Px"                ,&piplus_Px              , "piplus_Px[Npips]/D"    );
    MergedTree->Branch("piplus_Py"                ,&piplus_Py              , "piplus_Py[Npips]/D"    );
    MergedTree->Branch("piplus_Pz"                ,&piplus_Pz              , "piplus_Pz[Npips]/D"    );
    MergedTree->Branch("piplus_E"                 ,&piplus_E               , "piplus_E[Npips]/D"    );
    MergedTree->Branch("Vpiplus_X"                ,&Vpiplus_X              , "Vpiplus_X[Npips]/D"    );
    MergedTree->Branch("Vpiplus_Y"                ,&Vpiplus_Y              , "Vpiplus_Y[Npips]/D"    );
    MergedTree->Branch("Vpiplus_Z"                ,&Vpiplus_Z              , "Vpiplus_Z[Npips]/D"    );
    
    
    MergedTree->Branch("Npims"                     ,&Npims                   , "Npims/I"    );
    MergedTree->Branch("piminus_Px"                ,&piminus_Px              , "piminus_Px[Npims]/D"    );
    MergedTree->Branch("piminus_Py"                ,&piminus_Py              , "piminus_Py[Npims]/D"    );
    MergedTree->Branch("piminus_Pz"                ,&piminus_Pz              , "piminus_Pz[Npims]/D"    );
    MergedTree->Branch("piminus_E"                 ,&piminus_E               , "piminus_E[Npims]/D"    );
    MergedTree->Branch("Vpiminus_X"                ,&Vpiminus_X              , "Vpiminus_X[Npims]/D"    );
    MergedTree->Branch("Vpiminus_Y"                ,&Vpiminus_Y              , "Vpiminus_Y[Npims]/D"    );
    MergedTree->Branch("Vpiminus_Z"                ,&Vpiminus_Z              , "Vpiminus_Z[Npims]/D"    );
    
    if (pionCharge=="pi+") {
        MergedTree->Branch("Z"                      ,Zpips                  );
    } else if (pionCharge=="pi-") {
        MergedTree->Branch("Z"                      ,Zpims                  );
    }
    if (fdebug>3) {
        std::cout << "done SetInputAndOutputTTrees()" << std::endl;
    }
}



// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void GetBANDData(int MergedEvtId){
    BANDTree -> GetEntry( BANDEventIndicesToMerge[MergedEvtId] );
    
    bandhit* this_nHit = (bandhit*)nHits->At(nleadindex);
    Pn_Vect = this_nHit->getMomentumN();
    Pn.SetVectM( Pn_Vect , aux.Mn );
    
    // get first electron from BAND TTree to compare with our electron
    clashit* this_eHit = (clashit*)eHit;
    Band_e_Vect.SetMagThetaPhi ( this_eHit->getMomentum(), this_eHit->getTheta(), this_eHit->getPhi() );
    Band_data_e.SetVectM( Band_e_Vect , aux.Me );
    
    // neutron ToF
    n_ToF = this_nHit->getTof();

    // neutron hit position in BAND
    n_HitPos = TVector3(this_nHit->getX(),this_nHit->getY(),this_nHit->getZ());
    
    if (fdebug>1) {
        std::cout
        << "p(n): "   << Pn.P() << " GeV/c" << std::endl
        << "E(n): "   << Pn.E() << " GeV"   << std::endl
        << "ToF(n): " << n_ToF  << " ns"    << std::endl
        << "neutron hit position in BAND: ("
        << std::setprecision(3)
        << n_HitPos.X()  << "," << n_HitPos.Y()  << "," << n_HitPos.Z()  << ") cm"    << std::endl;

    if (fdebug>3) {
        std::cout << "BANDTree->GetEntry("<<BANDEventIndicesToMerge[MergedEvtId]<<")" << std::endl;
        
        std::cout
        << "GetBANDData(" << BANDeventID << "," << MergedEvtId  << ")"
        << std::endl
        << "BANDEventIndicesToMerge["<<MergedEvtId<<"]: "   << BANDEventIndicesToMerge[MergedEvtId] << ","
        << "BANDeventID: "                                  << BANDeventID << ","
        << "Ebeam: "                                        << Ebeam << ","
        << "eepipsPastCutsInEvent: "                        << eepipsPastCutsInEvent    << ","
        << "goodneutron: "                                  << goodneutron              << ","
        << "E(n): "                                         << Pn.E()                   << " GeV,"
        << "this_nHit->getMomentumN()->Mag(): "             << this_nHit->getMomentumN().Mag() << " GeV/c,"
        << "p(n): "                                         << Pn.P()                   << " GeV/c,"
        << "Pe: "                                           << e->P()                   << " GeV/c,"
        << std::endl;
    }
    }
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void GetSIDISData( int MergedEvtId ){
    if (fdebug>2) { std::cout << "GetSIDISData(" << MergedEvtId  << ")" << std::endl;}
    
    SIDISTree -> GetEntry( SIDISEventIndicesToMerge[MergedEvtId] );
    
    eepiPastCutsInEvent = false;
    Npions              = 0;
    eventnumber         = SIDISeventID;
    
    if (pionCharge=="pi+") {
        eepiPastCutsInEvent = eepipsPastCutsInEvent;
        Npions = Npips;
    } else if (pionCharge=="pi-") {
        eepiPastCutsInEvent = eepimsPastCutsInEvent;
        Npions = Npims;
    }
    for (int pipsIdx=0; pipsIdx<Npips; pipsIdx++){
        piplus.at(pipsIdx)   = TLorentzVector(piplus_Px[pipsIdx],
                                              piplus_Py[pipsIdx],
                                              piplus_Pz[pipsIdx],
                                              piplus_E[pipsIdx]) ;
        
        Vpiplus.at(pipsIdx)  = TVector3(Vpiplus_X[pipsIdx],
                                        Vpiplus_Y[pipsIdx],
                                        Vpiplus_Z[pipsIdx]) ;
        if (fdebug>4) {
            std::cout
            << "piplus.at(pipsIdx) = TLorentzVector("
            << piplus_Px[pipsIdx]   <<","
            << piplus_Py[pipsIdx]   <<","
            << piplus_Pz[pipsIdx]   <<","
            << piplus_E[pipsIdx]    <<")); "
            << std::endl
            << "piplus_qFrame_pT: " << piplus_qFrame_pT[pipsIdx]    <<","
            << "piplus_qFrame_pL: " << piplus_qFrame_pL[pipsIdx]    <<","
            << std::endl
            << "Z: "                << Zpips[pipsIdx]               <<","
            << "Z_LC: "             << Zpips_LC[pipsIdx]            <<","
            << std::endl;
        }
    }
    for (int pimsIdx=0; pimsIdx<Npims; pimsIdx++){
        piminus.at(pimsIdx)  = TLorentzVector(piminus_Px[pimsIdx],
                                              piminus_Py[pimsIdx],
                                              piminus_Pz[pimsIdx],
                                              piminus_E[pimsIdx]) ;
        Vpiminus.at(pimsIdx) = TVector3(Vpiminus_X[pimsIdx],
                                        Vpiminus_Y[pimsIdx],
                                        Vpiminus_Z[pimsIdx]) ;
    }
    if (fdebug>2) {
        std::cout
        << "Merging event " << EventNumbersToMerge[MergedEvtId]
        << " ("  << MergedEvtId << "/" << NeventsToMerge << ")" << std::endl;
        
        std::cout
        << "SIDISEventIndicesToMerge["<<MergedEvtId<<"]: "  << SIDISEventIndicesToMerge[MergedEvtId] << ","
        << "SIDISeventID: "                                 << SIDISeventID << ","
        << std::endl
        << "eepipsPastCutsInEvent: "                        << eepipsPastCutsInEvent << ","
        << "Nelectrons: "                                   << Ne       << ","
        << "Npions: "                                       << Npions   << ","
        << "Npips: "                                        << Npips    << ","
        << "Npims: "                                        << Npims    << ","
        << std::endl
        << "Pe: "                                           << e->P() << " GeV/c,"
        << "q: "                                            << q->P() << " GeV/c,"
        << "omega: "                                        << q->E() << " GeV,"
        << std::endl;
        if (Npips>0){
            std::cout << "1st pi+: Ppi: "                  << piplus[0].P() << " GeV/c,"
            << std::endl;
        }
        if (Npims>0){
            std::cout << "1st pi-: Ppi: "                  << piminus[0].P() << " GeV/c,"
            << std::endl;
        }
    }
}


//// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
//void ComputeElectronKinematics(){
//    // compute kinematics
//    // SIDISc12rSkimmer.C already computes few of the kinematical variables:
//    //     xB,  Q2, omega, W, Z, y, Z_LC
//}

//// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
//void ComputeNeutronAndProtonKinematics(){
//
//}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ComputeKinematics(TLorentzVector pi){
    
    // SIDISc12rSkimmer.C already computes few of the kinematical variables:
    //     xB,  Q2, omega, W, y
    
    // detected neutron
    En      = Pn.E();
    Es      = Pn.E();
    Ps      = Pn.P();
    theta_sq= Pn.Angle( q->Vect() );
    
    // rotate neutron to q-Frame
    TVector3 Pn_3Vector = RotateVectorTo_qFrame( Pn.Vect() );
    Pn_qFrame.SetVectM( Pn_3Vector, aux.Mn );
    alpha_n  = aux.ComputeLightConeFraction( Pn_qFrame );
    
    // initial state of the virtual moving proton in the nucleus
    E_init  = aux.Md - Pn.E();
    p_init  = d_rest - Pn;
    v_init  = p_init.Mag2() - aux.Mp2;
    // rotate (initial) virtual moving proton to the q-Frame
    TVector3 p_init_3Vector = RotateVectorTo_qFrame( p_init.Vect() );
    p_init_qFrame.SetVectM( p_init_3Vector, E_init );
    alpha_init  = aux.ComputeLightConeFraction( p_init_qFrame );
    
    // pion
    // move to q-Frame
    pi_qFrame.SetVectM( RotateVectorTo_qFrame( pi.Vect() ), aux.Mpi );
    
    // pion energy fraction    
    Zpi          = pi.E()/omega;
    Zpi_LC       = (pi_qFrame.E() - pi_qFrame.Pz()) / (q->E() - q->P());
    Zpi_LC_Prime = Zpi_LC / alpha_n;
    
    // Kinematics assuming scattering off a proton at rest
    // W is read off the SIDIS TTree
    // Kinematics for the virtual moving proton
    //      W = ( p_rest + q ).Mag();
    W_Prime   = ( p_init + *q ).Mag();
        
    M_x       = ( p_rest + *q - pi ).Mag();
    M_x_Prime = ( p_init + *q - pi ).Mag();
        
    xF        = 2. * (pi.Vect().Dot(q->Vect())) / (q->P() * W);
    xF_Prime  = 2. * (pi.Vect().Dot(q->Vect())) / (q->P() * W_Prime);
    
    eta_pi    = 0.5 * log((pi_qFrame.E()+pi_qFrame.Pz()) /
                          (pi_qFrame.E()-pi_qFrame.Pz()));
        
    // xB is read off the SIDIS TTree
    //     xB = Q2 / (2. * aux.Mp * omega);
    xB_Prime  = Q2 / (W_Prime*W_Prime - aux.Mp2 + Q2);
    
    
    if (fdebug>3) {
        std::cout
        << "ComputeElectronKinematics()"
        << std::endl
        << "Pe: "       << e->P()     << " GeV/c,"
        << "xB: "       << xB         << ","
        << "xB': "      << xB_Prime   << ","
        << std::endl
        << "(p_init + *q).Px(): " << (p_init + *q).Px() << ","
        << "(p_init + *q).Py(): " << (p_init + *q).Py() << ","
        << "(p_init + *q).Pz(): " << (p_init + *q).Pz() << ","
        << "(p_init + *q).E(): "  << (p_init + *q).E()  << ","
        << "(p_init + *q).Mag(): "<< (p_init + *q).Mag()  << ","
        << "(p_init + *q).Mag2(): "<< (p_init + *q).Mag2()  << ","
        << std::endl
        << "-Pn: "      << (-Pn).P() << " GeV/c,"
        << "W: "        << W        << " GeV/c2,"
        << "W': "       << W_Prime   << " GeV/c2,"
        << std::endl;
    }
}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void PrintDone(){
    std::cout << "Done merging SIDIS and BAND events from run " <<  SIDISrunID << ", execution time: "
    << double(clock() - tStart) / (double)CLOCKS_PER_SEC
    << " sec "<< std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Stream_e_pi_n_line_to_CSV(int piIdx,
                               bool passed_cuts_e_pi_kinematics,
                               bool passed_cuts_n){
    status = 0;
    TLorentzVector       * pi;
//    TLorentzVector  pi_qFrame;
    TVector3            * Vpi;
    double       pi_DC_sector;
    double       pi_qFrame_pT;
    double       pi_qFrame_pL;
    double    pi_qFrame_Theta;
    double      pi_qFrame_Phi;
    
    if (pionCharge=="pi+") {
        pi              = &piplus               [piIdx];
        Vpi             = &Vpiplus              [piIdx];
        pi_DC_sector    = pips_DC_sector        [piIdx];
        pi_qFrame_pT    = piplus_qFrame_pT      [piIdx];
        pi_qFrame_pL    = piplus_qFrame_pL      [piIdx];
        pi_qFrame_Theta = piplus_qFrame_Theta   [piIdx];
        pi_qFrame_Phi   = piplus_qFrame_Phi     [piIdx];
    }
    else if (pionCharge=="pi-") {
        pi              = &piminus              [piIdx];
        Vpi             = &Vpiminus             [piIdx];
        pi_DC_sector    = pims_DC_sector        [piIdx];
        pi_qFrame_pT    = piminus_qFrame_pT     [piIdx];
        pi_qFrame_pL    = piminus_qFrame_pL     [piIdx];
        pi_qFrame_Theta = piminus_qFrame_Theta  [piIdx];
        pi_qFrame_Phi   = piminus_qFrame_Phi    [piIdx];
    }
    else {
        std::cout << "bad pion charge in Stream_e_pi_line_to_CSV()!" << std::endl;
        return;
    }
        
    ComputeKinematics( *pi );

    
    std::vector<double> variables =
    {   (double)status, (double)SIDISrunID, (double)eventnumber,    (double)beam_helicity,
        e->P(),         e->Theta(),         e->Phi(),               Ve->Z(),
        pi->P(),        pi->Theta(),        pi->Phi(),              Vpi->Z(),
        Pn.P(),         Pn.Theta(),         Pn.Phi(),               Ve->Z(),
        n_HitPos.X(),   n_HitPos.Y(),       n_HitPos.Z(),
        Pn.E(),         alpha_n,
        Q2,             xB,                 omega,                  y,
        (double)e_DC_sector,                (double)pi_DC_sector,
        pi_qFrame_Theta,                    pi_qFrame_Phi,
        pi_qFrame_pT,                       pi_qFrame_pL,
        Pn_qFrame.Theta(),                  Pn_qFrame.Phi(),
        Pn_qFrame.Pt(),                     Pn_qFrame.Pz(),
        E_init,         p_init.P(),         alpha_init,             v_init,
        Zpi,            Zpi_LC,
        Zpi_LC_Prime,
        W,              M_x,
        xF,             eta_pi,
        W_Prime,        M_x_Prime,
        xF_Prime,
        xB_Prime,
    };
    
    if (fdebug>2){
        std::cout
        << "Stream_e_pi_n_line_to_CSV(piIdx "         << piIdx  << "), "
        << std::endl;
        if (passed_cuts_e_pi_kinematics)
            std::cout << "passed (e,e'π) kinematical cuts" << std::endl;
        if (passed_cuts_n)
            std::cout << "passed (e,e'πn) cuts" << std::endl;
        
        std::cout
        << "E(π): "         << pi->E()          << " GeV, "
        << "Z(π): "         << Zpi              << ", "
        << "Z(π)_LC: "      << Zpi_LC           << ", "
        << "Z'(π)_LC: "     << Zpi_LC_Prime     << ", "
        << std::endl;
        
        if ((pi->P() < 1.25) && (passed_cuts_e_pi_kinematics)){
            std::cout
            << " pi->P() = " << pi->P() << " < 1.25 GeV/c!"
            << std::endl
            << "How the F**K is this possible??"
            << std::endl
            << "SIDISrunID: "                   << SIDISrunID
            << ", event number "                << eventnumber
            << ", piIdx "                       << piIdx
            << ", passed_cuts_e_pi_kinematics " << passed_cuts_e_pi_kinematics
            << std::endl;
        }
    }
    
    if (passed_cuts_e_pi_kinematics){
        aux.StreamToCSVfile( CSVfile_e_pi_n, variables, csvprecisions );
    } else {
        if (fdebug>2) {
            std::cout << "Did not pass cuts, not writing to CSV file" << std::endl;
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void InitializeVariables(){
    Pn_Vect     = TVector3();
    Band_e_Vect = TVector3();
    e           = new TLorentzVector( 0, 0, 0, aux.Me );
    Band_data_e = TLorentzVector( 0, 0, 0, aux.Me );
    
    xB          = Q2        = omega     = -9999;
    y                                   = -9999;
    M_x         = M_x_Prime             = -9999;
    W           = W_Prime               = -9999;
    xF          = xF_Prime              = -9999;
    Ve                                  = new TVector3();
    n_ToF                               = -9999;
    n_HitPos                            = TVector3();
    
    piplus      .clear();
    piminus     .clear();
    Vpiplus     .clear();
    Vpiminus    .clear();
    piplus         .clear();
    piplus_qFrame  .clear();
    piminus        .clear();
    piminus_qFrame .clear();
    
    for (int piIdx=0; piIdx<NMAXPIONS; piIdx++) {
        piplus  .push_back( TLorentzVector(0,0,0,aux.Mpips) );
        Vpiplus .push_back( TVector3() );
        eepipsPastKinematicalCuts[piIdx]            = false;
        piminus .push_back( TLorentzVector(0,0,0,aux.Mpims) );
        Vpiminus.push_back( TVector3() );
        pimsPastSelectionCuts[piIdx]                = false;
        eepimsPastKinematicalCuts[piIdx]            = false;
        piplus_Px[piIdx]    = piplus_Py[piIdx]  = -9999;
        piplus_Pz[piIdx]    = piplus_E[piIdx]   = -9999;
        piminus_Px[piIdx]   = piminus_Py[piIdx] = -9999;
        piminus_Pz[piIdx]   = piminus_E[piIdx]  = -9999;
        Vpiplus_X[piIdx]    = Vpiplus_Y[piIdx]  = -9999;
        Vpiplus_Z[piIdx]                        = -9999;
        Vpiminus_X[piIdx]   = Vpiminus_Y[piIdx] = -9999;
        Vpiminus_Z[piIdx]                       = -9999;
        piminus_qFrame_pT[piIdx] = piminus_qFrame_pL[piIdx] = -9999;
        piplus_qFrame_pT[piIdx]  = piplus_qFrame_pL[piIdx]  = -9999;
        
        piplus        .push_back( TLorentzVector(0,0,0,aux.Mpips) );
        piplus_qFrame .push_back( TLorentzVector(0,0,0,aux.Mpips) );
        piminus       .push_back( TLorentzVector(0,0,0,aux.Mpims) );
        piminus_qFrame.push_back( TLorentzVector(0,0,0,aux.Mpims) );
        Pn_qFrame     = TLorentzVector(0,0,0,aux.Mn) ;
        p_init_qFrame = TLorentzVector(0,0,0,aux.Mp) ;
    }
    status                                          = 1; // 0 is good...
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MergeEventData(){
    if (fdebug>2){
        std::cout
        << "MergeEventData( event " << BANDeventID << " from run " << BANDrunID << ")"
        << std::endl        << std::endl;
    }

    
    // fill output TTree and CSV file only if the neutron is a "good-neutron"
    // Note, that we want to fill the ROOT TTree once,
    // and not multiple times as we do for the CSV file (for each pi-index)
    if (goodneutron){ MergedTree  -> Fill();}
    
    for (int piIdx=0; piIdx<Npions; piIdx++) {
        bool eepiPastKinematicalCuts = false;
        if (pionCharge=="pi+") {
            eepiPastKinematicalCuts = eepipsPastKinematicalCuts[piIdx];
        } else if (pionCharge=="pi-") {
            eepiPastKinematicalCuts = eepimsPastKinematicalCuts[piIdx];
        }
        
        if (goodneutron){
            Stream_e_pi_n_line_to_CSV(piIdx,
                                      eepiPastKinematicalCuts, goodneutron);
            
        }
    }
    if (fdebug>2) { std::cout << "-----------------------------------------------------------------" << std::endl; }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MoveTo_qFrame(){
    if (fdebug>1){
        std::cout << "Moving to q-Frame" <<std::endl;
    }
    // Move to the "q-frame" and define the pion momentum in this frame
    // This is not boost, only 3D rotation
    // q-frame is defined as follows:
    // z axis is defined by the q - parallel to q
    // x axis is defined by the e' - such that p(e') resides in the x-z plane

    // (1) define q-angles
    q_phi   = q->Phi();
    q_theta = q->Theta();

    // (2) rotate Pe and q according to q angles
    TVector3 Pe = e->Vect();
    TVector3 Pq = q->Vect();
    
    Pe     .RotateZ(-q_phi);
    Pe     .RotateY(-q_theta);
    Pe_phi = Pe.Phi();
    Pe     .RotateZ(-Pe_phi);
    Pq     .RotateZ(-q_phi);
    Pq     .RotateY(-q_theta);

    
    // (3) verify on q and Pe that the frame-change is done correctly
    //    RotateVectorTo_qFrame( &Pe );
    e_qFrame.SetVectM( Pe, aux.Me );
    //    RotateVectorTo_qFrame( &Pq );
    q_qFrame.SetVectM( Pq, q->M() );
    
    if (fdebug>2){
        aux.Print4Vector( *e, "e" );
        aux.Print4Vector( e_qFrame, "e in q-Frame" );
        aux.Print4Vector( *q , "q");
        aux.Print4Vector( q_qFrame, "q in q-Frame" );
    }
    
    
    // (4) rotate pions to this q-frame
    for (int piIdx=0; piIdx<Npips; piIdx++) {
        TVector3 Ppiplus = RotateVectorTo_qFrame( piplus.at(piIdx).Vect() );
        piplus_qFrame.at(piIdx).SetVectM( Ppiplus, aux.Mpi  );
        // Double_t variables like piplus_qFrame_pT are read-off from SIDIS TTree
    }
    for (int piIdx=0; piIdx<Npims; piIdx++) {
        TVector3 Ppiminus = RotateVectorTo_qFrame( piminus.at(piIdx).Vect() );
        piminus_qFrame.at(piIdx).SetVectM( Ppiminus, aux.Mpi );
        // Double_t variables like piminus_qFrame_pT are read-off from SIDIS TTree
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


//
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//void Print4Vector( TLorentzVector v, std::string label ){
//    std::cout << label << " 4-vector:"<<std::endl;
//    std::cout
//    << "(Px,Py,Pz,E) = (" << v.Px() << "," << v.Py() << "," << v.Pz() << "," << v.E()
//    << "), M = " << v.Mag()
//    << std::endl;
//}
