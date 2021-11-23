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



#define NMAXEVENTS 5000000
#define NMAXPIONS 20 // maximal allowed number of pions
#define r2d 180./3.1415 // radians to degrees

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
// Globals
// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
TString DataPath = "/volatile/clas12/users/ecohen/BAND/";
TString   skimmedBANDFilename;
TString  skimmedSIDISFilename;
TString            pionCharge; // "pi+" or "pi-"
TString               pionStr;
TString            csvheader = ((TString)"status,runnum,evnum,beam_helicity,"
                                +(TString)"e_P,e_Theta,e_Phi,e_Vz,"
                                +(TString)"pi_P,pi_Theta,pi_Phi,pi_Vz,"
                                +(TString)"n_P,n_Theta,n_Phi,n_Vz,"
                                +(TString)"Q2,W,xB,Zpi,"
                                +(TString)"omega,xF,y,"
                                +(TString)"M_X_ee_pi,M_X_ee_pi_n,xPrime1,xPrime2,"
                                +(TString)"theta_sq,WPrime,");

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
Double_t       Me  = 0.00051099895; // GeV/c2
Double_t       Mpims  = 0.13957039; // GeV/c2
Double_t       Mpips  = 0.13957039; // GeV/c2
Double_t            Mp  = 0.938272; // GeV/c2
Double_t             Mn = 0.939565; // GeV/c2
Double_t                Md = 1.875; // GeV/c2
Double_t             Mp2 = Mp * Mp;
Double_t                        xB; // Bjorken x
Double_t                        xF; // Feynman x
Double_t                        Q2;
Double_t                     omega;
Double_t                        w2; // omega^2
Double_t          Zpips[NMAXPIONS]; // energy fraction rest frame
Double_t          Zpims[NMAXPIONS]; // energy fraction rest frame
Double_t                     W, W2; // energy of the hadronic system
Double_t                   alpha_s; // light cone fraction of momentum of the recoil neutron
Double_t           WPrime, W2prime;  // moving proton
Double_t                   xPrime1; // moving proton defined as
// x' = Q2 / (2. * ((Md - Es) * omega + Pn_Vect*q->Vect() ));
Double_t                   xPrime2; // moving proton defined as
// x' = Q2 / (W'^2 - mN^2 + Q2)

Double_t                        Es; // spectator energy
Double_t                        Ps; // spectator momentum
Double_t                  theta_sq; // spectator angle with respect to momentum transfer
Double_t     M_X_ee_pi,M_X_ee_pi_n;
Double_t                         y;



Int_t                       NeventsToMerge;
Int_t               BANDrunID, BANDeventID;
Int_t             SIDISrunID, SIDISeventID;
Int_t          EventIDsToMerge[NMAXEVENTS];
Int_t                 Npions, Npips, Npims;
Int_t                  Ne, Np, Nn, Ngammas;
Int_t                               status;
std::vector<Int_t>                   BANDeventIDs;
std::vector<Int_t>                  SIDISeventIDs;
std::vector<Int_t>               BANDeventIndices;
std::vector<Int_t>              SIDISeventIndices;
std::vector<Int_t>            EventNumbersToMerge;
std::vector<std::size_t>  BANDEventIndicesToMerge;
std::vector<std::size_t> SIDISEventIndicesToMerge;

// SIDIS Tree
TLorentzVector     *target = new TLorentzVector(0, 0, 0, Md );
TLorentzVector     *Beam=0;
TLorentzVector        *e=0;
TLorentzVector        *q=0;
TLorentzVector          Pn; // neutron momentum
TLorentzVector Band_data_e; // electron information in BAND TTree
// reconstructed vertex position
TVector3             *Ve=0;
TVector3           Pn_Vect; // neutron 3-momentum
TVector3       Band_e_Vect;

std::vector<TVector3>               Vpiplus;
std::vector<TVector3>              Vpiminus;
std::vector<TLorentzVector>          piplus;
std::vector<TLorentzVector>         piminus;

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

double           piplus_Px[NMAXPIONS];
double           piplus_Py[NMAXPIONS];
double           piplus_Pz[NMAXPIONS];
double            piplus_E[NMAXPIONS];
double           Vpiplus_X[NMAXPIONS];
double           Vpiplus_Y[NMAXPIONS];
double           Vpiplus_Z[NMAXPIONS];

double          piminus_Px[NMAXPIONS];
double          piminus_Py[NMAXPIONS];
double          piminus_Pz[NMAXPIONS];
double           piminus_E[NMAXPIONS];
double          Vpiminus_X[NMAXPIONS];
double          Vpiminus_Y[NMAXPIONS];
double          Vpiminus_Z[NMAXPIONS];


clashit *eHit             = new clashit;
//TClonesArray      * eHit  = new TClonesArray("clashit"); // CLAS12 electrons in BAND analysis
TClonesArray      * nHits = new TClonesArray("bandhit"); // BAND neutrons in BAND analysis
TClonesArray    * mcParts = new TClonesArray("genpart");
// TClonesArray  &save_e_Hit = *eHit;
// TClonesArray     &saveHit = *nHits;
// TClonesArray      &saveMC = *mcParts;

void             OpenInputFiles (TString RunStr);
void            OpenOutputFiles (TString RunStr);
void            StreamToCSVfile (std::vector<Double_t> observables);
void            CloseInputFiles ();
void           CloseOutputFiles (TString OutDataPath);
void    MergeSIDISandBANDevents (int NMAXeventsToMerge=10,
                                 int PrintProgress=5000);
Int_t CreateListOfEventsToMerge (TTree * BANDTree,
                                 TTree * SIDISTree,
                                 int NMAXeventsToMerge=-1);
void  Stream_e_pi_n_line_to_CSV (int piIdx,
                                 bool passed_cuts_e_pi_kinematics,
                                 bool passed_cuts_n);


void        InitializeVariables ();
void               GetSIDISData ( int SIDISeventID, int MergedEvtId );
void                GetBANDData ( int BANDeventID, int MergedEvtId );
void             MergeEventData ();
void    SetInputAndOutputTTrees ();
void          ComputeKinematics ();
void                  PrintDone ();
void          PrintMonitorHello ();
void              SetPionCharge ( TString fpionCharge ) {
    pionCharge = fpionCharge;
    if (pionCharge=="pi+") {
        pionStr = "_e_piplus";
    } else if (pionCharge=="pi-") {
        pionStr = "_e_piminus";
    } else {
        pionStr = "_no_pion_charge_info";
    }
};
void                SetDataPath ( TString fDataPath )   {DataPath = fDataPath + "/";};
void               SetVerbosity ( int ffdebug )         {fdebug = ffdebug;};
void                  PrintTime ( TString prefix ){
    std::cout << prefix << ", after "
    << double(clock() - tStart) / (double)CLOCKS_PER_SEC
    << " sec "<< std::endl;
}



// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
// Main functionality
// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void MergeSIDISandBANDSkimmers(int RunNumber=6420,
                               TString fpionCharge="pi+", // "pi+" or "pi-"
                               int NMAXeventsToMerge=-1,
                               int ffdebug=1,
                               int PrintProgress=1000,
                               TString fDataPath="/volatile/clas12/users/ecohen/BAND/"){
    
    SetPionCharge    ( fpionCharge );
    SetVerbosity     ( ffdebug );
    SetDataPath      ( fDataPath );
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
        
    if (fdebug>2)
        PrintTime((TString)"Created list of (" + (TString)std::to_string(NeventsToMerge) + (TString)") events to merge");
    
    // assign TTree branches to variables
    SetInputAndOutputTTrees ();
        
    // step over list of events-to-merge and merge them...
    for (int MergedEvtId=0; MergedEvtId < NeventsToMerge; MergedEvtId++) {
        
        // initialize
        InitializeVariables ();

        // grab electron and pion information from SIDIS TTree
        GetSIDISData( SIDISeventID, MergedEvtId );
        
        // grab neturon information from BAND
        GetBANDData( BANDeventID, MergedEvtId );
        
        // compute kinematical variables, also for the neutron
        ComputeKinematics   ();
        
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
                            + "skimmed_SIDIS_inc_"  + RunStr + pionStr + ".root");
    if (fdebug>2) std::cout << "Opening " << skimmedSIDISFilename << std::endl;
    SIDISFile                     = new TFile( skimmedSIDISFilename );
    SIDISTree                     = (TTree*)SIDISFile->Get("tree");
}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void OpenOutputFiles (TString RunStr){
    
    TString skimmedMergedFilename = (DataPath + "merged_SIDIS_and_BAND_skimming/"
                                     + "skimmed_SIDIS_and_BAND_inc_"  + RunStr + pionStr + "_n" );
    
    if (fdebug>2) std::cout << "Opening output file: " << skimmedMergedFilename  << ".root/csv " << std::endl;
    
    
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


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void StreamToCSVfile (std::vector<Double_t> observables){
    for (auto v:observables) CSVfile_e_pi_n << v << ",";
    CSVfile_e_pi_n << std::endl;
    if (fdebug>3) {
        std::cout << "StreamToCSVfile()" << std::endl;
        std::cout << csvheader << std::endl;
        for (auto v:observables) std::cout << v << ",";
        std::cout << std::endl;
    }
}


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
    
    
    // First, define two vectors that containt the event IDs in each TTree
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



// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
Double_t ComputeLightConeFraction( TLorentzVector p ){
    // compute light-cone momentum fraction
    Double_t m = p.Mag();
    Double_t alpha = (p.E() - p.Z())/m;
    return alpha;
}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void SetInputAndOutputTTrees (){
    if (fdebug>3) {
        std::cout << "SetInputAndOutputTTrees()" << std::endl;
    }
    
    // SIDIS input tree branches
    SIDISTree  -> SetBranchAddress("eventnumber"                ,&SIDISeventID          );
    SIDISTree  -> SetBranchAddress("runnum"                     ,&SIDISrunID            );
    SIDISTree  -> SetBranchAddress("e"                          ,&e                     );
    SIDISTree  -> SetBranchAddress("Ve"                         ,&Ve                    );
    SIDISTree  -> SetBranchAddress("Beam"                       ,&Beam                  );
    SIDISTree  -> SetBranchAddress("beam_helicity"              ,&beam_helicity         );
    SIDISTree  -> SetBranchAddress("q"                          ,&q                     );
    SIDISTree  -> SetBranchAddress("omega"                     ,&omega                  );
    SIDISTree  -> SetBranchAddress("Z"                         ,&Zpips                  );
    SIDISTree  -> SetBranchAddress("Nelectrons"                ,&Ne                     );
    SIDISTree  -> SetBranchAddress("Ngammas"                   ,&Ngammas                );
    SIDISTree  -> SetBranchAddress("Nprotons"                  ,&Np                     );
    SIDISTree  -> SetBranchAddress("Nneutrons"                 ,&Nn                     );
    SIDISTree  -> SetBranchAddress("Npips"                     ,&Npips                  );
    SIDISTree  -> SetBranchAddress("Npims"                     ,&Npims                  );
    
    
    
    //    // branches that depend on pion charge
    if (pionCharge=="pi+") {
        SIDISTree  -> SetBranchAddress("eepipsPastCutsInEvent"     ,&eepipsPastCutsInEvent    );
        SIDISTree  -> SetBranchAddress("eepipsPastKinematicalCuts" ,&eepipsPastKinematicalCuts);
        SIDISTree -> SetBranchAddress("piplus_Px"                  ,&piplus_Px                );
        SIDISTree -> SetBranchAddress("piplus_Py"                  ,&piplus_Py                );
        SIDISTree -> SetBranchAddress("piplus_Pz"                  ,&piplus_Pz                );
        SIDISTree -> SetBranchAddress("piplus_E"                   ,&piplus_E                 );
        SIDISTree -> SetBranchAddress("Vpiplus_X"                  ,&Vpiplus_X                );
        SIDISTree -> SetBranchAddress("Vpiplus_Y"                  ,&Vpiplus_Y                );
        SIDISTree -> SetBranchAddress("Vpiplus_Z"                  ,&Vpiplus_Z                );
    } else if (pionCharge=="pi-") {
        SIDISTree  -> SetBranchAddress("eepimsPastCutsInEvent"     ,&eepimsPastCutsInEvent      );
        SIDISTree  -> SetBranchAddress("eepimsPastKinematicalCuts" ,&eepimsPastKinematicalCuts  );
        SIDISTree -> SetBranchAddress("piminus_Px"                  ,&piminus_Px                );
        SIDISTree -> SetBranchAddress("piminus_Py"                  ,&piminus_Py                );
        SIDISTree -> SetBranchAddress("piminus_Pz"                  ,&piminus_Pz                );
        SIDISTree -> SetBranchAddress("piminus_E"                   ,&piminus_E                 );
        SIDISTree -> SetBranchAddress("Vpiminus_X"                  ,&Vpiminus_X                );
        SIDISTree -> SetBranchAddress("Vpiminus_Y"                  ,&Vpiminus_Y                );
        SIDISTree -> SetBranchAddress("Vpiminus_Z"                  ,&Vpiminus_Z                );
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
    
    MergedTree->Branch("piplus_Px"                ,&piplus_Px              , "piplus_Px[20]/D"    );
    MergedTree->Branch("piplus_Py"                ,&piplus_Py              , "piplus_Py[20]/D"    );
    MergedTree->Branch("piplus_Pz"                ,&piplus_Pz              , "piplus_Pz[20]/D"    );
    MergedTree->Branch("piplus_E"                 ,&piplus_E               , "piplus_E[20]/D"    );
    MergedTree->Branch("Vpiplus_X"                ,&Vpiplus_X              , "Vpiplus_X[20]/D"    );
    MergedTree->Branch("Vpiplus_Y"                ,&Vpiplus_Y              , "Vpiplus_Y[20]/D"    );
    MergedTree->Branch("Vpiplus_Z"                ,&Vpiplus_Z              , "Vpiplus_Z[20]/D"    );
    
    
    MergedTree->Branch("piminus_Px"                ,&piminus_Px              , "piminus_Px[20]/D"    );
    MergedTree->Branch("piminus_Py"                ,&piminus_Py              , "piminus_Py[20]/D"    );
    MergedTree->Branch("piminus_Pz"                ,&piminus_Pz              , "piminus_Pz[20]/D"    );
    MergedTree->Branch("piminus_E"                 ,&piminus_E               , "piminus_E[20]/D"    );
    MergedTree->Branch("Vpiminus_X"                ,&Vpiminus_X              , "Vpiminus_X[20]/D"    );
    MergedTree->Branch("Vpiminus_Y"                ,&Vpiminus_Y              , "Vpiminus_Y[20]/D"    );
    MergedTree->Branch("Vpiminus_Z"                ,&Vpiminus_Z              , "Vpiminus_Z[20]/D"    );
    
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
void PrintMonitorHello(){
    std::cout << "Hello..." << std::endl;
    std::cout << "Is it me you're looking for?..." << std::endl;
    std::cout << "I can see it in your eyes..." << std::endl;
    std::cout << "I can see it in your smile..." << std::endl;
    std::cout << "You're all I've ever wanted, and your arms are open wide..." << std::endl;
    std::cout << "Cause you know just what to say, and you know just what to do..." << std::endl;
    std::cout << "And I want to tell you so much, I love you" << std::endl;
    std::cout << "I long to see the sunlight in your hair" << std::endl;
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void GetBANDData(int BANDeventID, int MergedEvtId){
    BANDTree -> GetEntry( BANDEventIndicesToMerge[MergedEvtId] );
    
    bandhit* this_nHit = (bandhit*)nHits->At(nleadindex);
    Pn_Vect = this_nHit->getMomentumN();
    Pn.SetVectM( Pn_Vect , Mn );
    
    // get first electron from BAND TTree to compare with our electron
    clashit* this_eHit = (clashit*)eHit;
    Band_e_Vect.SetMagThetaPhi ( this_eHit->getMomentum(), this_eHit->getTheta(), this_eHit->getPhi() );
    Band_data_e.SetVectM( Band_e_Vect , Me );

    if (fdebug>1) {
        std::cout << "p(n): " << Pn.P() << std::endl;
    }
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

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void GetSIDISData( int SIDISeventID, int MergedEvtId ){
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
        piplus.push_back( TLorentzVector(piplus_Px[pipsIdx], piplus_Py[pipsIdx], piplus_Pz[pipsIdx],
                                         piplus_E[pipsIdx]) );
        Vpiplus.push_back( TVector3(Vpiplus_X[pipsIdx], Vpiplus_Y[pipsIdx], Vpiplus_Z[pipsIdx]) );
    }
    for (int pimsIdx=0; pimsIdx<Npims; pimsIdx++){
        piminus.push_back( TLorentzVector(piminus_Px[pimsIdx], piminus_Py[pimsIdx], piminus_Pz[pimsIdx],
                                         piminus_E[pimsIdx]) );
        Vpiminus.push_back( TVector3(Vpiminus_X[pimsIdx], Vpiminus_Y[pimsIdx], Vpiminus_Z[pimsIdx]) );
    }
    if (fdebug>3) {
        std::cout
        << "Merging event " << EventNumbersToMerge[MergedEvtId]
        << " ("  << MergedEvtId+1 << "/" << NeventsToMerge << ")" << std::endl;
        
        std::cout
        << "GetSIDISData(" << SIDISeventID << "," << MergedEvtId  << ")"
        << std::endl
        << "SIDISEventIndicesToMerge["<<MergedEvtId<<"]: "  << SIDISEventIndicesToMerge[MergedEvtId] << ","
        << "SIDISeventID: "                                 << SIDISeventID << ","
        << "E(electron): "                                  << e->E() << " GeV,"
        << "eepipsPastCutsInEvent: "                        << eepipsPastCutsInEvent << ","
        << "Npions: "                                       << Npions   << ","
        << "Npips: "                                        << Npips    << ","
        << std::endl
        << "Pe: "                                           << e->P() << " GeV/c,"
        << "E(e): "                                         << e->E() << " GeV,"
        << "q: "                                            << q->P() << " GeV/c,"
        << "omega: "                                        << q->E() << " GeV,"
        << std::endl;
        if (Npips>0 && pionCharge=="pi+"){
            std::cout << "piplus_Px[0]: "                   << piplus_Px[0] <<  " GeV/c,"
            << std::endl;
        }
    }
}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void ComputeKinematics(){
    // compute kinematics
    // SIDISc12rSkimmer.C already computes few of the kinematical variables:
    //     xB,  Q2, omega, W, Z, y
    
    Es      = Pn.E();
    Ps      = Pn.P();
    omega   = q->E();
    w2      = omega * omega;
    y       = omega / Ebeam;
    theta_sq= Pn.Angle( q->Vect() );
    
    Q2      = -q->Mag2();
    omega   = q->E();
    
    W2      = Mp2 - Q2 + 2. * omega * Mp;
    W       = sqrt(W2);
    
    // W' from [E12-11-003A proposal, p. 13]
    W2prime = Mp2 - Q2 + 2. * omega * (Md - Es) + 2. * Ps * sqrt(Q2 + w2) * cos( theta_sq );
    WPrime  = sqrt(W2prime);

    xB      = Q2 / (2. * Mp * omega);
    xPrime1 = Q2 / (2. * ((Md - Es) * omega + Pn_Vect*q->Vect() ));
    xPrime2 = Q2 / (W2prime - Mp2 + Q2);
    
    // move to q-frame
    // compute light-cone fraction of momentum
    // alpha_s = ComputeLightConeFraction( Pn_qFrame );
    
    
    if (fdebug>3) {
        std::cout
        << "ComputeKinematics()"
        << std::endl
        << "Pe: "       << e->P()   << " GeV/c,"
        << "x: "        << xB       << ","
        << "x'(1): "       << xPrime1   << ","
        << "x'(2): "       << xPrime2   << ","
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
    
    TLorentzVector  * pi;
    TVector3        * Vpi;
    double          Zpi;
    if (pionCharge=="pi+") {
        pi  = &piplus[piIdx];
        Vpi = &Vpiplus[piIdx];
        Zpi = Zpips  [piIdx];
    }
    else if (pionCharge=="pi-") {
        pi  = &piminus[piIdx];
        Vpi = &Vpiminus[piIdx];
        Zpi = Zpims   [piIdx];
    }
    else {
        std::cout << "pion charge ill defined at Stream_e_pi_line_to_CSV(), returning " << std::endl;
        return;
    }
    if (fdebug>3) {
        std::cout
        << "Stream_e_pi_n_line_to_CSV(" << piIdx  << "," << passed_cuts_e_pi_kinematics << "," << passed_cuts_n << ")"
        << std::endl
        << "piIdx: "    << piIdx << ","
        << "passed_cuts_e_pi_kinematics: " << passed_cuts_e_pi_kinematics << ","
        << "passed_cuts_n: " << passed_cuts_n << ","
        << std::endl;
        std::cout
        << "*pi.X(): "  << pi->X()      << ","
        << "Pe: "       << e->P()       << " GeV/c,"
        << std::endl;
    }
    TLorentzVector pion = *pi;

    // write a (e,e'pi) event-line to CSV file

    // from Harut A., Aug-2, 2021:
    //    1: status, 2: runnum, 3: evnum, 4: helicity, 5: e_p, 6: e_theta, 7: e_phi, 8: vz_e, 9: pi_p, 10: pi_theta, 11: pi_phi, 12: vz_pi, 13: P_p, 14: P_theta, 15: P_phi, 16: vz_P, 17: Q2, 18: W, 19: x, 20: y, 21: z_pi, 22: z_P, 23: Mx(e:pi:P:X), 24: Mx(e:pi:X), 25: Mx(e:P:X), 26: zeta, 27: Mh, 28: pT_pi, 29: pT_P, 30: pTpT, 31: xF_pi, 32: xF_P, 33: eta_pi, 34: eta_P, 35: Delta_eta, 36: phi_pi (gamma*N COM), 37: phi_P (gamma*N COM), 38: Delta_phi.
    //
    //    1: status is a number indicating the quality of the event with the non-0 number indicating something was not not good (ex. out of fiducial region, not within the final cuts on energies of particles, missing or invarian masses....) The final observables will be done using status==0, while sensitivity of the observable to different cuts could be studied for various values of status>0
    //    2-3: run number and event number to identify the event
    //    4: helicity of the electron +1 along the beam and -1 opposite to it
    //    5,6,7,8 electron momentum,theta,phi_Lab, and z-vertex
    //    9,10,11,12 the same for pi+
    //    All other columns could be calculated from the first 16, but are included for cross check and minimizing the work in production of final observables. Some of them simple, like x,Q^2,W,y,  z=E\pion/nu, zome less trivial like Breit frame rapidities of pion (eta_pi) and proton eta_P, or corresponding values for X_Feynman variable xF_pi, xF_P.
    //    Some, like azimuthal angles of pion (phi_pi) and proton (phi_p) in the CM frame may also be confusing.
    //    In addition to the first 16 columns (mandatory) you can add as many columns as you are comfortable to fill, and consider relevant for your process.

    // ------------------------------------------------------------------------------------------------
    // compute kinematics that also relies on pion information
    // ------------------------------------------------------------------------------------------------
    xF          = 2. * (pion.Dot(*q)) / (q->Mag() * W);
    M_X_ee_pi   = ( (*Beam + *target) - (*e + *pi) ).Mag(); // missing mass of (e,e'pi)
    M_X_ee_pi_n = ( (*Beam + *target) - (*e + *pi + Pn) ).Mag(); // missing mass of (e,e'pi n)
    status      = 0;

    
    // now stream data to CSV file
    std::vector<double> variables =
    {   (double)status, (double)SIDISrunID,   (double)eventnumber,    (double)beam_helicity,
        e->P(),         e->Theta(),         e->Phi(),             Ve->Z(),
        pi->P(),        pi->Theta(),        pi->Phi(),            Vpi->Z(),
        Pn.P(),         Pn.Theta(),         Pn.Phi(),             Ve->Z(),
        Q2,             W,                  xB,                   Zpi,
        omega,          xF,                 y,
        M_X_ee_pi,      M_X_ee_pi_n,        xPrime1,              xPrime2,
        theta_sq,
        WPrime,
    };
    StreamToCSVfile(variables);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void InitializeVariables(){
    Pn_Vect     = TVector3();
    Band_e_Vect = TVector3();
    e           = new TLorentzVector( 0, 0, 0, Me );
    Band_data_e = TLorentzVector( 0, 0, 0, Me );
    
    xB          = Q2        = omega     = -9999;
    xF          = y                     = -9999;
    M_X_ee_pi   = M_X_ee_pi_n           = -9999;
    Ve                                  = new TVector3();
    piplus      .clear();
    piminus     .clear();
    Vpiplus     .clear();
    Vpiminus    .clear();
    for (int piIdx=0; piIdx<NMAXPIONS; piIdx++) {
        piplus  .push_back( TLorentzVector(0,0,0,Mpips) );
        Vpiplus .push_back( TVector3() );
        eepipsPastKinematicalCuts[piIdx]            = false;
        piminus .push_back( TLorentzVector(0,0,0,Mpims) );
        Vpiminus.push_back( TVector3() );
        pimsPastSelectionCuts[piIdx]                = false;
        eepimsPastKinematicalCuts[piIdx]            = false;
        piplus_Px[piIdx]    = piplus_Py[piIdx]  = piplus_Pz[piIdx]  = piplus_E[piIdx]   = -9999;
        piminus_Px[piIdx]   = piminus_Py[piIdx] = piminus_Pz[piIdx] = piminus_E[piIdx]  = -9999;
        Vpiplus_X[piIdx]    = Vpiplus_Y[piIdx]  = Vpiplus_Z[piIdx]  = -9999;
        Vpiminus_X[piIdx]   = Vpiminus_Y[piIdx] = Vpiminus_Z[piIdx] = -9999;
         
    }
    status                                          = 1; // 0 is good...
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MergeEventData(){
    if (fdebug>2){
        std::cout
        << "MergeEventData()"
        << std::endl
        << "Merging event " << BANDeventID << " from run " << BANDrunID
        << std::endl        << std::endl;
    }

    // fill output TTree and CSV file
    MergedTree  -> Fill();
    
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
    if (fdebug>3) { std::cout << "------------------------------------" << std::endl; }
}
