// last edit Sep-19, 2021

#ifdef __CINT__
#pragma link C++ class std::vector<TVector3>+;
#pragma link C++ class std::vector<TLorentzVector>+;
#endif

#include <cstdlib>
#include <iostream>
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
#include<time.h>

#include "/u/home/cohen/SIDIS_at_BAND/Auxiliary/bank.h"
#include "/u/home/cohen/SIDIS_at_BAND/Auxiliary/BBand.h"
#include "/u/home/cohen/SIDIS_at_BAND/Auxiliary/BEvent.h"
#include "/u/home/cohen/SIDIS_at_BAND/Auxiliary/constants.h"

#define NMAXEVENTS 5000000
#define NMAXPIONS 20 // maximal allowed number of pions
#define r2d 180./3.1415 // radians to degrees

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
// Globals
// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
TString DataPath = "/volatile/clas12/users/ecohen/BAND/";



// Input root files and trees
TFile * SIDISFile, * BANDFile;
TTree * SIDISTree, * BANDTree;

// Output root file and tree
TFile * MergedFile;
TTree * MergedTree;

// Output CSV file
std::ofstream   CSVfile;

// time
clock_t tStart = clock();

// globals
TString                         pionCharge; // "pi+" or "pi-"
Int_t               BANDrunID, BANDeventID;
Int_t             SIDISrunID, SIDISeventID;
Int_t          EventIDsToMerge[NMAXEVENTS];
Int_t  BANDEventIndicesToMerge[NMAXEVENTS];
Int_t SIDISEventIndicesToMerge[NMAXEVENTS];
Int_t             BANDeventIDs[NMAXEVENTS];
Int_t            SIDISeventIDs[NMAXEVENTS];
Int_t                 Npions, Npips, Npims;
Int_t                  Ne, Np, Nn, Ngammas;

// kinematics and observables
Double_t Mp      = 0.938272; // GeV/c2
Double_t Mn      = 0.939565; // GeV/c2
Double_t Md      = 1.875;    // GeV/c2
Double_t Mp2     = Mp * Mp;

// SIDIS Tree
TLorentzVector     *Beam=0;
TLorentzVector        *e=0;
TLorentzVector        *q=0;
TLorentzVector       *Pn=0; // neutron momentum
// reconstructed vertex position
TVector3             *Ve=0;
std::vector<TVector3>          *Vpiplus=0;
std::vector<TVector3>         *Vpiminus=0;
std::vector<TLorentzVector>     *piplus=0;
std::vector<TLorentzVector>    *piminus=0;

bool                     eepiPastCutsInEvent;
bool                     goodneutron = false;
bool    eepipsPastKinematicalCuts[NMAXPIONS];
bool    eepimsPastKinematicalCuts[NMAXPIONS];

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
Double_t                        xB;
Double_t                        Q2;
Double_t                     omega;
Double_t                        w2; // omega^2
Double_t          Zpips[NMAXPIONS]; // energy fraction rest frame
Double_t                         W; // energy of the hadronic system
Double_t                   alpha_s; // light cone fraction of momentum of the recoil neutron
Double_t                    WPrime; // moving proton
Double_t                    xPrime; // moving proton
Double_t                        Es; // spectator energy
Double_t                        Ps; // spectator momentum
Double_t                  theta_sq; // spectator angle with respect to momentum transfer

TClonesArray   * nHits = new TClonesArray("bandhit");
TClonesArray  &saveHit = *nHits;
TClonesArray * mcParts = new TClonesArray("genpart");
TClonesArray   &saveMC = *mcParts;

void             OpenInputFiles (TString RunStr);
void            OpenOutputFiles (TString RunStr);
void            StreamToCSVfile (std::vector<Double_t> observables, int fdebug=0);
void            CloseInputFiles ();
void           CloseOutputFiles (TString OutDataPath);
void    MergeSIDISandBANDevents (int NeventsToMerge=10,
                                 int fdebug=2,
                                 int PrintProgress=5000);
Int_t CreateListOfEventsToMerge (TTree * BANDTree,
                                 TTree * SIDISTree,
                                 int NeventsToMerge=-1,
                                 int fdebug=0);

void Stream_e_pi_n_line_to_CSV (int piIdx,
                                bool passed_cuts_e_pi_kinematics,
                                bool passed_cuts_n,
                                int fdebug );
void           StreamToCSVfile (TString pionCharge, // "pi+" or "pi-"
                                std::vector<Double_t> observables,
                                bool passed_cuts_e_pi_kinematics,
                                bool passed_cuts_n,
                                int fdebug);

void    SetInputAndOutputTTrees ();
void          ComputeKinematics ();
void                  PrintDone ();
void          PrintMonitorHello ();
void              SetPionCharge ( TString fpionCharge ){ pionCharge = fpionCharge; };

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
// Main functionality
// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void MergeSIDISandBANDSkimmers(int RunNumber=6420,
                               TString fpionCharge="pi+", // "pi+" or "pi-"
                               int NeventsToMerge=-1,
                               int fdebug=1,
                               int PrintProgress=5000){
    
    SetPionCharge    ( fpionCharge );
    char RunNumberStr[20];
    sprintf( RunNumberStr, "00%d", RunNumber );
    OpenInputFiles   ( (TString)RunNumberStr );
    OpenOutputFiles  ( (TString)RunNumberStr );
    MergeSIDISandBANDevents( NeventsToMerge, fdebug, PrintProgress );
    CloseOutputFiles (DataPath + "merged_SIDIS_and_BAND_skimming/");
    CloseInputFiles  ();
    PrintDone();
}



// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void MergeSIDISandBANDevents (int NeventsToMerge=10,
                              int fdebug=2,
                              int PrintProgress=5000){
        
    Int_t   NeventsBAND  = BANDTree->GetEntries();
    Int_t   NeventsSIDIS = SIDISTree->GetEntries();
    // Create a list of events to merge
    // this takes the most resources, and the largest amount of time.
    // typically, per 1 merged event, it takes about 14-40 ms
    // and we typically merge 1M events = 1e4 sec
    if (fdebug>1) {
        std::cout << "Create a list of events to merge" << std::endl
        << "stepping over "
        << NeventsBAND << " BAND and "
        << NeventsSIDIS << " SIDIS events"
        << std::endl
        << "Take some coffee, this takes some 20-40 ms per event to merge."
        << std::endl;
    }
    Int_t Nevents2Merge = CreateListOfEventsToMerge(BANDTree, SIDISTree, NeventsToMerge, fdebug);
    if (fdebug>1) {
        std::cout << std::endl;
        std::cout << "Done creating merge list of " << Nevents2Merge << " events, after "
        << double(clock() - tStart) / (double)CLOCKS_PER_SEC
        << " sec "<< std::endl;
        
        if (fdebug>2) {
            std::cout << "Merging events (BAND/SIDIS):"  << std::endl;
            for (int i=0; i<Nevents2Merge; i++) {
                auto         eventID = EventIDsToMerge[i];
                auto  BANDEventIndex = BANDEventIndicesToMerge[i];
                auto SIDISEventIndex = SIDISEventIndicesToMerge[i];
                std::cout << eventID << "(" << BANDEventIndex << "/" << SIDISEventIndex << "),";
            }
        }
    }
    
    SetInputAndOutputTTrees ();
    
    for (int MergedEvtId=0; MergedEvtId<Nevents2Merge; MergedEvtId++) {
        
        BANDTree -> GetEntry( BANDEventIndicesToMerge[MergedEvtId] );
                
        if (eepiPastCutsInEvent) {
            
            // compute kinematical variables
            ComputeKinematics ();
            
            // grab neturon information from BAND
            SIDISTree   -> GetEntry( SIDISEventIndicesToMerge[MergedEvtId] );
            
            // fill output TTree and CSV file
            MergedTree  -> Fill();
            
            if (pionCharge=="pi+") {
                Npions = Npips;
            } else if (pionCharge=="pi-") {
                Npions = Npims;
            }
            
            for (int piIdx=0; piIdx<Npions; piIdx++) {
                bool eepiPastKinematicalCuts = false;
                if (pionCharge=="pi+") {
                    eepiPastKinematicalCuts = eepipsPastKinematicalCuts[piIdx];
                } else if (pionCharge=="pi-") {
                    eepiPastKinematicalCuts = eepimsPastKinematicalCuts[piIdx];
                }
                
                
                Stream_e_pi_n_line_to_CSV(piIdx,
                                          eepiPastKinematicalCuts, goodneutron,
                                          fdebug );
            }
            
            
            if (fdebug>2){
                std::cout
                << "merging event " << BANDeventID << " from run " << BANDrunID
                << std::endl;
            }
        }
    } // end merged event loop
    
    
    if (fdebug>2){
        std::cout << "merged " << Nevents2Merge << " SIDIS and BAND events." << std::endl;
        PrintMonitorHello();
    }
}



// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void OpenInputFiles (TString RunStr){
    
    std::cout << "Opening " << DataPath + "neutron_skimming/"
    + "skimmed_neutrons_inc_"  + RunStr + ".root" << std::endl;
    
    TString   skimmedBANDFilename = (DataPath + "neutron_skimming/"
                                     + "skimmed_neutrons_inc_"  + RunStr + ".root");
    BANDFile                      = new TFile( skimmedBANDFilename );
    BANDTree                      = (TTree*)BANDFile->Get("neutrons");
    
    
    
    std::cout << "Opening " << DataPath + "SIDIS_skimming/"
    + "skimmed_SIDIS_inc_"  + RunStr + ".root" << std::endl;
    
    TString  skimmedSIDISFilename = (DataPath + "SIDIS_skimming/"
                                     + "skimmed_SIDIS_inc_"  + RunStr + ".root");
    SIDISFile                     = new TFile( skimmedSIDISFilename );
    SIDISTree                     = (TTree*)SIDISFile->Get("sidis");
    
    
}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void OpenOutputFiles (TString RunStr){
    TString csvheader = ("runID,eventID,"
                         +(TString)"livetime,current,"
                         +(TString)"xB,Q2,"
                         +(TString)"Ebeam,z,"
                         +(TString)"W,alpha_s,"
                         +(TString)"WPrime,xPrime,",
                         +(TString)"e_Px,e_Py,"
                         +(TString)"e_Pz,e_E,"
                         +(TString)"q_Px,q_Py,"
                         +(TString)"q_Pz,q_E,"
                         +(TString)"Ve_z,Vpiplus_z,"
                         +(TString)"goodneutron,");
    
    TString skimmedMergedFilename = (DataPath + "merged_SIDIS_and_BAND_skimming/"
                                     + "skimmed_SIDIS_and_BAND_inc_"  + RunStr );
    
    
    
    // Create output tree
    MergedFile = new TFile( skimmedMergedFilename + ".root" ,"RECREATE");
    MergedTree = new TTree( "T" , "Event information from merged SIDIS and BAND skimmers");
    
    // Create output csv files
    CSVfile.open( skimmedMergedFilename + ".csv" );
    CSVfile << csvheader << std::endl;
    
}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void CloseInputFiles (){
    SIDISFile->Close();
    BANDFile->Close();
}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void CloseOutputFiles (TString OutDataPath){
    
    // close output CSV
    CSVfile.close();
    
    int Nevents = MergedTree->GetEntries();
    
    // close output ROOT
    MergedFile->cd();
    MergedTree->Write();
    MergedFile->Close();
    
    std::cout << "output file ready in " << std::endl << OutDataPath << std::endl
    << "merged " << Nevents << " SIDIS and BAND (neutron) events" << std::endl;
}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void StreamToCSVfile (std::vector<Double_t> observables, int fdebug){
    if (fdebug>2) std::cout << "StreamToCSVfile()" << std::endl;
    for (auto v:observables) {
        CSVfile << v << ",";
    }
    CSVfile << std::endl;
}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
Int_t CreateListOfEventsToMerge(TTree * BANDTree,
                                TTree * SIDISTree,
                                int NeventsToMerge,
                                int fdebug){
    
    // Decide which event-indices to merge from the two TTrees
    // Since this functionallity takes the most time,
    // we need to be very creative to make it less time-consuming
    if (fdebug>2) {
        std::cout << "CreateListOfEventsToMerge()" << std::endl;
    }
    
    // first define two vectors that containt the event IDs in each TTree
    // BAND
    Int_t     fBANDeventID, NeventsBAND = BANDTree->GetEntries();
    BANDTree   -> SetBranchAddress("eventnumber"  ,&fBANDeventID);
    for (int BANDevent=0; BANDevent < NeventsBAND ; BANDevent++){
        BANDTree -> GetEntry(BANDevent);
        BANDeventIDs[BANDevent] = fBANDeventID;
    }
    Int_t     fSIDISeventID, NeventsSIDIS = SIDISTree->GetEntries();
    SIDISTree  -> SetBranchAddress("eventnumber"  ,&fSIDISeventID);
    for (int SIDISevent=0; SIDISevent < NeventsSIDIS ; SIDISevent++){
        SIDISTree -> GetEntry(SIDISevent);
        SIDISeventIDs[SIDISevent] = fSIDISeventID;
    }
    
    // now, we merge the events
    Int_t         NmergedEvents = 0;
    int      SIDISeventIndexMin = 0;
    SIDISEventIndicesToMerge[0] = 0; // initialize for loop efficiency
    for (int BANDeventIndex=0;
         BANDeventIndex < NeventsBAND ;
         BANDeventIndex++){
        
        
        for (int SIDISeventIndex=SIDISEventIndicesToMerge[NmergedEvents];
             SIDISeventIndex < NeventsSIDIS ;
             SIDISeventIndex++){
            
            if ( SIDISeventIDs[SIDISeventIndex] > BANDeventIDs[BANDeventIndex] ) {
                break;
            }
            
            if ( BANDeventIDs[BANDeventIndex] == SIDISeventIDs[SIDISeventIndex] ){
                
                EventIDsToMerge[NmergedEvents] = BANDeventID;
                BANDEventIndicesToMerge[NmergedEvents] = BANDeventIndex;
                SIDISEventIndicesToMerge[NmergedEvents] = SIDISeventIndex;
                NmergedEvents ++ ;
                
                if ((NeventsToMerge>0) && (NmergedEvents >= NeventsToMerge)){
                    return NmergedEvents;
                }
            }
            
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
    
    // SIDIS input tree branches
    
    SIDISTree  -> SetBranchAddress("eventnumber"        ,&SIDISeventID                  );
    SIDISTree  -> SetBranchAddress("runnum"             ,&SIDISrunID                    );
    SIDISTree  -> SetBranchAddress("e"                          ,&e                     );
    SIDISTree  -> SetBranchAddress("pi"                         ,&piplus                );
    SIDISTree  -> SetBranchAddress("Ve"                         ,&Ve                    );
    SIDISTree  -> SetBranchAddress("Vpi"                        ,&Vpiplus               );
    SIDISTree  -> SetBranchAddress("Beam"                       ,&Beam                  );
    SIDISTree  -> SetBranchAddress("beam_helicity"              ,&beam_helicity         );
    SIDISTree  -> SetBranchAddress("q"                          ,&q                     );
    SIDISTree  -> SetBranchAddress("xB"                        ,&xB                     );
    SIDISTree  -> SetBranchAddress("Q2"                        ,&Q2                     );
    SIDISTree  -> SetBranchAddress("omega"                     ,&omega                  );
    SIDISTree  -> SetBranchAddress("Z"                         ,&Zpips                  );
    SIDISTree  -> SetBranchAddress("W"                         ,&W                      );
    SIDISTree  -> SetBranchAddress("eepiPastCutsInEvent"       ,&eepiPastCutsInEvent    );
    SIDISTree  -> SetBranchAddress("eepipsPastKinematicalCuts" ,&eepipsPastKinematicalCuts);
    SIDISTree  -> SetBranchAddress("eepimsPastKinematicalCuts" ,&eepimsPastKinematicalCuts);
    SIDISTree  -> SetBranchAddress("Npips"                     ,&Npips                 );
    SIDISTree  -> SetBranchAddress("Npims"                     ,&Npims                 );
    SIDISTree  -> SetBranchAddress("Nelectrons"                ,&Ne                    );
    SIDISTree  -> SetBranchAddress("Ngammas"                   ,&Ngammas               );
    SIDISTree  -> SetBranchAddress("Nprotons"                  ,&Np                    );
    SIDISTree  -> SetBranchAddress("Nneutrons"                 ,&Nn                    );

    
    // BAND input tree branches
    BANDTree   -> SetBranchAddress("eventnumber"  ,&BANDeventID);
    BANDTree   -> SetBranchAddress("Runno"        ,&BANDrunID);
    BANDTree   -> SetBranchAddress("Ebeam"        ,&Ebeam);
    BANDTree   -> SetBranchAddress("gated_charge" ,&gated_charge);
    BANDTree   -> SetBranchAddress("livetime"     ,&livetime);
    BANDTree   -> SetBranchAddress("starttime"    ,&starttime);
    BANDTree   -> SetBranchAddress("current"      ,&current);
    BANDTree   -> SetBranchAddress("weight"       ,&weight);
    BANDTree   -> SetBranchAddress("nMult"        ,&nMult);
    BANDTree   -> SetBranchAddress("nHits"        ,&nHits);
    BANDTree   -> SetBranchAddress("goodneutron"  ,&goodneutron);
    BANDTree   -> SetBranchAddress("nleadindex"   ,&nleadindex);
    BANDTree   -> SetBranchAddress("genMult"      ,&genMult);
    BANDTree   -> SetBranchAddress("mcParts"      ,&mcParts);
    
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
    MergedTree->Branch("nHits"                  ,&nHits                 );
    MergedTree->Branch("goodneutron"            ,&goodneutron           );
    MergedTree->Branch("nleadindex"             ,&nleadindex            );
    MergedTree->Branch("e"                      ,&e                     );
    MergedTree->Branch("piplus"                 ,&piplus                );
    MergedTree->Branch("Ve"                     ,&Ve                    );
    MergedTree->Branch("Vpiplus"                ,&Vpiplus               );
    MergedTree->Branch("Beam"                   ,&Beam                  );
    MergedTree->Branch("q"                      ,&q                     );
    MergedTree->Branch("xB"                     ,&xB                    );
    MergedTree->Branch("Q2"                     ,&Q2                    );
    MergedTree->Branch("omega"                  ,&omega                 );
    MergedTree->Branch("Z"                      ,&Zpips                 );
    MergedTree->Branch("eepiPastCutsInEvent"    ,&eepiPastCutsInEvent );
    
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
void ComputeKinematics(){
    // compute kinematics
    // SIDISc12rSkimmer.C already computes few of the kinematical variables:
    //     xB,  Q2, omega, W, Z
    Es      = Pn->E();
    Ps      = Pn->P();
    w2      = omega * omega;
    theta_sq= Pn->Angle( q );
    xPrime  = Q2 / (2. * ((Md - Es) * omega + Pn->Vect()*q->Vect() ));
    W       = sqrt(Mp2 - Q2 + 2. * omega * Mp);
    WPrime  = sqrt(Mp2 - Q2 + 2. * omega * (Md - Es) + 2. * Ps * sqrt(Q2 + w2) * cos( theta_sq ));
    // move to q-frame
    // compute light-cone fraction of momentum
    // alpha_s = ComputeLightConeFraction( Pn_qFrame );
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void PrintDone(){
    std::cout << "Done. " << "Execution time: "
    << double(clock() - tStart) / (double)CLOCKS_PER_SEC
    << " sec "<< std::endl;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Stream_e_pi_n_line_to_CSV(int piIdx,
                               bool passed_cuts_e_pi_kinematics,
                               bool passed_cuts_n,
                               int fdebug ){
    TLorentzVector  pi;
    TVector3        Vpi;
    double          Zpi;
    if (pionCharge=="pi+") {
        pi  = piplus [piIdx];
        Vpi = Vpiplus[piIdx];
        Zpi = Zpips  [piIdx];
    }
    else if (pionCharge=="pi-") {
        pi  = piminus [piIdx];
        Vpi = Vpiminus[piIdx];
        Zpi = Zpims   [piIdx];
   }
    else {
        std::cout << "pion charge ill defined at Stream_e_pi_line_to_CSV(), returning " << std::endl;
        return;
    }
    
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
    xF  = 2. * (pi.Dot(q)) / (q.Mag() * W);
    M_X = ( Beam + target - e - pi - Pn ).Mag();
    
    // now stream data to CSV file
    std::vector<double> variables =
    {   (double)status, (double)runnum,     (double)evnum,      (double)beam_helicity,
        e.P(),          e.Theta(),          e.Phi(),            Ve.Z(),
        pi.P(),         pi.Theta(),         pi.Phi(),           Vpi.Z(),
        Pn.P(),         Pn.Theta(),         Pn.Phi(),           Vn.Z(),
        Q2,             W,                  xB,                 Zpi,
        omega,          xF,                 y,                  M_X,
    };
    StreamToCSVfile(pionCharge,
                    variables ,
                    passed_cuts_e_pi_kinematics, passed_cuts_n,
                    fdebug );
}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void StreamToCSVfile (std::vector<Double_t> observables,
                      bool passed_cuts_e_pi_kinematics,
                      bool passed_cuts_n,
                      int fdebug){
    if (fdebug>1) {
        std::cout << "streaming to CSVfile" << std::endl;
    }
    for (auto v:observables) CSVfile_e_pi_n << v << ",";
    CSVfile_e_pi_n << std::endl;
}
