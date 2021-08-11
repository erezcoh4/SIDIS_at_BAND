// last edit July-8, 2021 (EOC, mbp), see README

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
Int_t          EventIDsToMerge[NMAXEVENTS];
Int_t  BANDEventIndicesToMerge[NMAXEVENTS];
Int_t SIDISEventIndicesToMerge[NMAXEVENTS];
Int_t             BANDeventIDs[NMAXEVENTS];
Int_t            SIDISeventIDs[NMAXEVENTS];

// kinematics and observables
Double_t Mp      = 0.938272; // GeV/c2
Double_t Mn      = 0.939565; // GeV/c2
Double_t Md      = 1.875;    // GeV/c2
Double_t Mp2     = Mp * Mp;

// SIDIS Tree
TLorentzVector        *e=0;
TLorentzVector   *piplus=0;
TLorentzVector     *Beam=0;
TLorentzVector        *q=0;
// reconstructed vertex position
TVector3             *Ve=0;
TVector3        *Vpiplus=0;

// kinematics
Double_t             xB;
Double_t             Q2;
Double_t          omega;
Double_t              z; // energy fraction rest frame
Double_t              W; // energy of the hadronic system
Double_t        alpha_s; // light cone fraction of momentum of the recoil neutron

// moving proton
Double_t         WPrime;
Double_t         xPrime;



void             OpenInputFiles (TString RunStr);
void            OpenOutputFiles (TString RunStr, TString header);
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





// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
// Main functionality
// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void MergeSIDISandBANDSkimmers(int RunNumber=6420,
                               int NeventsToMerge=-1,
                               int fdebug=1,
                               int PrintProgress=5000){
    
    char RunNumberStr[20];
    sprintf( RunNumberStr, "00%d", RunNumber );
    OpenInputFiles   ( (TString)RunNumberStr );
    OpenOutputFiles  ( (TString)RunNumberStr,
                      ("runID,eventID,"
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
                       +(TString)"goodneutron,"));
    
    MergeSIDISandBANDevents( NeventsToMerge, fdebug, PrintProgress );
    
    CloseOutputFiles (DataPath + "merged_SIDIS_and_BAND_skimming/");
    CloseInputFiles  ();
    
    std::cout << "execution time: "
    << double(clock() - tStart) / (double)CLOCKS_PER_SEC
    << " sec "<< std::endl;
}



// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void    MergeSIDISandBANDevents (int NeventsToMerge=10, int fdebug=2, int PrintProgress=5000){

    Int_t   BANDrunID, BANDeventID, SIDISrunID, SIDISeventID;
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
    Int_t Nevents2Merge = CreateListOfEventsToMerge(BANDTree,
                                                    SIDISTree,
                                                    NeventsToMerge,
                                                    fdebug);
    if (fdebug>0) {
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




    

    bool      ePastSelectionCuts = false;
    bool piplusPastSelectionCuts = false;

    int     DC_layer;

    double  pips_PCAL_W, pips_PCAL_V, pips_PCAL_x, pips_PCAL_y,       pips_PCAL_z;
    double  pips_chi2PID;
    double  pips_PCAL_sector, pips_DC_sector, pips_Chi2N;
    double  pips_DC_x[3], pips_DC_y[3];

    double  e_E_PCAL, e_E_ECIN, e_E_ECOUT; // electron energy deposit in ECAL_out [GeV]
    double  e_PCAL_W,    e_PCAL_V;
    double  e_PCAL_x,    e_PCAL_y,  e_PCAL_z;
    double  e_PCAL_sector;
    double  e_DC_sector,            e_DC_Chi2N;
    double  e_DC_x[3],   e_DC_y[3];

    double  pips_E_PCAL, pips_E_ECIN,       pips_E_ECOUT;
    SIDISTree  -> SetBranchAddress("eventnumber"        ,&SIDISeventID                  );
    SIDISTree  -> SetBranchAddress("runnum"             ,&SIDISrunID                    );

    // output tree branches
    SIDISTree  -> SetBranchAddress("e_E_PCAL"          ,&e_E_PCAL              );
    SIDISTree  -> SetBranchAddress("e_E_ECIN"          ,&e_E_ECIN              );
    SIDISTree  -> SetBranchAddress("e_E_ECOUT"         ,&e_E_ECOUT             );
    SIDISTree  -> SetBranchAddress("pips_chi2PID"      ,&pips_chi2PID          );

    SIDISTree  -> SetBranchAddress("e_PCAL_W"          ,&e_PCAL_W              );
    SIDISTree  -> SetBranchAddress("e_PCAL_V"          ,&e_PCAL_V              );
    SIDISTree  -> SetBranchAddress("pips_PCAL_x"       ,&pips_PCAL_x           );
    SIDISTree  -> SetBranchAddress("pips_PCAL_y"       ,&pips_PCAL_y           );
    SIDISTree  -> SetBranchAddress("pips_PCAL_z"       ,&pips_PCAL_z           );
    SIDISTree  -> SetBranchAddress("e_PCAL_x"          ,&e_PCAL_x              );
    SIDISTree  -> SetBranchAddress("e_PCAL_y"          ,&e_PCAL_y              );
    SIDISTree  -> SetBranchAddress("e_PCAL_z"          ,&e_PCAL_z              );

    SIDISTree  -> SetBranchAddress("e_PCAL_sector"     ,&e_PCAL_sector                  );
    SIDISTree  -> SetBranchAddress("e_DC_sector"       ,&e_DC_sector                    );
    SIDISTree  -> SetBranchAddress("e_DC_Chi2N"        ,&e_DC_Chi2N                     );
    //    SIDISTree  -> SetBranchAddress("e_DC_x"            ,&e_DC_x                );
    //    SIDISTree  -> SetBranchAddress("e_DC_y"            ,&e_DC_y                );

    SIDISTree  -> SetBranchAddress("pips_PCAL_sector"          ,&pips_PCAL_sector       );
    SIDISTree  -> SetBranchAddress("pips_DC_sector"            ,&pips_DC_sector         );
    SIDISTree  -> SetBranchAddress("pips_Chi2N"                ,&pips_Chi2N             );
    //    SIDISTree  -> SetBranchAddress("pips_DC_x"                 ,&pips_DC_x             );
    //    SIDISTree  -> SetBranchAddress("pips_DC_y"                 ,&pips_DC_y             );

    SIDISTree  -> SetBranchAddress("pips_E_PCAL"               ,&pips_E_PCAL            );
    SIDISTree  -> SetBranchAddress("pips_E_ECIN"               ,&pips_E_ECIN            );

    SIDISTree  -> SetBranchAddress("pips_E_ECIN"               ,&pips_E_ECIN            );
    SIDISTree  -> SetBranchAddress("pips_E_ECOUT"              ,&pips_E_ECOUT           );
    SIDISTree  -> SetBranchAddress("DC_layer"                  ,&DC_layer               );

    SIDISTree  -> SetBranchAddress("e"                         ,&e                      );
    SIDISTree  -> SetBranchAddress("piplus"                    ,&piplus                 );
    SIDISTree  -> SetBranchAddress("Ve"                        ,&Ve                     );
    SIDISTree  -> SetBranchAddress("Vpiplus"                   ,&Vpiplus                );
    SIDISTree  -> SetBranchAddress("Beam"                      ,&Beam                   );
    SIDISTree  -> SetBranchAddress("q"                         ,&q                      );

    SIDISTree  -> SetBranchAddress("xB"                        ,&xB                     );
    SIDISTree  -> SetBranchAddress("Q2"                        ,&Q2                     );
    SIDISTree  -> SetBranchAddress("omega"                     ,&omega                  );
    SIDISTree  -> SetBranchAddress("z"                         ,&z                      );
    SIDISTree  -> SetBranchAddress("xPrime"                    ,&xPrime                 );
    SIDISTree  -> SetBranchAddress("W"                         ,&W                      );
    SIDISTree  -> SetBranchAddress("WPrime"                    ,&WPrime                 );
    SIDISTree  -> SetBranchAddress("alpha_s"                   ,&alpha_s                );
    
    
    
    SIDISTree  -> SetBranchAddress("ePastSelectionCuts"        ,&ePastSelectionCuts    );
    SIDISTree  -> SetBranchAddress("piplusPastSelectionCuts"   ,&piplusPastSelectionCuts);






    // BAND Tree
    double           Ebeam = 0;
    double    gated_charge = 0;
    double        livetime = 0;
    double       starttime = 0;
    double         current = 0;
    int        eventnumber = 0;
    bool       goodneutron = false;
    int         nleadindex = -1;
    double          weight = 0;
    //     Neutron info:
    int              nMult = 0;
    TClonesArray   * nHits = new TClonesArray("bandhit");
    TClonesArray  &saveHit = *nHits;
    //    MC info:
    int            genMult = 0;
    TClonesArray * mcParts = new TClonesArray("genpart");
    TClonesArray   &saveMC = *mcParts;



    BANDTree   -> SetBranchAddress("eventnumber"  ,&BANDeventID);
    BANDTree   -> SetBranchAddress("Runno"        ,&BANDrunID);

    BANDTree   -> SetBranchAddress("Ebeam"        ,&Ebeam);
    BANDTree   -> SetBranchAddress("gated_charge" ,&gated_charge);
    BANDTree   -> SetBranchAddress("livetime"     ,&livetime);
    BANDTree   -> SetBranchAddress("starttime"    ,&starttime);
    BANDTree   -> SetBranchAddress("current"      ,&current);
    BANDTree   -> SetBranchAddress("weight"       ,&weight);
    //    Neutron branches:
    BANDTree   -> SetBranchAddress("nMult"        ,&nMult);
    BANDTree   -> SetBranchAddress("nHits"        ,&nHits);
    //Branches to store if good Neutron event and leadindex
    BANDTree   -> SetBranchAddress("goodneutron"  ,&goodneutron);
    BANDTree   -> SetBranchAddress("nleadindex"   ,&nleadindex);
    //    MC branches:
    BANDTree   -> SetBranchAddress("genMult"      ,&genMult);
    BANDTree   -> SetBranchAddress("mcParts"      ,&mcParts);

    // Merged Tree - containing all variables...
    // run and event number (ID) have to be consistent in the merged tree,
    // so it does not matter from where we take them...
    MergedTree->Branch("Runno"              ,&SIDISrunID    );
    MergedTree->Branch("eventnumber"        ,&SIDISeventID  );
    // BAND branches
    MergedTree->Branch("Ebeam"              ,&Ebeam         );
    MergedTree->Branch("gated_charge"       ,&gated_charge  );
    MergedTree->Branch("livetime"           ,&livetime      );
    MergedTree->Branch("starttime"          ,&starttime     );
    MergedTree->Branch("current"            ,&current       );
    MergedTree->Branch("weight"             ,&weight        );
    MergedTree->Branch("nMult"              ,&nMult         );
    MergedTree->Branch("nHits"              ,&nHits         );
    MergedTree->Branch("goodneutron"        ,&goodneutron   );
    MergedTree->Branch("nleadindex"         ,&nleadindex    );
    // SIDIS branches
    MergedTree->Branch("e_E_PCAL"           ,&e_E_PCAL              );
    MergedTree->Branch("e_E_ECIN"           ,&e_E_ECIN              );
    MergedTree->Branch("e_E_ECOUT"          ,&e_E_ECOUT             );
    MergedTree->Branch("pips_chi2PID"       ,&pips_chi2PID          );

    MergedTree->Branch("e_PCAL_W"           ,&e_PCAL_W              );
    MergedTree->Branch("e_PCAL_V"           ,&e_PCAL_V              );
    MergedTree->Branch("pips_PCAL_x"        ,&pips_PCAL_x           );
    MergedTree->Branch("pips_PCAL_y"        ,&pips_PCAL_y           );
    MergedTree->Branch("pips_PCAL_z"        ,&pips_PCAL_z           );
    MergedTree->Branch("e_PCAL_x"           ,&e_PCAL_x              );
    MergedTree->Branch("e_PCAL_y"           ,&e_PCAL_y              );
    MergedTree->Branch("e_PCAL_z"           ,&e_PCAL_z              );

    MergedTree->Branch("e_PCAL_sector"      ,&e_PCAL_sector         );
    MergedTree->Branch("e_DC_sector"        ,&e_DC_sector           );
    MergedTree->Branch("e_DC_Chi2N"         ,&e_DC_Chi2N            );
    //    MergedTree->Branch("e_DC_x"             ,&e_DC_x                );
    //    MergedTree->Branch("e_DC_y"             ,&e_DC_y                );

    MergedTree->Branch("pips_PCAL_sector"   ,&pips_PCAL_sector      );
    MergedTree->Branch("pips_DC_sector"     ,&pips_DC_sector        );
    MergedTree->Branch("pips_Chi2N"         ,&pips_Chi2N            );
    //    MergedTree->Branch("pips_DC_x"          ,&pips_DC_x             );
    //    MergedTree->Branch("pips_DC_y"          ,&pips_DC_y             );

    MergedTree->Branch("pips_E_PCAL"        ,&pips_E_PCAL           );
    MergedTree->Branch("pips_E_ECIN"        ,&pips_E_ECIN           );

    MergedTree->Branch("pips_E_ECIN"        ,&pips_E_ECIN           );
    MergedTree->Branch("pips_E_ECOUT"       ,&pips_E_ECOUT          );
    MergedTree->Branch("DC_layer"           ,&DC_layer              );

    MergedTree->Branch("e"                  ,&e                     );
    MergedTree->Branch("piplus"             ,&piplus                );
    MergedTree->Branch("Ve"                 ,&Ve                    );
    MergedTree->Branch("Vpiplus"            ,&Vpiplus               );
    MergedTree->Branch("Beam"               ,&Beam                  );
    MergedTree->Branch("q"                  ,&q                     );

    MergedTree->Branch("xB"                 ,&xB                    );
    MergedTree->Branch("Q2"                 ,&Q2                    );
    MergedTree->Branch("omega"              ,&omega                 );
    MergedTree->Branch("z"                  ,&z                     );
    MergedTree->Branch("ePastSelectionCuts" ,&ePastSelectionCuts    );
    MergedTree->Branch("piplusPastSelectionCuts",&piplusPastSelectionCuts);


    // compute kinematics
    Q2      = -q.Mag2();
    omega   = q.E();
    Es      = Pn.E();
    Ps      = Pn.P();
    w2      = omega * omega;
    theta_sq= Pn.Angle();
    xB      = Q2 / (2. * Mp * omega);
    xPrime  = Q2 / (2. * ((Md - Es) * omega + Pn.Vect()*q.Vect() ));
    W       = sqrt(Mp2 - Q2 + 2. * omega * Mp);
    WPrime  = sqrt(Mp2 - Q2 + 2. * omega * (Md - Es) + 2. * Ps * sqrt(Q2 + w2) * cos( theta_sq ));
    
    
    alpha_s = ComputeLightConeFraction( Pn_qFrame );
     


    //    int NmergedEvents = 0;
    for (int MergedEvtId=0; MergedEvtId<Nevents2Merge; MergedEvtId++) {

        BANDTree -> GetEntry( BANDEventIndicesToMerge[MergedEvtId] );
        SIDISTree -> GetEntry( SIDISEventIndicesToMerge[MergedEvtId] );
        MergedTree -> Fill();
        StreamToCSVfile ({  (double)BANDrunID,  (double)BANDeventID,
            livetime,           current,
            xB,                 Q2,
            Ebeam,              z,
            W,                  alpha_s,
            WPrime,             xPrime,
            e->Px(),            e->Py(),
            e->Pz(),            e->E(),
            q->Px(),            q->Py(),
            q->Pz(),            q->E(),
            Ve->z(),            Vpiplus->z(),
            (double)goodneutron,
        },fdebug);

        if (fdebug>2){
            std::cout
            << "merging event " << BANDeventID << " from run " << BANDrunID
            << std::endl;
        }
    } // end merged event loop

    if (fdebug>2){
        std::cout << "merged " << Nevents2Merge << " SIDIS and BAND events." << std::endl;
    }
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
void OpenOutputFiles (TString RunStr, TString csvheader){
    
    
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
    
    Int_t     BANDeventID, NeventsBAND = BANDTree->GetEntries();
    BANDTree   -> SetBranchAddress("eventnumber"  ,&BANDeventID);
    for (int BANDevent=0; BANDevent < NeventsBAND ; BANDevent++){
        BANDTree -> GetEntry(BANDevent);
        BANDeventIDs[BANDevent] = BANDeventID;
    }
    Int_t     SIDISeventID, NeventsSIDIS = SIDISTree->GetEntries();
    SIDISTree  -> SetBranchAddress("eventnumber"  ,&SIDISeventID);
    for (int SIDISevent=0; SIDISevent < NeventsSIDIS ; SIDISevent++){
        SIDISTree -> GetEntry(SIDISevent);
        SIDISeventIDs[SIDISevent] = SIDISeventID;
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

    
    //    int SIDISeventIndexMin = 0;
    //
    //    for (int BANDevent=0; BANDevent < NeventsBAND ; BANDevent++){
    //
    //        BANDTree -> GetEntry(BANDevent);
    //
    //        // after filling the first evnet we do not have to go all the way back to
    //        // the first event in the SIDIS TTree
    //        SIDISeventIndexMin = SIDISEventIndicesToMerge[NmergedEvents];
    //
    //
    //        for (int SIDISevent=SIDISeventIndexMin; SIDISevent < NeventsSIDIS ; SIDISevent++){
    //
    //            SIDISTree -> GetEntry(SIDISevent);
    //
    //            if ( BANDeventID == SIDISeventID ){
    //
    ////                if (fdebug>3){
    ////                    std::cout
    ////                    << "merged event "   << BANDeventID         << " from run " << BANDrunID
    ////                    << " (in total "     << (NmergedEvents+1)   << " merges)"
    ////                    << std::endl;
    ////                }
    //                EventIDsToMerge[NmergedEvents] = BANDeventID;
    //                BANDEventIndicesToMerge[NmergedEvents] = BANDevent;
    //                SIDISEventIndicesToMerge[NmergedEvents] = SIDISevent;
    //                NmergedEvents ++ ;
    //                if ((NeventsToMerge>0) && (NmergedEvents >= NeventsToMerge)){
    //                    return NmergedEvents;
    //                }
    //            }
    //
    //            // in case the run number is identical in SIDIS and BAND trees,
    //            // we can cut the loop shorter by breaking if we passed the BAND event ID
    //            if ( SIDISeventID > BANDeventID ) {
    //                // if (fdebug>4){ std::cout << "SIDISeventID > BANDeventID, breaking..."<<std::endl ;}
    //                break;
    //            }
    //        }
    //    }
}



// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
Double_t ComputeLightConeFraction( TLorentzVector p ){
    // compute light-cone momentum fraction
    Double_t m = p.Mag();
    Double_t alpha = (p.E() - p.Z())/m;
    return alpha;
}
