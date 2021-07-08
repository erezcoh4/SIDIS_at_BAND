// last edit July-5, 2021 (EOC, mbp), see README


//#include "/u/home/cohen/BAND_analysis/clas12root/Erez_analysis/Auxiliary/reader.h"
#include "Auxiliary/bank.h"

#include "/u/home/cohen/BAND_analysis/clas12root/Erez_analysis/Auxiliary/BBand.h"
#include "/u/home/cohen/BAND_analysis/clas12root/Erez_analysis/Auxiliary/BEvent.h"

#include "/u/home/cohen/BAND_analysis/clas12root/Erez_analysis/Auxiliary/RCDB/Connection.h"

#include "/u/home/cohen/BAND_analysis/clas12root/Erez_analysis/Auxiliary/constants.h"
#include "/u/home/cohen/BAND_analysis/clas12root/Erez_analysis/Auxiliary/readhipo_helper.h"


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



void          OpenInputFiles (TString RunStr);
void         OpenOutputFiles (TString RunStr, TString header);

void         CloseInputFiles ();
void        CloseOutputFiles ();
void MergeSIDISandBANDevents (int NeventsToMerge=5,
                              int fdebug=2,
                              int PrintProgress=5000);


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
// Main functionality
// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void MergeSIDISandBANDSkimmers(int RunNumber=6420,
                               int NeventsToMerge=5,
                               int fdebug=2,
                               int PrintProgress=5000){
    
    char RunNumberStr[20];
    sprintf( RunNumberStr, "00%d", RunNumber );
    OpenInputFiles   ( (TString)RunNumberStr );
    OpenOutputFiles  ( (TString)RunNumberStr,
                      "variables" );
    
    MergeSIDISandBANDevents( NeventsToMerge, fdebug, PrintProgress );
    
    
    CloseOutputFiles ();
    CloseInputFiles  ();
    
}



// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void MergeSIDISandBANDevents(int NeventsToMerge, int fdebug, int PrintProgress){
    
    Int_t   BANDrunID, BANDeventID, SIDISrunID, SIDISeventID;
    Int_t   NeventsBAND  = BANDTree->GetEntries();
    Int_t   NeventsSIDIS = SIDISTree->GetEntries();
    
    SIDISTree  -> SetBranchAddress("eventnumber"  ,&SIDISeventID);
    SIDISTree  -> SetBranchAddress("runnum"       ,&SIDISrunID);
    
    // BAND Tree
    double         Ebeam = 0;
    double  gated_charge = 0;
    double      livetime = 0;
    double     starttime = 0;
    double       current = 0;
    int      eventnumber = 0;
    bool     goodneutron = false;
    int       nleadindex = -1;
    double        weight = 0;
    //     Neutron info:
    int            nMult = 0;
    TClonesArray  * nHits = new TClonesArray("bandhit");
    TClonesArray &saveHit = *nHits;
    
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
    
    if (fdebug>1) {
        std::cout
        << "stepping over "
        << NeventsBAND << " BAND and "
        << NeventsSIDIS << "SIDIS events"
        << std::endl;
    }
    int NmergedEvents = 0;
    for (int BANDevent=0; BANDevent < NeventsBAND ; BANDevent++){
        
        BANDTree -> GetEntry(BANDevent);
        
        for (int SIDISevent=0; SIDISevent < NeventsSIDIS ; SIDISevent++){
            
            SIDISTree -> GetEntry(SIDISevent);
            
            if (fdebug>2){
                std::cout
                << "BAND run "  << BANDrunID
                << ", event "   << BANDeventID
                << ", "
                << "SIDIS run " << SIDISrunID
                << ", event "   << SIDISeventID
                << std::endl;
            }
            
            if ( (BANDrunID == SIDISrunID) && (BANDeventID == SIDISeventID)){
                // Can merge the event...
                if (fdebug>1){
                    std::cout
                    << "merged event "   << BANDeventID     << " from run " << BANDrunID
                    << " which is the merged event number " << (NmergedEvents+1)
                    << std::endl;
                }
                
                MergedTree -> Fill();
                // record event
                NmergedEvents ++;
                
                if (NmergedEvents >= NeventsToMerge){
                    std::cout << "merged " << NmergedEvents << " events, breaking." << std::endl;
                    return;
                }
            }
            
        } // end SIDIS event loop
        
    } // end BAND event loop
    std::cout << "merged " << NmergedEvents << " SIDIS and BAND events." << std::endl;
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
void OpenOutputFiles (TString RunStr, TString header){
    
    
    TString skimmedMergedFilename = (DataPath + "merged_SIDIS_and_BAND_skimming/"
                                     + "skimmed_SIDIS_and_BAND_inc_"  + RunStr );
    
    
    
    // Create output tree
    MergedFile = new TFile( skimmedMergedFilename + ".root" ,"RECREATE");
    MergedTree = new TTree( "T" , "Event information from merged SIDIS and BAND skimmers");
    
    // Create output csv files
    CSVfile.open( skimmedMergedFilename + ".csv" );
    CSVfile << header << std::endl;
    
}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void CloseInputFiles (){
    SIDISFile->Close();
    BANDFile->Close();
}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void CloseOutputFiles (){
    
    // close output CSV
    CSVfile.close();
    
    // close output ROOT
    MergedFile->cd();
    MergedTree->Write();
    MergedFile->Close();
    
}

