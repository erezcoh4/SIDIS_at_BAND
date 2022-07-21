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

// Output CSV file
std::ofstream   CSVfile_e_e_pips, CSVfile_e_e_pims;

// time
clock_t tStart = clock();

// globals
int               Nevents_e_e_pips;
int               Nevents_e_e_pims;
int                         fdebug;
TString                   DataPath;
TString                   FileName;
TString                CSVFilename;
TString            csvheader = ((TString)"status,runnum,evnum,"
                                +(TString)"e_P,e_Theta,e_Phi,e_Vz,"
                                +(TString)"pi_P,pi_Theta,pi_Phi,pi_Vz,"
                                +(TString)"Q2,W,xB,Zpi,"
                                +(TString)"omega,xF,y,"
                                +(TString)"M_X_ee_pi,"
                                +(TString)"theta_sq,WPrime,"
                                +(TString)"e_DC_sector,pi_DC_sector,");

TFile                      RGAFile;
TTree                      RGATree;

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




// methods
void                   PrintDone ();
void               OpenInputFile ();
void             OpenOutputFiles ();
void       InitializeFileReading ();
//void            StreamToCSVfile (std::vector<Double_t> observables, bool passed_cuts_e_pi_kinematics);
void                SetVerbosity (int ffdebug )          {fdebug = ffdebug;} ;
void                 SetDataPath (TString fDataPath )    {DataPath = fDataPath;} ;
void                 SetFileName (TString fFileName )    {FileName = fFileName;} ;



// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
// Main functionality
// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void ReadIgorRGAFile(TString fFileName="ntupleNew",
                     int NMAXeventsToMerge=-1,
                     int ffdebug=1,
                     int PrintProgress=1000,
                     TString fDataPath="/Users/erezcohen/Desktop/data/BAND/RGA_Free_proton/"
                     ){
    
    SetVerbosity          ( ffdebug );
    SetDataPath           ( fDataPath );
    SetFileName           ( fFileName );
    OpenInputFile         ();
    OpenOutputFiles       ();
    InitializeFileReading ();
    CloseOutputFile       ();
    CloseInputFile        ();
    PrintDone             ();
}




// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void PrintDone(){
    std::cout << "Done reading RGA file " << FileName << ", execution time: "
    << double(clock() - tStart) / (double)CLOCKS_PER_SEC
    << " sec "<< std::endl;
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void OpenInputFile (){
    
    TString InputROOTFilename = DataPath + Filename + ".root";
    if (fdebug>2) std::cout << "Opening " << InputROOTFilename << std::endl;
    RGAFile = new TFile( InputROOTFilename );
    
    RGATree = (TTree*)RGAFile->Get("T");
    if (fdebug>2) RGATree->Print();
}

//CONTINUE HERE!
// ALso write a macro script that exexcutes this code

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void OpenOutputFiles (){
    
    CSVFilename_e_e_pips = DataPath + Filename + "_e_e_pips.csv";
    CSVFilename_e_e_pims = DataPath + Filename + "_e_e_pims.csv";
    
    if (fdebug>2) {
        std::cout << "Opening output files: "
        << CSVFilename_e_e_pips  << ".csv "
        << " and "
        << CSVFilename_e_e_pims  << ".csv "
        << std::endl;
    }
    
    // Write csv header output csv files
    CSVfile_e_e_pips.open( CSVFilename + ".csv" );
    CSVfile_e_e_pips << csvheader << std::endl;
    
    if (fdebug>2) std::cout << csvheader << std::endl;
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void InitializeFileReading (){
    Nevents_e_e_pips = Nevents_e_e_pims = 0;
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void CloseInputFile (){
    RGAFile->Close();
    
}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void CloseOutputFiles (TString OutDataPath){
    
    // close output CSV
    CSVfile_e_pi.close();
    
    std::cout << "Done producing output file. They are ready in " << std::endl << DataPath << std::endl;
    std::cout << "Wrote "
    << Nevents_e_e_pips << " (e,e'π+) events"
    << " and "
    << Nevents_e_e_pims << " (e,e'π-) events." << std::endl;
}

