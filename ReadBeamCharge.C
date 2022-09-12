// Read beam current from HIPO files


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
#include "Auxiliary/SIDISatBAND_auxiliary.cpp"

using namespace clas12;
SIDISatBAND_auxiliary  aux;
TString DataPath, prefix;
TString csvheader = ( (TString)"runnum,beam_charge,");
std::ofstream csvfile;

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void SetDataPath (TString fDataPath) {
    if (DataPath=="" || DataPath=="sidisdvcs" || DataPath=="sidis dvcs"){
        // sidis-dvcs train files, used since July 2022
        // (the 'usual' train files)
        DataPath = "/cache/clas12/rg-b/production/recon/spring2019/torus-1/pass1/v0/dst/train/sidisdvcs/";
        prefix   = "sidisdvcs_";
    }
    else if (DataPath=="inclusive" || DataPath=="inc"){
        // inclusive train files, used until July 2022
        // (inclusive train files were only generated in the beginning of RGB without any backup)
        DataPath = "/volatile/clas12/rg-b/production/recon/spring2019/torus-1/pass1/v0/dst/train_20200610/inc/";
        prefix   = "inc_";
    }
}



// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void ReadBeamCharge( int RunNumber=6420, int fdebug=0, TString fDataPath = "sidisdvcs" ){
    
    TString RunNumberStr = aux.GetRunNumberSTR ( RunNumber );
    SetDataPath(fDataPath);
    
    TString inputFile   = DataPath + prefix + RunNumberStr + ".hipo";
    TString outfilename = "/volatile/clas12/users/ecohen/BAND/metaData/beam_charge_"+RunNumberStr+".csv";
    csvfile.open( outfilename );
        csvfile << csvheader << std::endl;
    
    TChain fake("hipo");
    fake.Add(inputFile.Data());
    //get the hipo data
    auto files = fake.GetListOfFiles();
    
    // step over events and extract information....
    for(Int_t i=0;i<files->GetEntries();i++){
        
        //create the event reader
        if (fdebug) std::cout << "reading file " << i << std::endl << files->At(i)->GetTitle() << std::endl;
        
        clas12reader c12(files->At(i)->GetTitle(),{0});
        
        // process the run
        auto scal          = c12.scalerReader();
        auto RunBeamCharge = c12.getRunBeamCharge();
        
        // go to next events to ask for run number
        c12.next();
        auto run           = c12.runconfig()->getRun();
        if (fdebug)
            std::cout
            << "run " << RunNumber << ", beam charge: "
            << RunBeamCharge << std::endl;
        
        csvfile << RunNumber << "," << RunBeamCharge << "," << std::endl;
        
    } // end file loop
    csvfile.close();
}



