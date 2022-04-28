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


TString DataPath = "/volatile/clas12/rg-b/production/recon/spring2019/torus-1/pass1/v0/dst/train_20200610/inc/";



SIDISatBAND_auxiliary  aux;



void ReadBeamCharge( int RunNumber=6420, int fdebug=0 ){
    
//    TString outfilepath = "/volatile/clas12/users/ecohen/BAND/SIDIS_skimming/";
    TString RunNumberStr = aux.GetRunNumberSTR ( RunNumber );
    TString outfilename = "skimmed_SIDIS_inc_" + RunNumberStr;
    

    TString inputFile = DataPath + "inc_" + RunNumberStr + ".hipo";
    TChain fake("hipo");
    fake.Add(inputFile.Data());
    //get the hipo data
    auto files = fake.GetListOfFiles();
    
    // step over events and extract information....
    for(Int_t i=0;i<files->GetEntries();i++){
        
        //create the event reader
        if (fdebug) std::cout << "reading file " << i << std::endl;
        clas12reader c12(files->At(i)->GetTitle(),{0});
        
        // process the events...
        auto run           = c12.runconfig()->getRun();
        auto RunBeamCharge = c12.getRunBeamCharge();
        if (fdebug)
            std::cout << "run " << run        
            << ", beam charge: " << RunBeamCharge << std::endl;
        
    } // end file loop
}



