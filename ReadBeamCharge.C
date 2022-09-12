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


//TString DataPath = "/volatile/clas12/rg-b/production/recon/spring2019/torus-1/pass1/v0/dst/train_20200610/inc/";
TString DataPath = "/cache/clas12/rg-b/production/recon/spring2019/torus-1/pass1/v0/dst/train/sidisdvcs/";



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
        if (fdebug) std::cout << "reading file " << i << std::endl << files->At(i)->GetTitle() << std::endl;
        
        clas12reader c12(files->At(i)->GetTitle(),{0});
        
        // process the run
        auto scal          = c12.scalerReader();
        auto RunBeamCharge = c12.getRunBeamCharge();
        if (fdebug){
            std::cout << "beam charge: " << RunBeamCharge << std::endl;
        }
        
        // go to next events to ask for run number
        c12.next();
        auto run           = c12.runconfig()->getRun();
        if (fdebug) std::cout << "run " << run  << std::endl;
        
        
    } // end file loop
}



