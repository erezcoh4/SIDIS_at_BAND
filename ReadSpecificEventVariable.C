// Read specific event variable
// last update May-4 2022


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
#define NMAXPIONS 20 // maximal allowed number of pions
#include "clas12reader.h"
#include "Auxiliary/SIDISatBAND_auxiliary.cpp"


TString  indatapath = "/volatile/clas12/rg-b/production/recon/spring2019/torus-1/pass1/v0/dst/train_20200610/inc/";
TString outdatapath = "/volatile/clas12/users/ecohen/BAND/SIDIS_skimming/";


SIDISatBAND_auxiliary aux;
int                          Npims;
std::vector<region_part_ptr> electrons, neutrons, protons, pipluses, piminuses;
std::vector<TLorentzVector>  piplus, piminus;
std::ofstream                csvfile;
int              NeventsMaxToProcess;

void ReadSpecificEventVariable(int RunNumber=6420,
                    int NeventsMax=-1,
                    TString variable="theta_pims",
                    int fdebug=0 ){
    
    
    TString RunNumberStr = aux.GetRunNumberSTR ( RunNumber );
    TString   infilename = indatapath + "inc_" + RunNumberStr + ".hipo";
    TString  outfilename = outdatapath + "skimmed_SIDIS_inc_" + RunNumberStr + "_" + variable + ".csv";
    aux.OpenCSVfile( csvfile, outfilename, variable.Data() );
    
    
    TChain fake("hipo");
    fake.Add(infilename.Data());
    auto files = fake.GetListOfFiles();
    
    // step over events and extract information....
    for(Int_t i=0;i<files->GetEntries();i++){
        
        clas12reader c12(files->At(i)->GetTitle(),{0});
        NeventsMaxToProcess = NeventsMax;
        if (NeventsMax<0) NeventsMaxToProcess = c12.getReader().getEntries();
        int event = 0;
        
        while((c12.next()==true) && (event < NeventsMaxToProcess)){
            
            piminus.clear(); piminuses = c12.getByID( -211 ); Npims = piminuses .size();
            
            if (Npims>0){
                std::cout << "Npims: " << Npims << std::endl;
                for (int pimsIdx=0; pimsIdx < Npims; pimsIdx++) {
                    std::cout << "pimsIdx: " << pimsIdx << std::endl;
                    piminus[pimsIdx] = TLorentzVector(0,0,0,0.139570);
                    aux.SetParticle4Momentum( piminus[pimsIdx]  ,piminuses[pimsIdx]);
                    std::cout << "piminus[pimsIdx].Theta(): " << piminus[pimsIdx].Theta() << std::endl;
                    
                    if (variable == "theta_pims"){
                        aux.StreamToCSVfile (csvfile, {piminus[pimsIdx].Theta()} );
                    }
                }
            }
        }
    } // end file loop
    std::cout << "Done. see " << std::endl << outfilename << std::endl;
}

