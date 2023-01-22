// Read a HIPO file and convert a few of the variables to CSV
// clas12root -q ConvertHIPOFileToCSV.C
//


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
#include <string>
#include <iomanip>
#include <utility>
#include <iostream>
#include <stdexcept>
 
#include "clas12reader.h"
#include "Auxiliary/SIDISatBAND_auxiliary.cpp"

using namespace clas12;
SIDISatBAND_auxiliary  aux;
TString DataPath, prefix;
TString csvheader = ( (TString)"RunNumber,xB,omega,W,Q2,");
std::ofstream csvfile;

// auxiliary
// DCfid_SIDIS dcfid;
std::vector<region_part_ptr>  electrons;
int              RunNumber, EventNumber;
TLorentzVector       Beam, target, e, q;
double                 Q2, omega, W, xB;

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void  SetLorentzVector (TLorentzVector &p4, clas12::region_part_ptr rp){
    p4.SetXYZM(rp->par()->getPx(),
               rp->par()->getPy(),
               rp->par()->getPz(),
               p4.M());
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ComputeElectronKinematics(int fdebug){
    // compute event kinematics (from e-only information)
    q       = Beam - e;
    Q2      = -q.Mag2();
    omega   = q.E();
    xB      = Q2/(2. * aux.Mp * q.E());
    W       = 0; //sqrt((p_rest + q).Mag2());
}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void ConvertHIPOFileToCSV(TString   fInFilepath = "/work/cebaf24gev/sidis/reconstructed/rgc-unp-deut-22gev/hipo/",
                          TString  fOutFilepath = "/volatile/clas12/users/ecohen/BAND/Simulations/JLAB22GeV_Harut_Jan2023/",
                          TString     fFilename = "0000",
                          double         fEbeam = 22.0,
                          int            fdebug = 0){
    
    TString inputFile   = fInFilepath + fFilename + ".hipo";
    TString outfilename = fInFilepath + fFilename + ".csv";
    csvfile.open( outfilename );
    csvfile << csvheader << std::endl;
    
    TChain fake("hipo");
    fake.Add(inputFile.Data());
    //get the hipo data
    auto files  = fake.GetListOfFiles();
    RunNumber   = std::stoi(fFilename);
    EventNumber = 0;
    
    // step over events and extract information....
    for(Int_t i=0;i<files->GetEntries();i++){
        
        //create the event reader
        if (fdebug) std::cout << "reading file " << i << std::endl << files->At(i)->GetTitle() << std::endl;
        
        clas12reader c12(files->At(i)->GetTitle(),{0});
        
        // process the run
        auto scal          = c12.scalerReader();
        
        // go to next event
        while((c12.next()==true)){
            auto run           = c12.runconfig()->getRun();
            
            Beam.SetPxPyPzE (0, 0, fEbeam, fEbeam );
            
            // Get Particles By Type
            electrons          = c12.getByID( 11   );
            int             Ne = electrons.size();
            if (Ne==0) continue;
            double  leading_e_E;
            int     leading_e_index = 0;
            SetLorentzVector(e,electrons[0]);
            TLorentzVector e_tmp(0,0,0,0.000511);
            for (int eIdx=0; eIdx < Ne; eIdx++) {
                SetLorentzVector(e_tmp  ,electrons[eIdx]);
                double Ee = e_tmp.E();
                if (Ee > leading_e_E) {
                    leading_e_index = eIdx;
                    leading_e_E     = Ee;
                }
            }
            // set leading electron 4-momentum
            SetLorentzVector(e , electrons[leading_e_index]);
            
            
            if (fdebug){
                std::cout
                << "run " << RunNumber
                << std::endl;
            }
            
            csvfile
            << RunNumber    << ","
            << EventNumber  << ","
            << xB           << ","
            << omega        << ","
            << Q2           << ","
            << W            << ","
            << std::endl;
            
            EventNumber++;
        }
        
    } // end file loop
    csvfile.close();
    
}




