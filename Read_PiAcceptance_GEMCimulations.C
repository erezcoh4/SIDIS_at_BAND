
// clas12root Read_PiAcceptance_GEMCimulations.C\(\10, 1\)

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
#include "clas12reader.h"
#include "Auxiliary/DCfid_SIDIS.cpp"
#include "Auxiliary/csv_reader.h"
#define NMAXPIONS 20 // maximal allowed number of pions
#define r2d 180./3.1415 // radians to degrees
using namespace clas12;


// globals
auto db = TDatabasePDG::Instance();
TString DataPath, FileLabel, PiCharge;


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
// start clock
auto start = std::chrono::high_resolution_clock::now();
// declare methods
TVector3                GetParticleVertex (clas12::region_part_ptr rp);
void                     SetLorentzVector (TLorentzVector &p4, clas12::region_part_ptr rp);
void                      OpenOutputFiles (TString filelabel, TString header);
void                     CloseOutputFiles ();
void                      StreamToCSVfile (std::vector<Double_t> observables,
                                           int fdebug);
void                InitializeFileReading (int NeventsMax,int c12Nentries, int fdebug);
void                  InitializeVariables ();
void                      OpenResultFiles ();
void           ExtractElectronInformation (int fdebug);
void              ExtractPionsInformation (int fdebug);
void               ExtractPipsInformation (int pipsIdx, int fdebug );
void               ExtractPimsInformation (int pimsIdx, int fdebug );
void                    ComputeKinematics ();
void                   WriteEventToOutput (int fdebug);
void                        FinishProgram ();
void                   GetParticlesByType (int evnum, int fdebug );
void              Stream_e_pi_line_to_CSV (int piIdx, int fdebug );
void                          SetDataPath ( TString fDataPath )  {DataPath = fDataPath;}
void                         SetFileLabel ( TString fFileLabel ) {FileLabel = fFileLabel;}
void                          SetPiCharge ( TString fPiCharge )  {PiCharge = fPiCharge;}
bool   CheckIfElectronPassedSelectionCuts (Double_t e_PCAL_x, Double_t e_PCAL_y,
                                           Double_t e_PCAL_W,Double_t e_PCAL_V,
                                           Double_t e_E_PCAL,
                                           Double_t e_E_ECIN, Double_t e_E_ECOUT,
                                           TLorentzVector e,
                                           TVector3 Ve,
                                           Double_t e_DC_sector,
                                           Double_t e_DC_x[3],
                                           Double_t e_DC_y[3],
                                           Double_t e_DC_z[3],
                                           int torusBending);

bool       CheckIfPionPassedSelectionCuts (TString pionCharge, // "pi+" or "pi-"
                                           Double_t DC_sector,
                                           Double_t DC_x[3], Double_t DC_y[3], Double_t DC_z[3],
                                           Double_t chi2PID, Double_t p,
                                           TVector3 Ve,      TVector3 Vpi,
                                           int pipsIdx,
                                           int fdebug);
Double_t          Chi2PID_pion_lowerBound (Double_t p, Double_t C=0.88); // C(pi+)=0.88, C(pi-)=0.93
Double_t          Chi2PID_pion_upperBound (Double_t p, Double_t C=0.88); // C(pi+)=0.88, C(pi-)=0.93

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.


// meta-data
int           torusBending = -1; // -1 for In-bending, +1 for Out-bending
int    DC_layers[3] = {6,18,36};// Region 1 is denoted at DC detector 6, Region 2 is denoted 18, Region 3 - as 36
int                    DC_layer;
int                      runnum;
int                       evnum;
int               beam_helicity; // helicity of the electron +1 along the beam and -1 opposite to it
int                      status;
int         NeventsMaxToProcess;
int           Nevents_processed;

// number of particles per event
int         Ne, Nn, Np, Npips, Npims, Ngammas;
int                          Nd; // number of detected deuterons

// variables
double          Mp = 0.938;
double       Mp2 = Mp * Mp;
double       Md = 1.875612; // NIST
// leading electron
double            e_E_PCAL; // electron energy deposit in PCAL [GeV]
double            e_E_ECIN; // electron energy deposit in ECAL_in [GeV]
double           e_E_ECOUT; // electron energy deposit in ECAL_out [GeV]
double            e_PCAL_W;
double            e_PCAL_V;
double            e_PCAL_x;
double            e_PCAL_y;
double            e_PCAL_z;
double       e_PCAL_sector;
double         e_DC_sector;
double          e_DC_Chi2N;
double           e_DC_x[3];
double           e_DC_y[3];
double           e_DC_z[3];

bool             e_reconstructed;
bool            pi_reconstructed;
bool               e_passed_cuts;
bool              pi_passed_cuts;
bool     pi_passed_fiducial_cuts;
bool          pi_passed_PID_cuts;


// positive pions
bool         pipsReconstructed[NMAXPIONS];
bool     pipsPastSelectionCuts[NMAXPIONS];
bool      pipsPastFiducialCuts[NMAXPIONS];
bool           pipsPastPIDCuts[NMAXPIONS];
bool eepipsPastKinematicalCuts[NMAXPIONS];
double        pips_chi2PID[NMAXPIONS];
double         pips_PCAL_W[NMAXPIONS];
double         pips_PCAL_V[NMAXPIONS];
double         pips_PCAL_x[NMAXPIONS];
double         pips_PCAL_y[NMAXPIONS];
double         pips_PCAL_z[NMAXPIONS];
double    pips_PCAL_sector[NMAXPIONS];
double      pips_DC_sector[NMAXPIONS];
double          pips_Chi2N[NMAXPIONS];
double        pips_DC_x[NMAXPIONS][3];
double        pips_DC_y[NMAXPIONS][3];
double        pips_DC_z[NMAXPIONS][3];
double         pips_E_PCAL[NMAXPIONS];
double         pips_E_ECIN[NMAXPIONS];
double        pips_E_ECOUT[NMAXPIONS];
double               Zpips[NMAXPIONS]; // hadron rest-frame energy

double           piplus_Px[NMAXPIONS];
double           piplus_Py[NMAXPIONS];
double           piplus_Pz[NMAXPIONS];
double            piplus_E[NMAXPIONS];
double           Vpiplus_X[NMAXPIONS];
double           Vpiplus_Y[NMAXPIONS];
double           Vpiplus_Z[NMAXPIONS];



// negative pions
bool         pimsReconstructed[NMAXPIONS];
bool     pimsPastSelectionCuts[NMAXPIONS];
bool      pimsPastFiducialCuts[NMAXPIONS];
bool           pimsPastPIDCuts[NMAXPIONS];
bool eepimsPastKinematicalCuts[NMAXPIONS];
double        pims_chi2PID[NMAXPIONS];
double         pims_PCAL_W[NMAXPIONS];
double         pims_PCAL_V[NMAXPIONS];
double         pims_PCAL_x[NMAXPIONS];
double         pims_PCAL_y[NMAXPIONS];
double         pims_PCAL_z[NMAXPIONS];
double    pims_PCAL_sector[NMAXPIONS];
double      pims_DC_sector[NMAXPIONS];
double          pims_Chi2N[NMAXPIONS];
double        pims_DC_x[NMAXPIONS][3];
double        pims_DC_y[NMAXPIONS][3];
double        pims_DC_z[NMAXPIONS][3];
double         pims_E_PCAL[NMAXPIONS];
double         pims_E_ECIN[NMAXPIONS];
double        pims_E_ECOUT[NMAXPIONS];
double               Zpims[NMAXPIONS]; // hadron rest-frame energy

double           piminus_Px[NMAXPIONS];
double           piminus_Py[NMAXPIONS];
double           piminus_Pz[NMAXPIONS];
double            piminus_E[NMAXPIONS];
double           Vpiminus_X[NMAXPIONS];
double           Vpiminus_Y[NMAXPIONS];
double           Vpiminus_Z[NMAXPIONS];

// Output CSV file
std::ofstream   CSVfile;

// vectors in lab-frame
TLorentzVector          Beam, target, e, q;
std::vector<TLorentzVector>         piplus; // positive pions
std::vector<TLorentzVector>        piminus; // negative pions
// reconstructed vertex position
TVector3                                Ve;
std::vector<TVector3>              Vpiplus;
std::vector<TVector3>             Vpiminus;

// kinematics
Double_t     Ebeam, xB, Q2, omega, W, W2, xF, y, M_X;

// MC information
TLorentzVector          P_mc_particle;
TLorentzVector          e_g,    pi_g;
TVector3                V_mc_particle;
TVector3                Ve_g,   Vpi_g;

// auxiliary
DCfid_SIDIS dcfid;
std::vector<region_part_ptr>  electrons, neutrons, protons, pipluses, piminuses, gammas, deuterons;


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
// Main functionality
// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void Read_PiAcceptance_GEMCimulations(TString fPiCharge = "pips",
                                      int  NeventsMax=-1,
                                      int  fdebug=1,
                                      int  PrintProgress=50000,
                                      TString fFileLabel = "p_uniform_distribution",
                                      TString fDataPath = "/volatile/clas12/users/ecohen/GEMC/hipo/10.2/AcceptanceCorrection/"
                                      ){
    if (fdebug>2) std::cout << "SetDataPath ( "<<fDataPath<<"  ); " << std::endl;
    SetDataPath ( fDataPath  );
    
    if (fdebug>2) std::cout << "SetFileLabel ( "<<fFileLabel<<"  ); " << std::endl;
    SetFileLabel( fFileLabel );
    
    if (fdebug>2) std::cout << "SetPiCharge ( "<<fPiCharge<<"  ); "  << std::endl;
    SetPiCharge ( fPiCharge  );
    // open result files
    if (fdebug>2) std::cout << "OpenResultFiles (); "  << std::endl;
    OpenResultFiles();
    
    TString inputFile = DataPath + "/" + PiCharge + "/ee" + PiCharge + "_" + FileLabel + "_reconstructed.hipo";
    if (fdebug>2) std::cout << "inputFile: " << inputFile << std::endl;
    
    TChain fake("hipo");
    fake.Add(inputFile.Data());
    //get the hipo data
    auto files = fake.GetListOfFiles();
    
    // step over events and extract information....
    for(Int_t i=0;i<files->GetEntries();i++){
            
        //create the event reader
        if (fdebug) std::cout << "reading file " << i << std::endl;
        clas12reader c12(files->At(i)->GetTitle(),{0});
        InitializeFileReading( NeventsMax, c12.getReader().getEntries(), fdebug );
        int event = 0;

        // process the events...
        while((c12.next()==true) && (event < NeventsMaxToProcess)){
            
            runnum = c12.runconfig()->getRun();
            evnum  = c12.runconfig()->getEvent();
            if (fdebug>2) std::cout << "begin analysis of event " << evnum << " (run " << runnum << ")"  << std::endl;
            
            
            InitializeVariables();
            // Get Particles By Type
            electrons   = c12.getByID( 11   );
            neutrons    = c12.getByID( 2112 );
            protons     = c12.getByID( 2212 );
            pipluses    = c12.getByID( 211  );
            piminuses   = c12.getByID(-211  );
            gammas      = c12.getByID( 22   );
            deuterons   = c12.getByID( 1000010020 );
            GetParticlesByType ( evnum, fdebug );
            
            
            
            // add truth-information,
            // i.e. generated electron and generated pion information
            auto mcpbank = c12.mcparts();
            const Int_t Ngen=mcpbank->getRows();
            if (fdebug>1) std::cout << "Grabbing truth-information of " << Ngen << " particles" << std::endl;
            
            for( Int_t i_mc =0; i_mc< Ngen ; i_mc++){
                mcpbank -> setEntry(i_mc);
                
                P_mc_particle.SetXYZM( mcpbank->getPx() , mcpbank->getPy() , mcpbank->getPz() , mcpbank->getMass() );
                V_mc_particle.SetXYZ( mcpbank->getVx() , mcpbank->getVy() , mcpbank->getVz() );
                auto pid = mcpbank->getPid();
                
                if (fdebug>2){
                    std::cout << "MC particle PDG code " << pid
                    << std::setprecision(4)
                    << ", p: "    << P_mc_particle.P()          << " GeV/c, "
                    << ", theta: "<< P_mc_particle.Theta()      << ", " << P_mc_particle.Theta()*r2d  << " deg, "
                    << ", phi: "  << P_mc_particle.Phi()        << ", "<< P_mc_particle.Phi()*r2d    << " deg, "
                    << ", V(z): " << V_mc_particle.Z()          << " cm"
                    << std::endl;
                }
                
                if ( pid==11 ) {
                    e_g = P_mc_particle;
                    Ve_g = V_mc_particle;
                }
                else if ( (pid==211) && (PiCharge=="pips") ) {
                    pi_g = P_mc_particle;
                    Vpi_g = V_mc_particle;
                }
                else if ( (pid==-211) && (PiCharge=="pims") ) {
                    pi_g = P_mc_particle;
                    Vpi_g = V_mc_particle;
                }
                else  {
                    if (fdebug>2){
                        std::cout << "MC particle PDG code " << pid << " do not match generated particles: e (" << 11 << ") + ";
                        if ( PiCharge=="pips")
                            std::cout << " pi+ (" << 211;
                        else if ( PiCharge=="pims")
                            std::cout << " pi- (" << -211;
                        std::cout << ")" << std::endl;
                    }
                }
                
            }
            
            if (Ne>0){
                // We do not filter, so we write event non-reconstructed events.
                // But here we extract information from electrons and pions
                e_reconstructed = true;
                ExtractElectronInformation  (fdebug);
                ExtractPionsInformation     (fdebug);
            }
            else {
                e_reconstructed = false;
                if (fdebug>1) std::cout << "no electrons in event " << event << std::endl;
            }
            WriteEventToOutput              (fdebug);
            
            if (fdebug>1) {
                std::cout << "done processing event " << evnum
                << " (" << event << "/" << NeventsMaxToProcess<< ") "
                << std::endl << "------------------------------------------------------------" << std::endl ;
            }
            event++; Nevents_processed++;
            if (fdebug && event%PrintProgress==0) std::cout << std::setprecision(1) << " event " << event << std::endl;
        } // end event loop
        
    } // end file loop
    
    
    FinishProgram();
}



// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
TVector3 GetParticleVertex(clas12::region_part_ptr rp){
    TVector3 V(rp->par()->getVx(),
               rp->par()->getVy(),
               rp->par()->getVz());
    return V;
}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void SetLorentzVector (TLorentzVector &p4,clas12::region_part_ptr rp){
    p4.SetXYZM(rp->par()->getPx(),
               rp->par()->getPy(),
               rp->par()->getPz(),
               p4.M());
}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void OpenOutputFiles (TString header){
        
    // Create output csv files
    CSVfile.open( DataPath + "/" + PiCharge + "/ee" + PiCharge + "_" + FileLabel + ".csv" );
    CSVfile << header << std::endl;
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void CloseOutputFiles (){
    // close output CSV
    CSVfile.close();
    
    
    std::cout
    << "Done processesing "  <<  Nevents_processed          << " events,"
    << std::endl;
    
    
    
    std::cout << "output files ready in csv formats in " << std::endl
    << std::endl
    << DataPath + "/" + PiCharge + "/" + "ee" + PiCharge + "_" + FileLabel + ".csv"
    << std::endl;
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void StreamToCSVfile (std::vector<Double_t> observables,
                      int fdebug){
    if (fdebug>1) {
        std::cout << "streaming to CSVfile" << std::endl;
    }
    for (auto v:observables) CSVfile << std::fixed << v << ",";
    CSVfile << std::endl;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
int GetBeamHelicity( event_ptr p_event, int runnum, int fdebug ){
    // deprecated as of Aug-11, 2021,
    //    since
    // however we keep it here for more data
    // get beam helicity (+1 along the beam and -1 opposite to it)
    // [Christopher Dilks <dilks@jlab.org>, email from Aug-5, 2021]
    // for more items
    // [https://github.com/JeffersonLab/clas12root/blob/master/AccesssingBankDataInCpp.txt]
    
    //// helFlip: if true, REC::Event.helicity has opposite sign from reality
    //def helFlip
    //if(RG=="RGA") helFlip = true
    //else if(RG=="RGB") {
    //  helFlip = true
    //  if(runnum>=11093 && runnum<=11283) helFlip = false // fall, 10.4 GeV period only
    //  else if(runnum>=11323 && runnum<=11571) helFlip = false // winter
    //};
    //else if(RG=="RGK") helFlip = false
    //else if(RG=="RGF") helFlip = true
    if (fdebug>3) std::cout << "beam_helicity = c12.event()->getHelicity()" << std::endl;
    
    beam_helicity = p_event->getHelicity();
    
    if (fdebug>3) std::cout << "check spin flip" << std::endl;
    // we are working here on RGB data
    bool helFlip = true;
    if      (runnum>=11093 && runnum<=11283)    helFlip = false; // falls, 10.4 GeV period only
    else if (runnum>=11323 && runnum<=11571)    helFlip = false; // winter
    
    if (helFlip) {
        beam_helicity = -1 * beam_helicity;
    }
    if (fdebug>3) std::cout << "done GetBeamHelicity() " << std::endl;
    return beam_helicity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
double GetBeamEnergy (int fdebug){

    if (fdebug>3) std::cout << "set beam energy" << std::endl;
    double Ebeam = 10.2; // [GeV] ( for Fall-2019 the enrgy was 10.4096)
    return Ebeam;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void InitializeFileReading(int NeventsMax, int c12Nentries, int fdebug){
    if (fdebug>1) {
        std::cout << "InitializeFileReading( " << NeventsMax << " , " << c12Nentries << " , " << fdebug << ")" << std::endl;
    }
    Ebeam   = GetBeamEnergy( fdebug );
    Beam    .SetPxPyPzE (0, 0, Ebeam, Ebeam );
    target  .SetXYZM    (0, 0, 0,     Md    );
    
    NeventsMaxToProcess = NeventsMax;
    if (NeventsMax<0) NeventsMaxToProcess = c12Nentries;
    Nevents_processed           = 0;
    if (fdebug>1) {
        std::cout << "NeventsMaxToProcess =  " << NeventsMaxToProcess << "" << std::endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void InitializeVariables(){
    e = TLorentzVector(0,0,0,db->GetParticle( 11   )->Mass());
    
    xB          = Q2        = omega     = -9999;
    xF          = y         = M_X       = -9999;
    e_E_ECIN    = e_E_ECOUT = e_E_PCAL  = -9999;
    e_PCAL_W    = e_PCAL_V              = -9999;
    e_PCAL_x    = e_PCAL_y  = e_PCAL_z  = -9999;
    e_PCAL_sector                       = -9999;
    e_DC_sector = e_DC_Chi2N            = -9999;
    for (int regionIdx=0; regionIdx<3; regionIdx++) {
        e_DC_x[regionIdx]               = -9999;
        e_DC_y[regionIdx]               = -9999;
        e_DC_z[regionIdx]               = -9999;
    }
    Ve                                  = TVector3();

    piplus      .clear();
    piminus     .clear();
    Vpiplus     .clear();
    Vpiminus    .clear();
    pipluses    .clear();
    piminuses   .clear();
    electrons   .clear();
    neutrons    .clear();
    protons     .clear();
    gammas      .clear();
    
    e_reconstructed = false;
    e_passed_cuts   = false;
    for (int piIdx=0; piIdx<NMAXPIONS; piIdx++) {
        pips_chi2PID[piIdx]                         = -9999;
        pips_DC_sector[piIdx]                       = -9999;
        pips_PCAL_sector[piIdx]                     = -9999;
        pips_PCAL_W[piIdx] = pips_PCAL_V[piIdx]     = -9999;
        pips_PCAL_x[piIdx] = pips_PCAL_y[piIdx]     = -9999;
        pips_PCAL_z[piIdx]                          = -9999;
        pips_E_PCAL[piIdx]                          = -9999;
        pips_E_ECIN[piIdx] = pips_E_ECOUT[piIdx]    = -9999;
        
        pims_chi2PID[piIdx]                         = -9999;
        pims_DC_sector[piIdx]                       = -9999;
        pims_PCAL_sector[piIdx]                     = -9999;
        pims_PCAL_W[piIdx] = pims_PCAL_V[piIdx]     = -9999;
        pims_PCAL_x[piIdx] = pims_PCAL_y[piIdx]     = -9999;
        pims_PCAL_z[piIdx]                          = -9999;
        pims_E_PCAL[piIdx]                          = -9999;
        pims_E_ECIN[piIdx] = pims_E_ECOUT[piIdx]    = -9999;
        for (int regionIdx=0; regionIdx<3; regionIdx++) {
            pips_DC_x[piIdx][regionIdx]= pips_DC_y[piIdx][regionIdx]    = -9999;
            pips_DC_z[piIdx][regionIdx]                                 = -9999;
            pims_DC_x[piIdx][regionIdx]= pims_DC_y[piIdx][regionIdx]    = -9999;
            pims_DC_z[piIdx][regionIdx]                                 = -9999;
        }
        piplus  .push_back( TLorentzVector(0,0,0,db->GetParticle( 211 )->Mass()) );
        Vpiplus .push_back( TVector3() );
        
        piminus .push_back( TLorentzVector(0,0,0,db->GetParticle( -211 )->Mass()) );
        Vpiminus.push_back( TVector3() );
        
        piplus_Px[piIdx]    = piplus_Py[piIdx]  = piplus_Pz[piIdx]  = piplus_E[piIdx]   = -9999;
        piminus_Px[piIdx]   = piminus_Py[piIdx] = piminus_Pz[piIdx] = piminus_E[piIdx]  = -9999;
        Vpiplus_X[piIdx]    = Vpiplus_Y[piIdx]  = Vpiplus_Z[piIdx]  = -9999;
        Vpiminus_X[piIdx]   = Vpiminus_Y[piIdx] = Vpiminus_Z[piIdx] = -9999;
        
        pipsPastSelectionCuts[piIdx]                = false;
        pipsPastFiducialCuts[piIdx]                 = false;
        pipsPastPIDCuts[piIdx]                      = false;
        
        pimsPastSelectionCuts[piIdx]                = false;
        pimsPastFiducialCuts[piIdx]                 = false;
        pimsPastPIDCuts[piIdx]                      = false;

    }
    DC_layer                                        = -9999;
    status                                          = 1; // 0 is good...
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TString GetRunNumberSTR(int RunNumber, int fdebug){
    char RunNumberStr[20];
    sprintf( RunNumberStr, "00%d", RunNumber );
    if (fdebug>1) std::cout << "(SIDIS) skimming run " << RunNumberStr << std::endl;
    return (TString)RunNumberStr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void OpenResultFiles(){
    OpenOutputFiles(( (TString)"e_P,e_Theta,e_Phi,e_Vz,"
                     +(TString)"pi_P,pi_Theta,pi_Phi,pi_Vz,"
                     +(TString)"Npips,Npims,Nelectrons,Ngammas,Nprotons,Nneutrons,Ndeuterons,"
                     +(TString)"e_P_g,e_Theta_g,e_Phi_g,e_Vz_g,"
                     +(TString)"pi_P_g,pi_Theta_g,pi_Phi_g,pi_Vz_g,"
                     +(TString)"pi_reconstructed,pi_passed_cuts,pi_passed_fiducial_cuts,pi_passed_PID_cuts,"
                     +(TString)"e_reconstructed,e_passed_cuts,"
                     +(TString)"e_DC_sector,pi_DC_sector")
                    );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExtractElectronInformation(int fdebug){
    // ------------------------------------------------------------------------------------------------
    // extract electron information
    // ------------------------------------------------------------------------------------------------
    // find leading electron as the one with highest energy
    double  leading_e_E;
    int     leading_e_index = 0;
    SetLorentzVector(e,electrons[0]);
    TLorentzVector e_tmp(0,0,0,db->GetParticle(11)->Mass());
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
    // set leading electron vertex
    Ve              = GetParticleVertex( electrons[leading_e_index] );
    
    // detector information on electron
    auto e_PCAL_info= electrons[leading_e_index]->cal(PCAL);
    e_E_PCAL        = e_PCAL_info->getEnergy();
    e_PCAL_sector   = e_PCAL_info->getSector();
    e_PCAL_V        = e_PCAL_info->getLv();
    e_PCAL_W        = e_PCAL_info->getLw();
    e_E_ECIN        = electrons[leading_e_index]->cal(ECIN)->getEnergy();
    e_E_ECOUT       = electrons[leading_e_index]->cal(ECOUT)->getEnergy();
    
    // hit position in PCAL
    e_PCAL_x        = e_PCAL_info->getX();
    e_PCAL_y        = e_PCAL_info->getY();
    e_PCAL_z        = e_PCAL_info->getZ();
    
    // Drift Chamber tracking system
    auto e_DC_info  = electrons[leading_e_index]->trk(DC);
    e_DC_sector     = e_DC_info->getSector(); // tracking sector
    e_DC_Chi2N      = e_DC_info->getChi2N();  // tracking chi^2/NDF
    
    for (int regionIdx=0; regionIdx<3; regionIdx++) {
        int DC_layer = DC_layers[regionIdx];
        e_DC_x[regionIdx] = electrons[leading_e_index]->traj(DC,DC_layer)->getX();
        e_DC_y[regionIdx] = electrons[leading_e_index]->traj(DC,DC_layer)->getY();
        e_DC_z[regionIdx] = electrons[leading_e_index]->traj(DC,DC_layer)->getZ();
    }
    
    e_passed_cuts = CheckIfElectronPassedSelectionCuts(e_PCAL_x, e_PCAL_y,
                                                              e_PCAL_W, e_PCAL_V,
                                                              e_E_PCAL, e_E_ECIN,
                                                              e_E_ECOUT,
                                                              e, Ve,
                                                              e_PCAL_sector, // e_PCAL_sector should be consistent with e_DC_sector
                                                              e_DC_x, e_DC_y, e_DC_z,
                                                              torusBending );

    if (fdebug > 2) std::cout << "extracted electron information and computed kinematics" << std::endl;
    
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExtractPionsInformation(int fdebug){
    
    // positive pions
//    if (evnum>67 && evnum<74)
//        std::cout << "event " << evnum << " in ExtractPionsInformation, Npips= " << Npips << std::endl;

    for (int pipsIdx=0; pipsIdx < Npips; pipsIdx++) {
//        std::cout << "event " << evnum << " pipsIdx= " << pipsIdx << std::endl;
        ExtractPipsInformation( pipsIdx, fdebug );
    }
    // negative pions
    for (int pimsIdx=0; pimsIdx < Npims; pimsIdx++) {
        ExtractPimsInformation( pimsIdx, fdebug );
    }
    
    // done
    if (fdebug > 2) std::cout << "done extracting pion information" << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void WriteEventToOutput(int fdebug){
    if (fdebug>3) std::cout << "Writing (e,e'pi) event" << std::endl;
    
    // if we simulated a positive pion
    // stream to csv all positive pions as seperate rows
    if (PiCharge=="pips") {
        if (Npips==0) {
            // this means that no pions were reconstructed,
            // and we just want to record the event that we generated
            Stream_e_pi_line_to_CSV( 0, fdebug );
        }
        else {
            for (int pipsIdx=0; pipsIdx<Npips; pipsIdx++) {
                Stream_e_pi_line_to_CSV( pipsIdx, fdebug );
            }
        }
    }
    // if we simulated a negative pion
    // stream to csv all negative pions as seperate rows
    else if (PiCharge=="pims") {
        if (Npims==0) {
            // this means that no pions were reconstructed,
            // and we just want to record the event that we generated
            Stream_e_pi_line_to_CSV( 0, fdebug );
        }
        else {
        for (int pimsIdx=0; pimsIdx<Npims; pimsIdx++) {
            Stream_e_pi_line_to_CSV( pimsIdx, fdebug );
        }
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void FinishProgram(){
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    
    std::cout << "Done. Elapsed time: " << elapsed.count() << std::endl;
    CloseOutputFiles();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExtractPipsInformation( int pipsIdx, int fdebug ){
    if (fdebug>2)
        std::cout << "ExtractPipsInformation( pipsIdx=" << pipsIdx << ", fdebug=" << fdebug << " )" << std::endl;
    
    
    // extract positive pion information
    SetLorentzVector(piplus[pipsIdx]  ,pipluses[pipsIdx]);
    Zpips[pipsIdx]              = piplus[pipsIdx].E() / omega;
    Vpiplus[pipsIdx]            = GetParticleVertex( pipluses[pipsIdx] );
    pips_chi2PID[pipsIdx]       = pipluses[pipsIdx]->par()->getChi2Pid();
    
    // EC in and out
    pips_E_ECIN[pipsIdx]        = pipluses[pipsIdx]->cal(ECIN)->getEnergy();
    pips_E_ECOUT[pipsIdx]       = pipluses[pipsIdx]->cal(ECOUT)->getEnergy();
    // PCAL
    auto pips_PCAL_info         = pipluses[pipsIdx]->cal(PCAL);
    pips_E_PCAL[pipsIdx]        = pips_PCAL_info->getEnergy();
    pips_PCAL_sector[pipsIdx]   = pips_PCAL_info->getSector();
    pips_PCAL_V[pipsIdx]        = pips_PCAL_info->getLv();
    pips_PCAL_W[pipsIdx]        = pips_PCAL_info->getLw();
    pips_PCAL_x[pipsIdx]        = pips_PCAL_info->getX();
    pips_PCAL_y[pipsIdx]        = pips_PCAL_info->getY();
    pips_PCAL_z[pipsIdx]        = pips_PCAL_info->getZ();
    // DC
    auto pips_DC_info           = pipluses[pipsIdx]->trk(DC);
    pips_DC_sector[pipsIdx]     = pips_DC_info->getSector(); // tracking sector
    pips_Chi2N[pipsIdx]         = pips_DC_info->getChi2N();  // tracking chi^2/NDF
    for (int regionIdx=0; regionIdx<3; regionIdx++) {
        DC_layer = DC_layers[regionIdx];
        pips_DC_x[pipsIdx][regionIdx] = pipluses[pipsIdx]->traj(DC,DC_layer)->getX();
        pips_DC_y[pipsIdx][regionIdx] = pipluses[pipsIdx]->traj(DC,DC_layer)->getY();
        pips_DC_z[pipsIdx][regionIdx] = pipluses[pipsIdx]->traj(DC,DC_layer)->getZ();
        if (fdebug>3) {
            std::cout
            << "pips_DC_sector[pipsIdx="<<pipsIdx<<"]="
            << pips_DC_sector[pipsIdx]
            << ", DC_layer = " << DC_layer
            << std::endl
            << "pips_DC_x[pipsIdx="<<pipsIdx<<"][regionIdx="<<regionIdx<<"]="
            << pips_DC_x[pipsIdx][regionIdx]
            << std::endl
            << "pips_DC_y[pipsIdx="<<pipsIdx<<"][regionIdx="<<regionIdx<<"]="
            << pips_DC_y[pipsIdx][regionIdx]
            << std::endl
            << "pips_DC_z[pipsIdx="<<pipsIdx<<"][regionIdx="<<regionIdx<<"]="
            << pips_DC_z[pipsIdx][regionIdx]
            << std::endl;
        }
    }
    
    piplus_Px[pipsIdx]          = piplus[pipsIdx].Px();
    piplus_Py[pipsIdx]          = piplus[pipsIdx].Py();
    piplus_Pz[pipsIdx]          = piplus[pipsIdx].Pz();
    piplus_E[pipsIdx]           = piplus[pipsIdx].E();
    Vpiplus_X[pipsIdx]          = Vpiplus[pipsIdx].X();
    Vpiplus_Y[pipsIdx]          = Vpiplus[pipsIdx].Y();
    Vpiplus_Z[pipsIdx]          = Vpiplus[pipsIdx].Z();
    
//    if (evnum>67 && evnum<74)
//        std::cout << "event " << evnum << " before CheckIfPionPassedSelectionCuts( pipsIdx=" << pipsIdx << ")" << std::endl;

    
    pipsPastSelectionCuts[pipsIdx] = CheckIfPionPassedSelectionCuts("pi+",
                                                                     pips_DC_sector[pipsIdx],
                                                                     pips_DC_x[pipsIdx],
                                                                     pips_DC_y[pipsIdx],
                                                                     pips_DC_z[pipsIdx],
                                                                     pips_chi2PID[pipsIdx],  piplus[pipsIdx].P(),
                                                                     Ve,
                                                                     Vpiplus[pipsIdx],
                                                                    pipsIdx,
                                                                     fdebug);
//    if (evnum>67 && evnum<74)
//        std::cout << "event " << evnum << " after CheckIfPionPassedSelectionCuts( pipsIdx=" << pipsIdx << ")" << std::endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExtractPimsInformation( int pimsIdx, int fdebug ){
    // extract negative pion information
    SetLorentzVector(piminus[pimsIdx]  ,piminuses[pimsIdx]);
    Zpims[pimsIdx]              = piminus[pimsIdx].E() / omega;
    Vpiminus[pimsIdx]           = GetParticleVertex( piminuses[pimsIdx] );
    pims_chi2PID[pimsIdx]       = piminuses[pimsIdx]->par()->getChi2Pid();
    
    // EC in and out
    pims_E_ECIN[pimsIdx]        = piminuses[pimsIdx]->cal(ECIN)->getEnergy();
    pims_E_ECOUT[pimsIdx]       = piminuses[pimsIdx]->cal(ECOUT)->getEnergy();
    // PCAL
    auto pims_PCAL_info         = piminuses[pimsIdx]->cal(PCAL);
    pims_E_PCAL[pimsIdx]        = pims_PCAL_info->getEnergy();
    pims_PCAL_sector[pimsIdx]   = pims_PCAL_info->getSector();
    pims_PCAL_V[pimsIdx]        = pims_PCAL_info->getLv();
    pims_PCAL_W[pimsIdx]        = pims_PCAL_info->getLw();
    pims_PCAL_x[pimsIdx]        = pims_PCAL_info->getX();
    pims_PCAL_y[pimsIdx]        = pims_PCAL_info->getY();
    pims_PCAL_z[pimsIdx]        = pims_PCAL_info->getZ();
    // DC
    auto pims_DC_info           = piminuses[pimsIdx]->trk(DC);
    pims_DC_sector[pimsIdx]     = pims_DC_info->getSector(); // tracking sector
    pims_Chi2N[pimsIdx]         = pims_DC_info->getChi2N();  // tracking chi^2/NDF
    for (int regionIdx=0; regionIdx<3; regionIdx++) {
        DC_layer = DC_layers[regionIdx];
        pims_DC_x[pimsIdx][regionIdx] = piminuses[pimsIdx]->traj(DC,DC_layer)->getX();
        pims_DC_y[pimsIdx][regionIdx] = piminuses[pimsIdx]->traj(DC,DC_layer)->getY();
        pims_DC_z[pimsIdx][regionIdx] = piminuses[pimsIdx]->traj(DC,DC_layer)->getZ();
    }

    piminus_Px[pimsIdx]          = piminus[pimsIdx].Px();
    piminus_Py[pimsIdx]          = piminus[pimsIdx].Py();
    piminus_Pz[pimsIdx]          = piminus[pimsIdx].Pz();
    piminus_E[pimsIdx]           = piminus[pimsIdx].E();
    Vpiminus_X[pimsIdx]          = Vpiminus[pimsIdx].X();
    Vpiminus_Y[pimsIdx]          = Vpiminus[pimsIdx].Y();
    Vpiminus_Z[pimsIdx]          = Vpiminus[pimsIdx].Z();
    
    pimsPastSelectionCuts[pimsIdx] = CheckIfPionPassedSelectionCuts("pi-",
                                                                     pims_DC_sector[pimsIdx],
                                                                     pims_DC_x[pimsIdx],
                                                                     pims_DC_y[pimsIdx],
                                                                     pims_DC_z[pimsIdx],
                                                                     pims_chi2PID[pimsIdx],  piminus[pimsIdx].P(),
                                                                     Ve,
                                                                     Vpiminus[pimsIdx],
                                                                    pimsIdx,
                                                                     fdebug);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GetParticlesByType (int evnum, int fdebug){
    // get particles by type
    Ne      = electrons .size();
    Nn      = neutrons  .size();
    Np      = protons   .size();
    Npips   = pipluses  .size();
    Npims   = piminuses .size();
    Ngammas = gammas    .size();
    Nd      = deuterons.size();
    if (fdebug>2){
        std::cout
        << "particles in event "            << evnum        << " : "
        << "N(electrons): "                 << Ne           <<  ","
        << "N(protons): "                   << Np           <<  ","
        << "N(neutrons): "                  << Nn           <<  ","
        << "N(pi+): "                       << Npips        <<  ","
        << "N(pi-): "                       << Npims        <<  ","
        << "N(gammas): "                    << Ngammas      <<  ","
        << "N(deuterons): "                 << Nd           <<  ","
        << std::endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Stream_e_pi_line_to_CSV( int piIdx, int fdebug ){ // write a row of pion number piIdx to CSV file
    
    
    TLorentzVector  pi;
    TVector3        Vpi;
    double          Zpi;
    int    pi_DC_sector;
    
//    std::cout << "Stream_e_pi_line_to_CSV(evnum="<<evnum<<"): "<< "Npips: " << Npips << ", piIdx: " << piIdx << std::endl;
    if (PiCharge=="pips") {
        pi  = piplus [piIdx];
        Vpi = Vpiplus[piIdx];
        
        if ( Npips>piIdx )  pi_reconstructed = true;
        else                pi_reconstructed = false;
        
        pi_passed_cuts              = pipsPastSelectionCuts[piIdx];
        pi_passed_fiducial_cuts     = pipsPastFiducialCuts[piIdx];
//        std::cout << "if (PiCharge==pips) { "<< "pi_passed_fiducial_cuts: " << pi_passed_fiducial_cuts << std::endl;
        pi_passed_PID_cuts          = pipsPastPIDCuts[piIdx];
        pi_DC_sector                = pips_DC_sector[piIdx];
    }
    else if (PiCharge=="pims") {
        
        pi  = piminus [piIdx];
        Vpi = Vpiminus[piIdx];
        
        if ( Npims>piIdx )  pi_reconstructed = true;
        else                pi_reconstructed = false;
        
        pi_passed_cuts              = pimsPastSelectionCuts[piIdx];
        pi_passed_fiducial_cuts     = pimsPastFiducialCuts[piIdx];
        pi_passed_PID_cuts          = pimsPastPIDCuts[piIdx];
        pi_DC_sector                = pims_DC_sector[piIdx];
   }
    else {
        std::cout << "pion charge undefined at Stream_e_pi_line_to_CSV(), returning " << std::endl;
        return;
    }
    // write a (e,e'pi) event-line to CSV file
//    std::cout << "variables: "<< "pi_passed_fiducial_cuts: " << pi_passed_fiducial_cuts << std::endl;
    std::vector<double> variables =
    {   e.P(),              e.Theta(),          e.Phi(),                Ve.Z(),
        pi.P(),             pi.Theta(),         pi.Phi(),               Vpi.Z(),
        (double)Npips,      (double)Npims,      (double)Ne,             (double)Ngammas,
        (double)Np,         (double)Nn,         (double)Nd,
        e_g.P(),            e_g.Theta(),        e_g.Phi(),              Ve_g.Z(),
        pi_g.P(),           pi_g.Theta(),       pi_g.Phi(),             Vpi_g.Z(),
        (double)pi_reconstructed,               (double)pi_passed_cuts,
        (double)pi_passed_fiducial_cuts,        (double)pi_passed_PID_cuts,
        (double)e_reconstructed,                (double)e_passed_cuts,
        (double)e_DC_sector,                    (double)pi_DC_sector
    };
    StreamToCSVfile( variables, fdebug );
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
bool CheckIfPionPassedSelectionCuts(TString pionCharge, // "pi+" or "pi-"
                                    Double_t DC_sector,
                                    Double_t DC_x[3], Double_t DC_y[3], Double_t DC_z[3],
                                    Double_t chi2PID, Double_t p,
                                    TVector3 Ve,
                                    TVector3 Vpi,
                                    int piIdx,
                                    int fdebug){
    
    // decide if pion (pi+ or pi-) passed event selection cuts
    //
    // input:
    // --------
    // DC_x, DC_y   pi drift-chamber coordinates
    // chi2PID      pi chi2PID     (pips_chi2PID)
    // p            pi momentum    (pi.P())
    //
    // comments
    // ---------------
    // DC - fiducial cuts on DC
    if (fdebug>3) {
        std::cout << "CheckIfPionPassedSelectionCuts()" << std::endl;
    }
    if (DC_sector == 0) { if (fdebug>2){std::cout << "DC_sector=0 (funny...)" << std::endl;} return false;}
    
    int PDGcode;
    double    C;
    if (pionCharge=="pi+"){
        PDGcode = 211;
        C       = 0.88;
    } else if (pionCharge=="pi-") {
        PDGcode = -211;
        C       = 0.93;
    } else {
        std::cout << "pion charge is not defined in CheckIfPionPassedSelectionCuts(), returning false" << std::endl;
        if (pionCharge=="pi+"){
            pipsPastPIDCuts[piIdx] = false;
        } else {
            pimsPastPIDCuts[piIdx] = false;
        }
    }

    bool DCFidRegion[3];
    for (int regionIdx=0; regionIdx<3; regionIdx++) {
        // DC_e_fid:
        // sector:  1-6
        // layer:   1-3
        // bending: 0(out)/1(in)
        int bending  = 1 ? (torusBending==-1) : 0;
        // new version Aug-11,2021
        if (fdebug>3) {
            std::cout << "dcfid.DC_fid_th_ph_sidis(): "
            << DC_x[regionIdx] <<     ","
            << DC_y[regionIdx] <<     ","
            << DC_z[regionIdx] <<     ","
            << DC_sector       <<     ","
            << regionIdx+1     <<     ","
            << bending         <<     ","
            << std::endl;
            
        }
        bool DC_fid  = dcfid.DC_fid_th_ph_sidis(PDGcode,            // particle PID
                                                DC_x[regionIdx],    // x
                                                DC_y[regionIdx],    // y
                                                DC_z[regionIdx],    // z
                                                DC_sector,          // sector
                                                regionIdx+1,        // layer
                                                bending);           // torus bending
        DCFidRegion[regionIdx] = DC_fid;
    }
    if (pionCharge=="pi+"){
        pipsPastFiducialCuts[piIdx] = DCFidRegion[0] && DCFidRegion[1] && DCFidRegion[2];
//        pipsPastFiducialCuts[piIdx] = true;
//        if (evnum>67 && evnum<74)
//            std::cout << "evnum " << evnum << " piIdx " << piIdx << ", pipsPastFiducialCuts[piIdx]: " << pipsPastFiducialCuts[piIdx] << std::endl;
        
    } else {
        pimsPastFiducialCuts[piIdx] = DCFidRegion[0] && DCFidRegion[1] && DCFidRegion[2];
    }

//    if (evnum>67 && evnum<74)
//        std::cout << "after if in evnum " << evnum << " piIdx " << piIdx << ", pipsPastFiducialCuts[piIdx]: " << pipsPastFiducialCuts[piIdx] << std::endl;
    
    
    if (fdebug>3) {
        std::cout << "in CheckIfPionPassedSelectionCuts()"<< std::endl
        << "pion charge: "          << pionCharge               << ","
        << "DC_x[0]: "              << DC_x[0]                  << ","
        << "chi2PID:"               << chi2PID                  << ","
        << "Chi2PID_pion_lowerBound( p="<<p<<", C="<<C<<" ): "  << Chi2PID_pion_lowerBound( p, C ) << ","
        << "Chi2PID_pion_upperBound( p="<<p<<", C="<<C<<" ): "  << Chi2PID_pion_upperBound( p, C ) << ","
        << "fabs((Ve-Vpi).Z()): "   << fabs((Ve-Vpi).Z())       << ","
        << std::endl;
    }
    if(!
       // pion Identification Refinement - chi2PID vs. momentum
       (( Chi2PID_pion_lowerBound( p, C ) < chi2PID && chi2PID < Chi2PID_pion_upperBound( p , C ) )
       
       // Cut on the z-Vertex Difference Between Electrons and Hadrons.
       &&  ( fabs((Ve-Vpi).Z()) < 20. )
       )) {
        if (pionCharge=="pi+"){
            pipsPastPIDCuts[piIdx] = false;
        } else {
            pimsPastPIDCuts[piIdx] = false;
        }
    } else {
        if (pionCharge=="pi+"){
            pipsPastPIDCuts[piIdx] = true;
        } else {
            pimsPastPIDCuts[piIdx] = true;
        }
    }
    
    if (fdebug>3) {
        std::cout
        << "succesfully passed CheckIfPionPassedSelectionCuts(), return true"
        << std::endl;
    }
//    if (evnum>67 && evnum<74)
//        std::cout << "before return in evnum " << evnum << " piIdx " << piIdx << ", pipsPastFiducialCuts[piIdx]: " << pipsPastFiducialCuts[piIdx] << std::endl;

    
    if (pionCharge=="pi+"){
        
        pipsPastSelectionCuts[piIdx] = (pipsPastPIDCuts[piIdx] && pipsPastFiducialCuts[piIdx]);
        return pipsPastSelectionCuts[piIdx];
        
    } else {
        
        pimsPastSelectionCuts[piIdx] = (pimsPastPIDCuts[piIdx] && pimsPastFiducialCuts[piIdx]);
        return pimsPastSelectionCuts[piIdx];
        
    }
}





// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
bool CheckIfElectronPassedSelectionCuts(Double_t e_PCAL_x, Double_t e_PCAL_y,
                                     Double_t e_PCAL_W, Double_t e_PCAL_V,
                                     Double_t e_E_PCAL,
                                     Double_t e_E_ECIN, Double_t e_E_ECOUT,
                                     TLorentzVector e,
                                     TVector3 Ve,
                                     Double_t e_DC_sector,
                                     Double_t e_DC_x[3],
                                     Double_t e_DC_y[3],
                                     Double_t e_DC_z[3],
                                     int torusBending){
    
    // decide if electron in event passes event selection cuts
    
    // DC - fiducial cuts on DC
    // from bandsoft_tools/skimmers/electrons.cpp,
    // where eHit.getDC_x1() - x position in first region of the drift chamber
    // same for y1,x2,y2,...
    // eHit.getDC_sector() - sector
    // checking DC Fiducials
    // torusBending         torus magnet bending:   ( 1 = inbeding, -1 = outbending    )
    
    // sometimes the readout-sector is 0. This is funny
    // Justin B. Estee (June-21): I also had this issue. I am throwing away sector 0. The way you check is plot the (x,y) coordinates of the sector and you will not see any thing. Double check me but I think it is 0.
    if (e_DC_sector == 0) return false;

    for (int regionIdx=0; regionIdx<3; regionIdx++) {
        // DC_e_fid:
        // sector:  1-6
        // layer:   1-3
        // bending: 0(out)/1(in)
        // std::cout << "e_DC_sector: " << e_DC_sector << ", regionIdx: " << regionIdx << std::endl;
        int bending  = 1 ? (torusBending==-1) : 0;
        bool DC_fid  = dcfid.DC_fid_xy_sidis(11,                 // particle PID,
                                             e_DC_x[regionIdx],  // x
                                             e_DC_y[regionIdx],  // y
                                             e_DC_sector,        // sector
                                             regionIdx+1,        // layer
                                             bending);           // torus bending
        if (DC_fid == false) {
            return false;
        }
    }
    
    double              cutValue_e_PCAL_W = 19.0;
    double              cutValue_e_PCAL_V = 19.0;
    double              cutValue_e_E_PCAL = 0.07;
    double  cutValue_SamplingFraction_min = 0.17;
    double      cutValue_PCAL_ECIN_SF_min = 0.20;
    double      cutValue_Vz_min, cutValue_Vz_max;
    
    if (torusBending==-1){ // in-bending torus field
        // Spring 19 and Spring 2020 in-bending.
        cutValue_Vz_min = -13.0;
        cutValue_Vz_max = 12.0;
    } else if (torusBending==1){ // Out-bending torus field
        // Fall 2019 (without low-energy-run) was out-bending.
        cutValue_Vz_min = -18.0;
        cutValue_Vz_max = 10.0;
    }
    
    if(!(true
       // fiducial cuts on PCAL
       //fabs(e_PCAL_x)>0
       //&&  fabs(e_PCAL_y)>0
       &&  e_PCAL_W > cutValue_e_PCAL_W
       &&  e_PCAL_V > cutValue_e_PCAL_V
       
       // Electron Identification Refinement  - PCAL Minimum Energy Deposition Cut
       &&  e_E_PCAL > cutValue_e_E_PCAL
       
       // Sampling fraction cut
       && ((e_E_PCAL + e_E_ECIN + e_E_ECOUT)/e.P()) > cutValue_SamplingFraction_min
       && (e_E_ECIN/e.P() > cutValue_PCAL_ECIN_SF_min - e_E_PCAL/e.P()) // RGA AN puts "<" here mistakenly
       
       // Cut on z-vertex position: in-bending torus field -13.0 cm < Vz < +12.0 cm
       // Spring 19 and Spring 2020 in-bending.
       // Fall 2019 (without low-energy-run) was out-bending.
       &&  ((cutValue_Vz_min < Ve.Z()) && (Ve.Z() < cutValue_Vz_max))
       )) return false;
    
    return true;
}
// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
Double_t Chi2PID_pion_lowerBound( Double_t p, Double_t C){
    // compute lower bound for chi2PID for a pi+
    // based on RGA_Analysis_Overview_and_Procedures_Nov_4_2020-6245173-2020-12-09-v3.pdf
    // p. 75
    // "Strict cut"
    //
    // input:
    // -------
    // p        pion momentum
    // C        is the scaling factor for sigma (away from the mean value of the pion distribution)
    //
    
    return ( -C * 3 );
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
Double_t Chi2PID_pion_upperBound( Double_t p, Double_t C){
    // compute upper bound for chi2PID for a pi+
    // based on RGA_Analysis_Overview_and_Procedures_Nov_4_2020-6245173-2020-12-09-v3.pdf
    // p. 75
    // "Strict cut"
    //
    // input:
    // -------
    // p        pion momentum
    // C        is the scaling factor for sigma (away from the mean value of the pion distribution)
    //
    
    if (p<2.44)
        
        return C*3;
    
    else if (p<4.6)
        
        return C*( 0.00869 + 14.98587*exp(-p/1.18236)+1.81751*exp(-p/4.86394) ) ;
    
    else
        
        return C*( -1.14099 + 24.14992*exp(-p/1.36554) + 2.66876*exp(-p/6.80552) );
    
}
