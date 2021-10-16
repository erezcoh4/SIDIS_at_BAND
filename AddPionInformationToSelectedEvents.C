
#define NMAXPIONS 20 // maximal allowed number of pions
#define r2d 180./3.1415 // radians to degrees


#include <cstdlib>
#include <sstream>
#include <fstream>
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
#include "Auxiliary/DCfid_SIDIS.cpp"
#include "Auxiliary/SIDISatBAND_auxiliary.cpp"
#include "clas12reader.h"
using namespace clas12;

// start clock
auto start = std::chrono::high_resolution_clock::now();
struct event_id { int run_number, event_number ; } ;

TString             DataPath = "/volatile/clas12/rg-b/production/recon/spring2019/torus-1/pass1/v0/dst/train_20200610/inc/";
TString        EventListPath = "macros/eventlist/";
TString        EventListName = "example_eventlist";
TString          OutDataPath = ""; // or e.g. /Users/erezcohen/Desktop/Software/CLAS12/
TString          OutFilename = "NonameOutputFile";
TString      OutFullFilename = "";
TString            csvheader = "Npips,Npims,Nelectrons,Ngammas,Nprotons,Nneutrons,Ndeuterons,";

bool                      FoundEvent;
Int_t              NeventsInList = 0;
Int_t                  RunNumber = 0;
Int_t                EventNumber = 0;
Int_t                     fdebug = 0;
Int_t           NeventsProcessed = 0;
Int_t                 Ne, Nn, Np, Nd;
Int_t          Npips, Npims, Ngammas;
Double_t                Ebeam = 10.2;
TLorentzVector            Beam, e, q;
TVector3                          Ve;

std::vector<event_id>      eventlist;
std::ofstream             outcsvfile;

// a list of run numbers
std::vector<int>                run_numbers;
// and for each run, a list of event numbers in that run
std::vector<std::vector<int>>   event_numbers;


int          DC_layers[3] = {6,18,36};// Region 1 is denoted at DC detector 6, Region 2 is denoted 18, Region 3 - as 36
int                          DC_layer;

// positive pions
bool              pipsPastCutsInEvent;
bool pipsPastSelectionCuts[NMAXPIONS];
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

std::vector<TLorentzVector>         piplus;
std::vector<TVector3>              Vpiplus;


// negative pions
bool              pimsPastCutsInEvent;
bool pimsPastSelectionCuts[NMAXPIONS];
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
double               Zpims[NMAXPIONS]; // hadron rest-frame energy=

std::vector<TLorentzVector>        piminus;
std::vector<TVector3>             Vpiminus;


// auxiliary
DCfid_SIDIS             dcfid;
SIDISatBAND_auxiliary   aux;

auto                                 db = TDatabasePDG::Instance();
std::vector<region_part_ptr>        electrons, pipluses, piminuses;
std::vector<region_part_ptr>  neutrons, protons, gammas, deuterons;




// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
// Setters
void            SetDataPath (TString fDataPath)       { DataPath = fDataPath; };
void     SetOutFullFilename ()                        { OutFullFilename = OutDataPath + OutFilename + ".csv"; };
void           SetVerbosity (Int_t _fdebug_)          { fdebug = _fdebug_;    };
void               SetEbeam (Double_t _Ebeam_)        { Ebeam = _Ebeam_;    Beam.SetPxPyPzE(0, 0, Ebeam, Ebeam );};
void         SetOutDataPath (TString fOutDataPath)    {
    OutDataPath = fOutDataPath + "/";
    SetOutFullFilename();
};
void         SetOutFilename ()                        {
    OutFilename = EventListName + "_PionInformation";
    SetOutFullFilename();
};
void       SetEventListName (TString fEventListName)  {
    EventListName = fEventListName;
    SetOutFilename();
};


// methods implemented below
void                       PrintEventList ();
void                        ReadEventList ();
void                  AssignPionsToEvents (Int_t NeventsMax);
void                     CloseOutputFiles ();
void                        FinishProgram ();
void                  InitializeVariables ();
void                   WriteEventToOutput ();
void              ExtractPionsInformation ();
void               ExtractPipsInformation (int pipsIdx);
void               ExtractPimsInformation (int pimsIdx);
bool       CheckIfPionPassedSelectionCuts (TString pionCharge, // "pi+" or "pi-"
                                           Double_t DC_sector,
                                           Double_t DC_x[3], Double_t DC_y[3], Double_t DC_z[3],
                                           Double_t chi2PID, Double_t p,
                                           TVector3 Ve,      TVector3 Vpi);
void                      OpenOutputFiles ();
int              FindRunNumberIndexInList ( int run_number_to_find );
void          ExtractElectronInformation  ();

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void AddPionInformationToSelectedEvents(Int_t      NeventsMax=-1,
                                        TString   _EventListName_="example_eventlist",
                                        Int_t     _fdebug_=1,
                                        TString   _OutDataPath_="/volatile/clas12/users/ecohen/BAND/PionInformationInEventLists/"){
    
    SetVerbosity        (_fdebug_);
    aux.SetVerbosity    (_fdebug_);
    aux.loadCutValues   ("cutValues.csv");
    SetEventListName    (_EventListName_);
    SetOutDataPath      (_OutDataPath_);
    OpenOutputFiles     ();
    ReadEventList       ();
    AssignPionsToEvents ( NeventsMax );
    FinishProgram       ();
}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void AssignPionsToEvents(Int_t NeventsMax){
    if (fdebug>2) {std::cout << "AssignPionsToEvents("<<NeventsMax<<")"<< std::endl;}

    if (NeventsMax<0) NeventsMax = NeventsInList;
    
    int eventIdx = 0;
    for (int run_idx=0; run_idx<(int)run_numbers.size(); run_idx++){
        auto RunNumber = run_numbers.at(run_idx);
        if (fdebug>3) std::cout << "Processing run " <<  RunNumber << std::endl ;
        
        // every time that the run number changes, we open a new hipo file
        TString    inputFile = DataPath + "inc_" + aux.GetRunNumberSTR( RunNumber ) + ".hipo";
        TChain fake("hipo");
        fake.Add(inputFile.Data());
        
        //get the hipo data
        if (fdebug>3) std::cout << "Reading hipo file " << fake.GetListOfFiles()->At(0)->GetTitle() << std::endl;
        clas12reader c12(fake.GetListOfFiles()->At(0)->GetTitle(),{0});
        
        
        for (int evt_idx=0; evt_idx<(int)event_numbers.at(run_idx).size(); evt_idx++){
            auto EventNumber = event_numbers.at(run_idx).at(evt_idx);
            if (fdebug>3) std::cout << "Grabbing event " << EventNumber << " from run " << RunNumber << std::endl;
            
            InitializeVariables();
            while(c12.next()){
                
                if (fdebug>3) std::cout << "(c12.runconfig()->getEvent()==EventNumber): " << c12.runconfig()->getEvent()==EventNumber << std::endl;
                
                if (c12.runconfig()->getEvent()==EventNumber){
                    FoundEvent = true;
                    
                    if (fdebug>3) std::cout << "Found event " << EventNumber << " in run " << RunNumber << std::endl;
                    pipluses    = c12.getByID( 211  );          Npips   = pipluses  .size();
                    piminuses   = c12.getByID(-211  );          Npims   = piminuses .size();
                    electrons   = c12.getByID( 11   );          Ne      = electrons .size();

                    if (Ne>0 && (Npips>0 || Npims>0)){
                        if (fdebug>3) {
                            std::cout
                            << Ne << " e, "
                            << Npips << " pi+ & "
                            << Npims << " pi- in event "
                            << EventNumber << std::endl;
                        }
                        ExtractElectronInformation  ();
                        ExtractPionsInformation     ();
                        WriteEventToOutput          ();
                    } else {
                        if(fdebug>3) {
                            if (Ne==0) std::cout << "no electrons in event " << EventNumber << std::endl;
                            else std::cout << "no pions in event " << EventNumber << std::endl;
                        }
                    }
                    if ( eventIdx%10==0 ) std::cout << eventIdx << "/" <<  NeventsMax << std::endl;
                    eventIdx++;
                    break;
                }
            }
            if(FoundEvent == false ){
                if(fdebug>3) std::cout << "Did not find event " << EventNumber << " in run " << RunNumber << std::endl;
            }
            if (fdebug>3)  std::cout << "-------------------------------------------------------------------" << std::endl;
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void WriteEventToOutput(){

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void InitializeVariables(){
    FoundEvent = false;
    e  = TLorentzVector(0,0,0,db->GetParticle( 11 )->Mass());
    Ve = TVector3(-9999,-9999,-9999);

    piplus      .clear();
    piminus     .clear();
    Vpiplus     .clear();
    Vpiminus    .clear();
    pipluses    .clear();
    piminuses   .clear();
    pipsPastCutsInEvent = pimsPastCutsInEvent       = false;
    DC_layer                                        = -9999;
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
        pipsPastSelectionCuts[piIdx]                = false;

        piminus .push_back( TLorentzVector(0,0,0,db->GetParticle( -211 )->Mass()) );
        Vpiminus.push_back( TVector3() );
        pimsPastSelectionCuts[piIdx]                = false;

    }
    DC_layer                                        = -9999;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SetParticle4Momentum (TLorentzVector &p4,clas12::region_part_ptr rp){
    p4.SetXYZM(rp->par()->getPx(),
               rp->par()->getPy(),
               rp->par()->getPz(),
               p4.M());
}
void ExtractElectronInformation(){
    if (fdebug>3) std::cout << "ExtractElectronInformation()" << std::endl;
    // ------------------------------------------------------------------------------------------------
    // extract information from first electron
    // here only the reconstructed vertext position, as it affects the pion vertex cut position
    // ------------------------------------------------------------------------------------------------
    if (fdebug>3) std::cout << "electrons[0]->par()->getPx(): " << electrons[0]->par()->getPx() << " GeV/c" << std::endl;
    
    SetParticle4Momentum( e , electrons[0] );
    if (fdebug>3) std::cout << "e.E(): " << e.E() << " GeV/c" << std::endl;
    q = Beam - e;
    if (fdebug>3) std::cout << "q.E(): " << q.E() << " GeV/c" << std::endl;
    Ve = aux.GetParticleVertex( electrons[0] );
    if (fdebug>3) std::cout << "Ve.X(): " << Ve.X() << " cm" << std::endl;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExtractPionsInformation(){
    if (fdebug>3) std::cout << "ExtractPionsInformation()" << std::endl;

    // positive pions)
    for (int pipsIdx=0; pipsIdx < Npips; pipsIdx++) {
        ExtractPipsInformation( pipsIdx );
    }
    // negative pions
    for (int pimsIdx=0; pimsIdx < Npims; pimsIdx++) {
        ExtractPimsInformation( pimsIdx );
    }

    // done
    if (fdebug > 2) std::cout << "done extracting pion information" << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExtractPipsInformation( int pipsIdx ){
    if (fdebug>2)
        std::cout << "ExtractPipsInformation( pipsIdx=" << pipsIdx << ", fdebug=" << fdebug << " )" << std::endl;


    // extract positive pion information
    aux.SetParticle4Momentum(piplus[pipsIdx]  ,pipluses[pipsIdx]);
    Zpips[pipsIdx]              = piplus[pipsIdx].E() / q.E();
    Vpiplus[pipsIdx]            = aux.GetParticleVertex( pipluses[pipsIdx] );
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
    // ------------------------------------------------------------------------------------------------
    // now, check if pion passed event selection requirements
    // ------------------------------------------------------------------------------------------------
    pipsPastSelectionCuts[pipsIdx] = CheckIfPionPassedSelectionCuts("pi+",
                                                                    pips_DC_sector[pipsIdx],
                                                                    pips_DC_x[pipsIdx],
                                                                    pips_DC_y[pipsIdx],
                                                                    pips_DC_z[pipsIdx],
                                                                    pips_chi2PID[pipsIdx],  piplus[pipsIdx].P(),
                                                                    Ve,
                                                                    Vpiplus[pipsIdx]);
    if (pipsPastSelectionCuts[pipsIdx]) {
        pipsPastCutsInEvent = true;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExtractPimsInformation( int pimsIdx ){
    // extract negative pion information
    aux.SetParticle4Momentum(piminus[pimsIdx]  ,piminuses[pimsIdx]);
    Zpims[pimsIdx]              = piminus[pimsIdx].E() / q.E();
    Vpiminus[pimsIdx]           = aux.GetParticleVertex( piminuses[pimsIdx] );
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
    // ------------------------------------------------------------------------------------------------
    // now, check if pion passed event selection requirements
    // ------------------------------------------------------------------------------------------------
    pimsPastSelectionCuts[pimsIdx] = CheckIfPionPassedSelectionCuts("pi-",
                                                                    pims_DC_sector[pimsIdx],
                                                                    pims_DC_x[pimsIdx],
                                                                    pims_DC_y[pimsIdx],
                                                                    pims_DC_z[pimsIdx],
                                                                    pims_chi2PID[pimsIdx],  piminus[pimsIdx].P(),
                                                                    Ve,
                                                                    Vpiminus[pimsIdx]);
    if (pimsPastSelectionCuts[pimsIdx]) {
        pimsPastCutsInEvent = true;
    }
}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
bool CheckIfPionPassedSelectionCuts(TString pionCharge, // "pi+" or "pi-"
                                    Double_t DC_sector,
                                    Double_t DC_x[3], Double_t DC_y[3], Double_t DC_z[3],
                                    Double_t chi2PID, Double_t p,
                                    TVector3 Ve,
                                    TVector3 Vpi){

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
        return false;
    }

    for (int regionIdx=0; regionIdx<3; regionIdx++) {
        // DC_e_fid:
        // sector:  1-6
        // layer:   1-3
        // bending: 0(out)/1(in)
        int bending  = 1 ? (aux.torusBending==-1) : 0;
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
        if (DC_fid == false) {
            return false;
        }
    }

    if (fdebug>3) {
        std::cout << "in CheckIfPionPassedSelectionCuts()"<< std::endl
        << "pion charge: "          << pionCharge               << ","
        << "DC_x[0]: "              << DC_x[0]                  << ","
        << "chi2PID:"               << chi2PID                  << ","
        << "Chi2PID_pion_lowerBound( p="<<p<<", C="<<C<<" ): "  << aux.Chi2PID_pion_lowerBound( p, C ) << ","
        << "Chi2PID_pion_upperBound( p="<<p<<", C="<<C<<" ): "  << aux.Chi2PID_pion_upperBound( p, C ) << ","
        << "fabs((Ve-Vpi).Z()): "   << fabs((Ve-Vpi).Z())       << ","
        << std::endl;
    }
    if(
       // pi+ Identification Refinement - chi2PID vs. momentum
       ( aux.Chi2PID_pion_lowerBound( p, C ) < chi2PID && chi2PID < aux.Chi2PID_pion_upperBound( p , C ) )

       // Cut on the z-Vertex Difference Between Electrons and Hadrons.
       &&  ( fabs((Ve-Vpi).Z()) < aux.cutValue_Ve_Vpi_dz_max )
       ) {
           if (fdebug>3) { std::cout << "succesfully passed CheckIfPionPassedSelectionCuts(), return true" << std::endl; }

           return true;
       }
    return false;
}





// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
std::istream &operator >>(std::istream &ist, event_id &event) {
    char comma ;
    ist
    >> event.run_number   >> comma
    >> event.event_number >> comma;
    return ist ;
}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void ReadEventList(){
    std::ifstream data(EventListPath + EventListName + ".csv" );
    std::string   line;
    std::getline(data,line); // kill header line
    while(std::getline(data,line)) {
        event_id event_entry ;
        std::stringstream lineStream(line);
        lineStream >> event_entry ;
        eventlist.push_back( event_entry ) ;
        
        int   RunNumber = event_entry.run_number;
        int EventNumber = event_entry.event_number;
        int RunNumberIdxInList = FindRunNumberIndexInList( RunNumber );
        if (RunNumberIdxInList==-1) {
            // this means that the run number is not in our list - its a new run number
            run_numbers  .push_back( RunNumber );
            // initialize the event-number sublist of events in this run
            event_numbers.push_back( std::vector<int>{} );
            event_numbers.at( run_numbers.size()-1 ).push_back( EventNumber );
        } else {
            // this means that the run number is in our list,
            // and we can add the event number to the run sublist of events
            event_numbers.at( RunNumberIdxInList ).push_back( EventNumber );
        }
        NeventsInList++;
    }
    
    if (fdebug>2) PrintEventList();
};

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
int FindRunNumberIndexInList( int run_number_to_find ){
    // find the index of "run_number_to_find" in the list "run_numbers"
    for (int i=0; i<(int)run_numbers.size(); i++){
        int run_number = run_numbers.at(i);
        if ( run_number==run_number_to_find ) return i;
    }
    return -1;
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void PrintEventList(){
    if (fdebug>5){
        std::cout
        << "Event list to process"
        << std::endl
        << "run,event"
        << std::endl;
        
        for (auto event: eventlist){
            std::cout << event.run_number << "," << event.event_number << std::endl;
        }
    }
    
    std::cout
    << "Events to process in specific runs"
    << std::endl;
    for (int run_idx=0; run_idx<(int)run_numbers.size(); run_idx++){
        std::cout << "Run " << run_numbers.at(run_idx) << ": " ;
        
        for (int evt_idx=0; evt_idx<(int)event_numbers.at(run_idx).size(); evt_idx++){
            std::cout << event_numbers.at(run_idx).at(evt_idx) << "," ;
        }
        std::cout << std::endl;
    }
    
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void OpenOutputFiles (){
    if (fdebug>2) std::cout << "Opening output csvfile: " << OutFullFilename << std::endl;

    outcsvfile.open( OutFullFilename );
    if (fdebug>2) std::cout << "Writing to output csvfile: " << OutFullFilename << std::endl;

    outcsvfile << csvheader << std::endl;

    if (fdebug>2) std::cout << "Opened output csvfile: " << OutFullFilename << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void FinishProgram(){
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;

    if (fdebug>2) std::cout << "Done. Elapsed time: " << elapsed.count() << std::endl;
    CloseOutputFiles( );
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void CloseOutputFiles (){
    // close output CSV
    if (fdebug>2) std::cout << "Closing output csvfile " << OutFullFilename << std::endl;

    outcsvfile.close();
    if (fdebug>2){
        std::cout
        << "Done processesing "  <<  NeventsProcessed  << " events,"
        << std::endl
        << "See results in "
        << std::endl
        << OutFullFilename
        << std::endl;
    }
}
