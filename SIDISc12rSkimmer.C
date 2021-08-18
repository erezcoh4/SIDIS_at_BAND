// last edit Aug-18, 2021 (EOC, mbp), see README
// ToDo:
// (1) add multiple pions per event
// (2) Check beam helicity
// (3) Add pion DC fiducial cuts


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
#include "/u/home/cohen/SIDIS_at_BAND/Auxiliary/DCfid_SIDIS.cpp"
#include "/u/home/cohen/SIDIS_at_BAND/Auxiliary/csv_reader.h"
#define NMAXPIONS 20 // maximal allowed number of pions
using namespace clas12;



// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
// start clock
auto start = std::chrono::high_resolution_clock::now();
// declare methods
TVector3                GetParticleVertex (clas12::region_part_ptr rp);
void                     SetLorentzVector (TLorentzVector &p4, clas12::region_part_ptr rp);
void                      OpenOutputFiles (TString csvfilename, TString header);
void                     CloseOutputFiles (TString OutDataPath, TString outfilename);
void                      StreamToCSVfile (TString pionCharge, // "pi+" or "pi-"
                                           std::vector<Double_t> observables,
                                           bool IsSelectedEvent,
                                           int fdebug);
void                       printCutValues ();
void                        loadCutValues (TString cutValuesFilename = "cutValues.csv", int fdebug=0);
void                      SetOutputTTrees ();
double                       FindCutValue ( std::string cutName );
bool      ElectronPassedSelectionCriteria (Double_t e_PCAL_x, Double_t e_PCAL_y,
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
bool      PionPassedSelectionCutsCriteria (TString pionCharge, // "pi+" or "pi-"
                                           Double_t DC_sector,
                                           Double_t DC_x[3], Double_t DC_y[3], Double_t DC_z[3],
                                           Double_t chi2PID, Double_t p,
                                           TVector3 Ve,      TVector3 Vpi,
                                           int fdebug);
Double_t          Chi2PID_pion_lowerBound (Double_t p, Double_t C=0.88); // C(pi+)=0.88, C(pi-)=0.93
Double_t          Chi2PID_pion_upperBound (Double_t p, Double_t C=0.88); // C(pi+)=0.88, C(pi-)=0.93
int                       GetBeamHelicity (event_ptr p_event, int runnum, int fdebug);
double                      GetBeamEnergy (int fdebug);
TString                   GetRunNumberSTR (int RunNumber, int fdebug);
void                InitializeFileReading (int NeventsMax,int c12Nentries, int fdebug);
void                  InitializeVariables ();
void                      OpenResultFiles ( TString outfilepath, TString outfilename );
TObjArray *            OpenInputHipoFiles ( TString inputFile, int fdebug );
void           ExtractElectronInformation (int fdebug);
void              ExtractPionsInformation (int fdebug);
void               ExtractPipsInformation (int pipsIdx, int fdebug );
void               ExtractPimsInformation (int pimsIdx, int fdebug );
void                    ComputeKinematics ();
void                   WriteEventToOutput (int fdebug);
void                        FinishProgram (TString outfilepath, TString outfilename);
void                   GetParticlesByType (int event, int fdebug );
void              Stream_e_pi_line_to_CSV (TString pionCharge, int piIdx, bool IsSelectedEvent, int fdebug );
// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.

// globals
auto db = TDatabasePDG::Instance();
// cut values
std::vector<std::pair<std::string, double>> cutValues;
double               cutValue_Vz_min;
double               cutValue_Vz_max;
double             cutValue_e_PCAL_W;
double             cutValue_e_PCAL_V;
double             cutValue_e_E_PCAL;
double cutValue_SamplingFraction_min;
double        cutValue_Ve_Vpi_dz_max;


bool        ePastCutsInEvent = false;
bool     pipsPastCutsInEvent = false;
bool     pimsPastCutsInEvent = false;
bool         EventPassedCuts = false;

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
int       Nevents_passed_e_cuts;
int    Nevents_passed_pips_cuts;
int  Nevents_passed_e_pips_cuts;
int    Nevents_passed_pims_cuts;
int  Nevents_passed_e_pims_cuts;

// number of particles per event
int         Ne, Nn, Np, Npips, Npims, Ngammas;

// variables
double          Mp = 0.938;
double       Mp2 = Mp * Mp;
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

// positive pions
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
// negative pions
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
double               Zpims[NMAXPIONS]; // hadron rest-frame energy



// Output root file and tree
TFile * outFile_e_piplus, * outFile_e_piminus;
TTree * outTree_e_piplus, * outTree_e_piminus;
// Output CSV file
std::ofstream   CSVfile_e_piplus,  SelectedEventsCSVfile_e_piplus;
std::ofstream   CSVfile_e_piminus, SelectedEventsCSVfile_e_piminus;
// vectors in lab-frame
TLorentzVector                  Beam, e, q;
TLorentzVector           piplus[NMAXPIONS]; // leading positive pion
TLorentzVector          piminus[NMAXPIONS]; // leading negative pion
// reconstructed vertex position
TVector3                                Ve;
TVector3                Vpiplus[NMAXPIONS];
TVector3               Vpiminus[NMAXPIONS];

// kinematics
Double_t     Ebeam, xB, Q2, omega, W, W2;


// vectors in q-frame
//TLorentzVector       piplus_qFrame;
//Double_t        Ppips_t_q, Ppips_q;
// auxiliary
DCfid_SIDIS dcfid;
std::vector<region_part_ptr>  electrons, neutrons, protons, pipluses, piminuses, gammas;



// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
// Main functionality
// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void SIDISc12rSkimmer(int  RunNumber=6420,
                      int  NeventsMax=-1,
                      int  fdebug=1,
                      int  PrintProgress=5000,
                      int NpipsMin=1, // minimal number of pi+
                      TString DataPath = "/volatile/clas12/rg-b/production/recon/spring2019/torus-1/pass1/v0/dst/train_20200610/inc/" ){
    TString RunNumberStr = GetRunNumberSTR(RunNumber,fdebug);
    // read cut values
    loadCutValues("cutValues.csv",fdebug);
    
    // open result files
    TString outfilepath = "/volatile/clas12/users/ecohen/BAND/SIDIS_skimming/";
    TString outfilename = "skimmed_SIDIS_inc_" + RunNumberStr;
    OpenResultFiles( outfilepath, outfilename );
    auto files = OpenInputHipoFiles( DataPath + "inc_" + RunNumberStr + ".hipo", fdebug );
    if (fdebug) std::cout << files->GetEntries() << " files to analyze.. "  << std::endl;
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
            
            GetBeamHelicity    ( c12.event() , runnum, fdebug );
            InitializeVariables();
            // Get Particles By Type
            electrons   = c12.getByID( 11   );
            neutrons    = c12.getByID( 2112 );
            protons     = c12.getByID( 2212 );
            pipluses    = c12.getByID( 211  );
            piminuses   = c12.getByID(-211  );
            gammas      = c12.getByID( 22   );
            GetParticlesByType ( event, fdebug );
            
            
            // filter events, extract information, and compute event kinematics:
            // we keep only d(e,e’pi+)X and d(e,e’pi-)X events
            if(  0 < Ne
               &&
               ((0 < Npips && Npips < NMAXPIONS) || (0 < Npims && Npims < NMAXPIONS)) ){
                   
                ExtractElectronInformation  (fdebug);
                ComputeKinematics           ();
                ExtractPionsInformation     (fdebug);
                WriteEventToOutput          (fdebug);
                
            } else {
                if (fdebug>1) {
                    std::cout << "Skipped computations in this event as there are not enough particles: "
                    << "Ne = " << Ne << ",Npips = " << Npips << ",Npims = " << Npims << std::endl ;
                }
            }
            if (fdebug>1) {
                std::cout << "done processing event " << event
                << std::endl << "------------------------------------------------------------" << std::endl ;
            }
            event++;
            if (fdebug && event%PrintProgress==0) std::cout << std::setprecision(1) << " event " << event << std::endl;
        } // end event loop
        
        Nevents_processed = event;
    } // end file loop
    
    
    FinishProgram( outfilepath, outfilename);
}




// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
bool ElectronPassedSelectionCriteria(Double_t e_PCAL_x, Double_t e_PCAL_y,
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
    
    
    if(
       // fiducial cuts on PCAL
       fabs(e_PCAL_x)>0
       &&  fabs(e_PCAL_y)>0
       &&  e_PCAL_W > cutValue_e_PCAL_W
       &&  e_PCAL_V > cutValue_e_PCAL_V
       
       // Electron Identification Refinement  - PCAL Minimum Energy Deposition Cut
       &&  e_E_PCAL > cutValue_e_E_PCAL
       
       // Sampling fraction cut
       && ((e_E_PCAL + e_E_ECIN + e_E_ECOUT)/e.P()) > cutValue_SamplingFraction_min
       && (e_E_ECIN/e.P() > 0.2 - e_E_PCAL/e.P()) // RGA AN puts "<" here mistakenly
       
       // Cut on z-vertex position: in-bending torus field -13.0 cm < Vz < +12.0 cm
       // Spring 19 and Spring 2020 in-bending.
       // Fall 2019 (without low-energy-run) was out-bending.
       &&  ((cutValue_Vz_min < Ve.Z()) && (Ve.Z() < cutValue_Vz_max))
       ) return true;
    
    return false;
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
bool PionPassedSelectionCutsCriteria(TString pionCharge, // "pi+" or "pi-"
                                     Double_t DC_sector,
                                     Double_t DC_x[3], Double_t DC_y[3], Double_t DC_z[3],
                                     Double_t chi2PID, Double_t p,
                                     TVector3 Ve,
                                     TVector3 Vpi,
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
        std::cout << "PionPassedSelectionCutsCriteria()" << std::endl;
    }
    if (DC_sector == 0) return false;
    
    int PDGcode;
    double    C;
    if (pionCharge=="pi+"){
        PDGcode = 211;
        C       = 0.88;
    } else if (pionCharge=="pi-") {
        PDGcode = -211;
        C       = 0.93;
    } else {
        std::cout << "pion charge is not defined in PionPassedSelectionCutsCriteria(), returning false" << std::endl;
        return false;
    }
    
    for (int regionIdx=0; regionIdx<3; regionIdx++) {
        // DC_e_fid:
        // sector:  1-6
        // layer:   1-3
        // bending: 0(out)/1(in)
        int bending  = 1 ? (torusBending==-1) : 0;
        // new version Aug-11,2021
        if (fdebug>3) {
            std::cout << "dcfid.DC_fid_th_ph_sidis()"
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
        std::cout << "in PionPassedSelectionCutsCriteria()"<< std::endl
        << "pion charge: "          << pionCharge               << ","
        << "DC_x[0]: "              << DC_x[0]                  << ","
        << "chi2PID:"               << chi2PID                  << ","
        << "Chi2PID_pion_lowerBound( p="<<p<<", C="<<C<<" ): "  << Chi2PID_pion_lowerBound( p, C ) << ","
        << "Chi2PID_pion_upperBound( p="<<p<<", C="<<C<<" ): "  << Chi2PID_pion_upperBound( p, C ) << ","
        << "fabs((Ve-Vpi).Z()): "   << fabs((Ve-Vpi).Z())       << ","
        << std::endl;
    }
    if(
       // pi+ Identification Refinement - chi2PID vs. momentum
       ( Chi2PID_pion_lowerBound( p, C ) < chi2PID && chi2PID < Chi2PID_pion_upperBound( p , C ) )
       
       // Cut on the z-Vertex Difference Between Electrons and Hadrons.
       &&  ( fabs((Ve-Vpi).Z()) < cutValue_Ve_Vpi_dz_max )
       ) {
        if (fdebug>3) {
            std::cout << "done PionPassedSelectionCutsCriteria(), return false" << std::endl;
            
            return true;
        }
    }
    return false;
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
void OpenOutputFiles (TString outfilename,TString header){
    
    // Create output tree
    outFile_e_piplus  = new TFile( outfilename + "_e_piplus.root"  ,"RECREATE");
    outTree_e_piplus  = new TTree( "(e,e'pi+) events" , "Event information");
    outFile_e_piminus = new TFile( outfilename + "_e_piminus.root" ,"RECREATE");
    outTree_e_piminus = new TTree( "(e,e'pi-) events" , "Event information");
    
    // Create output csv files
    CSVfile_e_piplus.open( outfilename  + "_e_piplus.csv" );
    CSVfile_e_piplus << header << std::endl;
    CSVfile_e_piminus.open( outfilename + "_e_piminus.csv" );
    CSVfile_e_piminus << header << std::endl;
    
    SelectedEventsCSVfile_e_piplus.open( outfilename + "_e_piplus_selected_events.csv" );
    SelectedEventsCSVfile_e_piplus << header << std::endl;
    SelectedEventsCSVfile_e_piminus.open( outfilename + "_e_piminus_selected_events.csv" );
    SelectedEventsCSVfile_e_piminus << header << std::endl;

}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void CloseOutputFiles (TString OutDataPath, TString outfilename){
    // close output CSV
    CSVfile_e_piplus                .close();
    SelectedEventsCSVfile_e_piplus  .close();
    CSVfile_e_piminus               .close();
    SelectedEventsCSVfile_e_piminus .close();
    
    int Nentires_e_piplus  = outTree_e_piplus  -> GetEntries();
    int Nentires_e_piminus = outTree_e_piminus -> GetEntries();
    
    // close output ROOT
    outFile_e_piplus->cd();
    outTree_e_piplus->Write();
    outFile_e_piplus->Close();
    
    outFile_e_piminus->cd();
    outTree_e_piminus->Write();
    outFile_e_piminus->Close();
    
    
    std::cout
    << "Done processesing "  <<  Nevents_processed          << " events,"
    << std::endl
    << std::setprecision(3)
    << (float)Nevents_passed_e_cuts/Nevents_processed       << " events passed e cuts,"
    << std::endl
    << (float)Nevents_passed_pips_cuts/Nevents_processed    << " events passed pi+ cuts,"
    << std::endl
    << (float)Nevents_passed_e_pips_cuts/Nevents_processed  << " events passed (e,e'pi+) cuts,"
    << std::endl
    << (float)Nevents_passed_pims_cuts/Nevents_processed    << " events passed pi- cuts,"
    << std::endl
    << (float)Nevents_passed_e_pims_cuts/Nevents_processed  << " events passed (e,e'pi-) cuts,"
    << std::endl;
    
    
    
    std::cout << "output files ready in root/csv formats in " << std::endl
    << std::endl
    << "wrote "  << Nentires_e_piplus  << " to (e,e'pi+) root file, "
    << std::endl << outFile_e_piplus -> GetName()
    << std::endl << OutDataPath + outfilename + "_e_piplus_selected_events.csv"
    << std::endl
    << "and "    << Nentires_e_piminus << " to (e,e'pi-) root file. "
    << std::endl << outFile_e_piminus -> GetName()
    << std::endl << OutDataPath + outfilename + "_e_piminus_selected_events.csv"
    << std::endl;
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void StreamToCSVfile (TString pionCharge, // "pi+" or "pi-"
                      std::vector<Double_t> observables, bool IsSelectedEvent, int fdebug){
    if (fdebug>1) {
        std::cout << "streaming to CSVfile" << std::endl;
    }
    // decide which file to write...
    if (pionCharge=="pi+") {
        for (auto v:observables) CSVfile_e_piplus << v << ",";
        CSVfile_e_piplus << std::endl;
        
        if (IsSelectedEvent) {
            for (auto v:observables) SelectedEventsCSVfile_e_piplus << v << ",";
            SelectedEventsCSVfile_e_piplus << std::endl;
        }
        
    }
    else if (pionCharge=="pi-") {
        for (auto v:observables) CSVfile_e_piminus << v << ",";
        CSVfile_e_piminus << std::endl;
        
        if (IsSelectedEvent) {
            for (auto v:observables) SelectedEventsCSVfile_e_piminus << v << ",";
            SelectedEventsCSVfile_e_piminus << std::endl;
        }
    }
    else {
        std::cout << "pion charge ill-defined in StreamToCSVfile(), returning" << std::endl;
        return;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void loadCutValues(TString cutValuesFilename, int fdebug){
    
    // read cut values csv file
    csv_reader csvr;
    cutValues = csvr.read_csv("cutValues.csv");
    if (fdebug>0) { printCutValues(); }
    
    // assign specific cut values - to speed things up
    // by avoiding recalling FindCutValue() on every event
    
    // Cut on z-vertex position:
    // based on RGA_Analysis_Overview_and_Procedures_Nov_4_2020-6245173-2020-12-09-v3.pdf
    // p. 71
    if (torusBending==-1){ // in-bending torus field
        // Spring 19 and Spring 2020 in-bending.
        cutValue_Vz_min = FindCutValue("Vz_e_min_inbending");
        cutValue_Vz_max = FindCutValue("Vz_e_max_inbending");
    } else if (torusBending==1){ // Out-bending torus field
        // Fall 2019 (without low-energy-run) was out-bending.
        cutValue_Vz_min = FindCutValue("Vz_e_min_outbending");
        cutValue_Vz_max = FindCutValue("Vz_e_max_outbending");
        
    } else {
        std::cout
        << "Un-identified torus bending "
        << torusBending
        << ", return" << std::endl;
        return;
    }
        
    cutValue_e_PCAL_W               = FindCutValue("e_PCAL_W_min");
    cutValue_e_PCAL_V               = FindCutValue("e_PCAL_V_min");
    cutValue_e_E_PCAL               = FindCutValue("e_E_PCAL_min");
    cutValue_SamplingFraction_min   = FindCutValue("SamplingFraction_min");
    cutValue_Ve_Vpi_dz_max          = FindCutValue("(Ve-Vpi)_z_max");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void printCutValues(){
    std::cout << "Using cut values:" << std::endl;
    for (auto cut: cutValues) {
        std::cout << cut.first << ": " << cut.second << std::endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
double FindCutValue( std::string cutName ){
    for (auto cut: cutValues) {
        if (strcmp(cut.first.c_str(),cutName.c_str())==0){
            return cut.second;
        }
    }
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
void SetOutputTTrees(){
    // pi+
    outTree_e_piplus->Branch("eventnumber"          ,&evnum                 );
    outTree_e_piplus->Branch("runnum"               ,&runnum                );
    outTree_e_piplus->Branch("e_E_PCAL"             ,&e_E_PCAL              );
    outTree_e_piplus->Branch("e_E_ECIN"             ,&e_E_ECIN              );
    outTree_e_piplus->Branch("e_E_ECOUT"            ,&e_E_ECOUT             );
    outTree_e_piplus->Branch("pi_chi2PID"           ,&pips_chi2PID          );
    outTree_e_piplus->Branch("e_PCAL_W"             ,&e_PCAL_W              );
    outTree_e_piplus->Branch("e_PCAL_V"             ,&e_PCAL_V              );
    outTree_e_piplus->Branch("pi_PCAL_x"            ,&pips_PCAL_x           );
    outTree_e_piplus->Branch("pi_PCAL_y"            ,&pips_PCAL_y           );
    outTree_e_piplus->Branch("pi_PCAL_z"            ,&pips_PCAL_z           );
    outTree_e_piplus->Branch("e_PCAL_x"             ,&e_PCAL_x              );
    outTree_e_piplus->Branch("e_PCAL_y"             ,&e_PCAL_y              );
    outTree_e_piplus->Branch("e_PCAL_z"             ,&e_PCAL_z              );
    outTree_e_piplus->Branch("e_PCAL_sector"        ,&e_PCAL_sector         );
    outTree_e_piplus->Branch("e_DC_sector"          ,&e_DC_sector           );
    outTree_e_piplus->Branch("e_DC_Chi2N"           ,&e_DC_Chi2N            );
    outTree_e_piplus->Branch("e_DC_x"               ,&e_DC_x                );
    outTree_e_piplus->Branch("e_DC_y"               ,&e_DC_y                );
    outTree_e_piplus->Branch("e_DC_z"               ,&e_DC_z                );
    outTree_e_piplus->Branch("pi_PCAL_sector"       ,&pips_PCAL_sector      );
    outTree_e_piplus->Branch("pi_DC_sector"         ,&pips_DC_sector        );
    outTree_e_piplus->Branch("pi_Chi2N"             ,&pips_Chi2N            );
    outTree_e_piplus->Branch("pi_DC_x"              ,&pips_DC_x             );
    outTree_e_piplus->Branch("pi_DC_y"              ,&pips_DC_y             );
    outTree_e_piplus->Branch("pi_DC_z"              ,&pips_DC_z             );
    outTree_e_piplus->Branch("pi_E_PCAL"            ,&pips_E_PCAL           );
    outTree_e_piplus->Branch("pi_E_ECIN"            ,&pips_E_ECIN           );
    outTree_e_piplus->Branch("pi_E_ECIN"            ,&pips_E_ECIN           );
    outTree_e_piplus->Branch("pi_E_ECOUT"           ,&pips_E_ECOUT          );
    outTree_e_piplus->Branch("DC_layer"             ,&DC_layer              );
    outTree_e_piplus->Branch("e"                    ,&e                     );
    outTree_e_piplus->Branch("pi"                   ,&piplus                );
    outTree_e_piplus->Branch("Ve"                   ,&Ve                    );
    outTree_e_piplus->Branch("Vpi"                  ,&Vpiplus               );
    outTree_e_piplus->Branch("Beam"                 ,&Beam                  );
    outTree_e_piplus->Branch("beam_helicity"        ,&beam_helicity         );
    outTree_e_piplus->Branch("q"                    ,&q                     );
    outTree_e_piplus->Branch("Ebeam"                ,&Ebeam                 );
    outTree_e_piplus->Branch("xB"                   ,&xB                    );
    outTree_e_piplus->Branch("Q2"                   ,&Q2                    );
    outTree_e_piplus->Branch("omega"                ,&omega                 );
    outTree_e_piplus->Branch("W"                    ,&W                     );
    outTree_e_piplus->Branch("Z"                    ,&Zpips                 );

    outTree_e_piplus->Branch("EventPassedCuts"      ,&EventPassedCuts       );
    outTree_e_piplus->Branch("ePastCutsInEvent"     ,&ePastCutsInEvent      );
    outTree_e_piplus->Branch("piPastCutsInEvent"    ,&pipsPastCutsInEvent   );
    
    
    // pi-
    outTree_e_piminus->Branch("eventnumber"         ,&evnum                 );
    outTree_e_piminus->Branch("runnum"              ,&runnum                );
    outTree_e_piminus->Branch("e_E_PCAL"            ,&e_E_PCAL              );
    outTree_e_piminus->Branch("e_E_ECIN"            ,&e_E_ECIN              );
    outTree_e_piminus->Branch("e_E_ECOUT"           ,&e_E_ECOUT             );
    outTree_e_piminus->Branch("pi_chi2PID"          ,&pims_chi2PID          );
    outTree_e_piminus->Branch("e_PCAL_W"            ,&e_PCAL_W              );
    outTree_e_piminus->Branch("e_PCAL_V"            ,&e_PCAL_V              );
    outTree_e_piminus->Branch("pi_PCAL_x"           ,&pims_PCAL_x           );
    outTree_e_piminus->Branch("pi_PCAL_y"           ,&pims_PCAL_y           );
    outTree_e_piminus->Branch("pi_PCAL_z"           ,&pims_PCAL_z           );
    outTree_e_piminus->Branch("e_PCAL_x"            ,&e_PCAL_x              );
    outTree_e_piminus->Branch("e_PCAL_y"            ,&e_PCAL_y              );
    outTree_e_piminus->Branch("e_PCAL_z"            ,&e_PCAL_z              );
    outTree_e_piminus->Branch("e_PCAL_sector"       ,&e_PCAL_sector         );
    outTree_e_piminus->Branch("e_DC_sector"         ,&e_DC_sector           );
    outTree_e_piminus->Branch("e_DC_Chi2N"          ,&e_DC_Chi2N            );
    outTree_e_piminus->Branch("e_DC_x"              ,&e_DC_x                );
    outTree_e_piminus->Branch("e_DC_y"              ,&e_DC_y                );
    outTree_e_piminus->Branch("e_DC_z"              ,&e_DC_z                );
    outTree_e_piminus->Branch("pi_PCAL_sector"      ,&pims_PCAL_sector      );
    outTree_e_piminus->Branch("pi_DC_sector"        ,&pims_DC_sector        );
    outTree_e_piminus->Branch("pi_Chi2N"            ,&pims_Chi2N            );
    outTree_e_piminus->Branch("pi_DC_x"             ,&pims_DC_x             );
    outTree_e_piminus->Branch("pi_DC_y"             ,&pims_DC_y             );
    outTree_e_piminus->Branch("pi_DC_z"             ,&pims_DC_z             );
    outTree_e_piminus->Branch("pi_E_PCAL"           ,&pims_E_PCAL           );
    outTree_e_piminus->Branch("pi_E_ECIN"           ,&pims_E_ECIN           );
    outTree_e_piminus->Branch("pi_E_ECIN"           ,&pims_E_ECIN           );
    outTree_e_piminus->Branch("pi_E_ECOUT"          ,&pims_E_ECOUT          );
    outTree_e_piminus->Branch("DC_layer"            ,&DC_layer              );
    outTree_e_piminus->Branch("e"                   ,&e                     );
    outTree_e_piminus->Branch("pi"                  ,&piminus               );
    outTree_e_piminus->Branch("Ve"                  ,&Ve                    );
    outTree_e_piminus->Branch("Vpi"                 ,&Vpiminus              );
    outTree_e_piminus->Branch("Beam"                ,&Beam                  );
    outTree_e_piminus->Branch("beam_helicity"       ,&beam_helicity         );
    outTree_e_piminus->Branch("q"                   ,&q                     );
    outTree_e_piminus->Branch("Ebeam"               ,&Ebeam                 );
    outTree_e_piminus->Branch("xB"                  ,&xB                    );
    outTree_e_piminus->Branch("Q2"                  ,&Q2                    );
    outTree_e_piminus->Branch("omega"               ,&omega                 );
    outTree_e_piminus->Branch("W"                   ,&W                     );
    outTree_e_piplus->Branch("Z"                    ,&Zpims                 );

    outTree_e_piminus->Branch("EventPassedCuts"     ,&EventPassedCuts       );
    outTree_e_piminus->Branch("ePastCutsInEvent"     ,&ePastCutsInEvent      );
    outTree_e_piminus->Branch("piPastCutsInEvent"    ,&pimsPastCutsInEvent   );
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
double GetBeamEnergy (int fdebug){
    // ToDo:
    // make this automatic using GetBeamEnergy
    // rcdb crashes with
    /*
     *** Break *** segmentation violation



     ===========================================================
     There was a crash.
     This is the entire stack trace of all threads:
     ===========================================================
     #0  0x00007f1733bdf41c in waitpid () from /lib64/libc.so.6
     #1  0x00007f1733b5cf12 in do_system () from /lib64/libc.so.6
     #2  0x00007f17386cff95 in TUnixSystem::StackTrace() () from /site/12gev_phys/2.4/Linux_CentOS7.7.1908-gcc9.2.0/root/6.20.04/lib/libCore.so.6.20
     #3  0x00007f17386cd00c in TUnixSystem::DispatchSignals(ESignals) () from /site/12gev_phys/2.4/Linux_CentOS7.7.1908-gcc9.2.0/root/6.20.04/lib/libCore.so.6.20
     #4  <signal handler called>
     #5  0x00007f171f2b7584 in SIDISc12rSkimmer(int, int, int, bool, int, int, TString) () from /u/home/cohen/SIDIS_at_BAND/SIDISc12rSkimmer_C.so
     #6  0x00007f17391a308b in ?? ()
     #7  0x00007ffe3aadbba0 in ?? ()
     #8  0x00007ffe3aadbc90 in ?? ()
     #9  0x00007ffe3aadbc50 in ?? ()
     #10 0x00007ffe3aadbba0 in ?? ()
     #11 0x00007ffe00000001 in ?? ()
     #12 0x00000000019f8a20 in ?? ()
     #13 0x00007f1738abe220 in vtable for TString () from /site/12gev_phys/2.4/Linux_CentOS7.7.1908-gcc9.2.0/root/6.20.04/lib/libCore.so.6.20
     #14 0x0000005a00000061 in ?? ()
     #15 0x0000000008f9ae70 in ?? ()
     #16 0x00007ffe3aadc120 in ?? ()
     #17 0x0000000000861d80 in ?? ()
     #18 0x00007f172cb02df1 in cling::IncrementalExecutor::executeWrapper(llvm::StringRef, cling::Value*) const () from /site/12gev_phys/2.4/Linux_CentOS7.7.1908-gcc9.2.0/root/6.20.04/lib/libCling.so
     #19 0x00007f172ca8ef23 in cling::Interpreter::RunFunction(clang::FunctionDecl const*, cling::Value*) () from /site/12gev_phys/2.4/Linux_CentOS7.7.1908-gcc9.2.0/root/6.20.04/lib/libCling.so
     #20 0x00007f172ca90a5d in cling::Interpreter::EvaluateInternal(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const&, cling::CompilationOptions, cling::Value*, cling::Transaction**, unsigned long) () from /site/12gev_phys/2.4/Linux_CentOS7.7.1908-gcc9.2.0/root/6.20.04/lib/libCling.so
     #21 0x00007f172ca90d45 in cling::Interpreter::process(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const&, cling::Value*, cling::Transaction**, bool) () from /site/12gev_phys/2.4/Linux_CentOS7.7.1908-gcc9.2.0/root/6.20.04/lib/libCling.so
     #22 0x00007f172cb517ad in cling::MetaProcessor::process(llvm::StringRef, cling::Interpreter::CompilationResult&, cling::Value*, bool) () from /site/12gev_phys/2.4/Linux_CentOS7.7.1908-gcc9.2.0/root/6.20.04/lib/libCling.so
     #23 0x00007f172c9f711c in HandleInterpreterException(cling::MetaProcessor*, char const*, cling::Interpreter::CompilationResult&, cling::Value*) () from /site/12gev_phys/2.4/Linux_CentOS7.7.1908-gcc9.2.0/root/6.20.04/lib/libCling.so
     #24 0x00007f172ca0d65c in TCling::ProcessLine(char const*, TInterpreter::EErrorCode*) () from /site/12gev_phys/2.4/Linux_CentOS7.7.1908-gcc9.2.0/root/6.20.04/lib/libCling.so
     #25 0x00007f172ca0dae1 in TCling::ProcessLineSynch(char const*, TInterpreter::EErrorCode*) () from /site/12gev_phys/2.4/Linux_CentOS7.7.1908-gcc9.2.0/root/6.20.04/lib/libCling.so
     #26 0x00007f173858ca0a in TApplication::ExecuteFile(char const*, int*, bool) [clone .localalias] () from /site/12gev_phys/2.4/Linux_CentOS7.7.1908-gcc9.2.0/root/6.20.04/lib/libCore.so.6.20
     #27 0x00007f173858d6a7 in TApplication::ProcessLine(char const*, bool, int*) () from /site/12gev_phys/2.4/Linux_CentOS7.7.1908-gcc9.2.0/root/6.20.04/lib/libCore.so.6.20
     #28 0x00007f1735e4b462 in TRint::ProcessLineNr(char const*, char const*, int*) () from /site/12gev_phys/2.4/Linux_CentOS7.7.1908-gcc9.2.0/root/6.20.04/lib/libRint.so.6.20
     #29 0x00007f1735e4cb7b in TRint::Run(bool) () from /site/12gev_phys/2.4/Linux_CentOS7.7.1908-gcc9.2.0/root/6.20.04/lib/libRint.so.6.20
     #30 0x0000000000400c76 in main ()
     ===========================================================


     The lines below might hint at the cause of the crash.
     You may get help by asking at the ROOT forum http://root.cern.ch/forum
     Only if you are really convinced it is a bug in ROOT then please submit a
     report at http://root.cern.ch/bugs Please post the ENTIRE stack trace
     from above as an attachment in addition to anything else
     that might help us fixing this issue.
     ===========================================================
     #5  0x00007f171f2b7584 in SIDISc12rSkimmer(int, int, int, bool, int, int, TString) () from /u/home/cohen/SIDIS_at_BAND/SIDISc12rSkimmer_C.so
     #6  0x00007f17391a308b in ?? ()
     #7  0x00007ffe3aadbba0 in ?? ()
     #8  0x00007ffe3aadbc90 in ?? ()
     #9  0x00007ffe3aadbc50 in ?? ()
     #10 0x00007ffe3aadbba0 in ?? ()
     #11 0x00007ffe00000001 in ?? ()
     #12 0x00000000019f8a20 in ?? ()
     #13 0x00007f1738abe220 in vtable for TString () from /site/12gev_phys/2.4/Linux_CentOS7.7.1908-gcc9.2.0/root/6.20.04/lib/libCore.so.6.20
     #14 0x0000005a00000061 in ?? ()
     #15 0x0000000008f9ae70 in ?? ()
     #16 0x00007ffe3aadc120 in ?? ()
     #17 0x0000000000861d80 in ?? ()
     #18 0x00007f172cb02df1 in cling::IncrementalExecutor::executeWrapper(llvm::StringRef, cling::Value*) const () from /site/12gev_phys/2.4/Linux_CentOS7.7.1908-gcc9.2.0/root/6.20.04/lib/libCling.so
     #19 0x00007f172ca8ef23 in cling::Interpreter::RunFunction(clang::FunctionDecl const*, cling::Value*) () from /site/12gev_phys/2.4/Linux_CentOS7.7.1908-gcc9.2.0/root/6.20.04/lib/libCling.so
     #20 0x00007f172ca90a5d in cling::Interpreter::EvaluateInternal(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const&, cling::CompilationOptions, cling::Value*, cling::Transaction**, unsigned long) () from /site/12gev_phys/2.4/Linux_CentOS7.7.1908-gcc9.2.0/root/6.20.04/lib/libCling.so
     #21 0x00007f172ca90d45 in cling::Interpreter::process(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const&, cling::Value*, cling::Transaction**, bool) () from /site/12gev_phys/2.4/Linux_CentOS7.7.1908-gcc9.2.0/root/6.20.04/lib/libCling.so
     #22 0x00007f172cb517ad in cling::MetaProcessor::process(llvm::StringRef, cling::Interpreter::CompilationResult&, cling::Value*, bool) () from /site/12gev_phys/2.4/Linux_CentOS7.7.1908-gcc9.2.0/root/6.20.04/lib/libCling.so
     #23 0x00007f172c9f711c in HandleInterpreterException(cling::MetaProcessor*, char const*, cling::Interpreter::CompilationResult&, cling::Value*) () from /site/12gev_phys/2.4/Linux_CentOS7.7.1908-gcc9.2.0/root/6.20.04/lib/libCling.so
     ===========================================================

     */
    //rcdb info
    //        if (fdebug>3) std::cout << "reading RCDB info" << std::endl;
    //        auto& rcdbData = c12.rcdb()->current();//struct with all relevent rcdb values
    //
    //        // get beam energy
    //        if (fdebug>3) std::cout << "getting beam energy" << std::endl;
    //        Ebeam = rcdbData.beam_energy ;
    if (fdebug>3) std::cout << "set beam energy" << std::endl;
    double Ebeam = 10.2; // [GeV] ( for Fall-2019 the enrgy was 10.4096)
    return Ebeam;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void InitializeFileReading(int NeventsMax, int c12Nentries, int fdebug){
    Ebeam = GetBeamEnergy( fdebug );
    Beam.SetPxPyPzE(0,0,Ebeam,Ebeam);
    Int_t NeventsMaxToProcess = NeventsMax;
    if (NeventsMax<0) NeventsMaxToProcess = c12Nentries;
    Nevents_processed           = 0;
    Nevents_passed_e_cuts       = 0;
    Nevents_passed_pips_cuts    = 0;
    Nevents_passed_pims_cuts    = 0;
    Nevents_passed_e_pips_cuts  = 0;
    Nevents_passed_e_pims_cuts  = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void InitializeVariables(){
    e = TLorentzVector(0,0,0,db->GetParticle( 11   )->Mass());
    
    xB          = Q2        = omega     = -9999;
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
    ePastCutsInEvent                    = false;
    
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
        piplus[piIdx]                               = TLorentzVector(0,0,0,db->GetParticle( 211 )->Mass());
        Vpiplus[piIdx]                              = TVector3();
        pipsPastSelectionCuts[piIdx]                = false;
        
        piminus[piIdx]                              = TLorentzVector(0,0,0,db->GetParticle(-211 )->Mass());
        Vpiminus[piIdx]                             = TVector3();
        pimsPastSelectionCuts[piIdx]                = false;
    }
    DC_layer                                        = -9999;
    status                                          = 1; // 0 is good...
    
    pipsPastCutsInEvent                             = false;
    pimsPastCutsInEvent                             = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TString GetRunNumberSTR(int RunNumber, int fdebug){
    char RunNumberStr[20];
    sprintf( RunNumberStr, "00%d", RunNumber );
    if (fdebug>1) std::cout << "(SIDIS) skimming run " << RunNumberStr << std::endl;
    return (TString)RunNumberStr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void OpenResultFiles( TString outfilepath, TString outfilename ){
    OpenOutputFiles( outfilepath + outfilename,
                    ( (TString)"status, runnum,     evnum,      beam_helicity,"
                     +(TString)"e.P(),  e.Theta(),  e.Phi(),    Ve.Z(),"
                     +(TString)"pi.P(), pi.Theta(), pi.Phi(),   Vpi.Z(),"
                     +(TString)"Q2,     W,          x,          Zpi,"));
    // output tree branches
    SetOutputTTrees();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TObjArray * OpenInputHipoFiles( TString inputFile, int fdebug ){
    
    // ------------------------------------------------------------
    // open input file(s)
    // ------------------------------------------------------------
//    for(Int_t i=1;i<gApplication->Argc();i++){
//        TString opt=gApplication->Argv(i);
//        if((opt.Contains(".hipo"))){
//            inputFile=opt(5,opt.Sizeof());
//        }
//    }
    if(inputFile==TString())  {
        std::cout << " *** please provide a file name..." << std::endl;
        exit(0);
    }
    if (fdebug) std::cout << "Analysing hipo file " << inputFile << std::endl;
    // decide if torus is in-bending or out-bending
    if (torusBending==-1){
        std::cout << "using In-Bending torus data" << std::endl;
    }
    else {
        std::cout << "using Out-Bending torus data" << std::endl;
    }
    TChain fake("hipo");
    fake.Add(inputFile.Data());
    //get the hipo data
    auto files = fake.GetListOfFiles();
    gBenchmark->Start("timer");
    return files;
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
    if (fdebug > 2) std::cout << "extracted electron information and computed kinematics" << std::endl;
    
    // ------------------------------------------------------------------------------------------------
    // now, check if electron passed event selection requirements
    // ------------------------------------------------------------------------------------------------
    ePastCutsInEvent = ElectronPassedSelectionCriteria(e_PCAL_x, e_PCAL_y,
                                                              e_PCAL_W, e_PCAL_V,
                                                              e_E_PCAL, e_E_ECIN,
                                                              e_E_ECOUT,
                                                              e, Ve,
                                                              e_PCAL_sector, // e_PCAL_sector should be consistent with e_DC_sector
                                                              e_DC_x, e_DC_y, e_DC_z,
                                                              torusBending );
    if (ePastCutsInEvent)  Nevents_passed_e_cuts++ ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExtractPionsInformation(int fdebug){
    
    // positive pions)
    for (int pipsIdx=0; pipsIdx < Npips; pipsIdx++) {
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
    
    // (Maybe) write this event to "selected events csv-file"
    bool IsSelectedEvent = false;
        
    if (ePastCutsInEvent && pipsPastCutsInEvent) {
        IsSelectedEvent = true;
        outTree_e_piplus -> Fill();
        Nevents_passed_e_pips_cuts ++ ;
        for (int pipsIdx=0; pipsIdx<Npips; pipsIdx++) {
            Stream_e_pi_line_to_CSV( "pi+", pipsIdx, IsSelectedEvent, fdebug );
        }
    }
    
    if (ePastCutsInEvent && pimsPastCutsInEvent) {
        IsSelectedEvent = true;
        outTree_e_piminus -> Fill();
        Nevents_passed_e_pims_cuts ++ ;
        for (int pimsIdx=0; pimsIdx<Npims; pimsIdx++) {
            Stream_e_pi_line_to_CSV( "pi-", pimsIdx, IsSelectedEvent, fdebug );
        }
    }


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void FinishProgram(TString outfilepath, TString outfilename){
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    
    std::cout << "Done. Elapsed time: " << elapsed.count() << std::endl;
    CloseOutputFiles( outfilepath, outfilename );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ComputeKinematics(){
    // compute event kinematics
    q       = Beam - e;
    Q2      = -q.Mag2();
    omega   = q.E();
    xB      = Q2/(2. * Mp * q.E());
    W2      = Mp2 - Q2 + 2. * omega * Mp;
    W       = sqrt(W2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExtractPipsInformation( int pipsIdx, int fdebug ){
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
    }
    // ------------------------------------------------------------------------------------------------
    // now, check if pion passed event selection requirements
    // ------------------------------------------------------------------------------------------------
    pipsPastSelectionCuts[pipsIdx] = PionPassedSelectionCutsCriteria("pi+",
                                                                     pips_DC_sector[pipsIdx],
                                                                     pips_DC_x[pipsIdx],
                                                                     pips_DC_y[pipsIdx],
                                                                     pips_DC_z[pipsIdx],
                                                                     pips_chi2PID[pipsIdx],  piplus[pipsIdx].P(),
                                                                     Ve,
                                                                     Vpiplus[pipsIdx],
                                                                     fdebug);
    if (pipsPastSelectionCuts[pipsIdx]) {
        pimsPastCutsInEvent = true;
    }
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
    // ------------------------------------------------------------------------------------------------
    // now, check if pion passed event selection requirements
    // ------------------------------------------------------------------------------------------------
    pimsPastSelectionCuts[pimsIdx] = PionPassedSelectionCutsCriteria("pi-",
                                                                     pims_DC_sector[pimsIdx],
                                                                     pims_DC_x[pimsIdx],
                                                                     pims_DC_y[pimsIdx],
                                                                     pims_DC_z[pimsIdx],
                                                                     pims_chi2PID[pimsIdx],  piminus[pimsIdx].P(),
                                                                     Ve,
                                                                     Vpiminus[pimsIdx],
                                                                     fdebug);
    if (pimsPastSelectionCuts[pimsIdx]) {
        pimsPastCutsInEvent = true;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GetParticlesByType (int event, int fdebug){
    // get particles by type
    Ne      = electrons .size();
    Nn      = neutrons  .size();
    Np      = protons   .size();
    Npips   = pipluses  .size();
    Npims   = piminuses .size();
    Ngammas = gammas    .size();
    if (fdebug>2){
        std::cout
        << "particles in event "            << event        << " : "
        << "N(electrons): "                 << Ne           <<  ","
        << "N(protons): "                   << Np           <<  ","
        << "N(neutrons): "                  << Nn           <<  ","
        << "N(pi+): "                       << Npips        <<  ","
        << "N(pi-): "                       << Npims        <<  ","
        << "N(gammas): "                    << Ngammas      <<  ","
        << std::endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Stream_e_pi_line_to_CSV( TString pionCharge, int piIdx, bool IsSelectedEvent, int fdebug ){
    TLorentzVector  pi;
    TVector3        Vpi;
    double          Zpi;
    if (pionCharge=="pi+") {
        pi  = piplus [piIdx];
        Vpi = Vpiplus[piIdx];
        Zpi = Zpips  [piIdx];
    }
    else if (pionCharge=="pi-") {
        pi  = piminus [piIdx];
        Vpi = Vpiminus[piIdx];
        Zpi = Zpims   [piIdx];
   }
    else {
        std::cout << "pion charge ill defined at Stream_e_pi_line_to_CSV(), returning " << std::endl;
        return;
    }
    
    // write a (e,e'pi) event-line to CSV file

    // from Harut A., Aug-2, 2021:
    //    1: status, 2: runnum, 3: evnum, 4: helicity, 5: e_p, 6: e_theta, 7: e_phi, 8: vz_e, 9: pi_p, 10: pi_theta, 11: pi_phi, 12: vz_pi, 13: P_p, 14: P_theta, 15: P_phi, 16: vz_P, 17: Q2, 18: W, 19: x, 20: y, 21: z_pi, 22: z_P, 23: Mx(e:pi:P:X), 24: Mx(e:pi:X), 25: Mx(e:P:X), 26: zeta, 27: Mh, 28: pT_pi, 29: pT_P, 30: pTpT, 31: xF_pi, 32: xF_P, 33: eta_pi, 34: eta_P, 35: Delta_eta, 36: phi_pi (gamma*N COM), 37: phi_P (gamma*N COM), 38: Delta_phi.
    //
    //    1: status is a number indicating the quality of the event with the non-0 number indicating something was not not good (ex. out of fiducial region, not within the final cuts on energies of particles, missing or invarian masses....) The final observables will be done using status==0, while sensitivity of the observable to different cuts could be studied for various values of status>0
    //    2-3: run number and event number to identify the event
    //    4: helicity of the electron +1 along the beam and -1 opposite to it
    //    5,6,7,8 electron momentum,theta,phi_Lab, and z-vertex
    //    9,10,11,12 the same for pi+
    //    All other columns could be calculated from the first 16, but are included for cross check and minimizing the work in production of final observables. Some of them simple, like x,Q^2,W,y,  z=E\pion/nu, zome less trivial like Breit frame rapidities of pion (eta_pi) and proton eta_P, or corresponding values for X_Feynman variable xF_pi, xF_P.
    //    Some, like azimuthal angles of pion (phi_pi) and proton (phi_p) in the CM frame may also be confusing.
    //    In addition to the first 16 columns (mandatory) you can add as many columns as you are comfortable to fill, and consider relevant for your process.
    
    std::vector<double> variables =
    {   (double)status, (double)runnum,     (double)evnum,      (double)beam_helicity,
        e.P(),          e.Theta(),          e.Phi(),            Ve.Z(),
        pi.P(),         pi.Theta(),         pi.Phi(),           Vpi.Z(),
        Q2,             W,                  xB,                 Zpi,
    };
    StreamToCSVfile( pionCharge, variables , IsSelectedEvent, fdebug );
}





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Arxiv
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//    OpenOutputFiles ( outfilepath + outfilename,
//                     ( (TString)("status,runnum,evnum,")
//                      +(TString)("Ebeam,beam_helicity,")
//                      +(TString)("Ee,Pe,Pe_x,Pe_y,Pe_z,Pe_theta,Pe_phi,")
//                      +(TString)("Ngammas,Np,Nn,Npips,Npims,")
//                      +(TString)("omega,q,q_x,q_y,q_z,xB,Q2,z,W,")
//                      +(TString)("Epi,Ppi,Ppi_x,Ppi_y,Ppi_z,")
//                      +(TString)("Ppi_theta,Ppi_phi,")                              //                      +(TString)("Ppi_t_q,Ppips_q,")
//                      +(TString)("e_E_PCAL,e_E_ECIN,e_E_ECOUT,")
//                      +(TString)("Vx_e,Vy_e,Vz_e,Vx_pi,Vy_pi,Vz_pi,")
//                      +(TString)("pi_chi2PID,")
//                      +(TString)("e_PCAL_W,e_PCAL_V,pips_PCAL_W,pips_PCAL_V,")
//                      +(TString)("e_PCAL_x,e_PCAL_y,e_PCAL_z,e_PCAL_sector,")
//                      +(TString)("e_DC_sector,e_DC_Chi2N,")
//                      +(TString)("e_DC_x[region-1],e_DC_y[region-1],e_DC_z[region-1],")
//                      +(TString)("e_DC_x[region-2],e_DC_y[region-2],e_DC_z[region-2],")
//                      +(TString)("e_DC_x[region-3],e_DC_y[region-3],e_DC_z[region-3],")
//                      ));
//StreamToCSVfile({
//    (Double_t)status,
//    (Double_t)runnum,
//    (Double_t)event,
//    Ebeam,
//    (Double_t)beam_helicity,
//    e.E(),
//    e.P(),
//    e.Px(),
//    e.Py(),
//    e.Pz(),
//    e.Theta(),
//    e.Phi(),
//    (Double_t)Ngammas,
//    (Double_t)Np,
//    (Double_t)Nn,
//    (Double_t)Npips,
//    (Double_t)Npims,
//    q.E(),
//    q.P(),
//    q.Px(),
//    q.Py(),
//    q.Pz(),
//    xB,
//    Q2,
//    z,
//    W,
//    pi.E(),
//    pi.P(),
//    pi.Px(),
//    pi.Py(),
//    pi.Pz(),
//    pi.Theta(),
//    pi.Phi(),                                           // Ppips_t_q,          Ppips_q,
//    e_E_PCAL,
//    e_E_ECIN,
//    e_E_ECOUT,
//    Ve.X(),
//    Ve.Y(),
//    Ve.Z(),
//    Vpi.X(),
//    Vpi.Y(),
//    Vpi.Z(),
//    pi_chi2PID,
//    e_PCAL_W,
//    e_PCAL_V,
//    pi_PCAL_W,
//    pi_PCAL_V,
//    e_PCAL_x,
//    e_PCAL_y,
//    e_PCAL_z,
//    e_PCAL_sector,
//    e_DC_sector,
//    e_DC_Chi2N,
//    e_DC_x[0],
//    e_DC_y[0],
//    e_DC_z[0],
//    e_DC_x[1],
//    e_DC_y[1],
//    e_DC_z[1],
//    e_DC_x[2],
//    e_DC_y[2],
//    e_DC_z[2],
//}, IsSelectedEvent, fdebug);
// ------------------------------------------------------------------------------------------------

