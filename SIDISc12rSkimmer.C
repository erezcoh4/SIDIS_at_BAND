// last edit Aug-11, 2021 (EOC, mbp), see README
// ToDo:
//(1) Check addition of pi-
//(2) Check beam helicity
//(3) Add pion DC fiducial cuts


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
#include "/u/home/cohen/SIDIS_at_BAND/Auxiliary/DC_fiducial.cpp"
#include "/u/home/cohen/SIDIS_at_BAND/Auxiliary/csv_reader.h"
using namespace clas12;



// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
// declare methods
TVector3                GetParticleVertex (clas12::region_part_ptr rp);
void                     SetLorentzVector (TLorentzVector &p4, clas12::region_part_ptr rp);
void                      OpenOutputFiles (TString csvfilename, TString header);
void                     CloseOutputFiles (TString OutDataPath);
void                      StreamToCSVfile (std::vector<Double_t> observables,
                                           bool IsSelectedEvent,
                                           int fdebug);
void             StreamVariablesToCSVfile (std::ofstream CSVfile,
                                           std::ofstream SelectedEventsCSVfile,
                                           std::vector<Double_t> observables,
                                           bool IsSelectedEvent);
void                      ChangeAxesFrame (TString FrameName="q(z) frame");
void                        MoveTo_qFrame ();
void                       printCutValues ();
void                        loadCutValues (TString cutValuesFilename = "cutValues.csv", int fdebug=0);
void                      SetOutputTTrees ();
double                       FindCutValue ( std::string cutName );
bool EventPassedElectronSelectionCriteria (Double_t e_PCAL_x, Double_t e_PCAL_y,
                                           Double_t e_PCAL_W,Double_t e_PCAL_V,
                                           Double_t e_E_PCAL,
                                           Double_t e_E_ECIN, Double_t e_E_ECOUT,
                                           TLorentzVector e,
                                           TVector3 Ve,
                                           Double_t e_DC_sector,
                                           Double_t e_DC_x[3],
                                           Double_t e_DC_y[3],
                                           int torusBending);
bool EventPassedPionSelectionCutsCriteria (Double_t DC_x[3], Double_t DC_y[3],
                                           Double_t chi2PID, Double_t p,
                                           TVector3 Ve,      TVector3 Vpi );
Double_t          Chi2PID_pion_lowerBound (Double_t p, Double_t C=0.88); // C(pi+)=0.88, C(pi-)=0.93
Double_t          Chi2PID_pion_upperBound (Double_t p, Double_t C=0.88); // C(pi+)=0.88, C(pi-)=0.93
int                       GetBeamHelicity (clas12reader c12, int runnum );

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
// globals
// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
// cut values
std::vector<std::pair<std::string, double>> cutValues;
double               cutValue_Vz_min;
double               cutValue_Vz_max;
double             cutValue_e_PCAL_W;
double             cutValue_e_PCAL_V;
double             cutValue_e_E_PCAL;
double cutValue_SamplingFraction_min;
double        cutValue_Ve_Vpi_dz_max;


bool      ePastSelectionCuts = false;
bool     piPastSelectionCuts = false;
bool         EventPassedCuts = false;

// meta-data
int           torusBending = -1; // -1 for In-bending, +1 for Out-bending
int    DC_layers[3] = {6,18,36};// Region 1 is denoted at DC detector 6, Region 2 is denoted 18, Region 3 - as 36
int                    DC_layer;
int                      runnum;
int                       evnum;
int                  good_event;
int               beam_helicity; // helicity of the electron +1 along the beam and -1 opposite to it
int                      status;
//    int Fastest_pipsIdx = 0;
// define the leading jet as the one with greatest z
double            z_max_pi;
// variables
double                  Mp;
double                 Mp2;
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
// leading pion (pi+ or pi-)
double          pi_chi2PID;
double           pi_PCAL_W;
double           pi_PCAL_V;
double           pi_PCAL_x;
double           pi_PCAL_y;
double           pi_PCAL_z;
double      pi_PCAL_sector;
double        pi_DC_sector;
double            pi_Chi2N;
double          pi_DC_x[3];
double          pi_DC_y[3];
double           pi_E_PCAL;
double           pi_E_ECIN;
double          pi_E_ECOUT;
// leading pi+
double        pips_chi2PID;
double         pips_PCAL_W;
double         pips_PCAL_V;
double         pips_PCAL_x;
double         pips_PCAL_y;
double         pips_PCAL_z;
double    pips_PCAL_sector;
double      pips_DC_sector;
double          pips_Chi2N;
double        pips_DC_x[3];
double        pips_DC_y[3];
double         pips_E_PCAL;
double         pips_E_ECIN;
double        pips_E_ECOUT;
// leading pi-
double        pims_chi2PID;
double         pims_PCAL_W;
double         pims_PCAL_V;
double         pims_PCAL_x;
double         pims_PCAL_y;
double         pims_PCAL_z;
double    pims_PCAL_sector;
double      pims_DC_sector;
double          pims_Chi2N;
double        pims_DC_x[3];
double        pims_DC_y[3];
double         pims_E_PCAL;
double         pims_E_ECIN;
double        pims_E_ECOUT;



// Output root file and tree
TFile * outFile_e_piplus, * outFile_e_piminus;
TTree * outTree_e_piplus, * outTree_e_piminus;
// Output CSV file
std::ofstream   CSVfile_e_piplus,  SelectedEventsCSVfile_e_piplus;
std::ofstream   CSVfile_e_piminus, SelectedEventsCSVfile_e_piminus;
// vectors in lab-frame
TLorentzVector             Beam;
TLorentzVector                e;
TLorentzVector                q;
TLorentzVector               pi; // leading pion
TLorentzVector           piplus; // leading positive pion
TLorentzVector          piminus; // leading negative pion
TLorentzVector    piplus_leader; // leading positive pion
TLorentzVector   piminus_leader; // leading negative pion
TString       LeadingPionCharge; // "piplus" / "piminus"
// reconstructed vertex position
TVector3                     Ve;
TVector3                    Vpi;
TVector3                Vpiplus;
TVector3               Vpiminus;

// kinematics
Double_t           Ebeam = 10.2; // [GeV] ( for Fall-2019 the enrgy was 10.4096)
// beam energy was initiated here, but is later overridden with rcdbData information
Double_t                     xB;
Double_t                     Q2;
Double_t                  omega;
Double_t                      z; // energy fraction rest frame
Double_t                     W2; // squared hadronic c.m. energy
Double_t                      W; // hadronic c.m. energy


// vectors in q-frame
//TLorentzVector       piplus_qFrame;
//Double_t        Ppips_t_q, Ppips_q;
// auxiliary
DCFiducial dcfid;





// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
// Main functionality
// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void SIDISc12rSkimmer(int  RunNumber=6420,
                      int  NeventsMax=-1,
                      int  fdebug=1,
                      bool doApplySelectionCuts=true,
                      int  PrintProgress=5000,
                      int NpipsMin=1, // minimal number of pi+
                      TString DataPath = "/volatile/clas12/rg-b/production/recon/spring2019/torus-1/pass1/v0/dst/train_20200610/inc/" ){
    
    char RunNumberStr[20];
    sprintf( RunNumberStr, "00%d", RunNumber );
    if (fdebug>1) {
        std::cout << "(SIDIS) skimming run " << RunNumberStr << std::endl;
    }
    
    loadCutValues("cutValues.csv",fdebug);
    
    

    
    // defenitions
    auto db = TDatabasePDG::Instance();
    Mp      = db->GetParticle( 2212 )->Mass();
    Mp2     = Mp * Mp;
    e       = TLorentzVector(0,0,0,db->GetParticle( 11   )->Mass());
    pi      = TLorentzVector(0,0,0,db->GetParticle( 211  )->Mass());
    piplus  = TLorentzVector(0,0,0,db->GetParticle( 211  )->Mass());
    piminus = TLorentzVector(0,0,0,db->GetParticle(-211  )->Mass());
    Ve      = TVector3();
    Vpi     = TVector3();
    Vpiplus = TVector3();
    Vpiminus= TVector3();
    
    
    
    // open CSV file
    OpenOutputFiles ("/volatile/clas12/users/ecohen/BAND/SIDIS_skimming/skimmed_SIDIS_inc_"
                     + (TString)RunNumberStr ,
                     ( (TString)("status,runnum,evnum,")
                      +(TString)("Ebeam,beam_helicity,")
                      +(TString)("Ee,Pe,Pe_x,Pe_y,Pe_z,Pe_theta,Pe_phi,")
                      +(TString)("Ngammas,Np,Nn,Npips,Npims,")
                      +(TString)("omega,q,q_x,q_y,q_z,xB,Q2,z,W,")
                      +(TString)("Epi,Ppi,Ppi_x,Ppi_y,Ppi_z,")
                      +(TString)("Ppi_theta,Ppi_phi,")                              //                      +(TString)("Ppi_t_q,Ppips_q,")
                      +(TString)("e_E_PCAL,e_E_ECIN,e_E_ECOUT,")
                      +(TString)("Vx_e,Vy_e,Vz_e,Vx_pi,Vy_pi,Vz_pi,")
                      +(TString)("pi_chi2PID,")
                      +(TString)("e_PCAL_W,e_PCAL_V,pips_PCAL_W,pips_PCAL_V,")
                      +(TString)("e_PCAL_x,e_PCAL_y,e_PCAL_z,e_PCAL_sector,")
                      +(TString)("e_DC_sector,e_DC_Chi2N,")
                      +(TString)("e_DC_x[region-1],e_DC_y[region-1],")
                      +(TString)("e_DC_x[region-2],e_DC_y[region-2],")
                      +(TString)("e_DC_x[region-3],e_DC_y[region-3],")
                      ));
    
    // output tree branches
    SetOutputTTrees();
    
    // ------------------------------------------------------------
    // open input file(s)
    // ------------------------------------------------------------
    // Record start time
    auto start = std::chrono::high_resolution_clock::now();
    TString inputFile = (TString)DataPath + "inc_" + (TString)RunNumberStr + ".hipo";
    for(Int_t i=1;i<gApplication->Argc();i++){
        TString opt=gApplication->Argv(i);
        if((opt.Contains(".hipo"))){
            inputFile=opt(5,opt.Sizeof());
        }
    }
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
    // ------------------------------------------------------------

    
    
    
    
    // step over events and extract information....
    for(Int_t i=0;i<files->GetEntries();i++){
        
        //create the event reader
        clas12reader c12(files->At(i)->GetTitle(),{0});
        
        //rcdb info
        auto& rcdbData = c12.rcdb()->current();//struct with all relevent rcdb values
        
        // get beam energy
        Ebeam = rcdbData.beam_energy ;
        Beam.SetPxPyPzE(0,0,Ebeam,Ebeam);
        
        
        // first count number of events - I don't know how to read it from the HIPO file
        if (fdebug) std::cout << "reading file " << i << std::endl;
        
        Int_t NeventsMaxToProcess = NeventsMax;
        if (NeventsMax<0){
            NeventsMaxToProcess = c12.getReader().getEntries();
        }
        int event   = 0;
        int good_event  = 0;

            
        // now process the events from the first one...
        while((c12.next()==true) && (event < NeventsMaxToProcess)){
            
            runnum = c12.runconfig()->getRun();
            evnum  = c12.runconfig()->getEvent();
            beam_helicity = GetBeamHelicity(c12, runnum);
            
            
            if (fdebug>2) std::cout << "begin analysis of event " << event << std::endl;
            
            // initialize
            xB          = Q2        = omega     = -9999;
            e_E_ECIN    = e_E_ECOUT = e_E_PCAL  = -9999;
            pips_chi2PID                        = -9999;
            e_PCAL_W    = e_PCAL_V              = -9999;
            e_PCAL_x    = e_PCAL_y  = e_PCAL_z  = -9999;
            e_PCAL_sector                       = -9999;
            e_DC_sector = e_DC_Chi2N            = -9999;
            
            pi_DC_sector                        = -9999;
            pi_PCAL_sector                      = -9999;
            pi_PCAL_W   = pi_PCAL_V             = -9999;
            pi_PCAL_x   = pi_PCAL_y             = -9999;
            pi_PCAL_z                           = -9999;
            pi_E_PCAL                           = -9999;
            pi_E_ECIN   = pi_E_ECOUT            = -9999;
            
            pips_DC_sector                      = -9999;
            pips_PCAL_sector                    = -9999;
            pips_PCAL_W = pips_PCAL_V           = -9999;
            pips_PCAL_x = pips_PCAL_y           = -9999;
            pips_PCAL_z                         = -9999;
            pips_E_PCAL                         = -9999;
            pips_E_ECIN = pips_E_ECOUT          = -9999;
            
            pims_DC_sector                      = -9999;
            pims_PCAL_sector                    = -9999;
            pims_PCAL_W = pims_PCAL_V           = -9999;
            pims_PCAL_x = pims_PCAL_y           = -9999;
            pims_PCAL_z                         = -9999;
            pims_E_PCAL                         = -9999;
            pims_E_ECIN = pims_E_ECOUT          = -9999;

            DC_layer                            = -9999;
            Ve          = Vpi                   = TVector3();
            Vpiplus     = Vpiminus              = TVector3();
            
            for (int regionIdx=0; regionIdx<3; regionIdx++) {
                e_DC_x[regionIdx]   = e_DC_y[regionIdx]     = -9999;
                pips_DC_x[regionIdx]= pips_DC_y[regionIdx]  = -9999;
                pims_DC_x[regionIdx]= pims_DC_y[regionIdx]  = -9999;
                pi_DC_x[regionIdx]  = pi_DC_y[regionIdx]    = -9999;
            }
            ePastSelectionCuts  = piPastSelectionCuts       = false;
            LeadingPionCharge                   = "piplus"; // arbitrary initialisation
            status                              = 1; // 0 is good...
            
            
            //can get an estimate of the beam current to this event
            c12.getCurrApproxCharge(); //if called c12.scalerReader();
            
            // get particles by type
            auto  electrons = c12.getByID( 11   ); auto         Ne = electrons.size();
            auto   neutrons = c12.getByID( 2112 ); auto         Nn = neutrons.size();
            auto    protons = c12.getByID( 2212 ); auto         Np = protons.size();
            auto       pips = c12.getByID( 211  ); auto      Npips = pips.size();
            auto       pims = c12.getByID(-211  ); auto      Npims = pims.size();
            auto     gammas = c12.getByID( 22   ); auto    Ngammas = gammas.size();
            
            //Loop over all particles to see how to access detector info.
            if (fdebug>1){
                std::cout
                << "number of particles in event "  << event        << " : " << c12.getNParticles() << ", "
                << "N(electrons): "                 << Ne           <<  ","
                << "N(protons): "                   << Np           <<  ","
                << "N(neutrons): "                  << Nn           <<  ","
                << "N(pi+): "                       << Npips        <<  ","
                << "N(pi-): "                       << Npims        <<  ","
                << "N(gammas): "                    << Ngammas      <<  ","
                << std::endl;
            }
            
            
            // filter events, extract information, and compute event kinematics:
            // we keep only d(e,e’pi+)X and d(e,e’pi-)X events
            if(  Ne == 1
               &&
               (Npips > 0 || Npims > 0) ){
                   
                // ------------------------------------------------------------------------------------------------
                // extract electron information
                // ------------------------------------------------------------------------------------------------
                // set electron 4-momentum
                SetLorentzVector(e,electrons[0]);
                // set electron vertex
                Ve      = GetParticleVertex( electrons[0] );
                // compute event kinematics
                q       = Beam - e;
                Q2      = -q.Mag2();
                omega   = q.E();
                xB      = Q2/(2. * Mp * q.E());
                W2      = Mp2 - Q2 + 2. * omega * Mp;
                W       = sqrt(W2);
                
                
                // detector information on electron
                auto e_PCAL_info= electrons[0]->cal(PCAL);
                e_E_PCAL        = e_PCAL_info->getEnergy();
                e_PCAL_sector   = e_PCAL_info->getSector();
                e_PCAL_V        = e_PCAL_info->getLv();
                e_PCAL_W        = e_PCAL_info->getLw();
                e_E_ECIN        = electrons[0]->cal(ECIN)->getEnergy();
                e_E_ECOUT       = electrons[0]->cal(ECOUT)->getEnergy();
                
                // hit position in PCAL
                e_PCAL_x        = e_PCAL_info->getX();
                e_PCAL_y        = e_PCAL_info->getY();
                e_PCAL_z        = e_PCAL_info->getZ();
                
                // Drift Chamber tracking system
                auto e_DC_info  = electrons[0]->trk(DC);
                e_DC_sector     = e_DC_info->getSector(); // tracking sector
                e_DC_Chi2N      = e_DC_info->getChi2N();  // tracking chi^2/NDF
                
                for (int regionIdx=0; regionIdx<3; regionIdx++) {
                    int DC_layer = DC_layers[regionIdx];
                    e_DC_x[regionIdx] = electrons[0]->traj(DC,DC_layer)->getX();
                    e_DC_y[regionIdx] = electrons[0]->traj(DC,DC_layer)->getY();
                    
                }
                if (fdebug > 2) std::cout << "extracted electron information and computed kinematics" << std::endl;
                // ------------------------------------------------------------------------------------------------
                
                
                
                // ------------------------------------------------------------------------------------------------
                // extract pion information
                // ------------------------------------------------------------------------------------------------
                // select the fastest pion as the "leader"
                double     z_max_piminus = 0;
                double      z_max_piplus = 0;

                // pi+ (positive pions)
                if (Npips > 0) {
                    int  z_max_piplus_index = 0;
                    SetLorentzVector(piplus  ,pips[0]);
                    TLorentzVector piplus_tmp(0,0,0,db->GetParticle(211)->Mass());
                    for (int pipsIdx=0; pipsIdx < Npips; pipsIdx++) {
                        SetLorentzVector(piplus_tmp  ,pips[pipsIdx]);
                        z = piplus_tmp.E()/q.E();
                        if (z > z_max_piplus) {
                            z_max_piplus_index = pipsIdx;
                            z_max_piplus       = z;
                        }
                    }
                    SetLorentzVector(piplus , pips[z_max_piplus_index]);
                    Vpiplus = GetParticleVertex( pips[z_max_piplus_index] );
                    pips_chi2PID = pips[z_max_piplus_index]->par()->getChi2Pid();
                    
                    auto pips_PCAL_info= pips[z_max_piplus_index]->cal(PCAL);
                    pips_E_PCAL        = pips_PCAL_info->getEnergy();
                    pips_PCAL_sector   = pips_PCAL_info->getSector();
                    pips_PCAL_V        = pips_PCAL_info->getLv();
                    pips_PCAL_W        = pips_PCAL_info->getLw();
                    pips_E_ECIN        = pips[z_max_piplus_index]->cal(ECIN)->getEnergy();
                    pips_E_ECOUT       = pips[z_max_piplus_index]->cal(ECOUT)->getEnergy();
                    pips_PCAL_x        = pips_PCAL_info->getX();
                    pips_PCAL_y        = pips_PCAL_info->getY();
                    pips_PCAL_z        = pips_PCAL_info->getZ();
                    auto pips_DC_info  = pips[z_max_piplus_index]->trk(DC);
                    pips_DC_sector     = pips_DC_info->getSector(); // tracking sector
                    pips_Chi2N         = pips_DC_info->getChi2N();  // tracking chi^2/NDF
                    for (int regionIdx=0; regionIdx<3; regionIdx++) {
                        DC_layer = DC_layers[regionIdx];
                        pips_DC_x[regionIdx] = pips[z_max_piplus_index]->traj(DC,DC_layer)->getX();
                        pips_DC_y[regionIdx] = pips[z_max_piplus_index]->traj(DC,DC_layer)->getY();
                    }
                }
                
                // pi- (negative pions)
                if (Npims > 0) {
                    int  z_max_piminus_index = 0;
                    SetLorentzVector(piminus  ,pims[0]);
                    TLorentzVector piminus_tmp(0,0,0,db->GetParticle(211)->Mass());
                    for (int pimsIdx=0; pimsIdx < Npims; pimsIdx++) {
                        SetLorentzVector(piminus_tmp  ,pims[pimsIdx]);
                        z = piminus_tmp.E()/q.E();
                        if (z > z_max_piminus) {
                            z_max_piminus_index = pimsIdx;
                            z_max_piminus       = z;
                        }
                    }
                    SetLorentzVector(piminus , pims[z_max_piminus_index]);
                    Vpiminus = GetParticleVertex( pims[z_max_piminus_index] );
                    pims_chi2PID = pims[z_max_piminus_index]->par()->getChi2Pid();
                    
                    auto pims_PCAL_info= pims[z_max_piminus_index]->cal(PCAL);
                    pims_E_PCAL        = pims_PCAL_info->getEnergy();
                    pims_PCAL_sector   = pims_PCAL_info->getSector();
                    pims_PCAL_V        = pims_PCAL_info->getLv();
                    pims_PCAL_W        = pims_PCAL_info->getLw();
                    pims_E_ECIN        = pims[z_max_piminus_index]->cal(ECIN)->getEnergy();
                    pims_E_ECOUT       = pims[z_max_piminus_index]->cal(ECOUT)->getEnergy();
                    pims_PCAL_x        = pims_PCAL_info->getX();
                    pims_PCAL_y        = pims_PCAL_info->getY();
                    pims_PCAL_z        = pims_PCAL_info->getZ();
                    auto pims_DC_info  = pims[z_max_piminus_index]->trk(DC);
                    pims_DC_sector     = pims_DC_info->getSector(); // tracking sector
                    pims_Chi2N         = pims_DC_info->getChi2N();  // tracking chi^2/NDF
                    for (int regionIdx=0; regionIdx<3; regionIdx++) {
                        DC_layer = DC_layers[regionIdx];
                        pims_DC_x[regionIdx] = pims[z_max_piminus_index]->traj(DC,DC_layer)->getX();
                        pims_DC_y[regionIdx] = pims[z_max_piminus_index]->traj(DC,DC_layer)->getY();
                    }
                }
                
                // decide which is the leading pion (leading jet)
                // accoring to the z (fraction of momentum) it carries
                if (z_max_piplus >= z_max_piminus){
                    
                    // leading pi+
                    LeadingPionCharge = "piplus";
                    pi                = piplus;
                    Vpi               = Vpiplus;
                    
                    pi_chi2PID      = pips_chi2PID;
                    pi_E_PCAL       = pips_E_PCAL;
                    pi_PCAL_sector  = pips_PCAL_sector;
                    pi_PCAL_V       = pips_PCAL_V;
                    pi_PCAL_W       = pips_PCAL_W;
                    pi_E_ECIN       = pips_E_ECIN;
                    pi_E_ECOUT      = pips_E_ECOUT;
                    pi_PCAL_x       = pips_PCAL_x;
                    pi_PCAL_y       = pips_PCAL_y;
                    pi_PCAL_z       = pips_PCAL_z;
                    pi_DC_sector    = pips_DC_sector;
                    pi_Chi2N        = pips_Chi2N;
                    for (int regionIdx=0; regionIdx<3; regionIdx++) {
                        pi_DC_x[regionIdx] = pips_DC_x[regionIdx];
                        pi_DC_y[regionIdx] = pips_DC_y[regionIdx];
                    }
                } else {
                    
                    // leading pi-
                    LeadingPionCharge = "piminus";
                    pi                = piminus;
                    Vpi               = Vpiminus;
                    
                    pi_chi2PID      = pims_chi2PID;
                    pi_E_PCAL       = pims_E_PCAL;
                    pi_PCAL_sector  = pims_PCAL_sector;
                    pi_PCAL_V       = pims_PCAL_V;
                    pi_PCAL_W       = pims_PCAL_W;
                    pi_E_ECIN       = pims_E_ECIN;
                    pi_E_ECOUT      = pims_E_ECOUT;
                    pi_PCAL_x       = pims_PCAL_x;
                    pi_PCAL_y       = pims_PCAL_y;
                    pi_PCAL_z       = pims_PCAL_z;
                    pi_DC_sector    = pims_DC_sector;
                    pi_Chi2N        = pims_Chi2N;
                    for (int regionIdx=0; regionIdx<3; regionIdx++) {
                        pi_DC_x[regionIdx] = pims_DC_x[regionIdx];
                        pi_DC_y[regionIdx] = pims_DC_y[regionIdx];
                    }
                }
                
                // leading hadron rest frame energy fraction
                z = pi.E()/q.E();
                // done
                if (fdebug > 2){
                    std::cout << "selected the leading pion " << LeadingPionCharge << std::endl;
                }
                
                //                // temporarily fill pips 4-vector in q-frame
                //                piplus_qFrame = piplus;
                //                // move to the axes-frame where q is parallel to the z-direction,
                //                // and compute pion P(perp) and P(parallel) in this frame
                //                ChangeAxesFrame();
                //                Ppips_t_q   = piplus_qFrame.Pt();
                //                Ppips_q     = piplus_qFrame.Pz();
                //                if (fdebug > 2) Printf("Moved to q-Frame and computed pi+ momentum components");
                // ------------------------------------------------------------------------------------------------

                
                
                // ------------------------------------------------------------------------------------------------
                // decide if to write this event to "selected events csv-file"
                bool IsSelectedEvent = false;
                if ( doApplySelectionCuts ) {
                    
                    ePastSelectionCuts = EventPassedElectronSelectionCriteria(e_PCAL_x, e_PCAL_y,
                                                                              e_PCAL_W, e_PCAL_V,
                                                                              e_E_PCAL, e_E_ECIN,
                                                                              e_E_ECOUT,
                                                                              e, Ve,
                                                                              e_PCAL_sector, // e_PCAL_sector should be consistent with e_DC_sector
                                                                              e_DC_x, e_DC_y,
                                                                              torusBending );
                    
                    piPastSelectionCuts = EventPassedPionSelectionCutsCriteria(pi_DC_x,     pi_DC_y,
                                                                               pi_chi2PID,  pi.P(),
                                                                               Ve,          Vpi);
                    
                    EventPassedCuts = ( ePastSelectionCuts && piPastSelectionCuts );
                    if ( EventPassedCuts ) {
                        status          = 0;
                        IsSelectedEvent = true;
                        
                        if (LeadingPionCharge=="piplus") {
                            outTree_e_piplus -> Fill();
                        } else {
                            outTree_e_piminus -> Fill();
                        }
                    }
                }
                
                
                
                StreamToCSVfile({
                    (Double_t)status,
                    (Double_t)runnum,
                    (Double_t)event,
                    Ebeam,
                    (Double_t)beam_helicity,
                    e.E(),
                    e.P(),
                    e.Px(),
                    e.Py(),
                    e.Pz(),
                    e.Theta(),
                    e.Phi(),
                    (Double_t)Ngammas,
                    (Double_t)Np,
                    (Double_t)Nn,
                    (Double_t)Npips,
                    (Double_t)Npims,
                    q.E(),
                    q.P(),
                    q.Px(),
                    q.Py(),
                    q.Pz(),
                    xB,
                    Q2,
                    z,
                    W,
                    pi.E(),
                    pi.P(),
                    pi.Px(),
                    pi.Py(),
                    pi.Pz(),
                    pi.Theta(),
                    pi.Phi(),                                           // Ppips_t_q,          Ppips_q,
                    e_E_PCAL,
                    e_E_ECIN,
                    e_E_ECOUT,
                    Ve.X(),
                    Ve.Y(),
                    Ve.Z(),
                    Vpi.X(),
                    Vpi.Y(),
                    Vpi.Z(),
                    pi_chi2PID,
                    e_PCAL_W,
                    e_PCAL_V,
                    pi_PCAL_W,
                    pi_PCAL_V,
                    e_PCAL_x,
                    e_PCAL_y,
                    e_PCAL_z,
                    e_PCAL_sector,
                    e_DC_sector,
                    e_DC_Chi2N,
                    e_DC_x[0],
                    e_DC_y[0],
                    e_DC_x[1],
                    e_DC_y[1],
                    e_DC_x[2],
                    e_DC_y[2],
                }, IsSelectedEvent, fdebug);
                // ------------------------------------------------------------------------------------------------

                good_event ++ ;
                if (fdebug>2) {
                    std::cout
                    << "streamed to CSV file,"
                    << "Ee: "       << e.E()        << " GeV, "
                    << "E(pi): "    << pi.E()       << " GeV, "
                    << std::endl;
                }
            }
            
            if (fdebug>1) {
                std::cout << "done processing event " << event << std::endl << "------------------------------" << std::endl ;
            }
            event++;
            if (fdebug && event%PrintProgress==0) std::cout << std::setprecision(1) << " event " << event << std::endl;
        } // end event loop
    } // end file loop
    
    
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    
    std::cout << "Done. Elapsed time: " << elapsed.count() << std::endl;
    CloseOutputFiles("/volatile/clas12/users/ecohen/BAND/SIDIS_skimming/");
}




// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
bool EventPassedElectronSelectionCriteria(Double_t e_PCAL_x, Double_t e_PCAL_y,
                                          Double_t e_PCAL_W,Double_t e_PCAL_V,
                                          Double_t e_E_PCAL,
                                          Double_t e_E_ECIN, Double_t e_E_ECOUT,
                                          TLorentzVector e,
                                          TVector3 Ve,
                                          Double_t e_DC_sector,
                                          Double_t e_DC_x[3],
                                          Double_t e_DC_y[3],
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
        bool DC_fid  = dcfid.DC_e_fid(e_DC_x[regionIdx],
                                      e_DC_y[regionIdx],
                                      e_DC_sector,
                                      regionIdx+1,
                                      bending);
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
bool EventPassedPionSelectionCutsCriteria(Double_t DC_x[3], Double_t DC_y[3],
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
    // ToDo: Complete DC - fiducial cuts on DC! Need help from Dien with pion fiductial
    // cuts implementation in DC_fiducials.cpp
    //
    if(
       (DC_x[0]<-1e9 )
       ){
        return false;
    }
    
    double C;
    if (LeadingPionCharge=="piplus")
        C = 0.88;
    else if (LeadingPionCharge=="piminus")
        C = 0.93;
    
    if(
       // pi+ Identification Refinement - chi2PID vs. momentum
       ( Chi2PID_pion_lowerBound( p, C ) < chi2PID && chi2PID < Chi2PID_pion_upperBound( p , C ) )
       
       // Cut on the z-Vertex Difference Between Electrons and Hadrons.
       &&  ( fabs((Ve-Vpi).Z()) < cutValue_Ve_Vpi_dz_max )
       ) return true;
    
    return false;
}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
Double_t Chi2PID_pion_upperBound( Double_t p, Double_t C){
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
Double_t Chi2PID_pion_lowerBound( Double_t p, Double_t C){
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
    CSVfile_e_piplus.open( outfilename + "_e_piplus.csv" );
    CSVfile_e_piplus << header << std::endl;
    CSVfile_e_piminus.open( outfilename + "_e_piminus.csv" );
    CSVfile_e_piminus << header << std::endl;
    
    SelectedEventsCSVfile_e_piplus.open( outfilename + "_e_piplus_selected_events.csv" );
    SelectedEventsCSVfile_e_piplus << header << std::endl;
    SelectedEventsCSVfile_e_piminus.open( outfilename + "_e_pimius_selected_events.csv" );
    SelectedEventsCSVfile_e_piminus << header << std::endl;

}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void CloseOutputFiles (TString OutDataPath){
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
    
    std::cout << "output files ready in root/csv formats in " << std::endl
    << OutDataPath << std::endl
    << "wrote " << Nentires_e_piplus  << " to (e,e'pi+) root file, "
    << "and "   << Nentires_e_piminus << " to (e,e'pi-) root file. "
    << std::endl;
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void StreamToCSVfile (std::vector<Double_t> observables, bool IsSelectedEvent, int fdebug){
    if (fdebug>1) {
        std::cout << "streaming to CSVfile" << std::endl;
    }
    // decide which file to write...
    if (LeadingPionCharge=="piplus") {
        StreamVariablesToCSVfile (CSVfile_e_piplus,
                                  SelectedEventsCSVfile_e_piplus,
                                  observables,
                                  IsSelectedEvent);
    } else {
        StreamVariablesToCSVfile (CSVfile_e_piminus,
                                  SelectedEventsCSVfile_e_piminus,
                                  observables,
                                  IsSelectedEvent);
    }
}

void StreamVariablesToCSVfile (std::ofstream CSVfile,
                               std::ofstream SelectedEventsCSVfile,
                               std::vector<Double_t> observables, bool IsSelectedEvent){
    // write to file
    for (auto v:observables) {
        CSVfile << v << ",";
    }
    CSVfile << std::endl;
    
    if (IsSelectedEvent) {
        for (auto v:observables) {
            SelectedEventsCSVfile << v << ",";
        }
        SelectedEventsCSVfile << std::endl;
    }
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void ChangeAxesFrame(TString FrameName){
    if (FrameName == "q(z) frame")
        MoveTo_qFrame();
    return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MoveTo_qFrame(){
//    // move to a reference frame where
//    // q is the z axis
//    Double_t q_phi   = q.Phi();
//    Double_t q_theta = q.Theta();
//
//    q.RotateZ(-q_phi);
//    q.RotateY(-q_theta);
//
//    // rotate additional vectors to the same axes-frame
//    // and they reside on the x-z plane: v=(v_x,0,v_q)
//    piplus_qFrame.RotateZ(-q_phi);
//    piplus_qFrame.RotateY(-q_theta);
//    Double_t piplus_phi = piplus_qFrame.Phi();
//    piplus_qFrame.RotateZ(-piplus_phi);
    
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
int GetBeamHelicity( clas12reader c12, int runnum ){
    // get beam helicity (+1 along the beam and -1 opposite to it)
    // [Christopher Dilks <dilks@jlab.org>, email from Aug-5, 2021]
    // for more items [https://github.com/JeffersonLab/clas12root/blob/master/AccesssingBankDataInCpp.txt]
    
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
    
    beam_helicity = c12.event()->getHelicity();
    
    // we are working here on RGB data
    bool helFlip = true;
    if      (runnum>=11093 && runnum<=11283)    helFlip = false; // falls, 10.4 GeV period only
    else if (runnum>=11323 && runnum<=11571)    helFlip = false; // winter
    
    if (helFlip) {
        beam_helicity = -1 * beam_helicity;
    }
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
    outTree_e_piplus->Branch("pi_chi2PID"           ,&pi_chi2PID            );
    outTree_e_piplus->Branch("e_PCAL_W"             ,&e_PCAL_W              );
    outTree_e_piplus->Branch("e_PCAL_V"             ,&e_PCAL_V              );
    outTree_e_piplus->Branch("pi_PCAL_x"            ,&pi_PCAL_x             );
    outTree_e_piplus->Branch("pi_PCAL_y"            ,&pi_PCAL_y             );
    outTree_e_piplus->Branch("pi_PCAL_z"            ,&pi_PCAL_z             );
    outTree_e_piplus->Branch("e_PCAL_x"             ,&e_PCAL_x              );
    outTree_e_piplus->Branch("e_PCAL_y"             ,&e_PCAL_y              );
    outTree_e_piplus->Branch("e_PCAL_z"             ,&e_PCAL_z              );
    outTree_e_piplus->Branch("e_PCAL_sector"        ,&e_PCAL_sector         );
    outTree_e_piplus->Branch("e_DC_sector"          ,&e_DC_sector           );
    outTree_e_piplus->Branch("e_DC_Chi2N"           ,&e_DC_Chi2N            );
    outTree_e_piplus->Branch("e_DC_x"               ,&e_DC_x                );
    outTree_e_piplus->Branch("e_DC_y"               ,&e_DC_y                );
    outTree_e_piplus->Branch("pi_PCAL_sector"       ,&pi_PCAL_sector        );
    outTree_e_piplus->Branch("pi_DC_sector"         ,&pi_DC_sector          );
    outTree_e_piplus->Branch("pi_Chi2N"             ,&pi_Chi2N              );
    outTree_e_piplus->Branch("pi_DC_x"              ,&pi_DC_x               );
    outTree_e_piplus->Branch("pi_DC_y"              ,&pi_DC_y               );
    outTree_e_piplus->Branch("pi_E_PCAL"            ,&pi_E_PCAL             );
    outTree_e_piplus->Branch("pi_E_ECIN"            ,&pi_E_ECIN             );
    outTree_e_piplus->Branch("pi_E_ECIN"            ,&pi_E_ECIN             );
    outTree_e_piplus->Branch("pi_E_ECOUT"           ,&pi_E_ECOUT            );
    outTree_e_piplus->Branch("DC_layer"             ,&DC_layer              );
    outTree_e_piplus->Branch("e"                    ,&e                     );
    outTree_e_piplus->Branch("pi"                   ,&pi                    );
    outTree_e_piplus->Branch("Ve"                   ,&Ve                    );
    outTree_e_piplus->Branch("Vpi"                  ,&Vpi                   );
    outTree_e_piplus->Branch("Beam"                 ,&Beam                  );
    outTree_e_piplus->Branch("beam_helicity"        ,&beam_helicity         );
    outTree_e_piplus->Branch("q"                    ,&q                     );
    outTree_e_piplus->Branch("Ebeam"                ,&Ebeam                 );
    outTree_e_piplus->Branch("xB"                   ,&xB                    );
    outTree_e_piplus->Branch("Q2"                   ,&Q2                    );
    outTree_e_piplus->Branch("omega"                ,&omega                 );
    outTree_e_piplus->Branch("z"                    ,&z                     );
    outTree_e_piplus->Branch("W"                    ,&W                     );

    outTree_e_piplus->Branch("EventPassedCuts"      ,&EventPassedCuts       );
    outTree_e_piplus->Branch("ePastSelectionCuts"   ,&ePastSelectionCuts    );
    outTree_e_piplus->Branch("piPastSelectionCuts"  ,&piPastSelectionCuts   );
    
    
    // pi-
    outTree_e_piminus->Branch("eventnumber"         ,&evnum                 );
    outTree_e_piminus->Branch("runnum"              ,&runnum                );
    outTree_e_piminus->Branch("e_E_PCAL"            ,&e_E_PCAL              );
    outTree_e_piminus->Branch("e_E_ECIN"            ,&e_E_ECIN              );
    outTree_e_piminus->Branch("e_E_ECOUT"           ,&e_E_ECOUT             );
    outTree_e_piminus->Branch("pi_chi2PID"          ,&pi_chi2PID            );
    outTree_e_piminus->Branch("e_PCAL_W"            ,&e_PCAL_W              );
    outTree_e_piminus->Branch("e_PCAL_V"            ,&e_PCAL_V              );
    outTree_e_piminus->Branch("pi_PCAL_x"           ,&pi_PCAL_x             );
    outTree_e_piminus->Branch("pi_PCAL_y"           ,&pi_PCAL_y             );
    outTree_e_piminus->Branch("pi_PCAL_z"           ,&pi_PCAL_z             );
    outTree_e_piminus->Branch("e_PCAL_x"            ,&e_PCAL_x              );
    outTree_e_piminus->Branch("e_PCAL_y"            ,&e_PCAL_y              );
    outTree_e_piminus->Branch("e_PCAL_z"            ,&e_PCAL_z              );
    outTree_e_piminus->Branch("e_PCAL_sector"       ,&e_PCAL_sector         );
    outTree_e_piminus->Branch("e_DC_sector"         ,&e_DC_sector           );
    outTree_e_piminus->Branch("e_DC_Chi2N"          ,&e_DC_Chi2N            );
    outTree_e_piminus->Branch("e_DC_x"              ,&e_DC_x                );
    outTree_e_piminus->Branch("e_DC_y"              ,&e_DC_y                );
    outTree_e_piminus->Branch("pi_PCAL_sector"      ,&pi_PCAL_sector        );
    outTree_e_piminus->Branch("pi_DC_sector"        ,&pi_DC_sector          );
    outTree_e_piminus->Branch("pi_Chi2N"            ,&pi_Chi2N              );
    outTree_e_piminus->Branch("pi_DC_x"             ,&pi_DC_x               );
    outTree_e_piminus->Branch("pi_DC_y"             ,&pi_DC_y               );
    outTree_e_piminus->Branch("pi_E_PCAL"           ,&pi_E_PCAL             );
    outTree_e_piminus->Branch("pi_E_ECIN"           ,&pi_E_ECIN             );
    outTree_e_piminus->Branch("pi_E_ECIN"           ,&pi_E_ECIN             );
    outTree_e_piminus->Branch("pi_E_ECOUT"          ,&pi_E_ECOUT            );
    outTree_e_piminus->Branch("DC_layer"            ,&DC_layer              );
    outTree_e_piminus->Branch("e"                   ,&e                     );
    outTree_e_piminus->Branch("pi"                  ,&pi                    );
    outTree_e_piminus->Branch("Ve"                  ,&Ve                    );
    outTree_e_piminus->Branch("Vpi"                 ,&Vpi                   );
    outTree_e_piminus->Branch("Beam"                ,&Beam                  );
    outTree_e_piminus->Branch("beam_helicity"       ,&beam_helicity         );
    outTree_e_piminus->Branch("q"                   ,&q                     );
    outTree_e_piminus->Branch("Ebeam"               ,&Ebeam                 );
    outTree_e_piminus->Branch("xB"                  ,&xB                    );
    outTree_e_piminus->Branch("Q2"                  ,&Q2                    );
    outTree_e_piminus->Branch("omega"               ,&omega                 );
    outTree_e_piminus->Branch("z"                   ,&z                     );
    outTree_e_piminus->Branch("W"                   ,&W                     );

    outTree_e_piminus->Branch("EventPassedCuts"     ,&EventPassedCuts       );
    outTree_e_piminus->Branch("ePastSelectionCuts"  ,&ePastSelectionCuts    );
    outTree_e_piminus->Branch("piPastSelectionCuts" ,&piPastSelectionCuts   );
}
