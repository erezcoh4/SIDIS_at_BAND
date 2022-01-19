
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
void                     CloseOutputFiles (TString OutDataPath, TString outfilename);
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
void                        FinishProgram (TString outfilepath, TString outfilename);
void                   GetParticlesByType (int evnum, int fdebug );
void              Stream_e_pi_line_to_CSV (int fdebug );
void                          SetDataPath ( TString fDataPath )  {DataPath = fDataPath;}
void                         SetFileLabel ( TString fFileLabel ) {FileLabel = fFileLabel;}
void                          SetPiCharge ( TString fPiCharge )  {PiCharge = fPiCharge;}

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

// positive pions
bool     pipsPastSelectionCuts[NMAXPIONS];
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
bool     pimsPastSelectionCuts[NMAXPIONS];
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
TLorentzVector          e_g, pi_g;

// auxiliary
DCfid_SIDIS dcfid;
std::vector<region_part_ptr>  electrons, neutrons, protons, pipluses, piminuses, gammas, deuterons;


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
// Main functionality
// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void Read_PiAcceptance_GEMCimulations(int  NeventsMax=-1,
                                      int  fdebug=1,
                                      int  PrintProgress=50000,
                                      TString fDataPath = "/volatile/clas12/users/ecohen/GEMC/hipo/10.2/AcceptanceCorrection/",
                                      TString fFileLabel = "p_uniform_distribution",
                                      TString fPiCharge = "pips"
                                      ){
    
    SetDataPath ( fDataPath  );
    SetFileLabel( fFileLabel );
    SetPiCharge ( fPiCharge  );
    // open result files
    OpenResultFiles();
    
    TString inputFile = DataPath + "/" + PiCharge + "/ee" + PiCharge + "_" + FileLabel + ".hipo";
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
            
            // CONTINUE HERE:
            // add truth-information,
            // i.e. generated electron and generated pion information
            auto mcpbank = c12.mcparts();
            const Int_t Ngen=mcpbank->getRows();
            
            for( Int_t i_mc =0; i_mc< Ngen ; i_mc++){
              mcpbank -> setEntry(i_mc);
              
              P_mc_particle.SetXYZM( mcpbank->getPx() , mcpbank->getPy() , mcpbank->getPz() , mcpbank->getMass() );
              auto pid = mcpbank->getPid();

              cout <<
                "MC particle "  << i_mc     << " " << pid    <<
                " p4 = "        << p4.X()   << "," << p4.Y() << "," << p4.Z() << "," << p4.T()
                << ", mass "    << p4.M()
                << endl;
              
            }
            
            
            // do not filter any events here, extract information from all
            // we keep only d(e,e’pi+)X and d(e,e’pi-)X events
            
            ExtractElectronInformation  (fdebug);
            ExtractPionsInformation     (fdebug);
            WriteEventToOutput          (fdebug);
            
            if (fdebug>1) {
                std::cout << "done processing event " << evnum
                << " (" << event << "/" << NeventsMaxToProcess<< ") "
                << std::endl << "------------------------------------------------------------" << std::endl ;
            }
            event++; Nevents_processed++;
            if (fdebug && event%PrintProgress==0) std::cout << std::setprecision(1) << " event " << event << std::endl;
        } // end event loop
        
    } // end file loop
    
    
    FinishProgram( outfilepath, outfilename);
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
void CloseOutputFiles (TString OutDataPath, TString outfilename){
    // close output CSV
    CSVfile_e_piplus                .close();
    SelectedEventsCSVfile_e_piplus  .close();
    SelectedEventsCSVfile_e_piplus_kinematics  .close();
    CSVfile_e_piminus               .close();
    SelectedEventsCSVfile_e_piminus .close();
    SelectedEventsCSVfile_e_piminus_kinematics .close();

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
    << "\t" << (float)Nevents_passed_e_pips_cuts/Nevents_processed  << " passed (e,e'pi+) cuts,"
    << std::endl
    << "\t\t" << (float)Nevents_passed_e_pips_kinematics_cuts/Nevents_processed  << " also passed kinematical cuts,"
    << std::endl
    << (float)Nevents_passed_pims_cuts/Nevents_processed    << " events passed pi- cuts,"
    << std::endl
    << "\t" << (float)Nevents_passed_e_pims_cuts/Nevents_processed  << " passed (e,e'pi-) cuts,"
    << std::endl
    <<  "\t\t" << (float)Nevents_passed_e_pims_kinematics_cuts/Nevents_processed  << " also passed kinematical cuts,"
    << std::endl;
    
    
    
    std::cout << "output files ready in root/csv formats in " << std::endl
    << std::endl
    << "wrote "  << Nentires_e_piplus  << " to (e,e'pi+) root file, "
    << std::endl << outFile_e_piplus -> GetName()
    << std::endl << OutDataPath + outfilename + "_e_piplus_selected_*.csv"
    << std::endl
    << "and "    << Nentires_e_piminus << " to (e,e'pi-) root file. "
    << std::endl << outFile_e_piminus -> GetName()
    << std::endl << OutDataPath + outfilename + "_e_piminus_selected_*.csv"
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
    if (fdebug>1) {
        std::cout << "InitializeFileReading( " << NeventsMax << " , " << c12Nentries << " , " << fdebug << ")" << std::endl;
    }
    Ebeam   = GetBeamEnergy( fdebug );
    Beam    .SetPxPyPzE (0, 0, Ebeam, Ebeam );
    target  .SetXYZM    (0, 0, 0,     Md    );
    
    NeventsMaxToProcess = NeventsMax;
    if (NeventsMax<0) NeventsMaxToProcess = c12Nentries;
    Nevents_processed           = 0;
    Nevents_passed_e_cuts       = 0;
    Nevents_passed_pips_cuts    = 0;
    Nevents_passed_pims_cuts    = 0;
    Nevents_passed_e_pips_cuts  = 0;
    Nevents_passed_e_pims_cuts  = 0;
    Nevents_passed_e_pips_kinematics_cuts = 0;
    Nevents_passed_e_pims_kinematics_cuts = 0;
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
    ePastCutsInEvent                    = false;

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
        eepipsPastKinematicalCuts[piIdx]            = false;
        
        piminus .push_back( TLorentzVector(0,0,0,db->GetParticle( -211 )->Mass()) );
        Vpiminus.push_back( TVector3() );
        pimsPastSelectionCuts[piIdx]                = false;
        eepimsPastKinematicalCuts[piIdx]            = false;
        
        piplus_Px[piIdx]    = piplus_Py[piIdx]  = piplus_Pz[piIdx]  = piplus_E[piIdx]   = -9999;
        piminus_Px[piIdx]   = piminus_Py[piIdx] = piminus_Pz[piIdx] = piminus_E[piIdx]  = -9999;
        Vpiplus_X[piIdx]    = Vpiplus_Y[piIdx]  = Vpiplus_Z[piIdx]  = -9999;
        Vpiminus_X[piIdx]   = Vpiminus_Y[piIdx] = Vpiminus_Z[piIdx] = -9999;
         
    }
    DC_layer                                        = -9999;
    status                                          = 1; // 0 is good...
    
    pipsPastCutsInEvent                             = false;
    eepipsPastCutsInEvent                           = false;
    pimsPastCutsInEvent                             = false;
    eepimsPastCutsInEvent                           = false;
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
                     +(TString)"Npips,Npims,Nelectrons,Ngammas,Nprotons,Nneutrons,Ndeuterons,"));
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
    ePastCutsInEvent = CheckIfElectronPassedSelectionCuts(e_PCAL_x, e_PCAL_y,
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
    if (fdebug>3) std::cout << "Writing (e,e'pi) event" << std::endl;
    Stream_e_pi_line_to_CSV( fdebug );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void FinishProgram(TString outfilepath, TString outfilename){
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    
    std::cout << "Done. Elapsed time: " << elapsed.count() << std::endl;
    CloseOutputFiles( outfilepath, outfilename );
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
                                                                     Vpiplus[pipsIdx],
                                                                     fdebug);
    eepipsPastKinematicalCuts[pipsIdx] = eepiPassedKinematicalCriteria(piplus[pipsIdx],
                                                                       fdebug);
    if (pipsPastSelectionCuts[pipsIdx]) {
        pipsPastCutsInEvent = true;
        Nevents_passed_pips_cuts ++;
        if (eepipsPastKinematicalCuts[pipsIdx]) {
            eepipsPastCutsInEvent = true;
        }
    }
    
    piplus_Px[pipsIdx]          = piplus[pipsIdx].Px();
    piplus_Py[pipsIdx]          = piplus[pipsIdx].Py();
    piplus_Pz[pipsIdx]          = piplus[pipsIdx].Pz();
    piplus_E[pipsIdx]           = piplus[pipsIdx].E();
    Vpiplus_X[pipsIdx]          = Vpiplus[pipsIdx].X();
    Vpiplus_Y[pipsIdx]          = Vpiplus[pipsIdx].Y();
    Vpiplus_Z[pipsIdx]          = Vpiplus[pipsIdx].Z();
    
    
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
    pimsPastSelectionCuts[pimsIdx] = CheckIfPionPassedSelectionCuts("pi-",
                                                                     pims_DC_sector[pimsIdx],
                                                                     pims_DC_x[pimsIdx],
                                                                     pims_DC_y[pimsIdx],
                                                                     pims_DC_z[pimsIdx],
                                                                     pims_chi2PID[pimsIdx],  piminus[pimsIdx].P(),
                                                                     Ve,
                                                                     Vpiminus[pimsIdx],
                                                                     fdebug);
    eepimsPastKinematicalCuts[pimsIdx] = eepiPassedKinematicalCriteria(piminus[pimsIdx],
                                                                       fdebug);
    if (pimsPastSelectionCuts[pimsIdx]) {
        pimsPastCutsInEvent = true;
        Nevents_passed_pims_cuts ++;
        if (eepimsPastKinematicalCuts[pimsIdx]) {
            eepimsPastCutsInEvent = true;
        }
    }
    
    piminus_Px[pimsIdx]          = piminus[pimsIdx].Px();
    piminus_Py[pimsIdx]          = piminus[pimsIdx].Py();
    piminus_Pz[pimsIdx]          = piminus[pimsIdx].Pz();
    piminus_E[pimsIdx]           = piminus[pimsIdx].E();
    Vpiminus_X[pimsIdx]          = Vpiminus[pimsIdx].X();
    Vpiminus_Y[pimsIdx]          = Vpiminus[pimsIdx].Y();
    Vpiminus_Z[pimsIdx]          = Vpiminus[pimsIdx].Z();
    

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
void Stream_e_pi_line_to_CSV( int fdebug ){
    int piIdx = 0;
    TLorentzVector  pi;
    TVector3        Vpi;
    double          Zpi;
    if (piCharge=="pips") {
        pi  = piplus [piIdx];
        Vpi = Vpiplus[piIdx];
    }
    else if (piCharge=="pims") {
        pi  = piminus [piIdx];
        Vpi = Vpiminus[piIdx];
   }
    else {
        std::cout << "pion charge undefined at Stream_e_pi_line_to_CSV(), returning " << std::endl;
        return;
    }
    // write a (e,e'pi) event-line to CSV file
    std::vector<double> variables =
    {   e.P(),          e.Theta(),          e.Phi(),            Ve.Z(),
        pi.P(),         pi.Theta(),         pi.Phi(),           Vpi.Z(),
        (double)Npips, (double)Npims,       (double)Ne,         (double)Ngammas,
        (double)Np,    (double)Nn,          (double)Nd,
    };
    StreamToCSVfile( variables, fdebug );
}



