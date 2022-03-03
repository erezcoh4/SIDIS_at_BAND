
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
#include <TClonesArray.h>
#include <TH1.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include "clas12reader.h"
#include "Auxiliary/csv_reader.h"

#include "NPiEvent.hh"
#include "TaggedSIDISCuts.C"
#include "skimIO.C"
#include "skimHelper.C"

#define NMAXPIONS 20 // maximal allowed number of pions
#define r2d 180./3.1415 // radians to degrees
using namespace clas12;



// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
// start clock
auto start = std::chrono::high_resolution_clock::now();
// declare methods
void                      SetOutputTTrees ();
bool      CheckIfElectronPassedSelectionCuts (Double_t e_PCAL_x, Double_t e_PCAL_y,
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
bool      CheckIfPionPassedSelectionCuts (TString pionCharge, // "pi+" or "pi-"
                                           Double_t DC_sector,
                                           Double_t DC_x[3], Double_t DC_y[3], Double_t DC_z[3],
                                           Double_t chi2PID, Double_t p,
                                           TVector3 Ve,      TVector3 Vpi,
                                           int fdebug);
bool        eepiPassedKinematicalCriteria (TLorentzVector pi,
                                           int fdebug);
Double_t          Chi2PID_pion_lowerBound (Double_t p, Double_t C=0.88); // C(pi+)=0.88, C(pi-)=0.93
Double_t          Chi2PID_pion_upperBound (Double_t p, Double_t C=0.88); // C(pi+)=0.88, C(pi-)=0.93
int                       GetBeamHelicity (event_ptr p_event, int runnum, int fdebug);
double                      GetBeamEnergy (int fdebug);
TString                   GetRunNumberSTR (int RunNumber, int fdebug);
void                InitializeFileReading (int NeventsMax,int c12Nentries, int fdebug);
void           ExtractElectronInformation (int fdebug);
void              ExtractPionsInformation (int fdebug);
void               ExtractPipsInformation (int pipsIdx, int fdebug );
void               ExtractPimsInformation (int pimsIdx, int fdebug );
void                    ComputeKinematics ();
void                   WriteEventToOutput (int fdebug);
void                        FinishProgram (TString outfilepath, TString outfilename);
void                   GetParticlesByType (int evnum, int fdebug );
void              Stream_e_pi_line_to_CSV (TString pionCharge, int piIdx,
                                           bool passed_cuts_e_pi,
                                           bool passed_cuts_e_pi_kinematics,
                                           int fdebug );
// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.

// globals
auto db = TDatabasePDG::Instance();
// cut values

int         NeventsMaxToProcess;
int           Nevents_processed;
int       Nevents_passed_e_cuts;
int    Nevents_passed_pips_cuts;
int  Nevents_passed_e_pips_cuts;
int Nevents_passed_e_pips_kinematics_cuts;
int    Nevents_passed_pims_cuts;
int  Nevents_passed_e_pims_cuts;
int Nevents_passed_e_pims_kinematics_cuts;

NPiEvent npi;
SIDISIO io;
TaggedSIDISCuts cuts;

// vectors in q-frame
//TLorentzVector       piplus_qFrame;
//Double_t        Ppips_t_q, Ppips_q;
// auxiliary
std::vector<region_part_ptr>  electrons, neutrons, protons, pipluses, piminuses, gammas, deuterons;



// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
// Main functionality
// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void MCSIDISc12rSkimmer(int  NeventsMax=-1,
                      int  fdebug=1,
                      int  PrintProgress=50000,
                      int NpipsMin=1, // minimal number of pi+
                      TString inputfilepath = "/volatile/clas12/rg-b/production/recon/spring2019/torus-1/pass1/v0/dst/train_20200610/inc/",
                      TString outfilepath = "/volatile/clas12/users/akiral/BAND/SIDIS_skimming/",
                      TString outfilename = "skimmed_SIDIS_inc",
                      int setInclusive=0 ){
    // read cut values
    npi = new NPiEvent();
    io = new SIDISIO();
    cuts = new TaggedSIDISCuts();

    loadCutValues(cuts,"BANDcutValues.csv",fdebug);

    inclusive = setInclusive;
    if (inclusive == 1) std::cout << "Running as inclusive" << std::endl;

    // open result files
    io.OpenResultFiles( outfilepath, outfilename );

    TChain fake("hipo");
    fake.Add(inputfilepath.Data());
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
            
            GetBeamHelicity    ( c12.event() , runnum, fdebug );
            npi.InitializeVariables();
            // Get Particles By Type
            electrons   = c12.getByID( 11   );
            neutrons    = c12.getByID( 2112 );
            protons     = c12.getByID( 2212 );
            pipluses    = c12.getByID( 211  );
            piminuses   = c12.getByID(-211  );
            gammas      = c12.getByID( 22   );
            deuterons   = c12.getByID( 1000010020 );
            GetParticlesByType ( evnum, fdebug );

            /*for (int j = 0; j < 30; j++){
                c12.mcparts()->setEntry(j);
                //std::cout << j << " " << c12.mcparts()->getPid() << " " << c12.mcparts()->getPx() << " " << c12.mcparts()->getVx() << " " << c12.mcparts()->getMass() << std::endl;
                //std::cout << "a" << std::endl;
            }*/
            
            
            // filter events, extract information, and compute event kinematics:
            // we keep only d(e,e’pi+)X and d(e,e’pi-)X events
            if(  0 < Ne // after studying some MC and data, we need to kill events with more than 1 electron
               &&
               (inclusive == 1 || (0 < Npips && Npips < NMAXPIONS) || (0 < Npims && Npims < NMAXPIONS)) ){
                   
                ExtractElectronInformation  (fdebug);
                npi.ComputeKinematics             ();
                ExtractPionsInformation     (fdebug);
                WriteEventToOutput          (fdebug);
                
            } else {
                if (fdebug>1) {
                    std::cout << "Skipped computations in this event as there are not enough particles: "
                    << "Ne = " << Ne << ",Npips = " << Npips << ",Npims = " << Npims << std::endl ;
                }
            }
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


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
int GetBeamHelicity( NPiEvent npi, event_ptr p_event, int runnum, int fdebug ){
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

    npi.SetBeamHelicity(beam_helicity)

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
void WriteEventToOutput(int fdebug){
    // (Maybe) write this event to "selected events csv-file"
    bool            IsSelected_eepi = false;

    // Transfer pi vector data to TClonesArrays
    for (int i = 0; i < 20; i++) {
        new ((*piplusArray)[i]) TLorentzVector;
        (*piplusArray)[i] = &piplus[i];

        new ((*VpiplusArray)[i]) TVector3;
        (*VpiplusArray)[i] = &Vpiplus[i];

        new ((*piminusArray)[i]) TLorentzVector;
        (*piminusArray)[i] = &piminus[i];

        new ((*VpiminusArray)[i]) TVector3;
        (*VpiminusArray)[i] = &Vpiminus[i];
    }

    //ePastCutsInEvent = true;
    if (inclusive == 1) {
        pipsPastCutsInEvent = true;
        pimsPastCutsInEvent = true;
        eepipsPastCutsInEvent = true;
        eepimsPastCutsInEvent = true;
    }
    
    if (ePastCutsInEvent && pipsPastCutsInEvent) {
        IsSelected_eepi = true;
        outTree_e_piplus -> Fill();
        if (fdebug>3) std::cout << "Filling (e,e'pi+) TTree with this event!" << std::endl;
        
        Nevents_passed_e_pips_cuts ++ ;
        if (eepipsPastCutsInEvent) Nevents_passed_e_pips_kinematics_cuts ++;
        
        for (int pipsIdx=0; pipsIdx<Npips; pipsIdx++) {
            Stream_e_pi_line_to_CSV( "pi+", pipsIdx,
                                    pipsPastSelectionCuts[pipsIdx], eepipsPastKinematicalCuts[pipsIdx],
                                    fdebug );
        }
    }
    
    if (ePastCutsInEvent && pimsPastCutsInEvent) {
        IsSelected_eepi = true;
        outTree_e_piminus -> Fill();
        if (fdebug>3) std::cout << "Filling (e,e'pi-) TTree with this event!" << std::endl;
        Nevents_passed_e_pims_cuts ++ ;
        if (eepimsPastCutsInEvent) Nevents_passed_e_pims_kinematics_cuts ++;
        
        for (int pimsIdx=0; pimsIdx<Npims; pimsIdx++) {
            Stream_e_pi_line_to_CSV( "pi-", pimsIdx,
                                    pimsPastSelectionCuts[pimsIdx], eepimsPastKinematicalCuts[pimsIdx],
                                    fdebug );
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
void GetParticlesByType (NPiEvent npi, int evnum, int fdebug){
    // get particles by type
    npi.setNe      (electrons .size());
    npi.setNn      (neutrons  .size());
    npi.setNp      (protons   .size());
    npi.setNpips   (pipluses  .size());
    npi.setNpims   (piminuses .size());
    npi.setNgammas (gammas    .size());
    npi.setNd      (deuterons.size());
    if (fdebug>2){
        std::cout
        << "particles in event "            << evnum                << " : "
        << "N(electrons): "                 << electrons.size()     <<  ","
        << "N(protons): "                   << neutrons.size()      <<  ","
        << "N(neutrons): "                  << protons.size()       <<  ","
        << "N(pi+): "                       << pipluses.size()      <<  ","
        << "N(pi-): "                       << piminuses.size()     <<  ","
        << "N(gammas): "                    << gammas.size()        <<  ","
        << "N(deuterons): "                 << deuterons.size()     <<  ","
        << std::endl;
    }
}


    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    void ExtractElectronInformation(NPiEvent npi, int fdebug){
        // ------------------------------------------------------------------------------------------------
        // extract electron information
        // ------------------------------------------------------------------------------------------------
        // find leading electron as the one with highest energy
        double  leading_e_E;
        int     leading_e_index = 0;
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
        npi.setE(electrons[leading_e_index]);
        // set leading electron vertex
        np.setVe(GetParticleVertex( electrons[leading_e_index] ));
        
        // detector information on electron
        auto e_PCAL_info= electrons[leading_e_index]->cal(PCAL);
        npi.setE_E_PCAL        (e_PCAL_info->getEnergy());
        npi.setE_PCAL_sector   (e_PCAL_info->getSector());
        npi.setE_PCAL_V        (e_PCAL_info->getLv());
        npi.setE_PCAL_W        (e_PCAL_info->getLw());
        npi.setE_E_ECIN        (electrons[leading_e_index]->cal(ECIN)->getEnergy());
        npi.setE_E_ECOUT       (electrons[leading_e_index]->cal(ECOUT)->getEnergy());
        
        // hit position in PCAL
        npi.setE_PCAL_x        (e_PCAL_info->getX());
        npi.setE_PCAL_y        (e_PCAL_info->getY());
        npi.setE_PCAL_z        (e_PCAL_info->getZ());
        
        // Drift Chamber tracking system
        auto e_DC_info  = electrons[leading_e_index]->trk(DC);
        npi.setE_DC_sector     (e_DC_info->getSector()); // tracking sector
        npi.setE_DC_Chi2N      (e_DC_info->getChi2N());  // tracking chi^2/NDF
        
        for (int regionIdx=0; regionIdx<3; regionIdx++) {
            int DC_layer = DC_layers[regionIdx];
            npi.setE_DC_x   (electrons[leading_e_index]->traj(DC,DC_layer)->getX(), regionIdx);
            npi.setE_DC_y   (electrons[leading_e_index]->traj(DC,DC_layer)->getY(), regionIdx);
            npi.setE_DC_z   (electrons[leading_e_index]->traj(DC,DC_layer)->getZ(), regionIdx);
        }
        if (fdebug > 2) std::cout << "extracted electron information and computed kinematics" << std::endl;
        
        // ------------------------------------------------------------------------------------------------
        // now, check if electron passed event selection requirements
        // ------------------------------------------------------------------------------------------------
        npi.setEPastCutsInEvent(cuts.CheckIfElectronPassedSelectionCuts(npi.getE_PCAL_x(), npi.getE_PCAL_y(),
                                                                        npi.getE_PCAL_W(), npi.getE_PCAL_V(),
                                                                        npi.getE_E_PCAL(), npi.getE_E_ECIN(),
                                                                        npi.getE_E_ECOUT(),
                                                                        npi.getE(), npi.getVe(),
                                                                        npi.getE_PCAL_sector(), // e_PCAL_sector should be consistent with e_DC_sector
                                                                        npi.getE_DC_x(), npi.getE_DC_y(), npi.getE_DC_z(),
                                                                        torusBending ));
        if (npi.getEPastCutsInEvent)  Nevents_passed_e_cuts++ ;
    }

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
    void ExtractPipsInformation( int pipsIdx, int fdebug ){
        if (fdebug>2)
            std::cout << "ExtractPipsInformation( pipsIdx=" << pipsIdx << ", fdebug=" << fdebug << " )" << std::endl;
        
        
        // extract positive pion information
        npi.setPiplus               (pipluses[pipsIdx], pipsIdx);
        npi.setZpips                (piplus[pipsIdx].E() / omega, pipsIdx);
        npi.setVpiplus              (GetParticleVertex( pipluses[pipsIdx] ), pipsIdx);
        npi.setPips_chi2PID         (pipluses[pipsIdx]->par()->getChi2Pid(), pipsIdx);
        
        // EC in and out
        npi.setPips_E_ECIN          (pipluses[pipsIdx]->cal(ECIN)->getEnergy(), pipsIdx);
        npi.setPips_E_ECOUT         (pipluses[pipsIdx]->cal(ECOUT)->getEnergy(), pipsIdx);
        // PCAL
        auto pips_PCAL_info         = pipluses[pipsIdx]->cal(PCAL);
        npi.setPips_E_PCAL          (pips_PCAL_info->getEnergy(), pipsIdx);
        npi.setPips_PCAL_sector     (pips_PCAL_info->getSector(), pipsIdx);
        npi.setPips_PCAL_V          (pips_PCAL_info->getLv(), pipsIdx);
        npi.setPips_PCAL_W          (pips_PCAL_info->getLw(), pipsIdx);
        npi.setPips_PCAL_x          (pips_PCAL_info->getX(), pipsIdx);
        npi.setPips_PCAL_y          (pips_PCAL_info->getY(), pipsIdx);
        npi.setPips_PCAL_z          (pips_PCAL_info->getZ(), pipsIdx);
        // DC
        auto pips_DC_info           = pipluses[pipsIdx]->trk(DC);
        npi.setPips_DC_sector       (pips_DC_info->getSector(), pipsIdx); // tracking sector
        npi.setPips_Chi2N           (pips_DC_info->getChi2N(), pipsIdx);  // tracking chi^2/NDF
        for (int regionIdx=0; regionIdx<3; regionIdx++) {
            DC_layer = DC_layers[regionIdx];
            npi.setPips_DC_x        (pipluses[pipsIdx]->traj(DC,DC_layer)->getX(), pipsIdx, regionIdx);
            npi.setPips_DC_y        (pipluses[pipsIdx]->traj(DC,DC_layer)->getY(), pipsIdx, regionIdx);
            npi.setPips_DC_z        (pipluses[pipsIdx]->traj(DC,DC_layer)->getZ(), pipsIdx, regionIdx);
            if (fdebug>3) {
                std::cout
                << "pips_DC_sector[pipsIdx="<<pipsIdx<<"]="
                << npi.getPips_DC_sector(pipsIdx)
                << ", DC_layer = " << DC_layer
                << std::endl
                << "pips_DC_x[pipsIdx="<<pipsIdx<<"][regionIdx="<<regionIdx<<"]="
                << npi.getPips_DC_x(pipsIdx, regionIdx)
                << std::endl
                << "pips_DC_y[pipsIdx="<<pipsIdx<<"][regionIdx="<<regionIdx<<"]="
                << npi.getPips_DC_y(pipsIdx, regionIdx)
                << std::endl
                << "pips_DC_z[pipsIdx="<<pipsIdx<<"][regionIdx="<<regionIdx<<"]="
                << npi.getPips_DC_z(pipsIdx, regionIdx)
                << std::endl;
            }
        }
        // ------------------------------------------------------------------------------------------------
        // now, check if pion passed event selection requirements
        // ------------------------------------------------------------------------------------------------
        npi.setPipsPastSelectionCuts (cuts.CheckIfPionPassedSelectionCuts("pi+",
                                                                        npi.getPips_DC_sector(pipsIdx),
                                                                        npi.getPips_DC_x(pipsIdx),
                                                                        npi.getPips_DC_y(pipsIdx),
                                                                        npi.getPips_DC_z(pipsIdx),
                                                                        npi.getPips_chi2PID(pipsIdx),  npi.getPiplus(pipsIdx).P(),
                                                                        npi.getVe(),
                                                                        npi.getVpiplus(pipsIdx),
                                                                        fdebug), pipsIdx);
        if (npi.getPipsPastSelectionCuts[pipsIdx]) {
            setEepipsPastKinematicalCuts(eepiPassedKinematicalCriteria(piplus[pipsIdx], fdebug), pipsIdx);
        }
        
        if (npi.getPipsPastSelectionCuts(pipsIdx)) {
            setPipsPastCutsInEvent(true);
            Nevents_passed_pips_cuts ++;
            if (npi.getEepipsPastKinematicalCuts(pipsIdx)) {
                npi.setEepipsPastCutsInEvent(true);
            }
        }
        
        npi.setPiplus_Px    (piplus[pipsIdx].Px(), pipsIdx);
        npi.setPiplus_Py    (piplus[pipsIdx].Py(), pipsIdx);
        npi.setPiplus_Pz    (piplus[pipsIdx].Pz(), pipsIdx);
        npi.setPiplus_E     (piplus[pipsIdx].E(), pipsIdx);
        npi.setVpiplus_X    (Vpiplus[pipsIdx].X(), pipsIdx);
        npi.setVpiplus_Y    (Vpiplus[pipsIdx].Y(), pipsIdx);
        npi.setVpiplus_Z    (Vpiplus[pipsIdx].Z(), pipsIdx);
    }


    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    void ExtractPimsInformation( int pimsIdx, int fdebug ){
        if (fdebug>2)
            std::cout << "ExtractPimsInformation( pimsIdx=" << pimsIdx << ", fdebug=" << fdebug << " )" << std::endl;
        
        
        // extract negative pion information
        npi.setPiminus               (piminuses[pimsIdx], pimsIdx);
        npi.setZpims                (piminus[pimsIdx].E() / omega, pimsIdx);
        npi.setVpiminus              (GetParticleVertex( piminuses[pimsIdx] ), pimsIdx);
        npi.setPims_chi2PID         (piminuses[pimsIdx]->par()->getChi2Pid(), pimsIdx);
        
        // EC in and out
        npi.setPims_E_ECIN          (piminuses[pimsIdx]->cal(ECIN)->getEnergy(), pimsIdx);
        npi.setPims_E_ECOUT         (piminuses[pimsIdx]->cal(ECOUT)->getEnergy(), pimsIdx);
        // PCAL
        auto pims_PCAL_info         = piminuses[pimsIdx]->cal(PCAL);
        npi.setPims_E_PCAL          (pims_PCAL_info->getEnergy(), pimsIdx);
        npi.setPims_PCAL_sector     (pims_PCAL_info->getSector(), pimsIdx);
        npi.setPims_PCAL_V          (pims_PCAL_info->getLv(), pimsIdx);
        npi.setPims_PCAL_W          (pims_PCAL_info->getLw(), pimsIdx);
        npi.setPims_PCAL_x          (pims_PCAL_info->getX(), pimsIdx);
        npi.setPims_PCAL_y          (pims_PCAL_info->getY(), pimsIdx);
        npi.setPims_PCAL_z          (pims_PCAL_info->getZ(), pimsIdx);
        // DC
        auto pims_DC_info           = piminuses[pimsIdx]->trk(DC);
        npi.setPims_DC_sector       (pims_DC_info->getSector(), pimsIdx); // tracking sector
        npi.setPims_Chi2N           (pims_DC_info->getChi2N(), pimsIdx);  // tracking chi^2/NDF
        for (int regionIdx=0; regionIdx<3; regionIdx++) {
            DC_layer = DC_layers[regionIdx];
            npi.setPims_DC_x        (piminuses[pimsIdx]->traj(DC,DC_layer)->getX(), pimsIdx, regionIdx);
            npi.setPims_DC_y        (piminuses[pimsIdx]->traj(DC,DC_layer)->getY(), pimsIdx, regionIdx);
            npi.setPims_DC_z        (piminuses[pimsIdx]->traj(DC,DC_layer)->getZ(), pimsIdx, regionIdx);
            if (fdebug>3) {
                std::cout
                << "pims_DC_sector[pimsIdx="<<pimsIdx<<"]="
                << npi.getPims_DC_sector(pimsIdx)
                << ", DC_layer = " << DC_layer
                << std::endl
                << "pims_DC_x[pimsIdx="<<pimsIdx<<"][regionIdx="<<regionIdx<<"]="
                << npi.getPims_DC_x(pimsIdx, regionIdx)
                << std::endl
                << "pims_DC_y[pimsIdx="<<pimsIdx<<"][regionIdx="<<regionIdx<<"]="
                << npi.getPims_DC_y(pimsIdx, regionIdx)
                << std::endl
                << "pims_DC_z[pimsIdx="<<pimsIdx<<"][regionIdx="<<regionIdx<<"]="
                << npi.getPims_DC_z(pimsIdx, regionIdx)
                << std::endl;
            }
        }
        // ------------------------------------------------------------------------------------------------
        // now, check if pion passed event selection requirements
        // ------------------------------------------------------------------------------------------------
        npi.setPimsPastSelectionCuts (cuts.CheckIfPionPassedSelectionCuts("pi-",
                                                                        npi.getPims_DC_sector(pimsIdx),
                                                                        npi.getPims_DC_x(pimsIdx),
                                                                        npi.getPims_DC_y(pimsIdx),
                                                                        npi.getPims_DC_z(pimsIdx),
                                                                        npi.getPims_chi2PID(pimsIdx),  npi.getPiminus(pimsIdx).P(),
                                                                        npi.getVe(),
                                                                        npi.getVpiminus(pimsIdx),
                                                                        fdebug), pimsIdx);
        if (npi.getPimsPastSelectionCuts[pimsIdx]) {
            setEepimsPastKinematicalCuts(eepiPassedKinematicalCriteria(piminus[pimsIdx], fdebug), pimsIdx);
        }
        
        if (npi.getPimsPastSelectionCuts(pimsIdx)) {
            setPimsPastCutsInEvent(true);
            Nevents_passed_pims_cuts ++;
            if (npi.getEepimsPastKinematicalCuts(pimsIdx)) {
                npi.setEepimsPastCutsInEvent(true);
            }
        }
        
        npi.setPiminus_Px    (piminus[pimsIdx].Px(), pimsIdx);
        npi.setPiminus_Py    (piminus[pimsIdx].Py(), pimsIdx);
        npi.setPiminus_Pz    (piminus[pimsIdx].Pz(), pimsIdx);
        npi.setPiminus_E     (piminus[pimsIdx].E(), pimsIdx);
        npi.setVpiminus_X    (Vpiminus[pimsIdx].X(), pimsIdx);
        npi.setVpiminus_Y    (Vpiminus[pimsIdx].Y(), pimsIdx);
        npi.setVpiminus_Z    (Vpiminus[pimsIdx].Z(), pimsIdx);
    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Stream_e_pi_line_to_CSV( TString pionCharge, int piIdx,
                             bool passed_cuts_e_pi,
                             bool passed_cuts_e_pi_kinematics,
                             int fdebug ){
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
    
    // ------------------------------------------------------------------------------------------------
    // compute kinematics that also relies on pion information
    // ------------------------------------------------------------------------------------------------
    xF  = 2. * (pi.Dot(q)) / (q.Mag() * W);
    M_X = ( Beam + target - e - pi ).Mag();
    // now stream data to CSV file
    std::vector<double> variables =
    {   (double)status, (double)runnum,     (double)evnum,      (double)beam_helicity,
        e.P(),          e.Theta(),          e.Phi(),            Ve.Z(),
        pi.P(),         pi.Theta(),         pi.Phi(),           Vpi.Z(),
        Q2,             W,                  xB,                 Zpi,
        omega,          xF,                 y,                  M_X,
        (double)Npips, (double)Npims,       (double)Ne,         (double)Ngammas,
        (double)Np,    (double)Nn,          (double)Nd,
    };
    StreamToCSVfile( pionCharge, variables ,
                    passed_cuts_e_pi, passed_cuts_e_pi_kinematics,
                    fdebug );
}



