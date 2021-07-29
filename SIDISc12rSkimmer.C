// last edit July-29, 2021 (EOC, mbp), see README
// ToDo: Add pion DC fiducial cuts


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
void                      StreamToCSVfile (std::vector<Double_t> observables, bool IsSelectedEvent, int fdebug);
void                      ChangeAxesFrame (TString FrameName="q(z) frame");
void                        MoveTo_qFrame ();
void                       printCutValues ();
bool EventPassedElectronSelectionCriteria (Double_t e_PCAL_x, Double_t e_PCAL_y,
                                           Double_t e_PCAL_W,Double_t e_PCAL_V,
                                           Double_t E_PCAL_e,
                                           Double_t E_ECIN_e, Double_t E_ECOUT_e,
                                           TLorentzVector e,
                                           TVector3 Ve,
                                           Double_t e_DC_sector,
                                           Double_t e_DC_x[3],
                                           Double_t e_DC_y[3],
                                           int torusBending);
bool   EventPassedPiplusPastSelectionCutsCriteria (Double_t pips_DC_x[3], Double_t pips_DC_y[3],
                                           Double_t chi2PID, Double_t p,
                                           TVector3 Ve,
                                           TVector3 Vpiplus );
Double_t          Chi2PID_pips_lowerBound (Double_t p, Double_t C=0.88);
Double_t          Chi2PID_pips_upperBound (Double_t p, Double_t C=0.88);
double                       FindCutValue ( std::string cutName );


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
// globals
// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
// cut values
std::vector<std::pair<std::string, double>> cutValues;

// Output root file and tree
TFile * outFile;
TTree * outTree;
// Output CSV file
std::ofstream   CSVfile, SelectedEventsCSVfile;
// vectors in lab-frame
TLorentzVector        e;
TLorentzVector   piplus;
TLorentzVector     Beam;
TLorentzVector        q;
// reconstructed vertex position
TVector3             Ve;
TVector3        Vpiplus;

// kinematics
Double_t   Ebeam = 10.2; // GeV ( for Fall-2019 the enrgy was 10.4096 )
Double_t             xB;
Double_t             Q2;
Double_t          omega;
Double_t              z; // energy fraction rest frame

// meta-data
int        torusBending = -1; // -1 for In-bending, +1 for Out-bending
int        DC_layers[3] = {6,18,36};// Region 1 is denoted at DC detector 6, Region 2 is denoted 18, Region 3 - as 36

// vectors in q-frame
TLorentzVector       piplus_qFrame;
Double_t        Ppips_t_q, Ppips_q;

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
        std::cout << "processing " << RunNumberStr << std::endl;
    }
    
    // read cut values csv file
    csv_reader csvr;
    cutValues = csvr.read_csv("cutValues.csv");
    if (fdebug>0) { printCutValues(); }
    

    
    // defenitions
    auto db=TDatabasePDG::Instance();
    e       = TLorentzVector(0,0,0,db->GetParticle( 11   )->Mass());
    piplus  = TLorentzVector(0,0,0,db->GetParticle( 211  )->Mass());
    Ve      = TVector3();
    Vpiplus = TVector3();
    Beam.SetPxPyPzE(0,0,Ebeam,Ebeam);
    
    bool      ePastSelectionCuts = false;
    bool piplusPastSelectionCuts = false;
    bool         EventPassedCuts = false;
    
    int           event = 0;
    int      good_event = 0;
    int         Nevents = 0;
    int Fastest_pipsIdx = 0;
    
    double        Mp = db->GetParticle( 2212 )->Mass();
    double  E_PCAL_e = 0; // electron energy deposit in PCAL [GeV]
    double  E_ECIN_e = 0; // electron energy deposit in ECAL_in [GeV]
    double  E_ECOUT_e= 0; // electron energy deposit in ECAL_out [GeV]
    double  chi2PID_pips;
    double  e_PCAL_W,    e_PCAL_V;
    double  pips_PCAL_W, pips_PCAL_V;
    double  pips_PCAL_x, pips_PCAL_y,       pips_PCAL_z;
    double  e_PCAL_x,    e_PCAL_y,          e_PCAL_z;
    double  e_PCAL_sector;
    double  e_DC_sector,        e_DC_Chi2N;
    double  e_DC_x[3],   e_DC_y[3];
    double  pips_PCAL_sector;
    double  pips_DC_sector,     pips_Chi2N;
    double  pips_DC_x[3],pips_DC_y[3];
    
    double  E_PCAL_pips, E_ECIN_pips,       E_ECOUT_pips;
    int     DC_layer;
    int     runnum,      evnum;
    
    // open CSV file
    OpenOutputFiles ("/volatile/clas12/users/ecohen/BAND/SIDIS_skimming/skimmed_SIDIS_inc_"
                     + (TString)RunNumberStr ,
                     ((TString)("event,")
                      +(TString)("Ee,Pe,Pe_x,Pe_y,Pe_z,Ptheta_e,Pphi_e,")
                      +(TString)("Ngammas,Np,Nn,Npips,Npims,")
                      +(TString)("omega,q,q_x,q_y,q_z,xB,Q2,z,")
                      +(TString)("Epips,Ppips,Ppips_x,Ppips_y,Ppips_z,")
                      +(TString)("Ppips_t_q,Ppips_q,")
                      +(TString)("E_PCAL_e,E_ECIN_e,E_ECOUT_e,")
                      +(TString)("Vx_e,Vy_e,Vz_e,Vx_pips,Vy_pips,Vz_pips,")
                      +(TString)("chi2PID_pips,")
                      +(TString)("e_PCAL_W,e_PCAL_V,pips_PCAL_W,pips_PCAL_V,")
                      +(TString)("e_PCAL_x,e_PCAL_y,e_PCAL_z,e_PCAL_sector,")
                      +(TString)("e_DC_sector,e_DC_Chi2N,")
                      +(TString)("e_DC_x[region-1],e_DC_y[region-1],")
                      +(TString)("e_DC_x[region-2],e_DC_y[region-2],")
                      +(TString)("e_DC_x[region-3],e_DC_y[region-3],")
                      ));
    
    // output tree branches
    outTree->Branch("eventnumber"       ,&evnum                 );
    outTree->Branch("runnum"            ,&runnum                );
    outTree->Branch("E_PCAL_e"          ,&E_PCAL_e              );
    outTree->Branch("E_ECIN_e"          ,&E_ECIN_e              );
    outTree->Branch("E_ECOUT_e"         ,&E_ECOUT_e             );
    outTree->Branch("chi2PID_pips"      ,&chi2PID_pips          );
    
    outTree->Branch("e_PCAL_W"          ,&e_PCAL_W              );
    outTree->Branch("e_PCAL_V"          ,&e_PCAL_V              );
    outTree->Branch("pips_PCAL_x"       ,&pips_PCAL_x           );
    outTree->Branch("pips_PCAL_y"       ,&pips_PCAL_y           );
    outTree->Branch("pips_PCAL_z"       ,&pips_PCAL_z           );
    outTree->Branch("e_PCAL_x"          ,&e_PCAL_x              );
    outTree->Branch("e_PCAL_y"          ,&e_PCAL_y              );
    outTree->Branch("e_PCAL_z"          ,&e_PCAL_z              );
    
    outTree->Branch("e_PCAL_sector"     ,&e_PCAL_sector         );
    
    outTree->Branch("e_DC_sector"       ,&e_DC_sector           );
    outTree->Branch("e_DC_Chi2N"        ,&e_DC_Chi2N            );
    outTree->Branch("e_DC_x"            ,&e_DC_x                );
    outTree->Branch("e_DC_y"            ,&e_DC_y                );
    
    
    outTree->Branch("pips_PCAL_sector"  ,&pips_PCAL_sector      );
    outTree->Branch("pips_DC_sector"    ,&pips_DC_sector        );
    outTree->Branch("pips_Chi2N"        ,&pips_Chi2N            );
    outTree->Branch("pips_DC_x"         ,&pips_DC_x             );
    outTree->Branch("pips_DC_y"         ,&pips_DC_y             );
    
    outTree->Branch("E_PCAL_pips"       ,&E_PCAL_pips           );
    outTree->Branch("E_ECIN_pips"       ,&E_ECIN_pips           );
    
    outTree->Branch("E_ECIN_pips"       ,&E_ECIN_pips           );
    outTree->Branch("E_ECOUT_pips"      ,&E_ECOUT_pips          );
    outTree->Branch("DC_layer"          ,&DC_layer              );
    
    
    outTree->Branch("e"                 ,&e                     );
    outTree->Branch("piplus"            ,&piplus                );
    outTree->Branch("Ve"                ,&Ve                    );
    outTree->Branch("Vpiplus"           ,&Vpiplus               );
    outTree->Branch("Beam"              ,&Beam                  );
    outTree->Branch("q"                 ,&q                     );
    
    outTree->Branch("Ebeam"             ,&Ebeam                 );
    outTree->Branch("xB"                ,&xB                    );
    outTree->Branch("Q2"                ,&Q2                    );
    outTree->Branch("omega"             ,&omega                 );
    outTree->Branch("z"                 ,&z                     );
    outTree->Branch("EventPassedCuts"   ,&EventPassedCuts       );
    outTree->Branch("ePastSelectionCuts"        ,&ePastSelectionCuts            );
    outTree->Branch("piplusPastSelectionCuts"   ,&piplusPastSelectionCuts       );
    
    
    // Record start time
    auto start = std::chrono::high_resolution_clock::now();
    ///////////////////////////////////// input hipo file
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
    /////////////////////////////////////
    
    TChain fake("hipo");
    fake.Add(inputFile.Data());
    //get the hipo data
    auto files = fake.GetListOfFiles();
    gBenchmark->Start("timer");
        
    
    // step over events and extract information....
    for(Int_t i=0;i<files->GetEntries();i++){
        
        //create the event reader
        clas12reader c12(files->At(i)->GetTitle(),{0});
        
        // first count number of events - I don't know how to read it from the HIPO file
        if (fdebug) std::cout << "reading file " << i << std::endl;
        
        Int_t NeventsMaxToProcess = NeventsMax;
        if (NeventsMax<0){
            NeventsMaxToProcess = c12.getReader().getEntries();
        }
            
        // now process the events from the first one...
        while((c12.next()==true) && (event < NeventsMaxToProcess)){
            
            runnum = c12.runconfig()->getRun();
            evnum  = c12.runconfig()->getEvent();

            
            if (fdebug>2) std::cout << "begin analysis of event " << event << std::endl;
            Nevents++;
            
            // initialize
            xB          = Q2        = omega     = -9999;
            E_ECIN_e    = E_ECOUT_e = E_PCAL_e  = -9999;
            chi2PID_pips                        = -9999;
            e_PCAL_W    = e_PCAL_V              = -9999;
            e_PCAL_x    = e_PCAL_y  = e_PCAL_z  = -9999;
            e_PCAL_sector                       = -9999;
            e_DC_sector = e_DC_Chi2N            = -9999;
            pips_DC_sector                      = -9999;
            pips_PCAL_sector                    = -9999;
            pips_PCAL_W = pips_PCAL_V           = -9999;
            pips_PCAL_x = pips_PCAL_y           = -9999;
            pips_PCAL_z                         = -9999;
            E_PCAL_pips                         = -9999;
            E_ECIN_pips = E_ECOUT_pips          = -9999;
            DC_layer                            = -9999;
            Ve          = Vpiplus               = TVector3();
            for (int regionIdx=0; regionIdx<3; regionIdx++) {
                e_DC_x[regionIdx]   = e_DC_y[regionIdx]     = -9999;
                pips_DC_x[regionIdx]= pips_DC_y[regionIdx]  = -9999;
            }
            ePastSelectionCuts  = piplusPastSelectionCuts       = false;
            
            
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
            
            
            // filter D(e,e’pi+n)X reaction
            // and compute event kinematics
            if(  Ne == 1 && Npips >= NpipsMin ){
                
                // set electron 4-momentum
                SetLorentzVector(e,electrons[0]);
                // set electron vertex
                Ve      = GetParticleVertex( electrons[0] );
                // compute event kinematics
                q       = Beam - e;
                Q2      = -q.Mag2();
                omega   = q.E();
                xB      = Q2/(2. * Mp * q.E());
                
                // detector information on electron
                auto e_PCAL_info= electrons[0]->cal(PCAL);
                E_PCAL_e        = e_PCAL_info->getEnergy();
                e_PCAL_sector   = e_PCAL_info->getSector();
                e_PCAL_V        = e_PCAL_info->getLv();
                e_PCAL_W        = e_PCAL_info->getLw();
                E_ECIN_e        = electrons[0]->cal(ECIN)->getEnergy();
                E_ECOUT_e       = electrons[0]->cal(ECOUT)->getEnergy();
                
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
                
                // select the fastest pion as the "leader"
                Fastest_pipsIdx = 0;
                SetLorentzVector(piplus  ,pips[0]);
                TLorentzVector piplus_tmp(0,0,0,db->GetParticle(211)->Mass());
                for (int pipsIdx=0; pipsIdx < Npips; pipsIdx++) {
                    SetLorentzVector(piplus_tmp  ,pips[pipsIdx]);
                    if (piplus_tmp.E() > piplus.E()) {
                        Fastest_pipsIdx = pipsIdx;
                    }
                }
                SetLorentzVector(piplus , pips[Fastest_pipsIdx]);
                // set reconstructed vertex position
                Vpiplus = GetParticleVertex( pips[Fastest_pipsIdx] );
                // Chi2PID for the pion
                chi2PID_pips = pips[Fastest_pipsIdx]->par()->getChi2Pid();
                
                // detector information on electron
                auto pips_PCAL_info= pips[Fastest_pipsIdx]->cal(PCAL);
                E_PCAL_pips        = pips_PCAL_info->getEnergy();
                pips_PCAL_sector   = pips_PCAL_info->getSector();
                pips_PCAL_V        = pips_PCAL_info->getLv();
                pips_PCAL_W        = pips_PCAL_info->getLw();
                E_ECIN_pips        = pips[Fastest_pipsIdx]->cal(ECIN)->getEnergy();
                E_ECOUT_pips       = pips[Fastest_pipsIdx]->cal(ECOUT)->getEnergy();
                
                // hit position in PCAL
                pips_PCAL_x        = pips_PCAL_info->getX();
                pips_PCAL_y        = pips_PCAL_info->getY();
                pips_PCAL_z        = pips_PCAL_info->getZ();
                
                // Drift Chamber tracking system
                auto pips_DC_info  = pips[Fastest_pipsIdx]->trk(DC);
                pips_DC_sector     = pips_DC_info->getSector(); // tracking sector
                pips_Chi2N         = pips_DC_info->getChi2N();  // tracking chi^2/NDF
                
                for (int regionIdx=0; regionIdx<3; regionIdx++) {
                    DC_layer = DC_layers[regionIdx];
                    pips_DC_x[regionIdx] = pips[Fastest_pipsIdx]->traj(DC,DC_layer)->getX();
                    pips_DC_y[regionIdx] = pips[Fastest_pipsIdx]->traj(DC,DC_layer)->getY();
                    
                }
                // hadron rest frame energy fraction
                z = piplus.E()/q.E();
                if (fdebug > 2){
                    Printf("selected the fastest pion as the leader and extracted information");
                }
                
                // temporarily fill pips 4-vector in q-frame
                piplus_qFrame = piplus;
                // move to the axes-frame where q is parallel to the z-direction,
                // and compute pion P(perp) and P(parallel) in this frame
                ChangeAxesFrame();
                Ppips_t_q   = piplus_qFrame.Pt();
                Ppips_q     = piplus_qFrame.Pz();
                if (fdebug > 2) Printf("Moved to q-Frame and computed pi+ momentum components");
                
                
                
                // decide if to write this event to "selected events csv-file"
                bool IsSelectedEvent = false;
                if ( doApplySelectionCuts ) {
                    ePastSelectionCuts = EventPassedElectronSelectionCriteria(e_PCAL_x, e_PCAL_y, e_PCAL_W, e_PCAL_V,
                                                                      E_PCAL_e,  E_ECIN_e, E_ECOUT_e,
                                                                      e, Ve,
                                                                      e_PCAL_sector, // e_PCAL_sector should be consistent with e_DC_sector
                                                                      e_DC_x, e_DC_y,
                                                                      torusBending );
                    piplusPastSelectionCuts = EventPassedPiplusPastSelectionCutsCriteria(pips_DC_x,     pips_DC_y,
                                                                         chi2PID_pips,  piplus.P(),
                                                                         Ve, Vpiplus );
                    EventPassedCuts =( ePastSelectionCuts && piplusPastSelectionCuts );
                    if ( EventPassedCuts ) {
                        IsSelectedEvent = true;
                        outTree -> Fill();
                    }
                }
                
                
                
                StreamToCSVfile({(Double_t)event,
                    e.E(),              e.P(),          e.Px(),             e.Py(),
                    e.Pz(),             e.Theta(),      e.Phi(),
                    (Double_t)Ngammas,  (Double_t)Np,
                    (Double_t)Nn,       (Double_t)Npips,    (Double_t)Npims,
                    q.E(),              q.P(),          q.Px(),             q.Py(),             q.Pz(),         xB,             Q2,         z,
                    piplus.E(),         piplus.P(),     piplus.Px(),        piplus.Py(),        piplus.Pz(),
                    Ppips_t_q,          Ppips_q,
                    E_PCAL_e,           E_ECIN_e,       E_ECOUT_e,
                    Ve.X(),             Ve.Y(),         Ve.Z(),
                    Vpiplus.X(),        Vpiplus.Y(),    Vpiplus.Z(),
                    chi2PID_pips,
                    e_PCAL_W,           e_PCAL_V,       pips_PCAL_W,        pips_PCAL_V,
                    e_PCAL_x,           e_PCAL_y,       e_PCAL_z,           e_PCAL_sector,
                    e_DC_sector,        e_DC_Chi2N,
                    e_DC_x[0],          e_DC_y[0],
                    e_DC_x[1],          e_DC_y[1],
                    e_DC_x[2],          e_DC_y[2],
                }, IsSelectedEvent, fdebug);
                
                good_event ++ ;
                if (fdebug>2) {
                    std::cout
                    << "streamed to CSV file,"
                    << "Ee: "       << e.E()          << " GeV, "
                    << "E(pi+): "   << piplus  .E()   << " GeV, "
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
    std::cout << "Done. Elapsed time: " << elapsed.count()<< ", processed " << event << " events, " << good_event << " passed filter\n";
    CloseOutputFiles("/volatile/clas12/users/ecohen/BAND/SIDIS_skimming/");
}




// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
bool EventPassedElectronSelectionCriteria(Double_t e_PCAL_x, Double_t e_PCAL_y,
                                          Double_t e_PCAL_W,Double_t e_PCAL_V,
                                          Double_t E_PCAL_e,
                                          Double_t E_ECIN_e, Double_t E_ECOUT_e,
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
    
    
    // Cut on z-vertex position:
    // based on RGA_Analysis_Overview_and_Procedures_Nov_4_2020-6245173-2020-12-09-v3.pdf
    // p. 71
    Double_t Vz_min,Vz_max;
    if (torusBending==-1){ // in-bending torus field
        // Spring 19 and Spring 2020 in-bending.
        // -13.0 cm < Vz < +12.0 cm
        Vz_min = FindCutValue("Vz_e_min_inbending");
        Vz_max = FindCutValue("Vz_e_max_inbending");

        std::cout << "Vz min,max: " << Vz_min << ","<< Vz_min << std::endl;
//        Vz_max = 12.0;
        
    } else if (torusBending==1){ // Out-bending torus field
        // Fall 2019 (without low-energy-run) was out-bending.
        // -18.0 cm < Vz < +10.0 cm
        
        Vz_min = -18.0;
        Vz_max = 10.0;
        
    } else {
        std::cout
        << "Un-identified torus bending "
        << torusBending
        << ", return" << std::endl;
        return false;
    }
    
    if(
       // fiducial cuts on PCAL
       fabs(e_PCAL_x)>0
       &&  fabs(e_PCAL_y)>0
       &&  e_PCAL_W > 19
       &&  e_PCAL_V > 19
       
       // Electron Identification Refinement  - PCAL Minimum Energy Deposition Cut
       &&  E_PCAL_e > 0.07
       
       // Sampling fraction cut
       && ((E_PCAL_e + E_ECIN_e + E_ECOUT_e)/e.P()) > 0.17
       && (E_ECIN_e/e.P() > 0.2 - E_PCAL_e/e.P()) // RGA AN puts "<" here mistakenly
       
       // Cut on z-vertex position: in-bending torus field -13.0 cm < Vz < +12.0 cm
       // Spring 19 and Spring 2020 in-bending.
       // Fall 2019 (without low-energy-run) was out-bending.
       &&  ((Vz_min < Ve.Z()) && (Ve.Z() < Vz_max))
       ) return true;
    
    return false;
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
bool EventPassedPiplusPastSelectionCutsCriteria(Double_t pips_DC_x[3], Double_t pips_DC_y[3],
                                        Double_t chi2PID, Double_t p,
                                        TVector3 Ve,
                                        TVector3 Vpiplus ){
    // decide if pi+ selection cuts
    //
    // input:
    // --------
    // DC_x, DC_y   pi+ drift-chamber coordinates
    // chi2PID      pi+ chi2PID     (chi2PID_pips)
    // p            pi+ momentum    (Ppips.P())
    //
    
    // DC - fiducial cuts on DC
    // ToDo: Complete this! Need help from Dien with pion fiductial
    // cuts implementation in DC_fiducials.cpp
    if(
       (pips_DC_x[0]<-1e9 )
       ){
        return false;
    }
    
    
    if(
       // pi+ Identification Refinement - chi2PID vs. momentum
       ( Chi2PID_pips_lowerBound( p ) < chi2PID && chi2PID < Chi2PID_pips_upperBound( p ) )
       
       // Cut on the z-Vertex Difference Between Electrons and Hadrons.
       &&  ( fabs((Ve-Vpiplus).Z()) < 20.0)
       ) return true;
    
    return false;
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
Double_t Chi2PID_pips_lowerBound( Double_t p, Double_t C){
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
Double_t Chi2PID_pips_upperBound( Double_t p, Double_t C){
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
    outFile = new TFile( outfilename + ".root" ,"RECREATE");
    outTree = new TTree( "sidis" , "Event information");
    
    // Create output csv files
    CSVfile.open( outfilename + ".csv" );
    CSVfile << header << std::endl;
    
    SelectedEventsCSVfile.open( outfilename + "_selected_events.csv" );
    SelectedEventsCSVfile << header << std::endl;
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void CloseOutputFiles (TString OutDataPath){
    // close output CSV
    CSVfile.close();
    SelectedEventsCSVfile.close();
    
    int Nentires = outTree -> GetEntries();
    
    // close output ROOT
    outFile->cd();
    outTree->Write();
    outFile->Close();
    
    std::cout << "output files ready in root/csv formats in " << std::endl
    << OutDataPath << std::endl
    << "wrote " << Nentires << " to output root file"
    << std::endl;
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void StreamToCSVfile (std::vector<Double_t> observables, bool IsSelectedEvent, int fdebug){
    if (fdebug>1) {
        std::cout << "streaming to CSVfile" << std::endl;
    }
    
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
    // move to a reference frame where
    // q is the z axis
    Double_t q_phi   = q.Phi();
    Double_t q_theta = q.Theta();
    
    q.RotateZ(-q_phi);
    q.RotateY(-q_theta);
    
    // rotate additional vectors to the same axes-frame
    // and they reside on the x-z plane: v=(v_x,0,v_q)
    piplus_qFrame.RotateZ(-q_phi);
    piplus_qFrame.RotateY(-q_theta);
    Double_t piplus_phi = piplus_qFrame.Phi();
    piplus_qFrame.RotateZ(-piplus_phi);
    
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
