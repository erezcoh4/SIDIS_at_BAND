#include <vector>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TH1.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include <time.h>

#include "Auxiliary/bank.h"
#include "Auxiliary/BBand.h"
#include "Auxiliary/BEvent.h"
#include "Auxiliary/constants.h"
#include "Auxiliary/bandhit.cpp"
#include "Auxiliary/clashit.cpp"
#include "Auxiliary/genpart.cpp"

#define NMAXEVENTS 5000000
#define NMAXNML 20 // maximal number of XXX
#define NMAXPIONS 5 // maximal allowed number of pions
#define r2d 180./3.1415 // radians to degrees

// Output CSV file
std::ofstream   CSVfile_e_e_pips, CSVfile_e_e_pims;

// time
clock_t tStart = clock();

// globals
TString            csvheader = ((TString)"status,runnum,evnum,"
                                +(TString)"e_P,e_Theta,e_Phi,e_Vz,"
                                +(TString)"pi_P,pi_Theta,pi_Phi,pi_Vz,"
                                +(TString)"Q2,W,xB,Zpi,"
                                +(TString)"omega,xF,y,"
                                +(TString)"M_X_ee_pi,"
                                +(TString)"theta_sq,WPrime,"
                                +(TString)"e_DC_sector,pi_DC_sector,");
TString                   DataPath;
TString                   Filename;
TString       CSVFilename_e_e_pips;
TString       CSVFilename_e_e_pims;


int               Nevents_e_e_pips;
int               Nevents_e_e_pims;
int                         fdebug;
int               NeventsInRGATree;
Long64_t                    runnum;
Long64_t                     evnum;
int                     NMaxEvents;
int                            nml;
int                            nPi;



TFile                    * RGAFile;
TTree                    * RGATree;

Double_t       Me  = 0.00051099895; // GeV/c2
Double_t       Mpims  = 0.13957039; // GeV/c2
Double_t       Mpips  = 0.13957039; // GeV/c2
Double_t            Mp  = 0.938272; // GeV/c2
Double_t             Mn = 0.939565; // GeV/c2
Double_t                Md = 1.875; // GeV/c2
Double_t             Mp2 = Mp * Mp;
Double_t                        xB; // Bjorken x
Double_t                        xF; // Feynman x
Double_t                        Q2;
Double_t                     omega;
Double_t                        w2; // omega^2
Double_t          Zpips[NMAXPIONS]; // energy fraction rest frame
Double_t          Zpims[NMAXPIONS]; // energy fraction rest frame

// Input TTree variables
Float_t                Ep[NMAXNML]; // electron momentum
Float_t               Epx[NMAXNML]; // electron momentum
Float_t               Epy[NMAXNML]; // electron momentum
Float_t               Epz[NMAXNML]; // electron momentum
Float_t               Evz[NMAXNML]; // electron vertex in z-direction
Int_t                Epid[NMAXNML]; // electron PID
Int_t             Esector[NMAXNML]; // electron sector
Float_t            Etheta[NMAXNML];
Float_t              Ephi[NMAXNML];
Int_t          Efiducial1[NMAXNML];
Int_t          Efiducial2[NMAXNML];
Int_t          Efiducial3[NMAXNML];
Float_t      ePCAL_energy[NMAXNML];
Float_t           ePCAL_x[NMAXNML];
Float_t           ePCAL_y[NMAXNML];
Float_t           ePCAL_z[NMAXNML];
Float_t          ePCAL_lu[NMAXNML];
Float_t          ePCAL_lv[NMAXNML];
Float_t          ePCAL_lw[NMAXNML];
Float_t      eECIN_energy[NMAXNML];
Float_t     eECOUT_energy[NMAXNML];

Float_t                Nu[NMAXNML];
Float_t                qx[NMAXNML];
Float_t                qy[NMAXNML];
Float_t                qz[NMAXNML];




// methods
void                   PrintDone ();
void               OpenInputFile ();
void              CloseInputFile ();
void             OpenOutputFiles ();
void            CloseOutputFiles ();
void       InitializeFileReading ();
//void            StreamToCSVfile (std::vector<Double_t> observables, bool passed_cuts_e_pi_kinematics);
void                SetVerbosity (int ffdebug )          {fdebug = ffdebug;} ;
void                 SetDataPath (TString fDataPath )    {DataPath = fDataPath;} ;
void                 SetFileName (TString fFileName )    {Filename = fFileName;} ;
void               SetNMaxEVents (int fNMAXevents )      {NMaxEvents = fNMAXevents;};
void                  ReadEvents ();
void               SetInputTTree ();


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
// Main functionality
// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void ReadIgorRGAFile(TString fFileName="ntupleNew",
                     int fNMAXevents=-1,
                     int ffdebug=1,
                     int PrintProgress=1000,
                     TString fDataPath="/Users/erezcohen/Desktop/data/BAND/RGA_Free_proton/"
                     ){
    
    SetVerbosity          ( ffdebug );
    SetDataPath           ( fDataPath );
    SetFileName           ( fFileName );
    SetNMaxEVents         ( fNMAXevents );
    InitializeFileReading ();
    OpenInputFile         ();
    OpenOutputFiles       ();
    ReadEvents            ();
    CloseOutputFiles      ();
    CloseInputFile        ();
    PrintDone             ();
}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void ReadEvents (){
    if (fdebug>1) std::cout << "ReadEvents (NeventsInRGATree="<<NeventsInRGATree<<")" << std::endl;
    
    for (int event=0; event < NeventsInRGATree ; event++) {
        RGATree->GetEntry(event);
        
        if (fdebug>1){
            std::cout
            << "run "    << runnum      << ","
            << "event "  << evnum       << ","
            << "nml: "   << nml         << ","
            << "nPi: "   << nPi         << ","
            << std::endl;
            
            for (int e_idx=0; e_idx<nml; e_idx++) {
                std::cout
                << "Ep["<<e_idx<<"] "               << Ep[e_idx]            << ","
                << std::endl
                << "Epx["<<e_idx<<"] "              << Epx[e_idx]           << ","
                << "Epy["<<e_idx<<"] "              << Epy[e_idx]           << ","
                << "Epz["<<e_idx<<"] "              << Epz[e_idx]           << ","
                << std::endl
                << "Epid["<<e_idx<<"] "             << Epid[e_idx]          << ","
                << "Esector["<<e_idx<<"] "          << Esector[e_idx]       << ","
                << std::endl
                << "Etheta["<<e_idx<<"] "           << Etheta[e_idx]        << ","
                << "Ephi["<<e_idx<<"] "             << Ephi[e_idx]          << ","
                << std::endl
                << "Efiducial1["<<e_idx<<"] "       << Efiducial1[e_idx]    << ","
                << "Efiducial2["<<e_idx<<"] "       << Efiducial2[e_idx]    << ","
                << "Efiducial3["<<e_idx<<"] "       << Efiducial3[e_idx]    << ","
                << std::endl
                << "ePCAL_energy["<<e_idx<<"] "     << ePCAL_energy[e_idx]  << ","
                << "eECIN_energy["<<e_idx<<"] "     << eECIN_energy[e_idx]  << ","
                << "eECOUT_energy["<<e_idx<<"] "    << eECOUT_energy[e_idx] << ","
                << std::endl
                << "Nu["<<e_idx<<"] "               << Nu[e_idx]            << ","
                << std::endl
                << "qx["<<e_idx<<"] "              << qx[e_idx]           << ","
                << "qy["<<e_idx<<"] "              << qy[e_idx]           << ","
                << "qz["<<e_idx<<"] "              << qz[e_idx]           << ","
                << std::endl;
            }
        }
        

        if (fdebug>1) {
            std::cout
            << "event "
            << event
            << std::endl;
        }
        
        if ( (NMaxEvents>0) && (event>NMaxEvents) ){
            if (fdebug>1) { std::cout << std::endl << "Stop reading events at event " << event << std::endl;}
            break;
        }
    }
    if (fdebug>1) std::cout << "Done ReadEvents ()" << std::endl;
}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void SetInputTTree (){
    NeventsInRGATree = RGATree->GetEntries();
    RGATree  -> SetBranchAddress("runNum"           ,&runnum            );
    RGATree  -> SetBranchAddress("evNum"            ,&evnum             );
    RGATree  -> SetBranchAddress("nml"              ,&nml               );
    
    RGATree  -> SetBranchAddress("Epid"             ,Epid              );
    RGATree  -> SetBranchAddress("Ep"               ,Ep                );
    RGATree  -> SetBranchAddress("Epx"              ,Epx               );
    RGATree  -> SetBranchAddress("Epy"              ,Epy               );
    RGATree  -> SetBranchAddress("Epz"              ,Epz               );
    RGATree  -> SetBranchAddress("Evz"              ,Evz               );
    RGATree  -> SetBranchAddress("Etheta"           ,Etheta            );
    RGATree  -> SetBranchAddress("Ephi"             ,Ephi              );
    RGATree  -> SetBranchAddress("Efiducial1"       ,Efiducial1        );
    RGATree  -> SetBranchAddress("Efiducial2"       ,Efiducial2        );
    RGATree  -> SetBranchAddress("Efiducial3"       ,Efiducial3        );
    RGATree  -> SetBranchAddress("Efiducial3"       ,Efiducial3        );
    RGATree  -> SetBranchAddress("ePCAL_energy"     ,ePCAL_energy      );
    RGATree  -> SetBranchAddress("ePCAL_x"          ,ePCAL_x           );
    RGATree  -> SetBranchAddress("ePCAL_y"          ,ePCAL_y           );
    RGATree  -> SetBranchAddress("ePCAL_z"          ,ePCAL_z           );
    RGATree  -> SetBranchAddress("ePCAL_lu"         ,ePCAL_lu          );
    RGATree  -> SetBranchAddress("ePCAL_lv"         ,ePCAL_lv          );
    RGATree  -> SetBranchAddress("ePCAL_lw"         ,ePCAL_lw          );
    RGATree  -> SetBranchAddress("eECIN_energy"     ,eECIN_energy      );
    RGATree  -> SetBranchAddress("eECOUT_energy"    ,eECOUT_energy     );
    
    RGATree  -> SetBranchAddress("nPi"              ,&nPi              );
    
    RGATree  -> SetBranchAddress("Nu"               ,Nu                );
    RGATree  -> SetBranchAddress("qx"               ,qx                );
    RGATree  -> SetBranchAddress("qy"               ,qy                );
    RGATree  -> SetBranchAddress("qz"               ,qz                );
    
//    ******************************************************************************
//    *Tree    :T         : e'pion                                                 *
//    *Br    1 :runNum    : runNum/L                                               *
//    *Br    2 :evNum     : evNum/L                                                *
//    *Br    4 :Epid      : Epid[nml]/I                                            *
//    *Br    5 :Echarge   : Echarge[nml]/I                                         *
//    *Br    6 :Ep        : Ep[nml]/F                                              *
//    *Br    7 :Epx       : Epx[nml]/F                                             *
//    *Br    8 :Epy       : Epy[nml]/F                                             *
//    *Br    9 :Epz       : Epz[nml]/F                                             *
//    *Br   10 :Etheta    : Etheta[nml]/F                                          *
//    *Br   11 :Ephi      : Ephi[nml]/F                                            *
//    *Br   12 :Etime     : Etime[nml]/F                                           *
//    *Br   13 :Epath     : Epath[nml]/F                                           *
//    *Br   14 :Evz       : Evz[nml]/F                                             *
//    *Br   15 :Estat     : Estat[nml]/I                                           *
//    *Br   16 :EdE       : EdE[nml]/F                                             *
//    *Br   17 :Ebeta     : Ebeta[nml]/F                                           *
//    *Br   18 :Esector   : Esector[nml]/I                                         *
//    *Br   19 :Eregion   : Eregion[nml]/I                                         *
//    *Br   20 :Efiducial1 : Efiducial1[nml]/I                                     *
//    *Br   21 :Efiducial2 : Efiducial2[nml]/I                                     *
//    *Br   22 :Efiducial3 : Efiducial3[nml]/I                                     *
//    *Br   23 :ePCAL_time : ePCAL_time[nml]/F                                     *
//    *Br   24 :ePCAL_energy : ePCAL_energy[nml]/F                                 *
//    *Br   25 :ePCAL_path : ePCAL_path[nml]/F                                     *
//    *Br   26 :ePCAL_sector : ePCAL_sector[nml]/I                                 *
//    *Br   27 :ePCAL_x   : ePCAL_x[nml]/F                                         *
//    *Br   28 :ePCAL_y   : ePCAL_y[nml]/F                                         *
//    *Br   29 :ePCAL_z   : ePCAL_z[nml]/F                                         *
//    *Br   30 :ePCAL_lu  : ePCAL_lu[nml]/F                                        *
//    *Br   31 :ePCAL_lv  : ePCAL_lv[nml]/F                                        *
//    *Br   32 :ePCAL_lw  : ePCAL_lw[nml]/F                                        *
//    *Br   33 :ePCAL_status : ePCAL_status[nml]/I                                 *
//    *Br   34 :eECIN_time : eECIN_time[nml]/F                                     *
//    *Br   35 :eECIN_energy : eECIN_energy[nml]/F                                 *
//    *Br   36 :eECIN_path : eECIN_path[nml]/F                                     *
//    *Br   37 :eECIN_sector : eECIN_sector[nml]/I                                 *
//    *Br   38 :eECIN_x   : eECIN_x[nml]/F                                         *
//    *Br   39 :eECIN_y   : eECIN_y[nml]/F                                         *
//    *Br   40 :eECIN_z   : eECIN_z[nml]/F                                         *
//    *Br   41 :eECIN_lu  : eECIN_lu[nml]/F                                        *
//    *Br   42 :eECIN_lv  : eECIN_lv[nml]/F                                        *
//    *Br   43 :eECIN_lw  : eECIN_lw[nml]/F                                        *
//    *Br   44 :eECIN_status : eECIN_status[nml]/I                                 *
//    *Br   45 :eECOUT_time : eECOUT_time[nml]/F                                   *
//    *Br   46 :eECOUT_energy : eECOUT_energy[nml]/F                               *
//    *Br   47 :eECOUT_path : eECOUT_path[nml]/F                                   *
//    *Br   48 :eECOUT_sector : eECOUT_sector[nml]/I                               *
//    *Br   49 :eECOUT_x  : eECOUT_x[nml]/F                                        *
//    *Br   50 :eECOUT_y  : eECOUT_y[nml]/F                                        *
//    *Br   51 :eECOUT_z  : eECOUT_z[nml]/F                                        *
//    *Br   52 :eECOUT_lu : eECOUT_lu[nml]/F                                       *
//    *Br   53 :eECOUT_lv : eECOUT_lv[nml]/F                                       *
//    *Br   54 :eECOUT_lw : eECOUT_lw[nml]/F                                       *
//    *Br   55 :eECOUT_status : eECOUT_status[nml]/I                               *
//    *Br   56 :EdcX      : EdcX[nml][3]/F                                         *
//    *Br   57 :EdcY      : EdcY[nml][3]/F                                         *
//    *Br   58 :EdcZ      : EdcZ[nml][3]/F                                         *
//    *Br   59 :npiP      : npiP/I                                                 *
//    *Br   60 :piPpid    : piPpid[npiP]/I                                         *
//    *Br   61 :piPcharge : piPcharge[npiP]/I                                      *
//    *Br   62 :piPp      : piPp[npiP]/F                                           *
//    *Br   63 :piPpx     : piPpx[npiP]/F                                          *
//    *Br   64 :piPpy     : piPpy[npiP]/F                                          *
//    *Br   65 :piPpz     : piPpz[npiP]/F                                          *
//    *Br   66 :piPtheta  : piPtheta[npiP]/F                                       *
//    *Br   67 :piPphi    : piPphi[npiP]/F                                         *
//    *Br   68 :piPtime   : piPtime[npiP]/F                                        *
//    *Br   69 :piPpath   : piPpath[npiP]/F                                        *
//    *Br   70 :piPvz     : piPvz[npiP]/F                                          *
//    *Br   71 :piPstat   : piPstat[npiP]/I                                        *
//    *Br   72 :piPdE     : piPdE[npiP]/F                                          *
//    *Br   73 :piPbeta   : piPbeta[npiP]/F                                        *
//    *Br   74 :piPsector : piPsector[npiP]/I                                      *
//    *Br   75 :piPregion : piPregion[npiP]/I                                      *
//    *Br   76 :piPfiducial1 : piPfiducial1[npiP]/I                                *
//    *Br   77 :piPfiducial2 : piPfiducial2[npiP]/I                                *
//    *Br   78 :piPfiducial3 : piPfiducial3[npiP]/I                                *
//    *Br   79 :piPPCAL_time : piPPCAL_time[npiP]/F                                *
//    *Br   80 :piPPCAL_energy : piPPCAL_energy[npiP]/F                            *
//    *Br   81 :piPPCAL_path : piPPCAL_path[npiP]/F                                *
//    *Br   82 :piPPCAL_sector : piPPCAL_sector[npiP]/I                            *
//    *Br   83 :piPPCAL_x : piPPCAL_x[npiP]/F                                      *
//    *Br   84 :piPPCAL_y : piPPCAL_y[npiP]/F                                      *
//    *Br   85 :piPPCAL_z : piPPCAL_z[npiP]/F                                      *
//    *Br   86 :piPPCAL_lu : piPPCAL_lu[npiP]/F                                    *
//    *Br   87 :piPPCAL_lv : piPPCAL_lv[npiP]/F                                    *
//    *Br   88 :piPPCAL_lw : piPPCAL_lw[npiP]/F                                    *
//    *Br   89 :piPPCAL_status : piPPCAL_status[npiP]/I                            *
//    *Br   90 :piPECIN_time : piPECIN_time[npiP]/F                                *
//    *Br   91 :piPECIN_energy : piPECIN_energy[npiP]/F                            *
//    *Br   92 :piPECIN_path : piPECIN_path[npiP]/F                                *
//    *Br   93 :piPECIN_sector : piPECIN_sector[npiP]/I                            *
//    *Br   94 :piPECIN_x : piPECIN_x[npiP]/F                                      *
//    *Br   95 :piPECIN_y : piPECIN_y[npiP]/F                                      *
//    *Br   96 :piPECIN_z : piPECIN_z[npiP]/F                                      *
//    *Br   97 :piPECIN_lu : piPECIN_lu[npiP]/F                                    *
//    *Br   98 :piPECIN_lv : piPECIN_lv[npiP]/F                                    *
//    *Br   99 :piPECIN_lw : piPECIN_lw[npiP]/F                                    *
//    *Br  100 :piPECIN_status : piPECIN_status[npiP]/I                            *
//    *Br  101 :piPECOUT_time : piPECOUT_time[npiP]/F                              *
//    *Br  102 :piPECOUT_energy : piPECOUT_energy[npiP]/F                          *
//    *Br  103 :piPECOUT_path : piPECOUT_path[npiP]/F                              *
//    *Br  104 :piPECOUT_sector : piPECOUT_sector[npiP]/I                          *
//    *Br  105 :piPECOUT_x : piPECOUT_x[npiP]/F                                    *
//    *Br  106 :piPECOUT_y : piPECOUT_y[npiP]/F                                    *
//    *Br  107 :piPECOUT_z : piPECOUT_z[npiP]/F                                    *
//    *Br  108 :piPECOUT_lu : piPECOUT_lu[npiP]/F                                  *
//    *Br  109 :piPECOUT_lv : piPECOUT_lv[npiP]/F                                  *
//    *Br  110 :piPECOUT_lw : piPECOUT_lw[npiP]/F                                  *
//    *Br  111 :piPECOUT_status : piPECOUT_status[npiP]/I                          *
//    *Br  112 :piPdcX    : piPdcX[npiP][3]/F                                      *
//    *Br  113 :piPdcY    : piPdcY[npiP][3]/F                                      *
//    *Br  114 :piPdcZ    : piPdcZ[npiP][3]/F                                      *
//    *Br  115 :npiM      : npiM/I                                                 *
//    *Br  116 :piMpid    : piPpid[npiM]/I                                         *
//    *Br  117 :piMcharge : piPcharge[npiM]/I                                      *
//    *Br  118 :piMp      : piMp[npiM]/F                                           *
//    *Br  119 :piMpx     : piMpx[npiM]/F                                          *
//    *Br  120 :piMpy     : piMpy[npiM]/F                                          *
//    *Br  121 :piMpz     : piMpz[npiM]/F                                          *
//    *Br  122 :piMtheta  : piMtheta[npiM]/F                                       *
//    *Br  123 :piMphi    : piMphi[npiM]/F                                         *
//    *Br  124 :piMtime   : piMtime[npiM]/F                                        *
//    *Br  125 :piMpath   : piMpath[npiM]/F                                        *
//    *Br  126 :piMvz     : piMvz[npiM]/F                                          *
//    *Br  127 :piMstat   : piMstat[npiM]/I                                        *
//    *Br  128 :piMdE     : piMdE[npiM]/F                                          *
//    *Br  129 :piMbeta   : piMbeta[npiM]/F                                        *
//    *Br  130 :piMsector : piMsector[npiM]/I                                      *
//    *Br  131 :piMregion : piMregion[npiM]/I                                      *
//    *Br  132 :piMfiducial1 : piMfiducial1[npiM]/I                                *
//    *Br  133 :piMfiducial2 : piMfiducial2[npiM]/I                                *
//    *Br  134 :piMfiducial3 : piMfiducial3[npiM]/I                                *
//    *Br  135 :piMPCAL_time : piMPCAL_time[npiM]/F                                *
//    *Br  136 :piMPCAL_energy : piMPCAL_energy[npiM]/F                            *
//    *Br  137 :piMPCAL_path : piMPCAL_path[npiM]/F                                *
//    *Br  138 :piMPCAL_sector : piMPCAL_sector[npiM]/I                            *
//    *Br  139 :piMPCAL_x : piMPCAL_x[npiM]/F                                      *
//    *Br  140 :piMPCAL_y : piMPCAL_y[npiM]/F                                      *
//    *Br  141 :piMPCAL_z : piMPCAL_z[npiM]/F                                      *
//    *Br  142 :piMPCAL_lu : piMPCAL_lu[npiM]/F                                    *
//    *Br  143 :piMPCAL_lv : piMPCAL_lv[npiM]/F                                    *
//    *Br  144 :piMPCAL_lw : piMPCAL_lw[npiM]/F                                    *
//    *Br  145 :piMPCAL_status : piMPCAL_status[npiM]/I                            *
//    *Br  146 :piMECIN_time : piMECIN_time[npiM]/F                                *
//    *Br  147 :piMECIN_energy : piMECIN_energy[npiM]/F                            *
//    *Br  148 :piMECIN_path : piMECIN_path[npiM]/F                                *
//    *Br  149 :piMECIN_sector : piMECIN_sector[npiM]/I                            *
//    *Br  150 :piMECIN_x : piMECIN_x[npiM]/F                                      *
//    *Br  151 :piMECIN_y : piMECIN_y[npiM]/F                                      *
//    *Br  152 :piMECIN_z : piMECIN_z[npiM]/F                                      *
//    *Br  153 :piMECIN_lu : piMECIN_lu[npiM]/F                                    *
//    *Br  154 :piMECIN_lv : piMECIN_lv[npiM]/F                                    *
//    *Br  155 :piMECIN_lw : piMECIN_lw[npiM]/F                                    *
//    *Br  156 :piMECIN_status : piMECIN_status[npiM]/I                            *
//    *Br  157 :piMECOUT_time : piMECOUT_time[npiM]/F                              *
//    *Br  158 :piMECOUT_energy : piMECOUT_energy[npiM]/F                          *
//    *Br  159 :piMECOUT_path : piMECOUT_path[npiM]/F                              *
//    *Br  160 :piMECOUT_sector : piMECOUT_sector[npiM]/I                          *
//    *Br  161 :piMECOUT_x : piMECOUT_x[npiM]/F                                    *
//    *Br  162 :piMECOUT_y : piMECOUT_y[npiM]/F                                    *
//    *Br  163 :piMECOUT_z : piMECOUT_z[npiM]/F                                    *
//    *Br  164 :piMECOUT_lu : piMECOUT_lu[npiM]/F                                  *
//    *Br  165 :piMECOUT_lv : piMECOUT_lv[npiM]/F                                  *
//    *Br  166 :piMECOUT_lw : piMECOUT_lw[npiM]/F                                  *
//    *Br  167 :piMECOUT_status : piMECOUT_status[npiM]/I                          *
//    *Br  168 :piMdcX    : piMdcX[npiM][3]/F                                      *
//    *Br  169 :piMdcY    : piMdcY[npiM][3]/F                                      *
//    *Br  170 :piMdcZ    : piMdcZ[npiM][3]/F                                      *
//    *Br  171 :Q2        : Q2[nml]/F                                              *
//    *Br  172 :Nu        : Nu[nml]/F                                              *
//    *Br  173 :q         : q[nml]/F                                               *
//    *Br  174 :qx        : qx[nml]/F                                              *
//    *Br  175 :qy        : qy[nml]/F                                              *
//    *Br  176 :qz        : qz[nml]/F                                              *
//    *Br  177 :xB        : xB[nml]/F                                              *
//    *Br  178 :W2        : W2[nml]/F                                              *

}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void PrintDone(){
    std::cout << "Done reading RGA file " << Filename << ", execution time: "
    << double(clock() - tStart) / (double)CLOCKS_PER_SEC
    << " sec "<< std::endl;
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void OpenInputFile (){
    
    TString InputROOTFilename = DataPath + "/" + Filename + ".root";
    if (fdebug>2) std::cout << "Opening " << InputROOTFilename << std::endl;
    RGAFile = new TFile( InputROOTFilename );
    
    RGATree = (TTree*)RGAFile->Get("T");
    
    if (fdebug>2) RGATree->Print();
    
    NeventsInRGATree = RGATree->GetEntries();
    SetInputTTree();
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void OpenOutputFiles (){
    
    CSVFilename_e_e_pips = DataPath + "/" + Filename + "_e_e_pips.csv";
    CSVFilename_e_e_pims = DataPath + "/" + Filename + "_e_e_pims.csv";
    
    if (fdebug>2) {
        std::cout << "Opening output files: "
        << CSVFilename_e_e_pips
        << " and "
        << CSVFilename_e_e_pims
        << std::endl
        << csvheader
        << std::endl;
    }
    
    // Write csv header output csv files
    CSVfile_e_e_pips.open( CSVFilename_e_e_pips + ".csv" );
    CSVfile_e_e_pips << csvheader << std::endl;
    
    // Write csv header output csv files
    CSVfile_e_e_pims.open( CSVFilename_e_e_pims + ".csv" );
    CSVfile_e_e_pims << csvheader << std::endl;
        
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void InitializeFileReading (){
    Nevents_e_e_pips = Nevents_e_e_pims = 0;
    NeventsInRGATree = 0;
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void CloseInputFile (){
    RGAFile->Close();
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void CloseOutputFiles (){
    
    // close output CSV
    CSVfile_e_e_pips.close();
    CSVfile_e_e_pims.close();
    
    std::cout << "Done producing output file. They are ready in " << std::endl << DataPath << std::endl;
    std::cout << "Wrote "
    << Nevents_e_e_pips << " (e,e'π+) events"
    << " and "
    << Nevents_e_e_pims << " (e,e'π-) events." << std::endl;
}

