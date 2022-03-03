#ifndef NPIEVENT_H
#define NPIEVENT_H

#define NMAXPIONS 20 // maximal allowed number of pions

#include <utility>
#include <vector>
#include <string>

#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TVector3.h>

#include "skimHelper.C"

class NPiEvent {

public:
    NPiEvent();

    void SetOutputTTrees(TTree* outTree, bool piplus);
    void InitializeVariables();
    void ComputeKinematics();
    
    bool getEPastCutsInEvent () { return ePastCutsInEvent; }

    setBeamHelicity( double beam_helicity_in ) { beam_helicity = beam_helicity_in; }

    void setNe( int Ne_in ) { Ne = Ne_in; }
    void setNn( int Nn_in ) { Nn = Nn_in; }
    void setNp( int Np_in ) { Np = Np_in; }
    void setNpips( int Npips_in ) { Npips = Npips_in; }
    void setNpims( int Npims_in ) { Npims = Npims_in; }
    void setNgammas( int Ngammas_in ) { Ngammas = Ngammas_in; }
    void setNd( int Nd_in ) { Nd = Nd_in; }

    void setE ( clas12::region_part_ptr e_in ) { SetLorentzVector(e , e_in); }
    void setVe( TVector3 Ve_in ) { Ve = Ve_in; }

    double getE_PCAL_x() { return e_PCAL_x; }
    double getE_PCAL_y() { return e_PCAL_y; }
    double getE_PCAL_W() { return e_PCAL_W; }
    double getE_PCAL_V() { return e_PCAL_V; }
    double getE_E_PCAL() { return e_E_PCAL; }
    double getE_E_ECIN() { return e_E_ECIN; }
    double getE_E_ECOUT() { return e_E_ECOUT; }
    double getE_PCAL_sector() { return e_PCAL_sector; }

    double[3] getE_DC_x() { return e_DC_x; }
    double[3] getE_DC_y() { return e_DC_y; }
    double[3] getE_DC_z() { return e_DC_z; }

    TLorentzVector getE() { return e; } 
    TVector3 getVe() { return Ve; }

    void setE_E_PCAL ( double e_E_PCAL_in ) { e_E_PCAL = e_E_PCAL_in; }
    void setE_PCAL_sector ( double e_PCAL_sector_in ) { e_PCAL_sector = e_PCAL_sector_in; }
    void setE_PCAL_V ( double e_PCAL_V_in ) { e_PCAL_V = e_PCAL_V_in; }
    void setE_PCAL_W ( double e_PCAL_W_in ) { e_PCAL_W = e_PCAL_W_in; }
    void setE_E_ECIN ( double e_E_ECIN_in ) { e_E_ECIN = e_E_ECIN_in; }
    void setE_E_ECOUT ( double e_E_ECOUT_in ) { e_E_ECOUT = e_E_ECOUT_in; }
    void setE_PCAL_x ( double e_PCAL_x_in ) { e_PCAL_x = e_PCAL_x_in; }
    void setE_PCAL_y ( double e_PCAL_y_in ) { e_PCAL_y = e_PCAL_y_in; }
    void setE_PCAL_z ( double e_PCAL_z_in ) { e_PCAL_z = e_PCAL_z_in; }
    void setE_DC_sector ( double e_DC_sector_in ) { e_DC_sector = e_DC_sector_in; }
    void setE_DC_Chi2N ( double e_DC_Chi2N_in ) { e_DC_Chi2N = e_DC_Chi2N_in; }
    void setE_DC_x ( double e_DC_x_in, int regionIdx ) { e_DC_x[regionIdx] = e_DC_x_in; }
    void setE_DC_y ( double e_DC_y_in, int regionIdx ) { e_DC_y[regionIdx] = e_DC_y_in; }
    void setE_DC_z ( double e_DC_z_in, int regionIdx ) { e_DC_z[regionIdx] = e_DC_z_in; }

    void setEPastCutsInEvent ( bool ePastCutsInEvent_in ) { ePastCutsInEvent = ePastCutsInEvent_in; }


    
    void setPiplus ( clas12::region_part_ptr piplus_in, int pipsIdx ) { SetLorentzVector(piplus[pipsIdx] , piplus_in); }
        
    void setZpips ( double Zpips_in, int pipsIdx) { Zpips[pipsIdx] = Zpips_in; }
    void setVpiplus ( double Vpiplus_in, int pipsIdx) { Vpiplus[pipsIdx] = Vpiplus_in; }

    void setPips_chi2PID ( double pips_chi2PID_in, int pipsIdx ) { pips_chi2PID[pipsIdx] = pips_chi2PID_in; }
    void setPips_E_ECIN ( double pips_E_ECIN_in, int pipsIdx ) { pips_E_ECIN[pipsIdx] = pips_E_ECIN_in; }
    void setPips_E_ECOUT ( double pips_E_ECOUT_in, int pipsIdx ) { pips_E_ECOUT[pipsIdx] = pips_E_ECOUT_in; }
    void setPips_E_PCAL ( double pips_E_PCAL_in, int pipsIdx ) { pips_E_PCAL[pipsIdx] = pips_E_PCAL_in; }
    void setPips_PCAL_sector ( double pips_PCAL_sector_in, int pipsIdx ) { pips_PCAL_sector[pipsIdx] = pips_PCAL_sector_in; }
    void setPips_PCAL_V ( double pips_PCAL_V_in, int pipsIdx ) { pips_PCAL_V[pipsIdx] = pips_PCAL_V_in; }
    void setPips_PCAL_W ( double pips_PCAL_W_in, int pipsIdx ) { pips_PCAL_W[pipsIdx] = pips_PCAL_W_in; }
    void setPips_PCAL_x ( double pips_PCAL_x_in, int pipsIdx ) { pips_PCAL_x[pipsIdx] = pips_PCAL_x_in; }
    void setPips_PCAL_y ( double pips_PCAL_y_in, int pipsIdx ) { pips_PCAL_y[pipsIdx] = pips_PCAL_y_in; }
    void setPips_PCAL_z ( double pips_PCAL_z_in, int pipsIdx ) { pips_PCAL_z[pipsIdx] = pips_PCAL_z_in; }
    void setPips_DC_sector ( double pips_DC_sector_in, int pipsIdx ) { pips_DC_sector[pipsIdx] = pips_DC_sector_in; }
    void setPips_Chi2N ( double pips_Chi2N_in, int pipsIdx ) { pips_Chi2N[pipsIdx] = pips_Chi2N_in; }

    void setPips_DC_x ( double pips_DC_x_in, int pipsIdx, int regionIdx ) { pips_DC_x[pipsIdx][regionIdx] = pips_DC_x_in; }
    void setPips_DC_y ( double pips_DC_y_in, int pipsIdx, int regionIdx ) { pips_DC_y[pipsIdx][regionIdx] = pips_DC_y_in; }
    void setPips_DC_z ( double pips_DC_z_in, int pipsIdx, int regionIdx ) { pips_DC_z[pipsIdx][regionIdx] = pips_DC_z_in; }
        
    void setPipsPastSelectionCuts ( bool pipsPastSelectionCuts_in, int pipsIdx ) { pipsPastSelectionCuts = pipsPastSelectionCuts_in; }
    void setEepipsPastKinematicalCuts ( bool eepipsPastKinematicalCuts_in, int pipsIdx ) { eepipsPastKinematicalCuts = eepipsPastKinematicalCuts_in; }
    void setPipsPastCutsInEvent ( bool pipsPastCutsInEvent_in ) { pipsPastCutsInEvent = pipsPastCutsInEvent_in; }
    void setEepipsPastCutsInEvent ( bool eepipsPastCutsInEvent_in ) { eepipsPastCutsInEvent = eepipsPastCutsInEvent_in; }

        
    void setPiplus_Px ( double piplus_Px_in, int pipsIdx ) { piplus_Px[pipsIdx] = piplus_Px_in; }
    void setPiplus_Py ( double piplus_Py_in, int pipsIdx ) { piplus_Py[pipsIdx] = piplus_Py_in; }
    void setPiplus_Pz ( double piplus_Pz_in, int pipsIdx ) { piplus_Pz[pipsIdx] = piplus_Pz_in; }
    void setPiplus_E ( double piplus_E_in, int pipsIdx ) { piplus_E[pipsIdx] = piplus_E_in; }
    void setVpiplus_X ( double Vpiplus_X, int pipsIdx ) { Vpiplus_X[pipsIdx] = Vpiplus_X_in; }
    void setVpiplus_Y ( double Vpiplus_Y, int pipsIdx ) { Vpiplus_Y[pipsIdx] = Vpiplus_Y_in; }
    void setVpiplus_Z ( double Vpiplus_Z, int pipsIdx ) { Vpiplus_Z[pipsIdx] = Vpiplus_Z_in; }

    double getPips_DC_sector (int pipsIdx) { return pips_DC_sector[pipsIdx]; }
    double getPips_DC_x ( int pipsIdx, int regionIdx ) { return pips_DC_x[pipsIdx][regionsIdx]; }
    double getPips_DC_y ( int pipsIdx, int regionIdx ) { return pips_DC_y[pipsIdx][regionsIdx]; }
    double getPips_DC_z ( int pipsIdx, int regionIdx ) { return pips_DC_z[pipsIdx][regionsIdx]; }

    bool getPipsPastSelectionCuts() { return pipsPastSelectionCuts; }

    double getPips_chi2PID( int pipsIdx ) { return pips_chi2PID[pipsIdx]; }
    TLorentzVector getPiplus ( int pipsIdx ) { return piplus[pipsIdx]; }
    TVector3 getVpiplus( int pipsIdx ) { return Vpiplus[pipsIdx]; }

    bool getEepipsPastKinematicalCuts( int pipsIdx ) { return eepipsPastKinematicalCuts[pipsIdx]; }

//////////////////////////

    
    void setPiminus ( clas12::region_part_ptr piminus_in, int pimsIdx ) { SetLorentzVector(piminus[pimsIdx] , piminus_in); }
        
    void setZpims ( double Zpims_in, int pimsIdx) { Zpims[pimsIdx] = Zpims_in; }
    void setVpiminus ( double Vpiminus_in, int pimsIdx) { Vpiminus[pimsIdx] = Vpiminus_in; }

    void setPims_chi2PID ( double pims_chi2PID_in, int pimsIdx ) { pims_chi2PID[pimsIdx] = pims_chi2PID_in; }
    void setPims_E_ECIN ( double pims_E_ECIN_in, int pimsIdx ) { pims_E_ECIN[pimsIdx] = pims_E_ECIN_in; }
    void setPims_E_ECOUT ( double pims_E_ECOUT_in, int pimsIdx ) { pims_E_ECOUT[pimsIdx] = pims_E_ECOUT_in; }
    void setPims_E_PCAL ( double pims_E_PCAL_in, int pimsIdx ) { pims_E_PCAL[pimsIdx] = pims_E_PCAL_in; }
    void setPims_PCAL_sector ( double pims_PCAL_sector_in, int pimsIdx ) { pims_PCAL_sector[pimsIdx] = pims_PCAL_sector_in; }
    void setPims_PCAL_V ( double pims_PCAL_V_in, int pimsIdx ) { pims_PCAL_V[pimsIdx] = pims_PCAL_V_in; }
    void setPims_PCAL_W ( double pims_PCAL_W_in, int pimsIdx ) { pims_PCAL_W[pimsIdx] = pims_PCAL_W_in; }
    void setPims_PCAL_x ( double pims_PCAL_x_in, int pimsIdx ) { pims_PCAL_x[pimsIdx] = pims_PCAL_x_in; }
    void setPims_PCAL_y ( double pims_PCAL_y_in, int pimsIdx ) { pims_PCAL_y[pimsIdx] = pims_PCAL_y_in; }
    void setPims_PCAL_z ( double pims_PCAL_z_in, int pimsIdx ) { pims_PCAL_z[pimsIdx] = pims_PCAL_z_in; }
    void setPims_DC_sector ( double pims_DC_sector_in, int pimsIdx ) { pims_DC_sector[pimsIdx] = pims_DC_sector_in; }
    void setPims_Chi2N ( double pims_Chi2N_in, int pimsIdx ) { pims_Chi2N[pimsIdx] = pims_Chi2N_in; }

    void setPims_DC_x ( double pims_DC_x_in, int pimsIdx, int regionIdx ) { pims_DC_x[pimsIdx][regionIdx] = pims_DC_x_in; }
    void setPims_DC_y ( double pims_DC_y_in, int pimsIdx, int regionIdx ) { pims_DC_y[pimsIdx][regionIdx] = pims_DC_y_in; }
    void setPims_DC_z ( double pims_DC_z_in, int pimsIdx, int regionIdx ) { pims_DC_z[pimsIdx][regionIdx] = pims_DC_z_in; }
        
    void setPimsPastSelectionCuts ( bool pimsPastSelectionCuts_in, int pimsIdx ) { pimsPastSelectionCuts = pimsPastSelectionCuts_in; }
    void setEepimsPastKinematicalCuts ( bool eepimsPastKinematicalCuts_in, int pimsIdx ) { eepimsPastKinematicalCuts = eepimsPastKinematicalCuts_in; }
    void setPimsPastCutsInEvent ( bool pimsPastCutsInEvent_in ) { pimsPastCutsInEvent = pimsPastCutsInEvent_in; }
    void setEepimsPastCutsInEvent ( bool eepimsPastCutsInEvent_in ) { eepimsPastCutsInEvent = eepimsPastCutsInEvent_in; }

        
    void setPiminus_Px ( double piminus_Px_in, int pimsIdx ) { piminus_Px[pimsIdx] = piminus_Px_in; }
    void setPiminus_Py ( double piminus_Py_in, int pimsIdx ) { piminus_Py[pimsIdx] = piminus_Py_in; }
    void setPiminus_Pz ( double piminus_Pz_in, int pimsIdx ) { piminus_Pz[pimsIdx] = piminus_Pz_in; }
    void setPiminus_E ( double piminus_E_in, int pimsIdx ) { piminus_E[pimsIdx] = piminus_E_in; }
    void setVpiminus_X ( double Vpiminus_X, int pimsIdx ) { Vpiminus_X[pimsIdx] = Vpiminus_X_in; }
    void setVpiminus_Y ( double Vpiminus_Y, int pimsIdx ) { Vpiminus_Y[pimsIdx] = Vpiminus_Y_in; }
    void setVpiminus_Z ( double Vpiminus_Z, int pimsIdx ) { Vpiminus_Z[pimsIdx] = Vpiminus_Z_in; }

    double getPims_DC_sector (int pimsIdx) { return pims_DC_sector[pimsIdx]; }
    double getPims_DC_x ( int pimsIdx, int regionIdx ) { return pims_DC_x[pimsIdx][regionsIdx]; }
    double getPims_DC_y ( int pimsIdx, int regionIdx ) { return pims_DC_y[pimsIdx][regionsIdx]; }
    double getPims_DC_z ( int pimsIdx, int regionIdx ) { return pims_DC_z[pimsIdx][regionsIdx]; }

    bool getPimsPastSelectionCuts() { return pimsPastSelectionCuts; }

    double getPims_chi2PID( int pimsIdx ) { return pims_chi2PID[pimsIdx]; }
    TLorentzVector getPiminus ( int pimsIdx ) { return piminus[pimsIdx]; }
    TVector3 getVpiminus( int pimsIdx ) { return Vpiminus[pimsIdx]; }

    bool getEepimsPastKinematicalCuts( int pimsIdx ) { return eepimsPastKinematicalCuts[pimsIdx]; }

private:

    bool        ePastCutsInEvent = false;
    bool     pipsPastCutsInEvent = false;
    bool   eepipsPastCutsInEvent = false;
    bool     pimsPastCutsInEvent = false;
    bool         EventPassedCuts = false;
    bool   eepimsPastCutsInEvent = false;

    // meta-data
    int           torusBending = -1; // -1 for In-bending, +1 for Out-bending
    int    DC_layers[3] = {6,18,36};// Region 1 is denoted at DC detector 6, Region 2 is denoted 18, Region 3 - as 36
    int                    DC_layer;
    int                      runnum;
    int                       evnum;
    int               beam_helicity; // helicity of the electron +1 along the beam and -1 opposite to it
    int                      status;
    int                   inclusive; // tag to look at inclusive run

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

    // vectors in lab-frame
    TLorentzVector          Beam, target, e, q;
    std::vector<TLorentzVector>         piplus; // positive pions
    std::vector<TLorentzVector>        piminus; // negative pions
    TClonesArray *                 piplusArray; // positive pions
    TClonesArray *                 piminusArray; // negative pions
    // reconstructed vertex position
    TVector3                                Ve;
    std::vector<TVector3>              Vpiplus;
    std::vector<TVector3>             Vpiminus;
    TClonesArray *                VpiplusArray;
    TClonesArray *               VpiminusArray;

    // kinematics
    Double_t     Ebeam, xB, Q2, omega, W, W2, xF, y, M_X;

}

#endif