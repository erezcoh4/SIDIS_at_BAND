// Erez O. C., Oct-6, 2021
#ifndef __SIDISatBAND_auxiliary_H__
#define __SIDISatBAND_auxiliary_H__

#include "csv_reader.h"
#include <sstream>
#include <fstream>
#include <iostream>

#define    SHOW (x){ std::cout << (#x) << ": " << (x) << std::endl;};
#define SHOWArr (x){ for (auto _x:x) std::cout << (#_x) << ": " << (_x) << ", " << std::endl;};


//using namespace clas12;

class SIDISatBAND_auxiliary {
public:
    SIDISatBAND_auxiliary(int _fdebug_=0,int _torusBending_=-1);
    ~SIDISatBAND_auxiliary();

    Double_t Chi2PID_pion_upperBound (Double_t p, Double_t C);
    Double_t Chi2PID_pion_lowerBound (Double_t p, Double_t C);
    
    Double_t Chi2PID_Kaon_upperBound (Double_t p, Double_t C);
    Double_t Chi2PID_Kaon_lowerBound (Double_t p, Double_t C);
    
    bool eepiPassedKinematicalCriteria(Double_t Ebeam,
                                       Double_t omega,
                                       Double_t Q2,
                                       Double_t y,
                                       Double_t W,
                                       TLorentzVector pi,
                                       TLorentzVector e );
    bool eeKPassedKinematicalCriteria(Double_t Ebeam,
                                       Double_t omega,
                                       Double_t Q2,
                                       Double_t y,
                                       Double_t W,
                                       TLorentzVector K,
                                       TLorentzVector e );
    
    void                loadCutValues (std::string cutValuesFilename="cutValues.csv",
                                       int torusBending=1); //  -1 for In-bending, +1 for Out-bending
    void               printCutValues ();
    void            PrintMonitorHello ();
    void                 SetVerbosity (int _fdebug_)         {fdebug = _fdebug_;};
    void              SetTorusBending (int _torusBending_)   {torusBending = _torusBending_;};
    double               FindCutValue ( std::string cutName );
    TString           GetRunNumberSTR ( int RunNumber, TString fSkimming );
    Double_t    GetEbeamFromRunNumber ( Int_t RunNumber );
    void SetTorusBendingFromRunNumber ( Int_t RunNumber );
    void              StreamToCSVfile (std::ofstream&         csvfile,
                                       std::vector<Double_t>  observables,
                                       std::vector<int>       precisions);
    
    void                  OpenCSVfile (std::ofstream& csvfile,
                                       TString filename, std::string header);
    void                 Print4Vector ( TLorentzVector v, std::string label );
    double   ComputeLightConeFraction ( TLorentzVector p );
    double                  calcQStar ( TVector3 eP3, TVector3 piP3, double Ebeam );
    
    
    
    
    
    

    int                           fdebug;
    int                     torusBending; // -1 for In-bending, +1 for Out-bending

    
    std::vector<std::pair<std::string, double>> cutValues;
    double               cutValue_Vz_min;
    double               cutValue_Vz_max;
    double             cutValue_e_PCAL_W;
    double             cutValue_e_PCAL_V;
    double             cutValue_e_E_PCAL;
    double cutValue_SamplingFraction_min;
    double     cutValue_PCAL_ECIN_SF_min;
    double        cutValue_Ve_Vpi_dz_max;
    double               cutValue_Q2_min;
    double               cutValue_Q2_max;
    double                cutValue_W_min;
    double                cutValue_y_max;
    double          cutValue_e_theta_min;
    double          cutValue_e_theta_max;
    double         cutValue_pi_theta_min;
    double         cutValue_pi_theta_max;
    double              cutValue_Ppi_min;
    double              cutValue_Ppi_max;
    double               cutValue_Pe_min;
    double               cutValue_Pe_max;
    double              cutValue_Zpi_min;
    double              cutValue_Zpi_max;
    
    // Kaons
    double        cutValue_Ve_VK_dz_max;
    double         cutValue_K_theta_min;
    double         cutValue_K_theta_max;
    double              cutValue_PK_min; // Kaon momentum
    double              cutValue_PK_max; // Kaon momentum
    double              cutValue_ZK_min; // Kaon energy fraction
    double              cutValue_ZK_max; // Kaon energy fraction
    
    
    
    double         Me = 0.00051099895;// GeV/c2
    double        Mpi = 0.139570;     // GeV/c2
    double      Mpims = 0.13957039; // GeV/c2
    double      Mpips = 0.13957039; // GeV/c2
    double         MK = 0.493677;   // K GeV/c2
    double       MKms = 0.493677;   // K- GeV/c2
    double       MKps = 0.493677;   // K+ GeV/c2
    double         Mp = 0.938272;   // GeV/c2
    double         Mn = 0.939565;   // GeV/c2
    double         Md = 1.875;      // GeV/c2
    double        Mp2 = Mp * Mp;


protected:
    
    
private:
    
    
};

#endif
