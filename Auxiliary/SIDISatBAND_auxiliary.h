// Erez O. C., Oct-6, 2021
#ifndef __SIDISatBAND_auxiliary_H__
#define __SIDISatBAND_auxiliary_H__

#include "csv_reader.h"
#include <sstream>
#include <fstream>
#include <iostream>

#include "clas12reader.h"
using namespace clas12;

class SIDISatBAND_auxiliary {
public:
    SIDISatBAND_auxiliary(int _fdebug_=0,int _torusBending_=-1);
    ~SIDISatBAND_auxiliary();

    Double_t Chi2PID_pion_upperBound (Double_t p, Double_t C);
    Double_t Chi2PID_pion_lowerBound (Double_t p, Double_t C);
        
    TVector3        GetParticleVertex (clas12::region_part_ptr rp);
    void         SetParticle4Momentum (TLorentzVector &p4,clas12::region_part_ptr rp);
    void                loadCutValues (std::string cutValuesFilename="cutValues.csv");
    void               printCutValues ();
    void                 SetVerbosity (int _fdebug_)         {fdebug = _fdebug_;};
    void              SetTorusBending (int _torusBending_)   {torusBending = _torusBending_;};
    double               FindCutValue ( std::string cutName );
    TString           GetRunNumberSTR ( int RunNumber );
    Double_t    GetEbeamFromRunNumber ( Int_t RunNumber );
    void SetTorusBendingFromRunNumber ( Int_t RunNumber );
    void              StreamToCSVfile (std::ofstream& csvfile, std::vector<Double_t> observables);
    void                  OpenCSVfile (std::ofstream& csvfile, TString filename, std::string header);

    int                           fdebug;
    int                     torusBending; // -1 for In-bending, +1 for Out-bending

    
    std::vector<std::pair<std::string, double>> cutValues;
    double               cutValue_Vz_min;
    double               cutValue_Vz_max;
    double             cutValue_e_PCAL_W;
    double             cutValue_e_PCAL_V;
    double             cutValue_e_E_PCAL;
    double cutValue_SamplingFraction_min;
    double        cutValue_Ve_Vpi_dz_max;
    double               cutValue_Q2_min;
    double                cutValue_W_min;
    double                cutValue_y_max;
    double          cutValue_e_theta_min;
    double          cutValue_e_theta_max;
    double         cutValue_pi_theta_min;
    double         cutValue_pi_theta_max;
    double              cutValue_Ppi_min;
    double              cutValue_Ppi_max;
    double              cutValue_Zpi_min;
    double              cutValue_Zpi_max;
    
protected:
    
    
private:
    
    
};

#endif
