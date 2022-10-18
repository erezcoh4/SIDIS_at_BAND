// Erez O. C., Oct-6, 2021
#include "SIDISatBAND_auxiliary.h"
#define NMAXPIONS 5 // maximal allowed number of pions
#define r2d 180./3.1415 // radians to degrees



SIDISatBAND_auxiliary::SIDISatBAND_auxiliary(int _fdebug_, int _torusBending_){
    SetVerbosity    (_fdebug_);
    SetTorusBending (_torusBending_);
}

SIDISatBAND_auxiliary::~SIDISatBAND_auxiliary(){}



// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
Double_t SIDISatBAND_auxiliary::Chi2PID_pion_lowerBound( Double_t p, Double_t C){
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
Double_t SIDISatBAND_auxiliary::Chi2PID_pion_upperBound( Double_t p, Double_t C){
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

//// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
//TVector3 SIDISatBAND_auxiliary::GetParticleVertex(clas12::region_part_ptr rp){
//    TVector3 V(rp->par()->getVx(),
//               rp->par()->getVy(),
//               rp->par()->getVz());
//    return V;
//}

//// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
//void SIDISatBAND_auxiliary::SetParticle4Momentum (TLorentzVector &p4,clas12::region_part_ptr rp){
//    p4.SetXYZM(rp->par()->getPx(),
//               rp->par()->getPy(),
//               rp->par()->getPz(),
//               p4.M());
//}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SIDISatBAND_auxiliary::loadCutValues(std::string cutValuesFilename,
                                          int torusBending){
    if (fdebug>2) {
        std::cout << "SIDISatBAND_auxiliary::loadCutValues()" << std::endl;
    }
    // read cut values csv file
    csv_reader csvr;
    cutValues = csvr.read_csv(cutValuesFilename);
    
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
    cutValue_PCAL_ECIN_SF_min       = FindCutValue("PCAL_ECIN_SF_min");
    cutValue_Ve_Vpi_dz_max          = FindCutValue("(Ve-Vpi)_z_max");
    cutValue_Q2_min                 = FindCutValue("Q2_min");
    cutValue_Q2_max                 = FindCutValue("Q2_max");
    cutValue_W_min                  = FindCutValue("W_min");
    cutValue_y_max                  = FindCutValue("y_max");
    cutValue_e_theta_min            = FindCutValue("e_theta_min");
    cutValue_e_theta_max            = FindCutValue("e_theta_max");
    cutValue_pi_theta_min           = FindCutValue("pi_theta_min");
    cutValue_pi_theta_max           = FindCutValue("pi_theta_max");
    cutValue_Ppi_min                = FindCutValue("Ppi_min");
    cutValue_Ppi_max                = FindCutValue("Ppi_max");
    cutValue_Pe_min                 = FindCutValue("Pe_min");
    cutValue_Pe_max                 = FindCutValue("Pe_max");
    cutValue_Zpi_min                = FindCutValue("Zpi_min");
    cutValue_Zpi_max                = FindCutValue("Zpi_max");
    
    if (fdebug>2) { printCutValues(); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SIDISatBAND_auxiliary::printCutValues(){
    std::cout << "Using the following cut values:" << std::endl;
    std::cout <<
    "Vz_min: "                 << cutValue_Vz_min                 << ", " << std::endl <<
    "Vz_max: "                 << cutValue_Vz_max                 << ", " << std::endl <<
    "e_PCAL_W: "               << cutValue_e_PCAL_W               << ", " << std::endl <<
    "e_PCAL_V: "               << cutValue_e_PCAL_V               << ", " << std::endl <<
    "e_E_PCAL: "               << cutValue_e_E_PCAL               << ", " << std::endl <<
    "SamplingFraction_min: "   << cutValue_SamplingFraction_min   << ", " << std::endl <<
    "PCAL_ECIN_SF_min: "       << cutValue_PCAL_ECIN_SF_min       << ", " << std::endl <<
    "Ve_Vpi_dz_max: "          << cutValue_Ve_Vpi_dz_max          << ", " << std::endl <<
    "Q2_min: "                 << cutValue_Q2_min                 << ", " << std::endl <<
    "Q2_max: "                 << cutValue_Q2_max                 << ", " << std::endl <<
    "W_min: "                  << cutValue_W_min                  << ", " << std::endl <<
    "y_max: "                  << cutValue_y_max                  << ", " << std::endl <<
    "e_theta_min: "            << cutValue_e_theta_min            << ", " << std::endl <<
    "e_theta_max: "            << cutValue_e_theta_max            << ", " << std::endl <<
    "pi_theta_min: "           << cutValue_pi_theta_min           << ", " << std::endl <<
    "pi_theta_max: "           << cutValue_pi_theta_max           << ", " << std::endl <<
    "Ppi_min: "                << cutValue_Ppi_min                << ", " << std::endl <<
    "Ppi_max: "                << cutValue_Ppi_max                << ", " << std::endl <<
    "Pe_min: "                 << cutValue_Pe_min                 << ", " << std::endl <<
    "Pe_max: "                 << cutValue_Pe_max                 << ", " << std::endl <<
    "Zpi_min: "                << cutValue_Zpi_min                << ", " << std::endl <<
    "Zpi_max: "                << cutValue_Zpi_max               << ", " << std::endl <<
    std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
double SIDISatBAND_auxiliary::FindCutValue( std::string cutName ){
    for (auto cut: cutValues) {
        if (strcmp(cut.first.c_str(),cutName.c_str())==0){
            return cut.second;
        }
    }
    return -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TString SIDISatBAND_auxiliary::GetRunNumberSTR( int RunNumber ){
    char RunNumberStr[20];
    sprintf( RunNumberStr, "00%d", RunNumber );
    if (fdebug>1) std::cout << "(SIDIS) skimming run " << RunNumberStr << std::endl;
    return (TString)RunNumberStr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Double_t SIDISatBAND_auxiliary::GetEbeamFromRunNumber ( Int_t RunNumber ){
    if (6420 <= RunNumber && RunNumber <= 6598){
        return 10.2; // GeV
    }
    else if (11362 <= RunNumber && RunNumber <= 11571){
        return 10.4; // GeV
    }
    else if (6164 <= RunNumber && RunNumber <= 6399){
        return 10.6; // GeV
    }
    else{
        return 0;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SIDISatBAND_auxiliary::SetTorusBendingFromRunNumber ( Int_t RunNumber ){
    // -1 for In-bending, +1 for Out-bending
    // For BAND data
    // Spring 19 and Spring 2020 was in-bending
    // Fall 2019 (without low-energy-run) was out-bending
    
    if (6420 <= RunNumber && RunNumber <= 6598){
        this->torusBending = -1;
    }
    else if (11362 <= RunNumber && RunNumber <= 11571){
        this->torusBending = -1;
    }
    else if (6164 <= RunNumber && RunNumber <= 6399){
        this->torusBending = +1;
    }
    else{
        this->torusBending = 0;
    }
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void SIDISatBAND_auxiliary::StreamToCSVfile (std::ofstream&         csvfile,
                                             std::vector<Double_t>  observables,
                                             std::vector<int>       precisions){
    
    //    for (auto v:observables) csvfile << v << ",";
    //    csvfile << std::endl;
    
    for (int j=0; j < observables.size(); j++){
        int precision = 9;
        if (j < precisions.size()) precision = precisions.at(j);
        auto v = observables.at(j);
        csvfile << std::setprecision(precision) << std::fixed << v << ",";
    }
    csvfile << std::endl;
    
    
    if (fdebug>3) {
        std::cout << "StreamToEventCSVfile()" << std::endl;
        for (auto v:observables) std::cout << v << ",";
        std::cout << std::endl;
    }
}



// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void SIDISatBAND_auxiliary::OpenCSVfile (std::ofstream& csvfile,
                                         TString filename,
                                         std::string header){
    csvfile.open( filename );
    csvfile << header << "," << std::endl;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SIDISatBAND_auxiliary::Print4Vector( TLorentzVector v, std::string label ){
    std::cout << label << " 4-vector:"<<std::endl;
    std::cout
    << std::setprecision(2)
    << "(Px = " << v.Px()
    << ", Py = " << v.Py()
    << ", Pz = " << v.Pz()
    << ", E = " << v.E()
    << "), M = " << v.Mag()
    << " GeV"
    << std::endl;
}




// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
double SIDISatBAND_auxiliary::ComputeLightConeFraction( TLorentzVector p ){
    // compute light-cone momentum fraction
    double m = p.Mag();
    double alpha = (p.E() - p.Z())/m;
    return alpha;
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void SIDISatBAND_auxiliary::PrintMonitorHello(){
    std::cout << "Hello..." << std::endl;
    std::cout << "Is it me you're looking for?..." << std::endl;
    std::cout << "I can see it in your eyes..." << std::endl;
    std::cout << "I can see it in your smile..." << std::endl;
    std::cout << "You're all I've ever wanted, and your arms are open wide..." << std::endl;
    std::cout << "Cause you know just what to say, and you know just what to do..." << std::endl;
    std::cout << "And I want to tell you so much, I love you" << std::endl;
    std::cout << "I long to see the sunlight in your hair" << std::endl;
}





// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
bool SIDISatBAND_auxiliary::eepiPassedKinematicalCriteria(Double_t Ebeam,
                                                          Double_t omega,
                                                          Double_t Q2,
                                                          Double_t y,
                                                          Double_t W,
                                                          TLorentzVector pi,
                                                          TLorentzVector e){
    Double_t Zpi = pi.E()/omega;
    
    if (fdebug>2) {
        std::cout
        << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        << std::endl
        << "eepiPassedKinematicalCriteria()"
        << std::endl
        << std::setprecision(3)
        << "Q²: "       << Q2 << " (GeV/c)²,"
        << "W: "        << W << " GeV/c²,"
        << "y: "        << y << ","
        << std::endl
        << "𝜃(e): "    << e.Theta()*r2d    << "˚,"
        << "p(e): "    << e.P()            << " GeV/c,"
        << std::endl
        << "𝜃(pi): "   << pi.Theta()*r2d   << "˚,"
        << "p(π): "    << pi.P()           << " GeV/c,"
        << std::endl
        << "z(π): "    << Zpi              << ","
        << std::endl;
        
        if (fdebug>4) {
            std::cout
            << cutValue_Q2_min << " < " << "Q2="<< Q2 << " < " << cutValue_Q2_max
            << std::endl
            << cutValue_W_min  << " < " << "W="<<W
            << std::endl
            << "y="<<y               << " < " << cutValue_y_max
            << std::endl
            << cutValue_e_theta_min
            << " < "    << "𝜃(e) = "<< e.Theta()*r2d << "˚"
            << " < "    << cutValue_e_theta_max
            << std::endl
            << cutValue_pi_theta_min
            << " < "    << "𝜃(π) = " <<pi.Theta()*r2d << "˚"
            << " < "    << cutValue_pi_theta_max
            << std::endl
            << cutValue_Pe_min  << " < "  << "p(e)="<< e.P()  << " < " << Ebeam
            << std::endl
            << cutValue_Zpi_min << " < "  << "z(π)="<< Zpi    << " < " << cutValue_Zpi_max
            << std::endl;
            
            if ( cutValue_Q2_min < Q2  && Q2 < cutValue_Q2_max  )
                std::cout << "1 passed Q2 cut" << std::endl;
            else
                std::cout << "1 failed Q2 cut" << std::endl;
            if ( cutValue_W_min < W  )
                std::cout << "2 passed W cut" << std::endl;
            else
                std::cout << "2 failed W cut" << std::endl;
            if ( y < cutValue_y_max  )
                std::cout << "3 passed y cut" << std::endl;
            else
                std::cout << "3 failed y cut" << std::endl;
            if ( cutValue_e_theta_min < e.Theta()*r2d  &&  e.Theta()*r2d < cutValue_e_theta_max  )
                std::cout << "4 passed 𝜃(e) cut" << std::endl;
            else
                std::cout << "4 failed 𝜃(e) cut" << std::endl;
            if (cutValue_pi_theta_min < pi.Theta()*r2d && pi.Theta()*r2d < cutValue_pi_theta_max )
                std::cout << "5 passed 𝜃(π) cut" << std::endl;
            else
                std::cout << "5 failed 𝜃(π) cut" << std::endl;
            if (cutValue_Ppi_min < pi.P()         &&         pi.P() < cutValue_Ppi_max)
                std::cout << "6 passed p(π) cut" << std::endl;
            else
                std::cout << "6 failed p(π) cut" << std::endl;
            if (cutValue_Pe_min < e.P()          &&          e.P() < Ebeam)
                std::cout << "7 passed p(e) cut" << std::endl;
            else
                std::cout << "7 failed p(e) cut" << std::endl;
            if (cutValue_Zpi_min < Zpi            &&            Zpi < cutValue_Zpi_max)
                std::cout << "8 passed Z(π) cut" << std::endl;
            else
                std::cout << "8 failed Z(π) cut" << std::endl;
        }
        
    }
    
    if(   (      cutValue_Q2_min < Q2             &&             Q2 < cutValue_Q2_max       )
       && (       cutValue_W_min < W                                                        )
       && (                                                       y < cutValue_y_max        )
       && ( cutValue_e_theta_min < e.Theta()*r2d  &&  e.Theta()*r2d < cutValue_e_theta_max  )
       && (cutValue_pi_theta_min < pi.Theta()*r2d && pi.Theta()*r2d < cutValue_pi_theta_max )
       && (     cutValue_Ppi_min < pi.P()         &&         pi.P() < cutValue_Ppi_max      )
       && (      cutValue_Pe_min < e.P()          &&          e.P() < Ebeam                 )
       && (     cutValue_Zpi_min < Zpi            &&            Zpi < cutValue_Zpi_max      )
       ) {
        if (fdebug>2) {
            std::cout << "succesfully passed (e,e'π) kinematical cuts"
            << std::endl
            << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            << std::endl;
        }
        return true;
    }
    else {
        if (fdebug>2) {
            std::cout << "Did not pass (e,e'π) kinematical cuts"
            << std::endl
            << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            << std::endl;
        }
    }
    
    return false;
}



// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
double calcQStar(TVector3 eP3, TVector3 piP3, double Ebeam){
    //
    // Compute the novel SIDIS observable, q∗
    // designed to be maximally resilient against resolution
    // effects while delivering the same sensitivity to TMD dynamics as Pt
    //
    // [Gao et al., "A Better Angle on Hadron Transverse Momentum Distributions at the EIC"]
    // arXiv:2209.11211v1
    //
    // by Natalie Wright Oct-12, 2022
    //
    // input
    // ------
    // eP3      TVector3    electron momentum (in q-Frame)
    // piP3     TVector3    π momentum (in q-Frame)
    //
    // return
    // ------
    // qstar    double      coplanarity momentum transfer part
    //
    
    
    
    double tan_phi_acop = piP3.Y()/piP3.X(); // EIC Frame
    double       eta_pi = TMath::ATanH(-piP3.Z()/piP3.Mag()); // - sign from EIC frame
    double        eta_e = TMath::ATanH(-eP3.Z()/eP3.Mag());
    double    delta_eta = eta_pi - eta_e;
    
    double        qstar =   2 * Ebeam * TMath::Exp(eta_pi) * tan_phi_acop
                            /
                            (1 + TMath::Exp(delta_eta));
    
    return qstar;    
}
