
#include <vector>
#include <utility>

#include "Auxiliary/DCfid_SIDIS.cpp"

#define r2d 180./3.1415 // radians to degrees

class TaggedSIDISCuts {
public:
    bool CheckIfElectronPassedSelectionCuts(Double_t e_PCAL_x, Double_t e_PCAL_y, Double_t e_PCAL_W, 
                                                Double_t e_PCAL_V, Double_t e_E_PCAL, Double_t e_E_ECIN, 
                                                Double_t e_E_ECOUT, TLorentzVector e, TVector3 Ve, Double_t e_DC_sector, 
                                                Double_t e_DC_x[3], Double_t e_DC_y[3], Double_t e_DC_z[3], int torusBending);
    bool CheckIfPionPassedSelectionCuts(TString pionCharge,  Double_t DC_sector, Double_t DC_x[3],  Double_t DC_y[3], 
                                                Double_t DC_z[3], Double_t chi2PID, Double_t p, TVector3 Ve,  TVector3 Vpi, 
                                                int fdebug);
    bool eepiPassedKinematicalCriteria(TLorentzVector pi, int fdebug);
    Double_t Chi2PID_pion_lowerBound( Double_t p, Double_t C);
    Double_t Chi2PID_pion_upperBound( Double_t p, Double_t C);

    void setCutValue_Vz_min(double cutValue_Vz_min_in) { cutValue_Vz_min = cutValue_Vz_min_in; }
    void setCutValue_Vz_max(double cutValue_Vz_max_in) { cutValue_Vz_max = cutValue_Vz_max_in; }
    void setCutValue_e_PCAL_W(double cutValue_e_PCAL_W_in) { cutValue_e_PCAL_W = cutValue_e_PCAL_W_in; }
    void setCutValue_e_PCAL_V(double cutValue_e_PCAL_V_in) { cutValue_e_PCAL_V = cutValue_e_PCAL_V_in; }
    void setCutValue_e_E_PCAL(double cutValue_e_E_PCAL_in) { cutValue_e_E_PCAL = cutValue_e_E_PCAL_in; }
    void setCutValue_SamplingFraction_min(double cutValue_SamplingFraction_min_in) { cutValue_SamplingFraction_min = cutValue_SamplingFraction_min_in; }
    void setCutValue_PCAL_ECIN_SF_min(double cutValue_PCAL_ECIN_SF_min_in) { cutValue_PCAL_ECIN_SF_min = cutValue_PCAL_ECIN_SF_min_in; }
    void setCutValue_Ve_Vpi_dz_max(double cutValue_Ve_Vpi_dz_max_in) { cutValue_Ve_Vpi_dz_max = cutValue_Ve_Vpi_dz_max_in; }
    void setCutValue_Q2_min(double cutValue_Q2_min_in) { cutValue_Q2_min = cutValue_Q2_min_in; }
    void setCutValue_W_min(double cutValue_W_min_in) { cutValue_W_min = cutValue_W_min_in; }
    void setCutValue_y_max(double cutValue_y_max_in) { cutValue_y_max = cutValue_y_max_in; }
    void setCutValue_e_theta_min(double cutValue_e_theta_min_in) { cutValue_e_theta_min = cutValue_e_theta_min_in; }
    void setCutValue_e_theta_max(double cutValue_e_theta_max_in) { cutValue_e_theta_max = cutValue_e_theta_max_in; }
    void setCutValue_pi_theta_min(double cutValue_pi_theta_min_in) { cutValue_pi_theta_min = cutValue_pi_theta_min_in; }
    void setCutValue_pi_theta_max(double cutValue_pi_theta_max_in) { cutValue_pi_theta_max = cutValue_pi_theta_max_in; }
    void setCutValue_Ppi_min(double cutValue_Ppi_min_in) { cutValue_Ppi_min = cutValue_Ppi_min_in; }
    void setCutValue_Ppi_max(double cutValue_Ppi_max_in) { cutValue_Ppi_max = cutValue_Ppi_max_in; }
    void setCutValue_Zpi_min(double cutValue_Zpi_min_in) { cutValue_Zpi_min = cutValue_Zpi_min_in; }
    void setCutValue_Zpi_max(double cutValue_Zpi_max_in) { cutValue_Zpi_max = cutValue_Zpi_max_in; }

private:
    DCfid_SIDIS dcfid;

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

};

bool TaggedSIDISCuts::CheckIfElectronPassedSelectionCuts(Double_t e_PCAL_x, Double_t e_PCAL_y,
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
    
    
    if(!(true
       // fiducial cuts on PCAL
       //fabs(e_PCAL_x)>0
       //&&  fabs(e_PCAL_y)>0
       &&  e_PCAL_W > cutValue_e_PCAL_W
       &&  e_PCAL_V > cutValue_e_PCAL_V
       
       // Electron Identification Refinement  - PCAL Minimum Energy Deposition Cut
       &&  e_E_PCAL > cutValue_e_E_PCAL
       
       // Sampling fraction cut
       && ((e_E_PCAL + e_E_ECIN + e_E_ECOUT)/e.P()) > cutValue_SamplingFraction_min
       && (e_E_ECIN/e.P() > cutValue_PCAL_ECIN_SF_min - e_E_PCAL/e.P()) // RGA AN puts "<" here mistakenly
       
       // Cut on z-vertex position: in-bending torus field -13.0 cm < Vz < +12.0 cm
       // Spring 19 and Spring 2020 in-bending.
       // Fall 2019 (without low-energy-run) was out-bending.
       &&  ((cutValue_Vz_min < Ve.Z()) && (Ve.Z() < cutValue_Vz_max))
       )) return false;
    
    return true;
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
bool TaggedSIDISCuts::CheckIfPionPassedSelectionCuts(TString pionCharge, // "pi+" or "pi-"
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
    //std::cout << DC_sector << " " << DC_x[0] << " " << DC_y[0] << " " << DC_z[0] << std::endl;
    if (fdebug>3) {
        std::cout << "CheckIfPionPassedSelectionCuts()" << std::endl;
    }
    if (DC_sector == 0) { if (fdebug>2){std::cout << "DC_sector=0 (funny...)" << std::endl;} return false;}
    
    int PDGcode;
    double    C;
    if (pionCharge=="pi+"){
        PDGcode = 211;
        C       = 0.88;
    } else if (pionCharge=="pi-") {
        PDGcode = -211;
        C       = 0.93;
    } else {
        std::cout << "pion charge is not defined in CheckIfPionPassedSelectionCuts(), returning false" << std::endl;
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
            std::cout << "dcfid.DC_fid_th_ph_sidis(): "
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
//    return true;
    
    if (fdebug>3) {
        std::cout << "in CheckIfPionPassedSelectionCuts()"<< std::endl
        << "pion charge: "          << pionCharge               << ","
        << "DC_x[0]: "              << DC_x[0]                  << ","
        << "chi2PID:"               << chi2PID                  << ","
        << "Chi2PID_pion_lowerBound( p="<<p<<", C="<<C<<" ): "  << Chi2PID_pion_lowerBound( p, C ) << ","
        << "Chi2PID_pion_upperBound( p="<<p<<", C="<<C<<" ): "  << Chi2PID_pion_upperBound( p, C ) << ","
        << "fabs((Ve-Vpi).Z()): "   << fabs((Ve-Vpi).Z())       << ","
        << std::endl;
    }
    if(!
       // pi+ Identification Refinement - chi2PID vs. momentum
       (( Chi2PID_pion_lowerBound( p, C ) < chi2PID && chi2PID < Chi2PID_pion_upperBound( p , C ) )
       
       // Cut on the z-Vertex Difference Between Electrons and Hadrons.
       &&  ( fabs((Ve-Vpi).Z()) < cutValue_Ve_Vpi_dz_max )
       )) {
        return false;
    }
    if (fdebug>3) { std::cout << "succesfully passed CheckIfPionPassedSelectionCuts(), return true" << std::endl; }
    return true;
}

bool TaggedSIDISCuts::eepiPassedKinematicalCriteria(TLorentzVector pi, int fdebug){
    double Zpi = pi.E()/omega;
    if(   (      cutValue_Q2_min < Q2)
       && (       cutValue_W_min < W)
       && (                    y < cutValue_y_max )
       && ( cutValue_e_theta_min < e.Theta()*r2d  && e.Theta()*r2d  < cutValue_e_theta_max  )
       && (cutValue_pi_theta_min < pi.Theta()*r2d && pi.Theta()*r2d < cutValue_pi_theta_max )
       && (     cutValue_Ppi_min < pi.P()         &&         pi.P() < cutValue_Ppi_max      )
       && (     cutValue_Zpi_min < Zpi            &&            Zpi < cutValue_Zpi_max      )
       ) {
        if (fdebug>3) { std::cout << "succesfully passed (e,e'pi) kinematical cuts" << std::endl; }
        return true;
    }
    return false;
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
Double_t TaggedSIDISCuts::Chi2PID_pion_lowerBound( Double_t p, Double_t C){
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
Double_t TaggedSIDISCuts::Chi2PID_pion_upperBound( Double_t p, Double_t C){
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