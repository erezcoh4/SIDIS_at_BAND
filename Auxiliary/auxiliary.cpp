

#include "auxiliary.h"



// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
bool auxiliary::EventPassedElectronSelectionCriteria(Double_t e_PCAL_x, Double_t e_PCAL_y,
                                          Double_t e_PCAL_W,Double_t e_PCAL_V,
                                          Double_t E_PCAL_e,
                                          Double_t E_ECIN_e, Double_t E_ECOUT_e,
                                          TLorentzVector e,
                                          TVector3 Ve){
    
    // decide if electron in event passes event selection cuts
    
    // Cut on z-vertex position:
    // based on RGA_Analysis_Overview_and_Procedures_Nov_4_2020-6245173-2020-12-09-v3.pdf
    // p. 71
    Double_t Vz_min,Vz_max;
    if (torusBending==-1){
        
        // in-bending torus field -13.0 cm < Vz < +12.0 cm
        // Spring 19 and Spring 2020 in-bending.
        Vz_min = -13.0;
        Vz_max = 12.0;
        
    } else if (torusBending==-1){
        
        // Out-bending torus field -18.0 cm < Vz < +10.0 cm
        // Fall 2019 (without low-energy-run) was out-bending.
        Vz_min = -18.0;
        Vz_max = 10.0;
        
    } else {
        
        std::cout << "Un-identified torus bending, return";
        return false;
    }
    
    if(
       // fiducial cuts
       // fiducial cuts on PCAL
       fabs(e_PCAL_x)>0
       &&  fabs(e_PCAL_y)>0
       &&  e_PCAL_W > 19
       &&  e_PCAL_V > 19
       
       // fiducial cuts on DC
       
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
bool auxiliary::EventPassedPiPlusSelectionCriteria( Double_t DC_x, Double_t DC_y,
                                         Double_t chi2PID, Double_t p, TVector3 Ve, TVector3 Vpiplus ){
    // decide if pi+ selection cuts
    //
    // input:
    // --------
    // DC_x, DC_y   pi+ drift-chamber coordinates
    // chi2PID      pi+ chi2PID     (chi2PID_pips)
    // p            pi+ momentum    (Ppips.P())
    //
    
    if(
       // fiducial cuts on DC
       (DC_x>0 && DC_y>0 )
       
       // pi+ Identification Refinement - chi2PID vs. momentum
       && ( Chi2PID_pips_lowerBound( p ) < chi2PID && chi2PID < Chi2PID_pips_upperBound( p ) )
       
       // Cut on the z-Vertex Difference Between Electrons and Hadrons.
       &&  ( fabs((Ve-Vpiplus).Z()) < 20.0)
       ) return true;
    
    return false;
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
Double_t auxiliary::Chi2PID_pips_lowerBound( Double_t p, Double_t C){
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
Double_t auxiliary::Chi2PID_pips_upperBound( Double_t p, Double_t C){
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
TVector3 auxiliary::GetParticleVertex(clas12::region_part_ptr rp){
    TVector3 V(rp->par()->getVx(),
               rp->par()->getVy(),
               rp->par()->getVz());
    return V;
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void auxiliary::SetLorentzVector (TLorentzVector &p4,clas12::region_part_ptr rp){
    p4.SetXYZM(rp->par()->getPx(),
               rp->par()->getPy(),
               rp->par()->getPz(),
               p4.M());
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void auxiliary::OpenCSVfiles (TString csvfilename,TString header){
    CSVfile.open( csvfilename + ".csv" );
    CSVfile << header << std::endl;
    
    SelectedEventsCSVfile.open( csvfilename + "_selected_events.csv" );
    SelectedEventsCSVfile << header << std::endl;
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void auxiliary::CloseCSVfiles (){
    CSVfile.close();
    SelectedEventsCSVfile.close();
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void auxiliary::StreamToCSVfile (std::vector<Double_t> observables, bool IsSelectedEvent, int fdebug){
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
void auxiliary::ChangeAxesFrame(TString FrameName){
    if (FrameName == "q(z) frame")
        MoveTo_qFrame();
    return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void auxiliary::MoveTo_qFrame(TLorentzVector q){
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

    neutron_qFrame.RotateZ(-q_phi);
    neutron_qFrame.RotateY(-q_theta);
    Double_t neutron_phi = neutron_qFrame.Phi();
    neutron_qFrame.RotateZ(-neutron_phi);

}

