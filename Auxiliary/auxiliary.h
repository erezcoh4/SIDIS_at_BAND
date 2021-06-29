//
//  auxiliary.h
//
//
//  Created by Erez Cohen on 13/06/2021.
//

#ifndef tmp_h
#define tmp_h

#include <stdio.h>

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
// declare methods
TVector3                GetParticleVertex (clas12::region_part_ptr rp);
void                     SetLorentzVector (TLorentzVector &p4, clas12::region_part_ptr rp);
void                         OpenOutputFiles (TString csvfilename, TString header);
void                        CloseOutputFiles ();
void                      StreamToCSVfile (std::vector<Double_t> observables, bool IsSelectedEvent, int fdebug);
void                      ChangeAxesFrame (TString FrameName="q(z) frame");
void                        MoveTo_qFrame ();
bool EventPassedElectronSelectionCriteria (Double_t e_PCAL_x, Double_t e_PCAL_y,
                                           Double_t e_PCAL_W,Double_t e_PCAL_V,
                                           Double_t E_PCAL_e,
                                           Double_t E_ECIN_e, Double_t E_ECOUT_e,
                                           TLorentzVector e,
                                           TVector3 Ve);
bool   EventPassedPiPlusSelectionCriteria (Double_t DC_x, Double_t DC_y,
                                           Double_t chi2PID, Double_t p, TVector3 Ve, TVector3 Vpiplus );
Double_t          Chi2PID_pips_lowerBound (Double_t p, Double_t C=0.88);
Double_t          Chi2PID_pips_upperBound (Double_t p, Double_t C=0.88);


#endif /* auxiliary_h */
