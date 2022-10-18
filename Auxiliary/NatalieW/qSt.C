#include "TMath.h"
#include "TVector3.h"


double calcQStar(TVector3 eP3, TVector3 piP3, double Ebeam){

  double tan_phi_acop = piP3.Y()/piP3.X(); // EIC Frame
  double eta_pi = TMath::ATanH(-piP3.Z()/piP3.Mag()); // - sign from EIC frame
  double eta_e = TMath::ATanH(-eP3.Z()/eP3.Mag());
  double delta_eta = eta_pi - eta_e;

  double qstar = 2*Ebeam*TMath::Exp(eta_pi)*tan_phi_acop/(1 + TMath::Exp(delta_eta));

  return qstar;

}

void qSt(){
  TVector3 el(.4, .0003, 2);
  TVector3 pi(1,.4,.992);

  cout << "q* = " << calcQStar(el, pi, 4) << endl;

}