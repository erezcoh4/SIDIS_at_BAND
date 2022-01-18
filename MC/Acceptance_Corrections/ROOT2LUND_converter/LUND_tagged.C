#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TSystem.h"

#include <fstream>
#include <iostream>

using namespace std;


int LUND_tagged(TString inFileName, TString outFileName) {

	if(gSystem->AccessPathName(inFileName)){
        	cout << inFileName << " does not exist...skipping" << endl;
		return 0;
    	} 

	cout << "Converting file " << inFileName << endl;
	cout << "Saving as " << outFileName << endl;

	TFile* inFile = new TFile(inFileName);
	TTree* T = (TTree*)inFile->Get("T");

	const int nParticles = 2;
	TString particleType[nParticles] = {"e-", "n"};
        int particleID[nParticles] = {11, 2112};
        double particleMass[nParticles] = {0.511e-3, 0.93957};   // GeV

	int targA = 2;
	int targZ = 1;
	double targP = 0.; // polarization
	double beamP = 0.; // polarization
	int interactN = 1;
	int beamType = 11;
	double beamE = 10.2;	// GeV


	double pvec[nParticles][3];
	double weight, radweight, zvtx;
	T->SetBranchAddress("pe", pvec[0]);
	T->SetBranchAddress("pn", pvec[1]);
	T->SetBranchAddress("zvtx", &zvtx);
	T->SetBranchAddress("rad", &radweight);

	int nEvents = T->GetEntries();

	ofstream outfile;
	outfile.open(outFileName); 

	TString formatstring, outstring;

	for (int i = 0; i<nEvents; i++) {

		weight = 1.;

		T->GetEntry(i);

		weight = radweight;

		// LUND header for the event:
		formatstring = "%i \t %i \t %i \t %.3f \t %.3f \t %i \t %.1f \t %i \t %i \t %.3f \n";
		outstring = Form(formatstring, nParticles, targA, targZ, targP, beamP, beamType, beamE, interactN, i, weight);

		outfile << outstring; 

		double vx = 0.;
		double vy = 0.;
		double vz = zvtx;	

		// LUND info for each particle in the event
		for (int j = 0; j<nParticles; j++) {

			double px = pvec[j][0];
			double py = pvec[j][1];
			double pz = pvec[j][2];

			double E = sqrt(px*px + py*py + pz*pz + particleMass[j]*particleMass[j] );			

			formatstring = "%i \t %.3f \t %i \t %i \t %i \t %i \t %.5f \t %.5f \t %.5f \t %.5f \t %.5f \t %.5f \t %.5f \t %.5f \n";
			outstring = Form(formatstring, j+1, 0.0, 1, particleID[j], 0, 0, px, py, pz, E, particleMass[j], vx, vy, vz);  
		
			outfile << outstring;
		
		}	

	}

	outfile.close();

	return 1;

}
