#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TSystem.h"

#include <fstream>
#include <iostream>

using namespace std;

int LUND_tagged(TString xsv, TString xpv, int fileNum);

void LUND_tagged_loop(TString fXSV, TString fXPV) {

	int nFiles = 1000;

	for(int i = 1; i < nFiles+1; i++) {
		LUND_tagged(fXSV, fXPV, i);
	}

}


int LUND_tagged(TString xsv, TString xpv, int fileNum) {

	TString genPath = Form("/volatile/clas12/users/ecohen/GEMC/generator/10.2/%s/%s",xsv.Data(),xpv.Data());
	TString inFileName = Form("%s/band_%s_%s_%05i.root", genPath.Data(), xsv.Data(), xpv.Data(), fileNum);

	if(gSystem->AccessPathName(inFileName)){
        	cout << inFileName << " does not exist...skipping" << endl;
		return 0;
    	} 

	cout << "Converting file " << inFileName << endl;
	TFile* inFile = new TFile(inFileName);
	TString lundPath = Form("/volatile/clas12/users/tkutz/GEMC/LUND/10.2/%s/%s", xsv.Data(), xpv.Data());
	TString outfilename = Form("%s/lund_%s_%s_%05i.dat", lundPath.Data(), xsv.Data(), xpv.Data(), fileNum);
	cout << "Saving as " << outfilename << endl;

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
	outfile.open(outfilename); 

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
