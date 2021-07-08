// last edit July-5, 2021 (EOC, mbp), see README




// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
// Globals
// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.

TString              DataPath = "/volatile/clas12/users/ecohen/BAND/";

// Input root files and trees
TFile * SIDISFile, * BANDFile;
TTree * SIDISTree, * BANDTree;

// Output root file and tree
TFile * outFile;
TTree * outTree;




void          OpenInputFiles (TString RunStr);
void         OpenOutputFiles (TString RunStr, TString header);

void         CloseInputFiles ();
void        CloseOutputFiles ();
void MergeSIDISandBANDevents (int NeventsToMerge=100,
                              int fdebug=1,
                              int PrintProgress=5000


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
// Main functionality
// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void MergeSIDISandBANDSkimmers(int RunNumber=6420,
                               int NeventsToMerge=100,
                               int fdebug=1,
                               int PrintProgress=5000){
    
    char RunNumberStr[20];
    sprintf( RunNumberStr, "00%d", RunNumber );
    OpenInputFiles   ( (TString)RunNumberStr );
    OpenOutputFiles  ( (TString)RunNumberStr );
    
    MergeSIDISandBANDevents();
    
    
    CloseOutputFiles ();
    CloseInputFiles  ();
    
}



// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void MergeSIDISandBANDevents(int NeventsToMerge, int fdebug, int PrintProgress){
    
    Int_t   BANDrunID, BANDeventID, SIDISrunID, SIDISeventID;
    
    BANDTree  -> SetBranchAddress("eventnumber"  ,&BANDeventID);
    BANDTree  -> SetBranchAddress("runnum"       ,&BANDrunID);
    SIDISTree -> SetBranchAddress("eventnumber"  ,&SIDISeventID);
    SIDISTree -> SetBranchAddress("Runno"        ,&SIDISrunID);

    int NmergedEvents = 0;
    for (int BANDevent=0; BANDevent < NeventsBAND ; BANDevent++){
        
        BANDTree -> GetEntry(BANDevent);
        
        
        for (int SIDISevent=0; SIDISevent < NeventsSIDIS ; SIDISevent++){
            
            SIDISTree -> GetEntry(SIDISevent);
            
            if (fdebug>1){
                std::cout
                << "BAND run "  << BANDrunID
                << ", event "   << BANDeventID
                << ", "
                << "SIDIS run " << SIDISrunID
                << ", event "   << SIDISeventID
            }
            
            if ( (BANDrunID == SIDISrunID) && (BANDeventID == SIDISeventID)){
                // merge...
                // record event
                NmergedEvents ++;
                if (NmergedEvents >= NeventsToMerge){
                    std::cout << "merged " << NmergedEvents << " events, breaking." << std::cout;
                    return;
                }
            }
            
        } // end SIDIS event loop
        
    } // end BAND event loop
    std::cout << "merged " << NmergedEvents << " SIDIS and BAND events." << std::cout;
}



// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void OpenInputFiles (TString RunStr){
    
    TString   skimmedBANDFilename = (DataPath + "neutron_skimming/"
                                     + "skimmed_SIDIS_inc_"  + RunStr + ".root");
    BANDFile                      = new TFile( skimmedBANDFilename );
    BANDTree                      = (TTree*)BANDFile->Get("neutrons");
    
    
    
    TString  skimmedSIDISFilename = (DataPath + "SIDIS_skimming/"+ RunStr + ".root");
    SIDISFile                     = new TFile( skimmedBANDFilename );
    SIDISTree                     = (TTree*)SIDISFile->Get("sidis");
    
    
}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void OpenOutputFiles (TString RunStr, TString header){
    
    
    TString skimmedMergedFilename = (DataPath + "merged_SIDIS_and_BAND_skimming/"
                                     + "skimmed_SIDIS_and_BAND_inc_"  + RunStr );
    
    
    
    // Create output tree
    outFile = new TFile( skimmedMergedFilename + ".root" ,"RECREATE");
    outTree = new TTree( "T" , "Event information from merged SIDIS and BAND skimmers");
    
    // Create output csv files
    CSVfile.open( skimmedMergedFilename + ".csv" );
    CSVfile << header << std::endl;
    
}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void CloseInputFiles (){
    SIDISFile->Close();
    BANDFile->Close();
}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void CloseOutputFiles (){
    
    // close output CSV
    CSVfile.close();
    
    // close output ROOT
    outFile->cd();
    outTree->Write();
    outFile->Close();
    
}

