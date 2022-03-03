#include <TString.h>

#include "TaggedSIDISCuts.C"

void loadCutValues      (TString cutValuesFilename, NPiEvent event, int fdebug);
void StreamToCSVfile    (TString pionCharge, std::vector<Double_t> observables, bool passed_cuts_e_pi, 
                            bool passed_cuts_e_pi_kinematics, int fdebug);
void printCutValues     ();
double FindCutValue     ( std::string cutName );
void OpenResultFiles    ( TString outfilepath, TString outfilename );

class SIDISIO {
public:
    void OpenOutputFiles (TString outfilename,TString header);
    void CloseOutputFiles (TString OutDataPath, TString outfilename);

private:
    // Output root file and tree
    TFile * outFile_e_piplus, * outFile_e_piminus;
    TTree * outTree_e_piplus, * outTree_e_piminus;
    // Output CSV file
    std::ofstream   CSVfile_e_piplus,  SelectedEventsCSVfile_e_piplus,  SelectedEventsCSVfile_e_piplus_kinematics;
    std::ofstream   CSVfile_e_piminus, SelectedEventsCSVfile_e_piminus, SelectedEventsCSVfile_e_piminus_kinematics;
};

void loadCutValues(TaggedSIDISCuts cuts, TString cutValuesFilename, int fdebug){
    // read cut values csv file
    csv_reader csvr;
    cutValues = csvr.read_csv("BANDcutValues.csv");
    if (fdebug>2) { printCutValues(); }
    
    // assign specific cut values - to speed things up
    // by avoiding recalling FindCutValue() on every event
    
    // Cut on z-vertex position:
    // based on RGA_Analysis_Overview_and_Procedures_Nov_4_2020-6245173-2020-12-09-v3.pdf
    // p. 71
    if (torusBending==-1){ // in-bending torus field
        // Spring 19 and Spring 2020 in-bending.
        cut.setCutValue_Vz_min(FindCutValue("Vz_e_min_inbending"));
        cut.setCutValue_Vz_max(FindCutValue("Vz_e_max_inbending"));
    } else if (torusBending==1){ // Out-bending torus field
        // Fall 2019 (without low-energy-run) was out-bending.
        cut.setCutValue_Vz_min(FindCutValue("Vz_e_min_outbending"));
        cut.setCutValue_Vz_max(FindCutValue("Vz_e_max_outbending"));
        
    } else {
        std::cout
        << "Un-identified torus bending "
        << torusBending
        << ", return" << std::endl;
        return;
    }
        
    cut.setCutValue_e_PCAL_W               (FindCutValue("e_PCAL_W_min"));
    cut.setCutValue_e_PCAL_V               (FindCutValue("e_PCAL_V_min"));
    cut.setCutValue_e_E_PCAL               (FindCutValue("e_E_PCAL_min"));
    cut.setCutValue_SamplingFraction_min   (FindCutValue("SamplingFraction_min"));
    cut.setCutValue_PCAL_ECIN_SF_min       (FindCutValue("PCAL_ECIN_SF_min"));
    cut.setCutValue_Ve_Vpi_dz_max          (FindCutValue("(Ve-Vpi)_z_max"));
    cut.setCutValue_Q2_min                 (FindCutValue("Q2_min"));
    cut.setCutValue_W_min                  (FindCutValue("W_min"));
    cut.setCutValue_y_max                  (FindCutValue("y_max"));
    cut.setCutValue_e_theta_min            (FindCutValue("e_theta_min"));
    cut.setCutValue_e_theta_max            (FindCutValue("e_theta_max"));
    cut.setCutValue_pi_theta_min           (FindCutValue("pi_theta_min"));
    cut.setCutValue_pi_theta_max           (FindCutValue("pi_theta_max"));
    cut.setCutValue_Ppi_min                (FindCutValue("Ppi_min"));
    cut.setCutValue_Ppi_max                (FindCutValue("Ppi_max"));
    cut.setCutValue_Zpi_min                (FindCutValue("Zpi_min"));
    cut.setCutValue_Zpi_max                (FindCutValue("Zpi_max"));
}

void SIDISIO::OpenOutputFiles (TString outfilename,TString header){
    
    // Create output tree
    outFile_e_piplus  = new TFile( outfilename + "_e_piplus.root"  ,"RECREATE");
    outTree_e_piplus  = new TTree( "tree" , "(e,e'pi+) event information");
    outFile_e_piminus = new TFile( outfilename + "_e_piminus.root" ,"RECREATE");
    outTree_e_piminus = new TTree( "tree" , "(e,e'pi-) event  information");
    
    // Create output csv files
    CSVfile_e_piplus.open( outfilename  + "_e_piplus.csv" );
    CSVfile_e_piplus << header << std::endl;
    CSVfile_e_piminus.open( outfilename + "_e_piminus.csv" );
    CSVfile_e_piminus << header << std::endl;
    
    SelectedEventsCSVfile_e_piplus.open( outfilename + "_e_piplus_selected_eepi.csv" );
    SelectedEventsCSVfile_e_piplus << header << std::endl;
    SelectedEventsCSVfile_e_piminus.open( outfilename + "_e_piminus_selected_eepi.csv" );
    SelectedEventsCSVfile_e_piminus << header << std::endl;

    SelectedEventsCSVfile_e_piplus_kinematics.open( outfilename + "_e_piplus_selected_eepi_kinematics.csv" );
    SelectedEventsCSVfile_e_piplus_kinematics << header << std::endl;
    SelectedEventsCSVfile_e_piminus_kinematics.open( outfilename + "_e_piminus_selected_eepi_kinematics.csv" );
    SelectedEventsCSVfile_e_piminus_kinematics << header << std::endl;
}

void SIDISIO::CloseOutputFiles (TString OutDataPath, TString outfilename){
    // close output CSV
    CSVfile_e_piplus                .close();
    SelectedEventsCSVfile_e_piplus  .close();
    SelectedEventsCSVfile_e_piplus_kinematics  .close();
    CSVfile_e_piminus               .close();
    SelectedEventsCSVfile_e_piminus .close();
    SelectedEventsCSVfile_e_piminus_kinematics .close();

    int Nentires_e_piplus  = outTree_e_piplus  -> GetEntries();
    int Nentires_e_piminus = outTree_e_piminus -> GetEntries();
    
    // close output ROOT
    outFile_e_piplus->cd();
    outTree_e_piplus->Write();
    outFile_e_piplus->Close();
    
    outFile_e_piminus->cd();
    outTree_e_piminus->Write();
    outFile_e_piminus->Close();
    
    
    std::cout
    << "Done processesing "  <<  Nevents_processed          << " events,"
    << std::endl
    << std::setprecision(3)
    << (float)Nevents_passed_e_cuts/Nevents_processed       << " events passed e cuts,"
    << std::endl
    << (float)Nevents_passed_pips_cuts/Nevents_processed    << " events passed pi+ cuts,"
    << std::endl
    << "\t" << (float)Nevents_passed_e_pips_cuts/Nevents_processed  << " passed (e,e'pi+) cuts,"
    << std::endl
    << "\t\t" << (float)Nevents_passed_e_pips_kinematics_cuts/Nevents_processed  << " also passed kinematical cuts,"
    << std::endl
    << (float)Nevents_passed_pims_cuts/Nevents_processed    << " events passed pi- cuts,"
    << std::endl
    << "\t" << (float)Nevents_passed_e_pims_cuts/Nevents_processed  << " passed (e,e'pi-) cuts,"
    << std::endl
    <<  "\t\t" << (float)Nevents_passed_e_pims_kinematics_cuts/Nevents_processed  << " also passed kinematical cuts,"
    << std::endl;
    
    
    
    std::cout << "output files ready in root/csv formats in " << std::endl
    << std::endl
    << "wrote "  << Nentires_e_piplus  << " to (e,e'pi+) root file, "
    << std::endl << outFile_e_piplus -> GetName()
    << std::endl << OutDataPath + outfilename + "_e_piplus_selected_*.csv"
    << std::endl
    << "and "    << Nentires_e_piminus << " to (e,e'pi-) root file. "
    << std::endl << outFile_e_piminus -> GetName()
    << std::endl << OutDataPath + outfilename + "_e_piminus_selected_*.csv"
    << std::endl;
}

void StreamToCSVfile (TString pionCharge, // "pi+" or "pi-"
                      std::vector<Double_t> observables,
                      bool passed_cuts_e_pi,
                      bool passed_cuts_e_pi_kinematics,
                      int fdebug){
    if (fdebug>1) {
        std::cout << "streaming to CSVfile" << std::endl;
    }
    // decide which file to write...
    if (pionCharge=="pi+") {
        for (auto v:observables) CSVfile_e_piplus << std::fixed << v << ",";
        CSVfile_e_piplus << std::endl;
        
        if (passed_cuts_e_pi) {
            for (auto v:observables) SelectedEventsCSVfile_e_piplus << std::fixed << v << ",";
            SelectedEventsCSVfile_e_piplus << std::endl;
            
            if (passed_cuts_e_pi_kinematics){
                for (auto v:observables) SelectedEventsCSVfile_e_piplus_kinematics << std::fixed << v << ",";
                SelectedEventsCSVfile_e_piplus_kinematics << std::endl;
            }
        }
        
    }
    else if (pionCharge=="pi-") {
        for (auto v:observables) CSVfile_e_piminus << v << ",";
        CSVfile_e_piminus << std::endl;
        
        if (passed_cuts_e_pi) {
            for (auto v:observables) SelectedEventsCSVfile_e_piminus << v << ",";
            SelectedEventsCSVfile_e_piminus << std::endl;
            
            if (passed_cuts_e_pi_kinematics){
                for (auto v:observables) SelectedEventsCSVfile_e_piminus_kinematics << v << ",";
                SelectedEventsCSVfile_e_piminus_kinematics << std::endl;
            }
        }
    }
    else {
        std::cout << "pion charge ill-defined in StreamToCSVfile(), returning" << std::endl;
        return;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void printCutValues(){
    std::cout << "Using cut values:" << std::endl;
    for (auto cut: cutValues) {
        std::cout << cut.first << ": " << cut.second << std::endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
double FindCutValue( std::string cutName ){
    for (auto cut: cutValues) {
        if (strcmp(cut.first.c_str(),cutName.c_str())==0){
            return cut.second;
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void OpenResultFiles( TString outfilepath, TString outfilename ){
    OpenOutputFiles( outfilepath + outfilename,
                    ( (TString)"status,runnum,evnum,beam_helicity,"
                     +(TString)"e_P,e_Theta,e_Phi,e_Vz,"
                     +(TString)"pi_P,pi_Theta,pi_Phi,pi_Vz,"
                     +(TString)"Q2,W,xB,Zpi,omega,"
                     +(TString)"xF,y,M_X,"
                     +(TString)"Npips,Npims,Nelectrons,Ngammas,Nprotons,Nneutrons,Ndeuterons,"));
    // output tree branches
    SetOutputTTrees();
}