

# Skimming for Semi Inclusive DIS (SIDIS) on a deuteron target with $pi^+$ and $pi^-$ tagging

    Main open ToDo items:
    (1) Validate pion DC fiducial cuts
    
    
    
    
    
    
# Revisions

Sep-9, 2021    
-------------
1. Added number of particles from main species to output CSV file:
    "Npips,Npims,Nelectrons,Ngammas,Nprotons,Nneutrons,",
The reason for this is that we observe too many events with z > 0.9, which are featured by MX = 2 (namely, exclusive events).
We want to understand and perhaps throw away these unwanted events  


Sep-2, 2021    
-------------
1. Added "PrintProgress" to python script
2. Updated "Nruns" option in python script, such that '-1' is for all runs 


Aug-31, 2021    
-------------
1. Debugged pion DC information (DC_x/y/z seem funny in ROOT TTree results)
    There were two problems, one could be solved and the other is a general problem with the hipo files that we disregard for now.
    
    A. Sometimes the readout-sector is 0. This is funny, and, of-course, wrong. We should throw away these events.
    Justin B. Estee (June-21): "I also had this issue. I am throwing away sector 0. The way you check is plot the (x,y) coordinates of the sector and you will not see any thing."
    
    B. TTrree branch writing was using "/F" for "double" variables instead of "/D". This was a bug that was subsequently corrected to "/D".
    
    
    
    
Aug-27, 2021    
-------------
1. Added pion fiducial requirements in the DC using DCfid_SIDIS::DC_fid_th_ph_sidis()
            
    

Aug-26, 2021    
-------------
1. corrected TTree variables 
    A. ones that were wrongly filled as singlets, e.g. e_DC_x, e_DC_y, and e_DC_z - were D, should be F[3] for 3 layers
    B. piplus/piminus that were changed to std::vector<> type


2. added variables to CSV files: 
    "xF"           Feynman x = 2(p_pi \cdot q)/|q|W              
    "y"             omega / Ebeam
    "M_X"       missing mass of the reaction = || beam + target - e' - pi ||^2 
    
    
3. added "selected_eepi" CSV file to include events that pass electron and pion selection
    A. fill ROOT TTree with "selected_eepi" events
    
4. changed "selected_events" CSV file to include only events that pass electron and pion selection and kinematical cuts,
following [SIDIS_analysis_note_final-5721534-2021-01-11-v15]
    
    - 0.3 < z < 1
    - 1 (GeV/c)2 < Q2
    - 2 GeV < W
    - 5˚ < \theta_e < 35˚
    - 5˚ < \theta_\pi < 35˚
    - y < 0.75 (minimal electron momentum ~ 2.65 GeV/c)
    -  1.25 GeV/c < p_pi < 5 GeV/c
 



Aug-22, 2021    
-------------
1. Added variables to the output ROOT TTree 
        “Npips”         number of positive pions per event
        “Npims”         number of negative pions per event        
        “Nprotons”      number of protons
        “Nneutrons”     number of neutrons
        “Ngammas”       number of gammas
        “Nelectrons”    number of electrons
    2. DC_layer was always over-written by the last layer. This bug was fixed.
        DC layes is an array of dimensions [3], which equals {6,18,36}. 
        Region 1 is denoted at DC detector 6, Region 2 is denoted 18, Region 3 - as 36
    3. changed output TTree names from "(e,e'pi+) events" to "tree", which is easier to read-off
        
    
    
    
    
    
    
    
# Goals
    
    identify and analyse the following events:
 
        $$D(e,e’pi+)X  $$
        $$D(e,e’pi+n)X $$
        $$D(e,e’pi-)X  $$
        $$D(e,e’pi-n)X $$
        
        
    ## Execute:
    ---------------------
    ./macros/skim_multiple_runs --Nruns=3 --NeventsMax=10
    clas12root -q "SIDISc12rSkimmer.C+(6420,100)"
    root -l -q "MergeSIDISandBANDSkimmers.C(6420,100)"
    
    
    ## ToDo list:
    ---------------------    
    (1) add fiducial cuts for pions in the DC
    (2) Compare electron-only with BAND analysis results or RGA results (e.g. single-spin assymetries)
    (3) Verify that event matching between SIDIS and BAND skimming  is done correctly
    (4) Add neuton fiducial cuts
    (5) Verify we use the same values for the fiducial cuts as the general BAND group
    


    ## output

    Output files are saved in
    /volatile/clas12/users/ecohen/BAND/

        /volatile/clas12/users/ecohen/BAND/SIDIS_skimming/
        /volatile/clas12/users/ecohen/BAND/merged_SIDIS_and_BAND_skimming/



    ## MergeSIDISandBANDSkimmers.C
    ---------------------------------------------
    * merger of Semi-Inclusive DIS skimming and BAND-neutron skimming *





 ---------------------------------------------
 void MergeSIDISandBANDSkimmers (   int RunNumber=6420,
                                    int NeventsToMerge=-1, // all possible events
                                    int fdebug=1,
                                    int PrintProgress=5000 )
  ---------------------------------------------
Merge SIDIS skimming (based on SIDISc12rSkimmer.C) with BAND skimming (based on band-soft tools)


### Execute:
-----------------

    root -l MergeSIDISandBANDSkimmers.C
    root -l -q "MergeSIDISandBANDSkimmers.C(6420,-1)"










    ## SIDISc12rSkimmer.C
---------------------------------------------
* hipo skimmer for Semi-Inclusive DIS at CLAS12 with BAND *

 ---------------------------------------------
 void SIDISc12rSkimmer (   int RunNumber=6420,
                            int NeventsMax=-1, // -1 is "process all events"
                            int fdebug=1,
                            int PrintProgress=5000,
                            TString DataPath="/cache/clas12/rg-b/production/recon/spring2019/torus-1/pass1/v0/dst/train/sidisdvcs/")
  ---------------------------------------------
Read BAND data and write important observables to csv file
 
 
execute:
-----------------

   clas12root -q "SIDISc12rSkimmer.C+(6420,100)"
   clas12root SIDISc12rSkimmer.C
   


 input:
 -----------------
 RunNumber              Run number, e.g. 6420 (see "BAND good-runs file list")
 NeventsMax             Maximal number of events to process. (-1) is all events. Default is 100 events.
 fdebug                 Verbosity
 doApplySelectionCuts   If this flag is set to TRUE, events that pass the event selection cuts are recorded to a seperate csv file
 PrintProgress          Print event N events (default=5000)
 NpipsMin               Minimal number of pi+ required for defining an event
 DataPath
    examples:
        "/lustre19/expphy/cache/clas12/rg-b/production/recon/fall2019/torus+1/pass1/v1/dst/train/edeutcut/inc/"
        (inputFile = "edeutcut_" + RunNumberStr + ".hipo")
 
        "/cache/clas12/rg-b/production/recon/spring2019/torus-1/pass1/v0/dst/train/sidisdvcs/"
        (inputFile = "sidisdvcs_" + RunNumberStr + ".hipo")
        
 
 
 
 comments:
 -----------------
 DVCS train recommended by Florian for SIDIS skimming:
	"/cache/clas12/rg-b/production/recon/spring2019/torus-1/pass1/v0/dst/train/sidisdvcs/"
 Tyler: There is a semi-exclusive filter to the DVCS train, so use RGM data trains:
	"/volatile/clas12/rg-b/production/recon/spring2019/torus-1/pass1/v0/dst/train_20200610/"


 [Originally based on Ex1_CLAS12Reader.C]
