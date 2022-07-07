

# Semi Inclusive DIS (SIDIS) on a deuteron target with $pi^+$ and $pi^-$ tagging and a fast neutron recoiling from a SRC pair breakup 

This repository is responsible for 
(A) skimming for (e,e'$\pi^{\pm}$) events
(B) merging (e,e'$\pi^{\pm}$) events with fast neutron recoiling from BAND data
(C) adding pion information to an arbitrary event list    
    
    
    
## Revisions and release notes
------------------------------------------------------------------------------

July-7, 2022    
-------------
1. Added information about the neutron time of flight into MergeSIDISandBANDSkimmers()
    The csv file now has two more ouput columns:
    "n_E,n_ToF"


June-10, 2022    
-------------
1. Updated good run list 
from 
    *good_runs_10-2.txt*        (Florian H., Sep-9, 2021, 101 runs) 
to 
    *good_runs_10-2-final.txt*  (Efrain S., May-4, 2022,  86 runs )
    

2. Updated load_SIDIS_data() in PythonAnalysis to speed things up:
1 run (6420):
apply_further_selection_cuts_to_data() time: 1.48 sec
 
e_e_pi['piplus']:            167.2 MB  
e_e_pi_pass_cuts['piplus']:  106.5 MB  
e_e_pi['piminus']:           71.1  MB  
e_e_pi_pass_cuts['piminus']: 65.5  MB  

So, 100 runs may take only ~ 150 sec but will occupy about (167.2+106.5+71.1+65.5) x 100 = 41 GB, 
which is extremely large and impedes the ability of Python to handle the data

Reduced DataFrames occupies to about 1/5 of the sizes
reduced_e_e_pi['piplus']:    33.2  MB  
reduced_e_e_pi['piminus']:   14.1  MB  

which will result in a total memory usage of about 41 GB / 5 = 8 GB

The memory occupancy reduction was implemented in load_SIDIS_data() by adding these lines:
eepi = pd.read_csv(...
                    usecols=['runnum','evnum',...],
                    dtype={'runnum':int,'evnum': int...}) 


 *After change*
 1234060 events: 13.1 +/- sec (10 us/event)
 
 

May-4, 2022    
-------------
1. Added a script to read a specific variable e.g. \theta(pi^-) ReadSpecificVariable.C read_specific_variable macro

2. Updated good-run list from Efrain' email for 10.2 GeV: <good_runs_10-2-final.txt>
 
3. Updated GEMC for pion acceptance simulation, to use specific gcard and yaml cards for rg_spring2019
 

Apr-28, 2022    
-------------
1. Added a script to read beam charge ReadBeamCharge.C and read_beam_charge macro


Apr-7, 2022    
-------------
1. Added to SIDIS skimming results information of the pion direction with respect to q
    This is since the physics has no preference in azimuthal angle with respect to q, and so we want to use the data to apply acceptance correction (We've learned that at this point the MC is not ready for acceptance correction as the acceptance and efficiency don't reproduce the data)
    
    defined the "q-frame" as follows:
    z axis is defined by the q - parallel to q 
    x axis is defined by the e' - such that p(e') resides in the x-z plane
     

    The relevant variables that were added to the csv-file are \theta and \phi of the pion (and the electron) in the q-frame  
    
2. Added 
    *(Npips < NMAXPIONS) && (Npims < NMAXPIONS)* 
    condition to analyzing an event in SIDIS skimming, as event 131227479 (run 6420) included 22 pi+ and 3pi- in a manner that caused memory crash

Mar-26, 2022    
-------------
1. Added SIDIS kinematical cuts to uniform GEMC simulations, to verify that cuts do not alter $\phi$ distribution    


Mar-24, 2022    
-------------
1. Added M_X to the uniform GEMC simulation data, in order to apply acceptance correction using MC on which the same cuts were applied as the data (p-theta acceptance matching, and M_X )
    Where M_X = || Beam + target - e - pi  || is the (e,e'\pi) reaction missing mass
    
2. Added information about the pion and electron sector to tagged-SIDIS merging

3. Updated and cleaned Python analysis scheme, mostly in *acceptance_correction_tools*


Mar-10, 2022    
-------------
1. Added information about the pion and electron sector to SIDIS skimming 


Mar-1, 2022    
-------------
1. Added ROOT files with no cuts

2. Split simulation files to individual small files of 100k events

Feb-17, 2022    
-------------
1. Added a script to run GEMC simulation for pion acceptance maps production

2. Changed electron direction generation to uniform in acceptance map production, to avoid a "hole" in the acceptance map produced by artificial close-proximity between the electron and the pion   


Feb-10, 2022    
-------------
1. Added fiducial cuts a PID requirements to *Read_PiAcceptance_GEMCimulations.C*.
    Now, each "event" (simulated pion) has another two variables:
    pi_reconstructed , pi_passed_cuts  


Jan-18, 2022    
-------------
1. Commence working on pion acceptance corrections 

    A. In order to produce pion acceptance maps, a python script was compiled in:
    
    MC/Acceptance_Corrections/produce_LUND_file_for_GEMC.py
    
    This script generates multiple d(e,e'pi) events with the same electron, and a pion with momentum in which (p,cos\theta,\phi) are distributed uniformly
    The output of this script are LUND files for pi+ and pi-, that are processed on the ifarm with dedicated GEMC running commands
    
    B. Another script based on .C was compiled to read the GEMC result hipo files (after cooking), in 
    
        MC/Acceptance_Corrections/Read_PiAcceptance_GEMCimulations.C
        macros/read_gemc_pion_acceptance      


Nov-25, 2021    
-------------
1. Fixed a bug in writing a long event number to the SIDIS-BAND merged CSV file by using std::fixed

2. Fixed a bug in pion information streamline: 
Replaced *piplus.push_back(TLorentzVector(piplus_Px[pipsIdx], piplus_Py[pipsIdx], piplus_Pz[pipsIdx], piplus_E[pipsIdx]))*
By *piplus.at(pipsIdx) = TLorentzVector(piplus_Px[pipsIdx], piplus_Py[pipsIdx], piplus_Pz[pipsIdx],piplus_E[pipsIdx]) ;*
And similarly for Vpiplus, piminus, Vpiminus

3. Fixed a bug in merging event indices: The first rown in the CSV file was always evnum=0, which resulted from a bug in SIDISeventID assignment in *MergeSIDISandBANDevents()*

4. Fixed a bug that did not write z(π) correctly to the CSV file


Nov-23, 2021    
-------------
1. Moved to merging BAND events with SIDIS, such that these files are from neutron-BAND files of final skims, from */volatile/clas12/users/segarrae/BAND/v3.1/10.2/final/tagged*



Nov-11, 2021    
-------------
1. Expanded the definition of x' for a moving proton
    xPrime1 = Q2 / (2. x ((Md - Es) x omega + Pn_Vect x q->Vect() ));
    xPrime2 = Q2 / (W2prime - Mp2 + Q2);
    


Oct-19, 2021    
-------------
1. Added a python script that runs *MergeSIDISandBANDSkimmers.C*

2. Moved **StreamToCSVfile** to Auxiliary

3. Added a seperate csv file for event-data, that stores relevant information per event: The number of positive and negative pions, event number, Xb etc.

Oct-6, 2021    
-------------
1. Added a clas12root script to add pion information to selected list of events *MergeSIDISandBANDSkimmers.C* 

Oct-3, 2021    
-------------
1. removed Vn variable from "merging" script, as there is no information about neutron track and thus none about the neutron vertex

2. replaced std::vector<TLorentzVector> variable "piplus" and "piminus" in output TTree with float array branches, as it just cauases too much hustle in ROOT structure

3. Added requirement for "goodneutron", and "eepiPastCutsInEvent" to speed up CreateListOfEventsToMerge() by reducing the number of events we need to go through in the merging stage  


 
Sep-30, 2021    
-------------
1. Added float array branches to sidis skimming TTree, that will allow easier reading in the sidis-neutron merging stage
    piminus_Px[NMAXPIONS]
    piminus_Py[NMAXPIONS]
    piminus_Pz[NMAXPIONS]
    piminus_E[NMAXPIONS] 
    Vpiminus_X[NMAXPIONS]
    Vpiminus_Y[NMAXPIONS]
    Vpiminus_Z[NMAXPIONS]


Sep-22, 2021    
-------------
1. Continue polishing neutron and SIDIS merging
    1. Added a Python scipt to control neutron and SIDIS merging
    2. Updated neutron TTree name to "calib", for "ncalibration_newclass" skimmer     
    
2. Updated auxilary BAND classes from *bandsoft_tools* commit "cfb8879..3f0c9c8"

3. Added "y" variable to SIDIS TTree

4. Updated the name of the python script used for plotting the kinematical distributions from SIDIS data to *SIDISKinematicalDistributions.ipynb*

 


Sep-19, 2021    
-------------
1. Refreshed SIDIS and BAND Skimmers merging:
    1. Update MergeSIDISandBANDSkimmers.C script
        1. Deleted all unnecessary and masive detector features like DC position etc. for electrons and pions
        2. Updated pion momentum to an array of size [NMAXPIONS] (similarly for pion vertex, Z,...)
        
2. added eepipsPastKinematicalCuts and eepimsPastKinematicalCuts variables to SIDISc12rSkimmer.C output TTrees

3. updated python script file names and added a script for neutron skimming:
skim_sidis_multiple_runs, skim_neutron_multiple_runs

Sep-12, 2021    
-------------
1. Added a python script that merges CSV dataframes from multiple runs 
    "combine_sidis_skim_multiple_runs"
    

Sep-9, 2021    
-------------
1. Added a counter for the number of deuterons in the event
2. Added number of particles from main species to output CSV file:
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



# Merge SIDIS and BAND Skimmers

    This script assumes that the input trees from SIDISc12rSkimmer.C,
    have a boolean flag "eepipsPastCutsInEvent/O" that states if event passed (e,e'\pi) event selection criteria


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
    root -l -q 'MergeSIDISandBANDSkimmers.C(6420,"pi+",3,3)'
    root -l -q 'MergeSIDISandBANDSkimmers.C(6420,"pi+",10,3,1000,"/Users/erezcohen/Desktop/data/BAND")'
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
        
 
 
 ## Using BAND skimmer

 ./code [outputFile] [MC/DATA] [time sifts] [inputFile]

 /u/home/cohen/BAND_analysis/bandsoft_tools/bin/neutrons  /volatile/clas12/users/ecohen/BAND/neutron_skimming/skimmed_neutrons_inc_006420.root 1 1 /volatile/clas12/rg-b/production/recon/spring2019/torus-1/pass1/v0/dst/train_20200610/inc/inc_006420.hipo

 (For SIDIS train /u/home/cohen/BAND_analysis/bandsoft_tools/bin/neutrons  /volatile/clas12/users/ecohen/BAND/BAND_skimmed_neutrons_sidisdvcs_run6420.root DATA 1 /cache/clas12/rg-b/production/recon/spring2019/torus-1/pass1/v0/dst/train/sidisdvcs/sidisdvcs_006420.hipo)
 
 
 comments:
 -----------------
1. DVCS train recommended by Florian for SIDIS skimming:
	"/cache/clas12/rg-b/production/recon/spring2019/torus-1/pass1/v0/dst/train/sidisdvcs/"
2. Tyler: There is a semi-exclusive filter to the DVCS train, so use RGM data trains:
	"/volatile/clas12/rg-b/production/recon/spring2019/torus-1/pass1/v0/dst/train_20200610/"
3. The following objects are taken from *bandsoft_tools* for Auxiliary:

genpart/genpart.cpp
hipolib/dictionary.h
hipolib/dictionary.cpp
include/constants.h
clashit/clashit.cpp
hipolib/bank.h
hipolib/bank.cpp
include/bandhit.h
bandhit/bandhit.cpp
include/BScintillator.h
banklib/BScintillator.cpp
include/BScaler.h
banklib/BScaler.cpp
include/BParticle.h
banklib/BParticle.cpp
include/BEvent.h
banklib/BEvent.cpp
include/BConfig.h
banklib/BConfig.cpp
include/BCalorimeter.h
banklib/BCalorimeter.cpp
include/BBand.h
banklib/BBand.cpp



 [Originally based on Ex1_CLAS12Reader.C]
