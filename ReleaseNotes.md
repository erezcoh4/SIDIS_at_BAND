

# SRC-SIDIS revisions and release-notes 
-----------------------------------------------------------------

Mar-24, 2023 (commit_cfbc431)   
-------------

1. Moved to using only pions from the forward detector 
    The reason is that the central detector is not calibrated for our data-sets.
    The way that this was done is by updating *ExtractPipsInformation* and *ExtractPimsInformation* in **SIDISc12rSkimmer.C**
      
    if( pipluses[pipsIdx]->getRegion() != FD ) return;
    
2. Added a branch in TTree called “pi_region”, in order that in the event-selection cuts notebook, we show only events with Npi>0 where pi is from the FD
      
    where FD = 2000, see in Clas12Banks/clas12defs.h:  static const short FD = 2000;      
      
3. Updated Python analysis to step over event selection cuts one-by-one in **EventSelectionCuts.ipynb**




Mar-3, 2023 (commit_cfbc431)   
-------------

1. Added simulation to **SIDISc12rSkimmer.C** so that simulation will undergo the exact same analysis code as data


Feb-13, 2023    
-------------

1. Updated **Read_PiAcceptance_GEMCimulations.C** to read "white" GEMC simulations with the updated variables

2. Added *e_e_pi_GEMC* data dictionary to analysis results - "white" spectrum of (e,e'π) events, 
and the *do_e_e_pi_GEMC* option in **load_SIDIS_data()**


Feb-9, 2023    
-------------

1. Updated minimal Mx cut to be Mx > 1.7 in **apply_Mx_cut( df_dict=None, Mx_min = 1.7, Mx_max = 5 )**

2. Added *subdirname* to **load_SIDIS_ratio()**

3. Added *pT* and *phi* limits to **extract_SIDIS_ratio()** and to **compute_ratio_pips_to_pims()**, in order to study the dependence of the results on these variables and check hidden artificial Q2 dependence, where *phi* is the angle in the reaction frame (i.e. virtual-photon frame)

Feb-3, 2023    
-------------

1. Updated python analysis and cross-section extraction to include 10.4 and 10.6 GeV data
    **load_SIDIS_data()** 
        **load_SIDIS_data()**


Jan-31, 2023 (commit_2e77f15)    
-------------

1. Expanded the analysis to include, in addition to 10.2 GeV data, also 10.4 GeV and 10.6 data
        good_runs_10-2-final.txt [/cache/clas12/rg-b/production/recon/spring2019/torus-1/pass1/v0/dst/train/sidisdvcs/]
        good_runs_10-4.txt [/cache/clas12/rg-b/production/recon/spring2020/torus-1/pass1/v1/dst/train/sidisdvcs/]
        good_runs_10-6.txt [/cache/clas12/rg-b/production/recon/spring2019/torus-1/pass1/v0/dst/train/sidisdvcs/]
    
    (Good runs are taken from BAND analysis bandsoft_tools
    **bandsoft_tools/runlists**)
    

2. Added shell scripts to get the 10.4 and 10.6 GeV data from the cache **shell_scripts/jcache_get_files_sidisdvcs_10-4.csh** 
    **shell_scripts/jcache_get_files_sidisdvcs_10-6.csh** 


3. Corrected run-number TString format in SIDISatBAND_auxiliary::GetRunNumberSTR( int RunNumber )

*old version: sprintf( RunNumberStr, "00%d", RunNumber );*
*corrected version: sprintf( RunNumberStr, "%06d", RunNumber );*
    


Jan-26, 2023    
-------------

1. Updated ReadBeamCharge.C to read beam charge from RGA runs same as RGB runs

2. Added a shell script for RGA file extraction

Jan-17, 2023    
-------------

1. Updated python analysis code to make **p(e,e'π)** be analyzed similary to **d(e,e'π)** since we now extract the data from the farm in the same way and apply the same cuts on the data. 
    The only difference between **p(e,e'π)** and **d(e,e'π)** now is the charge-weights assigned on a run by run basis, that we need to complete for the free-proton data


Jan-12, 2023 (commit_ad849e5)    
-------------

1. Started analysis of free-p data *p(e,e'π)* from RGA with the same code for deuteron data *d(e,e'π)* from RGB

The reason for this change is that the small data files Igor produced were not enough and he is too busy to further produce more files from RGA
The RGA free-p data are taken from **nSIDIS** train located in:
    
    /cache/clas12/rg-a/production/recon/spring2019/torus-1/pass1/v1/dst/train/nSidis

For RGA use nSidis, they key difference is sidisdvcs has e_p > 1 GeV and nSidis has e_p > 2 GeV. [Harut A., Dec-2022] 
This is done by the exact same code as BAND skimming **SIDISc12rSkimmer.C**




Nov-24, 2022 (commit_5b2a758)    
-------------
1. Added *alpha_pi* and *zeta_pi* to d(e,e'πn) data in **MergeSIDISandBANDSkimmers.C**

    // pion light cone fraction in the lab-frame
    alpha_pi     = (pi_qFrame.E() - pi_qFrame.Pz()) / (Md/2);
    // LC momentum fraction of the active proton carried by the produced pion
    alpha_n      = (Pn_qFrame.E() - Pn_qFrame.Pz()) / (Md/2);
    zeta_pi      = alpha_pi / ( 2 - alpha_n );

2. Added a possibility to bin the tagged cross-section ratio as a function of *zeta_pi* in **Tagged_SIDIS_ratio.ipynb**
    (e,e π+ )
    3820 events after original cut (100.0 %)
    1733 events after p-theta cut (45.4 %)
    1680 events after Kinematical cut (44.0 %)
    (e,e π- )
    1231 events after original cut (100.0 %)
    811 events after p-theta cut (65.9 %)
    793 events after Kinematical cut (64.4 %) 




Nov-1, 2022    
-------------
1. Added *prefix* and *subdirname* to d(e,e'πn) data loading in **load_SIDIS_data()**



Oct-27, 2022  (sidisdvcs_27Oct2022_commit_2fe215f)   
-------------
1. Added qStar calculation also to tagged data in **MergeSIDISandBANDSkimmers.C**


Oct-18, 2022  (sidisdvcs_v18Oct2022_commit_574cf9b)   
-------------
1. Added qStar calculation to untagged data to **SIDISc12rSkimmer.C**, based on Natalie W. code, and output in SIDISc12rSkimmer, by copying *calcQStar()* into **SIDISatBAND_auxiliary.C**

2. Corrected a bug in the cache shell script to extract files using jcache

Oct-11, 2022    
-------------
1. Added option to limit *Q2* values in **extract_SIDIS_ratio()** and **compute_ratio_pips_to_pims()**


Sep-30, 2022    
-------------
1. Added uncertainties to *N(π)* in final results csv files in **extract_SIDIS_ratio()**


Sep-23, 2022    
-------------
1. Found and corrected a bug in **apply_p_theta_acceptance_cut_single_set()** that caused results to be unstable, see slides **BugFix_AcceptanceMatchingCut**



Sep-13, 2022    
-------------
1. Added weight_option to **extract_SIDIS_ratio()** with a default value to weight by 'beam-charge', in which case the number of events is measured in 1/nC

2. Added a column **'weight'** to **e_e_pi_pass_cuts** which is 1/beam-charge

Sep-6, 2022    
-------------
1. Added **W_d** back to output variables in the csv file

2. Updated **apply_Kinematical_cuts()** such that cut values are read off from *BANDcutValues.csv* file 


Aug-26, 2022    
-------------
1. Stabilized **ReadIgorRGAFile.C** to add free proton data from RGA to our analysis chain

2. Corrected the definition of **Z_LC** based on a conversation with Mark Strikman 


Aug-19, 2022    
-------------
1. Added a mechanism to reproduce (e,e'π+)/(e,e'π-) ratios in different M_x bins in **extract_SIDIS_ratio()**

2. Changed M_x cut from M_x>2.5 GeV which was wrong when we used a standing deuteron in M_x definition to M_x>1.3


Aug-4, 2022    
-------------
1. Replaced *ComputeKinematics* in SIDISc12rSkimmer by *ComputeElecctronKinematics* 

2. Added a computation of the following variables

Z_LC = (Ppi_q.E() + Ppi_q.Pz()) / (q.E() + q.P()); // z on the light-cone

3. Replaced the cut on W from 
   W_d = sqrt((d_rest + q).Mag2()) > 2.5 GeV/c2
   to
   W_p = sqrt((p_res + q).Mag2()) > 2.0 GeV/c2
   in *apply_Kinematical_cuts* and *SIDISc12rSkimmer*

4. Moved several generic routines from *SIDISc12rSkimmer* and *MergeSIDISandBANDSkimmers* to *SIDISatBAND_auxiliary*

5. Revisied all kinematical variables, now left with only Z_LC which definition seems odd

6. Added a varying precision to each variable that is streamed out to the CSV  

 


July-28, 2022    
-------------
1. Updated ReadIgorRGAFile.C script

2. Changed *W* variable in *SIDISc12rSkimmer* to  
 W_standing_d = sqrt((standing_d + q).Mag2());
 W_standing_p = sqrt((standing_p + q).Mag2());
 
 
3. Rotated neutron to q-Frame

July-25, 2022    
-------------
1. changed NMAXPIONS from 20 to 5 (we do not need to account for events with too many pions as we are interested in pions with high-z and in such events all pions are soft) 

2. Fixed a bug that wrote piplus_qFrame_pT and piminus_qFrame_pT wrongly as 0 to the TTree in SIDISc12rSkimmer()

3. Corrected W and W' from quantities calculated on the proton
    W       = sqrt( Mp2 - Q2 + 2. * omega * Mp );
    WPrime  = sqrt( Mp2 - Q2 + 2. * omega * (Md - Es) + 2. * Ps * sqrt(Q2 + w2) * cos( theta_sq ) );
to ones calculated on the deuteron
    W       = sqrt((P+q)2)   = sqrt( (target + q).Mag2()) 
    WPrime  = sqrt((p_i+q)2) = sqrt( (-p_n + q).Mag2() ) 


July-21, 2022    
-------------
1. Updated kinematical cut on W>2.5 instead of W>2 to match untagged data to the tagged one in *BANDcutValues.csv* and in the python method *apply_Kinematical_cuts()*
    and corrected to read *BANDcutValues.csv* instead of the wrong file *cutValues.csv* 

2. Replaced "inclusive" files with "sidisdvcs" files as the original "raw" data files for this analysis
   
    Until today we used the inclusive train files for RGB SIDIS analysis, from the directory.
    */volatile/clas12/rg-b/production/recon/spring2019/torus-1/pass1/v0/dst/train_20200610/inc/*
    On July-21, 2022 the inclusive train files that we were using for RGB SIDIS analysis were deleted. 
    Florian:  volatile is just storage which from time to time clears old files which are not used.
    The file system has no backup so the files are lost for now.
    Also I believe these inclusive train files were only generated in the beginning of RGB without any backup.
    The usual train files are at 
    */cache/clas12/rg-b/production/recon/spring2019/torus-1/pass1/v0/dst/train*
    These have a sidisdvcs skim
    However, we need first to get all files from long term storage mss back to cache via "jcache" command
    
    The files are different. 
    From cooking logbook https://github.com/zhaozhiwen/clas12_cooking_rgb/
    GitHub - zhaozhiwen/clas12_cooking_rgb


    sidisdvcs train:
    forward: 11:X+,X-:Xn electron + any particle in FD
    targetPDG: 2112 (not sure what this does to data filtering, probably only for cuts below)
    electron: Q2>0.95 && W>1.95 &&  p>1 && vz>-25 && vz<20 extra electron kinematic cut

    inclusive train:
    forward: 11:X+,X-:Xn electron + any particle in FD
    tagger: X+:X-:Xn any particle forward tagger
    central: X+:X-:Xn any particle in CD
    OR
    forward: -11:X+,X-:Xn positron + any particle in FD
    tagger: X+:X-:Xn any particle forward tagger
    central: X+:X-:Xn any particle in CD

    sidisdvcs skim has no positron filter and extra electron kinematic cuts + targetPDG.
    
    

3. Added transverse (pT) and longitudinal (pL) pion momentum with respect to the virtual photon as outputs in *SIDISc12rSkimmer.C* 



July-19, 2022    
-------------
1. Started to code a script to read off RGA data-file produced by Igor: *ReadIgorRGAFile.C*



July-14, 2022    
-------------
1. Added cuts to match BAND skimming:
    2 < Q2 < 10 (GeV/c)2
    3 < p(e) < 10.6 GeV/c 
     
    both in SIDISc12rSkimmer.C and in event_election_tools.py

2. Added a flad to check if (e,e'πn) passed event selection cuts in StreamToCSVfile()    
    This was done since before this correction, (e,e'πn) lines were written to CSV file even if they did not pass the event selection cuts, and e.g. pion momentum was too low (for a pion which was not the leading one)
    This does not impact the final results of the cross-section ratio, which were binned in z and minimally at z>0.3, but still relevant for e.g. kinematical distribution plots 



July-12, 2022    
-------------
1. Added a flag to cancel acceptance matching cut in p-theta plane to apply_further_selection_cuts_to_data()
  
2. Splitted apply_p_theta_acceptance_cut() to subsets of up to 500k events. 
 The reason is that it failed to apply the cut for N > 1M events. 
 This is a bug; you should open a ticket on github. 
 1M is the monotonic index size cutoff for hash tables, 
 and the logic beyond that is borked
 See [https://stackoverflow.com/questions/33814223/strange-error-in-pandas-indexing-with-range-when-length-1-000-000]

  
July-7, 2022    
-------------
1. Added information about the neutron time of flight into MergeSIDISandBANDSkimmers()
    The csv file now has two more ouput columns:
    "n_E,n_ToF"
    
2. Added a flag of applying p(n)>275 MeV/c cut to apply_further_selection_cuts_to_data()
    by removing this cut from load_SIDIS_data() and adding it as an independent function apply_minPn_cut() 

3. Updated compute_ratio_pips_to_pims() to return also the number of events with pi+ and the one with pi-


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
        
    
    
    
    
    
    
