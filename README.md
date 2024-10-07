

# SRC-SIDIS in CLAS12 + BAND
 
## p(e,e'π±), d(e,e'π±), and d(e,e'π±n) when tagging and a fast neutron recoiling from a SRC pair breakup 

This repository is responsible for 
(A) skimming data from CLAS12 for (e,e'π±) events
(B) merging (e,e'π±n) events with fast neutron recoiling from BAND data
(C) adding pion information to an arbitrary event list    
    
    
## Code flow 

### Untagged d(e,e'π) data
[SIDISc12rSkimmer.C (C++/ROOT) on ifarm] -> [Download csv files to local machine] -> [Untagged_SIDIS_ratio.ipynb (Python)]      

### Tagged d(e,e'πn) data
[SIDISc12rSkimmer.C (C++/ROOT) on ifarm] -> [MergeSIDISandBANDSkimmers.C (C++/ROOT) on ifarm] -> [Download csv files to local machine] -> [Tagged_SIDIS_ratio.ipynb (Python)]      
    
### Free proton p(e,e'π) data
[SIDISc12rSkimmer.C (C++/ROOT) on ifarm] -> [Download csv files to local machine] -> [FreeP_SIDIS_ratio.ipynb (Python)]
    
    
    
# Goals
    
    identify and analyse the following events:
 
        $$d(e,e’π+)X  $$
        $$d(e,e’π+n)X $$
        $$d(e,e’π-)X  $$
        $$d(e,e’π-n)X $$
        $$p(e,e’π+)X  $$
        $$p(e,e’π-)X  $$


# Prequisits
---------------------------------------
### Load modules on ifarm:

module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles
module purge
module load sqlite/dev
module load clas12



# Data runs used for analysis
---------------------------------------
1. BAND-RGB:

    Good runs are taken from BAND analysis 
    **bandsoft_tools/runlists**:
    
    10.2 GeV, spring 2019, 86 runs 6421 - 6597     
    good_runs_10-2-final.txt    
[/cache/clas12/rg-b/production/recon/spring2019/torus-1/pass1/v0/dst/train/sidisdvcs/]
    
    10.4 GeV, spring 2020, 145 runs 11362 - 11571
    good_runs_10-4.txt
[/cache/clas12/rg-b/production/recon/spring2020/torus-1/pass1/v1/dst/train/sidisdvcs/]
    
    10.6 GeV, spring 2019, 99 runs 6164 - 6399
    good_runs_10-6.txt
[/cache/clas12/rg-b/production/recon/spring2019/torus-1/pass1/v0/dst/train/sidisdvcs/]    

2. RGA free-proton 

    rga_nsidis_runs_10-6.txt compiled from nSidis train 


# Final SIDIS results
---------------------------------------
pandas.DataFrame(
{["$x_B$"],["$\Delta x_B$"],
 ['$N(\pi_{+})$'],['$N(\pi_{-})$'],
 ['$R$'],['$\Delta R_{+}$'],['$\Delta R_{-}$'],
  ['$\Delta N(\pi_{+})$'],['$\Delta N(\pi_{+})$']})
 
 
# Links and references
---------------------------------------
[SIDIS_analysis_note_final-5721534-2021-01-11-v15.pdf, p.1-2]
[THayward_thesis.pdf p.8]
BAND analysis note [https://www.overleaf.com/project/6159d43f5128373fedfddb47] 

rapidity: [https://www.jlab.org/conferences/radiative2016/talks/wednesday/weiss.pdf]
 
 
# Kinematical observables
---------------------------------------
 
SIDIS kinematic variables

    // QE y-variable
    y         = omega / Ebeam;
    
    // Feynman x 
    xF  = 2. * (pi.Vect().Dot(q.Vect())) / (q.P() * W)
    
    Zpi = pi.E()/omega;
    Zpi_LC = (pi_qFrame.E() + pi_qFrame.Pz()) / (q->E() + q->P());
    
        
    // barion rapidity  η_br = log(Pπ+/Pπ−) / 2
    eta_pi     = 0.5 * log( (pi_qFrame.E()+pi_qFrame.Pz())
                            /
                            (pi_qFrame.E()-pi_qFrame.Pz()))


Spectator energy-momentum and angle

    Es      = Pn.E();
    Ps      = Pn.P();
    omega   = q->E();
    theta_sq= Pn.Angle( q->Vect() );

proton initial state (in the nucleus)
    
    Einit            = Md - Pn.E();
    Pinit            = ( 0 - Pn.Vect(), Md - Pn.E() );
    Pinit_qFrame     = ( RotateVectorTo_qFrame( Pinit.Vect() ) , M_d - Pn.E() ); 
    Pinit_virtuality = Pinit_qFrame.Mag2() - Mp2; 


W - hadronic invariant mass

    // for a standing proton
    W       = (p_rest + q).Mag();
    W_Prime = ( p_init + *q ).Mag();   


Mx - missing mass of the emerging hadron

    M_x       = ( p_rest + q - pi ).Mag();    
    M_x_Prime = ( p_init + q - pi ).Mag();
    
    
alpha - Light cone momentum fraction    

    
    alpha_pi     = (pi_qFrame.E() - pi_qFrame.Pz()) / (Md/2);    
    alpha_n      = (Pn_qFrame.E() - Pn_qFrame.Pz()) / (Md/2);
    zeta_pi      = alpha_pi / ( 2 - alpha_n );

    
    
    
    
    
    
    
    
    

        
# SIDIS skimmer
        
    ## Execute:
    ---------------------
    ./macros/skim_multiple_runs --Nruns=3 --NeventsMax=10
    clas12root -q "SIDISc12rSkimmer.C+(6420,100)"
    root -l -q "MergeSIDISandBANDSkimmers.C(6420,100)"
    
    
    ## output

    Output files are saved in
    /volatile/clas12/users/ecohen/BAND/


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

    
    
    
    
    
    
    // write a (e,e'pi) event-line to CSV file

    // from Harut A., Aug-2, 2021:
    //    1: status, 2: runnum, 3: evnum, 4: helicity, 5: e_p, 6: e_theta, 7: e_phi, 8: vz_e, 9: pi_p, 10: pi_theta, 11: pi_phi, 12: vz_pi, 13: P_p, 14: P_theta, 15: P_phi, 16: vz_P, 17: Q2, 18: W, 19: x, 20: y, 21: z_pi, 22: z_P, 23: Mx(e:pi:P:X), 24: Mx(e:pi:X), 25: Mx(e:P:X), 26: zeta, 27: Mh, 28: pT_pi, 29: pT_P, 30: pTpT, 31: xF_pi, 32: xF_P, 33: eta_pi, 34: eta_P, 35: Delta_eta, 36: phi_pi (gamma*N COM), 37: phi_P (gamma*N COM), 38: Delta_phi.
    //
    //    1: status is a number indicating the quality of the event with the non-0 number indicating something was not not good (ex. out of fiducial region, not within the final cuts on energies of particles, missing or invarian masses....) The final observables will be done using status==0, while sensitivity of the observable to different cuts could be studied for various values of status>0
    //    2-3: run number and event number to identify the event
    //    4: helicity of the electron +1 along the beam and -1 opposite to it
    //    5,6,7,8 electron momentum,theta,phi_Lab, and z-vertex
    //    9,10,11,12 the same for pi+
    //    All other columns could be calculated from the first 16, but are included for cross check and minimizing the work in production of final observables. Some of them simple, like x,Q^2,W,y,  z=E\pion/nu, zome less trivial like Breit frame rapidities of pion (eta_pi) and proton eta_P, or corresponding values for X_Feynman variable xF_pi, xF_P.
    //    Some, like azimuthal angles of pion (phi_pi) and proton (phi_p) in the CM frame may also be confusing.
    //    In addition to the first 16 columns (mandatory) you can add as many columns as you are comfortable to fill, and consider relevant for your process.




 [Originally based on Ex1_CLAS12Reader.C]
