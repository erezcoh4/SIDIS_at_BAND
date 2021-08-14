

# Skimming for Semi Inclusive DIS (SIDIS) on a deuteron target with $pi^+$ and $pi^-$ tagging

    last edit: Aug-11, 2021 (EOC, mbp)
    
    ## We are interested in events:
 
        $$D(e,e’pi+)X  $$
        $$D(e,e’pi+n)X $$
        $$D(e,e’pi-)X  $$
        $$D(e,e’pi-n)X $$
        
        
    ## Execute:
    ---------------------
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
                            bool doApplySelectionCuts=true,
                            int PrintProgress=5000,
                            int NpipsMin=1,
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


 ToDo list:
 
 Priority 1:
 Add pi+ fiducial cuts on DC
 
 
 

 
 variable description:
 ---------------------------------
 event              event number
 Ee                 electron energy                             [GeV]
 Pe                 electron 3-momentum                         [GeV/c]
 Pe_x               electron x-momentum                         [GeV/c]
 Pe_y               electron y-momentum                         [GeV/c]
 Pe_z               electron z-momentum                         [GeV/c]
 Ptheta_e           electron theta                              [rad]
 Pphi_e             electron phi                                [rad]
 Ngammas            number of identified photons in event
 Np                 number of identified protons in event
 Nn                 number of identified neutrons in event
 Npips              number of identified pi+ in event
 Npims              number of identified pi- in event
 omega              energy transfer q.E()                       [GeV]
 q                  3-momentum transfer                         [GeV/c]
 q_x                x-momentum transfer                         [GeV/c]
 q_y                y-momentum transfer                         [GeV/c]
 q_z                z-momentum transfer                         [GeV/c]
 xB                 Q2/(2*Mp*q.E())                             [GeV/c]
 Q2                 4-momentum transfer -q.Mag2()               [(GeV/c)^2]
 z                  hadron rest frame energy Epi/q.E()
 Epips              pi+ energy                                  [GeV]
 Ppips              pi+ 3-momentum                              [GeV/c]
 Ppips_x            pi+ x-momentum                              [GeV/c]
 Ppips_y            pi+ y-momentum                              [GeV/c]
 Ppips_z            pi+ z-momentum                              [GeV/c]
 En                 neutron energy                              [GeV]
 Pn                 neutron 3-momentum                          [GeV/c]
 Pn_x               neutron x-momentum                          [GeV/c]
 Pn_y               neutron y-momentum                          [GeV/c]
 Pn_z               neutron z-momentum                          [GeV/c]
 Ppips_t_q          pi+ momentum transverse to q direction      [GeV/c]
 Ppips_q            pi+ momentum parallel to q direction        [GeV/c]
 e_E_PCAL           electron energy deposit in PCAL             [GeV]
 e_E_ECIN           electron energy deposit in ECAL(in)         [GeV]
 e_E_ECOUT          electron energy deposit in ECAL(out)        [GeV]
 Vx_e               electron vertex position in x-direction     [cm]
 Vy_e               electron vertex position in y-direction     [cm]
 Vz_e               electron vertex position in z-direction     [cm]
 Vx_pips            pi+ vertex position in x-direction          [cm]
 Vy_pips            pi+ vertex position in y-direction          [cm]
 Vz_pips            pi+ vertex position in z-direction          [cm]
 Vx_n               neutron vertex position in x-direction      [cm]
 Vy_n               neutron vertex position in y-direction      [cm]
 Vz_n               neutron vertex position in z-direction      [cm]
 pips_chi2PID       chi2 of PID for pi+
 chi2PID_n          chi2 of PID for n
 e_PCAL_W           W of the electron in PCAL
 e_PCAL_V           V of the electron in PCAL
 pips_PCAL_W        W of the pi+ in PCAL
 pips_PCAL_V        V of the pi+ in PCAL
 e_PCAL_x           x-hit position of the electron in PCAL
 e_PCAL_y           y-hit position of the electron in PCAL
 e_PCAL_z           z-hit position of the electron in PCAL
 e_PCAL_sector      sector
 e_DC_sector        electron DC sector
 e_DC_Chi2N         electron tracking chi2/Ndf
 
 
 [Originally based on Ex1_CLAS12Reader.C]
