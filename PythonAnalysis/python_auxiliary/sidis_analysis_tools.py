import numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib as mpl
software_path = '/Users/erezcohen/Desktop/Software/'
import sys; sys.path.insert(0, software_path + '/mySoftware/Python/');
import sys; sys.path.insert(0, software_path + '/CLAS12/BAND/SIDIS_at_BAND/PythonAnalysis/python_auxiliary/');
from my_tools               import *;
from plot_tools             import *;
from my_data_analysis_tools import *;
#from event_selection_tools  import *;


r2d = 180./np.pi

pi_charge_names  = ['piplus'   ,'piminus'  ]
pi_labels        = ['\pi^{+}'  ,'\pi^{-}'  ]
pi_prints        = ['π+'       ,'π-'       ]
pi_colors        = ['royalblue','salmon'   ]

NphiBins    = 180
phi_min     = -180 # deg.
phi_max     = 180  # deg.
phi_bins    = np.linspace(phi_min, phi_max,NphiBins+1)
phi_centers = (phi_bins[1:]+phi_bins[:-1])/2
Nphi_pts    = len(phi_centers)
phi_xticks  = [-180,-60,60,180]
phi_xlim    = [-180,180]

# (e,e'π) and (e,e'πn) events
global e_e_pi          , e_e_pi_n;
# (e,e'π) and (e,e'πn) events after all selection cuts
# (p-theta acceptance mathcing, Mx)
global e_e_pi_pass_cuts, e_e_pi_n_pass_cuts;
# (e,e'π) simulated events, generated from uniform electron direction and uniform pion direction and momoentum
global e_e_pi_GEMC     , e_e_pi_GEMC_pass_cuts;


e_e_pi, e_e_pi_n                                   = dict(),dict()
e_e_pi_pass_cuts, e_e_pi_n_pass_cuts               = dict(),dict()
e_e_pi_GEMC                                        = dict()
e_e_pi_GEMC_pass_cuts                              = dict()



# ------------------------------------------------------------------------------------------------ #
def load_SIDIS_ratio(
    # z_bins   = np.arange(0.3,0.8,0.1),
    #                            z_widths = 0.01*np.ones(5),
                               xlabel   = "Bjorken $x$",
                               # x_bins   = np.linspace(0.2,0.6,11),
                               prefix   = 'Untagged_SIDIS_ratio_',
                               suffix   = '',
                               doPlotResults=False,
                               data_path= '/Users/erezcohen/Desktop/data/BAND/Results/'):
    '''
    Load SIDIS ratio results
    last update July-19, 2022
    
    
    '''
    
    SIDIS_results = dict()
    
    z_min_arr, z_max_arr, Zavg_pips_arr, Zavg_pims_arr = [],[],[],[]
    print('Reading files from ' + data_path)
    filelist = os.listdir(data_path)
    for filename in filelist:
        if prefix in filename and suffix in filename:
            print( 'reading',filename )
            filenameparts = filename.split('_')
            z_min     = float(filenameparts[3][4:8])
            Zavg_pips = float(filenameparts[5][5:9])
            Zavg_pims = float(filenameparts[6][5:9])
            z_max     = float(filenameparts[7][4:8])
            filelabel = 'Zmin%.3f_Zmean_pips%.3f_pims%.3f_Zmax%.3f'%(z_min,Zavg_pips,Zavg_pims,z_max)
            # print( filelabel )
            
            df = pd.read_csv( data_path + '/' + filename )
            SIDIS_results[filelabel] = df
            z_min_arr.append(z_min)
            z_max_arr.append(z_max)
            Zavg_pips_arr.append(Zavg_pips)
            Zavg_pims_arr.append(Zavg_pims)

            # x_bins
            
    
    
    
#     for z_bin,z_width in zip(z_bins,z_widths):
#         z_min,z_max = z_bin-z_width,z_bin+z_width
#         filelabel = 'z_%.3f-%.3f'%(z_bin-z_width,z_bin+z_width)
#         filename  =  data_path + prefix + filelabel + suffix  + '.csv'
#         df = pd.read_csv(filename)
#         SIDIS_results[prefix + filelabel + suffix] = df
        
        
    if doPlotResults:#{
        x     = np.array(df["$x_B$"])
        x_err = np.array(df["$\Delta x_B$"])
        # x     = (x_bins[1:] + x_bins[:-1])/2
        # x_err = (x_bins[1:] - x_bins[:-1])/2


        fig = plt.figure(figsize=(9,6))
        ax  = fig.add_subplot(1,1,1)
        # for z_bin,z_width in zip(z_bins,z_widths):
        for z_min,z_max,Zavg_pips,Zavg_pims in zip( z_min_arr, z_max_arr, Zavg_pips_arr, Zavg_pims_arr ):
            Zavg = (Zavg_pips+Zavg_pims)/2.
            # z_min,z_max = z_bin-z_width,z_bin+z_width
            # filelabel = 'z_%.3f-%.3f'%(z_bin-z_width,z_bin+z_width)
            # filename  =  prefix + filelabel + suffix
            filelabel = 'Zmin%.3f_Zmean_pips%.3f_pims%.3f_Zmax%.3f'%(z_min,Zavg_pips,Zavg_pims,z_max)
            
            df = SIDIS_results[filelabel]
            y    = df['$R$']
            y_err= (df['$\Delta R_{+}$'],df['$\Delta R_{-}$'])
            # plot
            l=ax.errorbar(x=x, xerr=x_err,  y=y, yerr=y_err,
                        marker='o',markeredgecolor='k',
                        label='$%.3f<z<%.3f, \\bar{z}=%.3f$'%(z_min,z_max,Zavg))

        set_axes(ax,xlabel,"$N(e,e'\pi^+)/N(e,e'\pi^-)$",
                 title="$\pi^+/\pi^-$ ratio as a function of $x_B$ without a tagged neutron",
                 do_add_grid=True, do_add_legend=True, fontsize=18,
                );
        plt.legend(bbox_to_anchor=(1,1.05),loc='best',fontsize=18)
        
    print('Done.')
    return SIDIS_results
# ------------------------------------------------------------------------------------------------ #






# ------------------------------------------------------------------------------------------------ #
def compute_ratio_pips_to_pims(df_dict,
                               var='xB', bins=np.linspace(0,1,10),
                               weight_option = '',
                               z_min=0, z_max=1):#{
    '''
    last edit July-19, 2022
    
    weight_option: None, 'Acc. correction as f(phi)'
    
    return:
    -----------------
    R_pips_to_pims         np.array()   number of π+ events in each x-bin / number of π-
    R_pips_to_pims_errup   np.array()   err-up in number of π+ events in each x-bin / number of π-
    R_pips_to_pims_errdw   np.array()   err-dw in number of π+ events in each x-bin / number of π-
    N_pips                 np.array()   number of π+ events in each x-bin
    N_pims,                np.array()   number of π- events in each x-bin
    Zavg_pips              float        mean z-value in the range z_min < z < z_max for π+
    Zavg_pims              float        mean z-value in the range z_min < z < z_max for π-
    
    '''
    # z_min,z_max are z limits on the pion outgoing momentum
    df_pips = df_dict['piplus']
    df_pims = df_dict['piminus']
    
    # cut on z
    df_pips = df_pips[ (z_min < df_pips.Zpi) & (df_pips.Zpi < z_max) ]
    Zavg_pips = np.mean( np.array(df_pips.Zpi)  )
    
    df_pims = df_pims[ (z_min < df_pims.Zpi) & (df_pims.Zpi < z_max) ]
    Zavg_pims = np.mean( np.array(df_pims.Zpi)  )


    pips = df_pips[var]
    pims = df_pims[var]
    if weight_option == 'Acc. correction as f(phi)':#{
        phi_pips = np.array( df_pips.pi_Phi )*r2d
        phi_pims = np.array( df_pims.pi_Phi )*r2d
    #}
        
    R_pips_to_pims, R_pips_to_pims_err = [],[]
    N_pips, N_pims = [],[]
    for x_min,x_max in zip(bins[:-1],bins[1:]):#{
        

        if weight_option == 'Acc. correction as f(phi)':#{
            # each event is weighted by the acceptance correction weight
            
            phi_pips_in_bin = phi_pips[ (x_min < pips) & (pips < x_max) ]
            W_pips_in_bin   = [ Compute_acceptance_correction_weight( 'piplus' , phi ) for phi in phi_pips_in_bin ]
            Npips_in_bin    = np.sum( W_pips_in_bin )
            
            phi_pims_in_bin = phi_pims[ (x_min < pims) & (pims < x_max) ]
            W_pims_in_bin   = [ Compute_acceptance_correction_weight( 'piminus', phi ) for phi in phi_pims_in_bin ]
            Npims_in_bin    = np.sum( W_pims_in_bin )
            
        else:
            # no weight, no acceptance correction
            
            pips_in_bin      = pips[ (x_min < pips) & (pips < x_max) ]
            Npips_in_bin     = len(pips_in_bin)
            
            pims_in_bin      = pims[ (x_min < pims) & (pims < x_max) ]
            Npims_in_bin     = len(pims_in_bin)
        #}

        R     = Npips_in_bin / np.max([Npims_in_bin,1])
        R_err = R * np.sqrt( 1./np.max([1,Npips_in_bin]) + 1./np.max([1,Npims_in_bin]) )

        N_pips            .append(Npips_in_bin)
        N_pims            .append(Npims_in_bin)
        R_pips_to_pims    .append(R)
        R_pips_to_pims_err.append(R_err)
    #}
    R_pips_to_pims_errup,R_pips_to_pims_errdw = get_err_up_dw(R_pips_to_pims, R_pips_to_pims_err)
    
    return [np.array(R_pips_to_pims),
            np.array(R_pips_to_pims_errup),
            np.array(R_pips_to_pims_errdw),
            np.array(N_pips),
            np.array(N_pims),
            Zavg_pips,
            Zavg_pims]
#}
# ------------------------------------------------------------------------------------------------ #








# ------------------------------------------------------------------------------------------------ #
def get_err_up_dw(x, xerr,lim_dw = 0,lim_up = 10):
    '''
    last edit Apr-28, 2022
    '''
    errup=xerr
    errdw=xerr
    for i in range(len(x)):#{
        if (x[i]+errup[i]) > lim_up:   errup[i] = lim_up-x[i]
        if lim_dw > (x[i]-errdw[i]):   errdw[i] = x[i]-lim_dw
    #}
    return errup,errdw
#}
# ------------------------------------------------------------------------------------------------ #






# ------------------------------------------------------------------------------------------------ #
def load_SIDIS_data(runs_filename  = "good_runs_10-2-final.txt",
                    main_data_path = '/Users/erezcohen/Desktop/data/BAND/',
                    Nruns          = 1,
                    do_e_e_pi      = True,
                    do_e_e_pi_n    = True,
                    fdebug         = 2):#{
    '''
    Load SIDIS data, and fill e_e_pi and e_e_pi_n with data
    last update July-7, 2022
    
    input:
    -------------
    do_e_e_pi       flag to read (e,e'π) data - takes much time for a large number of runs
    do_e_e_pi_n     flag to read (e,e'πn) data - takes less time
    
    Comments:
    -------------
    e_e_pi & e_e_pi_n are dict(), of the following keys: ['piplus','piminus']
    e.g. :
    e_e_pi['piplus'] = pandas.DataFrame( (e,e'π) events data )
    
    '''
    global e_e_pi, e_e_pi_n;

    e_e_pi_data_path   = main_data_path + 'SIDIS_skimming/'
    e_e_pi_n_data_path = main_data_path + 'merged_SIDIS_and_BAND_skimming/'
    runfile_path = "/Users/erezcohen/Desktop/Software/CLAS12/BAND/SIDIS_at_BAND/macros/runlists/";
    
    # Using readlines()
    runs_file     = open( runfile_path + runs_filename, 'r')
    run_fileLines = runs_file.readlines()
    if Nruns==-1: Nruns = len(run_fileLines)


    runs = []
    for line in run_fileLines[0:Nruns]:#{
        run = int(line.strip())
        runs.append(run)
    #}
    runs = np.array(runs)

    
    for runnum,runIdx in zip(runs,range(len(runs))):
        if fdebug>1: print('Run number ',runnum,'(%d/%d runs)'%(runIdx+1,len(runs)))
        for pi_charge_name,pi_print in zip(pi_charge_names,pi_prints):
            if do_e_e_pi:#{
                eepi   = pd.read_csv(e_e_pi_data_path
                                     +'skimmed_SIDIS_inc_00%d_e_%s_selected_eepi_kinematics.csv'%
                                     (runnum,pi_charge_name),
                                     usecols=['runnum','evnum',
                                              'e_P','e_Theta','e_Phi',
                                              'pi_P', 'pi_Theta', 'pi_Phi',
                                              'Q2', 'W', 'xB', 'Zpi',
                                              'M_X', 'e_DC_sector', 'pi_DC_sector'],
                                     dtype={'runnum':int,'evnum': int,
                                            'e_DC_sector':int, 'pi_DC_sector':int,
                                            'e_P':np.half,'e_Theta':np.half,'e_Phi':np.half,
                                            'pi_P':np.half,'pi_Theta':np.half, 'pi_Phi':np.half,
                                            'Q2':np.half, 'W':np.half, 'xB':np.half, 'Zpi':np.half,
                                            'M_X':np.half})
                
                if runIdx==0: e_e_pi[pi_charge_name] = eepi
                else:         e_e_pi[pi_charge_name] = pd.concat([e_e_pi[pi_charge_name],eepi])

                
                if fdebug>1: print('Loaded',len(eepi)," (e,e'"+pi_print+") events")

            #}
            if do_e_e_pi_n:#{
                eepin = pd.read_csv(e_e_pi_n_data_path
                                +'skimmed_SIDIS_and_BAND_inc_00%d_e_%s_n.csv'%
                                (runnum,pi_charge_name))
                
                if fdebug>1: print('Loaded',len(eepin)," (e,e'"+pi_print+"n) events")

                if runIdx==0: e_e_pi_n[pi_charge_name] = eepin
                else:         e_e_pi_n[pi_charge_name] = pd.concat([e_e_pi_n[pi_charge_name],eepin])
            #}
    print('Done loading files.')
    if fdebug>0:
        print('')
        print('Total statistics:')
        for pi_charge_name,pi_print in zip(pi_charge_names,pi_prints):
            if do_e_e_pi:   print(len(e_e_pi[pi_charge_name])," (e,e'"+pi_print+") events")
            if do_e_e_pi_n: print(len(e_e_pi_n[pi_charge_name])," (e,e'"+pi_print+"n) events")
    #}
#}
# ------------------------------------------------------------------------------------------------ #












# ------------------------------------------------------------------------------------------------ #
def apply_further_selection_cuts_to_data(fdebug=2,
                                         NeventsMax=-1,
                                         NMaxPerSubset = 500000,
                                         doAcceptanceMatchingCut = True,
                                         doApply_minPn_cut       = True):#{
    '''
    e_e_pi_pass_cuts, e_e_pi_n_pass_cuts, e_e_pi_GEMC_pass_cuts = apply_further_selection_cuts_to_data(fdebug=2)
    last edit July-12, 2022
    
    Apply selection cuts not previously imposed
    
    The cuts applied for (e,e'π) events:
    1. pi+/pi- acceptance matching cut in p-theta plane
    2. Missing mass cut on (e,e'\pi) events
    
    The cuts applied for (e,e'πn) events:
    1. pi+/pi- acceptance matching cut in p-theta plane
    
    '''
    global e_e_pi, e_e_pi_n, e_e_pi_GEMC
    global e_e_pi_pass_cuts, e_e_pi_n_pass_cuts, e_e_pi_GEMC_pass_cuts
    
    
    # (e,e'\pi n) SIDIS data complete this -  need to add sector ID in the (e,e'\pi n) data
    print('Applying selection cuts not previously imposed')
    
    # print number of events retained on every cut
    if fdebug < 1: return
    Nevents      = dict()
    frac_Nevents = dict()
    
    # (1) Data
    if (e_e_pi=={}) is False:#{
        print('(1) DATA')
        print("(e,e'π)")
        
        # (e,e'\pi) SIDIS data
        if doAcceptanceMatchingCut:
            e_e_pi_after_p_theta_cut = apply_p_theta_acceptance_cut( e_e_pi,
                                                                    NeventsMax=NeventsMax,
                                                                    NMaxPerSubset=NMaxPerSubset,
                                                                    fdebug=fdebug )
        else:
            e_e_pi_after_p_theta_cut = dict()
            for pi_ch,pi_print in zip(pi_charge_names,pi_prints):#{
                df = e_e_pi[pi_ch]
                e_e_pi_after_p_theta_cut[pi_ch] = df[0:NeventsMax]
            #}
            
            
        e_e_pi_after_Mx_cut      = apply_Mx_cut( e_e_pi_after_p_theta_cut )
        
        e_e_pi_after_Kinematical_cuts = apply_Kinematical_cuts( e_e_pi_after_Mx_cut )
        
        e_e_pi_pass_cuts         = e_e_pi_after_Kinematical_cuts;

        for pi_ch,pi_print in zip(pi_charge_names,pi_prints):#{
            print('(e,e',pi_print,')')
            if NeventsMax < 0:
                Nevents[pi_ch + ' original'] = len(e_e_pi[pi_ch])
            else:
                Nevents[pi_ch + ' original'] = NeventsMax
            frac_Nevents[pi_ch + ' original'] = 1
            print(Nevents[pi_ch + ' original'],'events before cut')

            Nevents[pi_ch +' p-theta cut'] = len(e_e_pi_after_p_theta_cut[pi_ch])
            frac_Nevents[pi_ch + ' p-theta cut'] = float(Nevents[pi_ch +' p-theta cut'])/ Nevents[pi_ch + ' original']
            print(Nevents[pi_ch +' p-theta cut'],'events after p-theta cut (%.1f'%(100.*frac_Nevents[pi_ch + ' p-theta cut']),'%)')


            Nevents[pi_ch +' Mx cut'] = len(e_e_pi_after_Mx_cut[pi_ch])
            frac_Nevents[pi_ch + ' Mx cut'] = float(Nevents[pi_ch +' Mx cut'])/Nevents[pi_ch + ' original']
            print(Nevents[pi_ch +' Mx cut'],'events after M_X cut (%.1f'%(100.*frac_Nevents[pi_ch + ' Mx cut']),'%)')

            Nevents[pi_ch +' Kinematical cut'] = len(e_e_pi_after_Kinematical_cuts[pi_ch])
            frac_Nevents[pi_ch + ' Kinematical cut'] = float(Nevents[pi_ch +' Kinematical cut'])/Nevents[pi_ch + ' original']
            print(Nevents[pi_ch +' Kinematical cut'],'events after further Kinematical cut (%.1f'%(100.*frac_Nevents[pi_ch + ' Kinematical cut']),'%)')

        #}
        print(' ')
    #}
    if (e_e_pi_n=={}) is False:#{
        print("(e,e'πn)")
        
        # tagged (e,e'\pi) SIDIS data
        if doApply_minPn_cut:#{
            e_e_pi_n_after_minPn_cut   = apply_minPn_cut( e_e_pi_n )
        else:
            e_e_pi_n_after_minPn_cut = dict()
            for pi_ch,pi_print in zip(pi_charge_names,pi_prints):#{
                df = e_e_pi_n[pi_ch]
                e_e_pi_n_after_minPn_cut[pi_ch] = df[0:NeventsMax]
            #}
        #}

            
        if doAcceptanceMatchingCut:
            e_e_pi_n_after_p_theta_cut = apply_p_theta_acceptance_cut( e_e_pi_n_after_minPn_cut,
                                                                      NeventsMax=NeventsMax,
                                                                      NMaxPerSubset=NMaxPerSubset,
                                                                      fdebug=fdebug )
        else:
            e_e_pi_n_after_p_theta_cut = e_e_pi_n_after_minPn_cut
            
        e_e_pi_n_after_Kinematical_cuts = apply_Kinematical_cuts( e_e_pi_n_after_p_theta_cut )
        e_e_pi_n_pass_cuts  = e_e_pi_n_after_Kinematical_cuts;

        for pi_ch,pi_print in zip(pi_charge_names,pi_prints):#{
            print('(e,e',pi_print,')')
            if NeventsMax < 0:
                Nevents[pi_ch + 'n original'] = len(e_e_pi_n[pi_ch])
            else:
                Nevents[pi_ch + 'n original'] = NeventsMax
                
            frac_Nevents[pi_ch + 'n original'] = 1
            print(Nevents[pi_ch + 'n original'],'events before cut')

            Nevents[pi_ch +'n min p(n) cut'] = len(e_e_pi_n_after_minPn_cut[pi_ch])
            frac_Nevents[pi_ch + 'n min p(n) cut'] = float(Nevents[pi_ch +'n min p(n) cut'])/ Nevents[pi_ch + 'n original']
            print(Nevents[pi_ch +'n min p(n) cut'],'events after min p(n) cut (%.1f'%(100.*frac_Nevents[pi_ch + 'n min p(n) cut']),'%)')

            
            Nevents[pi_ch +'n p-theta cut'] = len(e_e_pi_n_after_p_theta_cut[pi_ch])
            frac_Nevents[pi_ch + 'n p-theta cut'] = float(Nevents[pi_ch +'n p-theta cut'])/ Nevents[pi_ch + 'n original']
            print(Nevents[pi_ch +'n p-theta cut'],'events after p-theta cut (%.1f'%(100.*frac_Nevents[pi_ch + 'n p-theta cut']),'%)')
            
            Nevents[pi_ch +'n Kinematical cut'] = len(e_e_pi_n_after_Kinematical_cuts[pi_ch])
            print(Nevents.keys())
            frac_Nevents[pi_ch + 'n Kinematical cut'] = float(Nevents[pi_ch +'n Kinematical cut'])/Nevents[pi_ch + 'n original']
            print(Nevents[pi_ch +'n Kinematical cut'],'events after further Kinematical cut (%.1f'%(100.*frac_Nevents[pi_ch + 'n Kinematical cut']),'%)')

        #}
        print(' ')
    #}
    
    # (2) MC
    if (e_e_pi_GEMC=={}) is False:#{
        print('(2) MC')
        
        # (e,e'\pi) - (uniform) MC for acceptance correction (uniform in e and \pi)
        e_e_pi_GEMC_after_eepi_cuts       = dict()

        # Apply (e,e'pi) SIDIS kinematical cuts while asking if pion was accepted,
        # externally (here, and not in the CLAS12ROOT script) since we
        # want to retain and record also the events that did not pass these cuts, in the simulation
        # whereas in data we just omit events that did not pass these cuts
        for pi_ch in pi_charge_names:#{
            e_e_pi_GEMC_after_eepi_cuts[pi_ch] = e_e_pi_GEMC[pi_ch][(e_e_pi_GEMC[pi_ch].pi_passed_cuts==1) & (e_e_pi_GEMC[pi_ch].eepiPastKinematicalCuts==1)];
        #}
        e_e_pi_GEMC_after_p_theta_cut = apply_p_theta_acceptance_cut( e_e_pi_GEMC_after_eepi_cuts )
        e_e_pi_GEMC_after_Mx_cut      = apply_Mx_cut(  e_e_pi_GEMC_after_p_theta_cut )
        e_e_pi_GEMC_pass_cuts         = e_e_pi_GEMC_after_Mx_cut;


        
        # print number of events retained on every cut in the uniform GEMC
        if fdebug<2: return
        for pi_ch,pi_print in zip(pi_charge_names,pi_prints):#{
            print('(e,e',pi_print,') in uniform GEMC simulation')

            Nevents[pi_ch + ' GEMC original'] = len(e_e_pi_GEMC[pi_ch])
            frac_Nevents[pi_ch + ' GEMC original'] = 1
            print(Nevents[pi_ch + ' GEMC original'],'events before cut')

            Nevents[pi_ch +' GEMC p-theta cut'] = len(e_e_pi_GEMC_after_p_theta_cut[pi_ch])
            frac_Nevents[pi_ch + ' GEMC p-theta cut'] = float(Nevents[pi_ch +' GEMC p-theta cut'])/ Nevents[pi_ch + ' GEMC original']
            print(Nevents[pi_ch +' GEMC p-theta cut'],'events after p-theta cut (%.1f'%(100.*frac_Nevents[pi_ch + ' GEMC p-theta cut']),'%)')


            Nevents[pi_ch +' GEMC Mx cut'] = len(e_e_pi_GEMC_after_Mx_cut[pi_ch])
            frac_Nevents[pi_ch + ' GEMC Mx cut'] = float(Nevents[pi_ch +' GEMC Mx cut'])/Nevents[pi_ch + ' GEMC original']
            print(Nevents[pi_ch +' GEMC Mx cut'],'events after M_X cut (%.1f'%(100.*frac_Nevents[pi_ch + ' GEMC Mx cut']),'%)')
        #}
        print(' ')
    #}
    print('Done applying selection cuts not previously imposed')
    return e_e_pi_pass_cuts, e_e_pi_n_pass_cuts, e_e_pi_GEMC_pass_cuts
#}
# ------------------------------------------------------------------------------------------------ #







# ------------------------------------------------------------------------------------------------ #
def apply_Kinematical_cuts( df_dict,
                           Q2_min=2,     Q2_max=10,
                           Pe_min=3,     Pe_max=10.6,
                           Ppi_min=1.25, Ppi_max=5 ,
                           W_min=2.5,    W_max=10    ):#{
    '''
        df_dict_after_cut = apply_Kinematical_cuts(df_dict)
        
        Apply kinematical cuts to match BAND neutron skimming
        last update July-14, 2022
        
    '''
    df_dict_after_cut = dict()
    for pi_charge_name in pi_charge_names:#{
        df = df_dict[pi_charge_name]
        df = df[ (W_min < df.W)      & (df.W < W_max)     ]
        df = df[ (Q2_min < df.Q2)    & (df.Q2 < Q2_max)   ]
        df = df[ (Pe_min < df.e_P)   & (df.e_P < Pe_max)  ]
        df = df[ (Ppi_min < df.pi_P) & (df.pi_P < Ppi_max)]
        df_dict_after_cut[pi_charge_name] = df
    #}
    return df_dict_after_cut
#}
# ------------------------------------------------------------------------------------------------ #






# ------------------------------------------------------------------------------------------------ #
def apply_p_theta_acceptance_cut( df_dict=None,
                                 NeventsMax=-1,
                                 NMaxPerSubset = 500000,
                                 fdebug=1):
    '''
        df_dict_after_cut = apply_p_theta_acceptance_cut(df_dict)
        
        Apply a π+/π- acceptance matching cut on the in p-\theta plane
        last update July-12, 2022
        
    '''
    import numpy as np
    print("Apply a π+/π- acceptance matching cut on the in p-theta plane")
    df_dict_after_cut = dict()
    
    
    for pi_ch in pi_charge_names:
        
        if NeventsMax > 0: NeventsMax = np.min( [NeventsMax, len(df_dict[pi_ch])] )
        else:              NeventsMax = len(df_dict[pi_ch])
        df = df_dict[pi_ch][0:NeventsMax]
        if fdebug: print('Applying p-theta on cut for '+pi_ch+' on %d events'%NeventsMax)

        # subdivide samples into sub-samples of no more than 1M events each
        Nsubsets  = np.max([1,(int)(NeventsMax/NMaxPerSubset)]);
        if fdebug: print('Subdividing into %d subsets up to %d events'%(Nsubsets,NMaxPerSubset))
        
        for subset_idx in range(Nsubsets):#{
            i_min_subset = NMaxPerSubset*subset_idx
            i_max_subset = np.min([NMaxPerSubset*(subset_idx+1), NeventsMax])
            subset_df = df[i_min_subset:i_max_subset]
            good_indices = np.array([])
            if fdebug>1:
                print('subset %d of index %d-%d'%(subset_idx,i_min_subset,i_max_subset-1))
            
            # focus on each sector and select the events that pass the p-theta cut
            for sector in range(1,7):#{
                df_in_sector   = subset_df[subset_df.pi_DC_sector == sector]
                pi_P           = np.array(df_in_sector.pi_P)
                theta_min_pi   = pi_min_theta_cut( pi_charge = 'any', sector=sector, p=pi_P )
                good_indices_in_sector = []
                good_indices_in_sector = df_in_sector[ df_in_sector.pi_Theta*r2d > theta_min_pi ].index
                good_indices = np.concatenate([good_indices,good_indices_in_sector])
            #}
            good_indices = np.unique(good_indices)
            subset_df_in_cut = subset_df.loc[good_indices]
            if subset_idx==0: df_after_cut = subset_df_in_cut
            else:             df_after_cut = pd.concat([df_after_cut , subset_df_in_cut])
        #}
        df_dict_after_cut[pi_ch] = df_after_cut
    return df_dict_after_cut
# ------------------------------------------------------------------------------------------------ #







# ------------------------------------------------------------------------------------------------ #
def apply_minPn_cut( df_dict=None, Pn_min = 0.275 ): #{
    '''
        df_dict_after_cut = apply_minPn_cut( df_dict, Pn_min = 0.275)
        
        Apply a cut on the minimal neutron momentum
        
    '''
    print("Apply a cut on the minimal neutron momentum p > %g GeV/c"%Pn_min)
    df_dict_after_cut = dict()
    for pi_charge_name in pi_charge_names:#{
        df = df_dict[pi_charge_name]
        df = df[ Pn_min < df.n_P ]
        df_dict_after_cut[pi_charge_name] = df
    #}
    return df_dict_after_cut
#}
# ------------------------------------------------------------------------------------------------ #





# ------------------------------------------------------------------------------------------------ #
def pi_min_theta_cut( pi_charge = 'pips', sector=1, p=2 ):
    '''
    p-theta cut optimized by Alex K., Feb-2022
    last update June-10, 2022
    
    input:
    --------
    pi_charge   'pips'/'pims'/'piplus'/'piminus'/'any'/'all'
    p           pion momentum in GeV/c
    sector      sector in the DC
    
    return:
    --------
    theta_min   minimal theta as a function of p in deg.
    
    '''
    pips_parameters = [[5.83055147, 1.16171268],
                       [5.78469108, 0.80349681],
                       [5.844136,   0.53329847],
                       [5.81600646, 0.62782215],
                       [5.75400247, 0.88392704],
                       [5.4041441,  2.12325929]]
    
    pims_parameters = [[ 6.72342364, 14.60912754],
                       [ 6.85593523, 14.2618346 ],
                       [ 6.66919406, 14.64725313],
                       [ 6.77248325, 14.33494692],
                       [ 6.63926263, 14.4344934 ],
                       [ 6.85481318, 14.3687773 ]]


        
    if pi_charge=='any' or pi_charge=='all':
        theta_min_pips = pips_parameters[sector-1][0] + pips_parameters[sector-1][1]/p
        theta_min_pims = pims_parameters[sector-1][0] + pims_parameters[sector-1][1]/p
        theta_min = [np.max([theta_min_pips[i],theta_min_pims[i]]) for i in range(len(theta_min_pips))]
        
    elif pi_charge=='pips' or pi_charge=='piplus':
        parameters = pips_parameters
        theta_min = parameters[sector-1][0] + parameters[sector-1][1]/p
            
    elif pi_charge=='pims' or pi_charge=='piminus':
        parameters = pims_parameters
        theta_min = parameters[sector-1][0] + parameters[sector-1][1]/p
    

       
    return theta_min
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def apply_Mx_cut( df_dict=None, Mx_min = 2.5, Mx_max = 5 ): #{
    '''
        df_dict_after_cut = apply_Mx_cut(df_dict)
        
        Apply a cut on the missing mass of a (e,e'\pi) reaction
        
    '''
    print("Apply a cut on the missing mass of a (e,e'π) reaction: %.1f<Mx<%.1f GeV/c2"%(Mx_min,Mx_max))
    df_dict_after_cut = dict()
    for pi_charge_name in pi_charge_names:#{
        df = df_dict[pi_charge_name]
        df = df[ (Mx_min < df.M_X) & (df.M_X < Mx_max)]
        df_dict_after_cut[pi_charge_name] = df
    #}
    return df_dict_after_cut
#}
# ------------------------------------------------------------------------------------------------ #


