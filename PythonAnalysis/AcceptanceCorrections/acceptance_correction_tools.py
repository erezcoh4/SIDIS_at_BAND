import numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib as mpl
import sys; sys.path.insert(0, '/Users/erezcohen/Desktop/Software/mySoftware/Python/'); 
from my_tools               import *; 
from plot_tools             import *;
from my_data_analysis_tools import *;


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
global e_e_pi, e_e_pi_n; 
# (e,e'π) and (e,e'πn) events after all selection cuts
# (p-theta acceptance mathcing, Mx)
global e_e_pi_pass_cuts, e_e_pi_n_pass_cuts;

# (e,e'π) simulated events, generated from uniform electron direction and uniform pion direction and momoentum
global e_e_pi_GEMC;
global e_e_pi_GEMC_pass_cuts; 
global h, h_err;


# (binary 0/1) flag of tight fiducial
# The areas in \phi that are "good", which we want to keep (0 or 1 for each of phi_centers)
global TightFiducialPhi;
global idx_good_phi # same as 'TightFiducialPhi' albeit the indices in phi that are "good"
global idx_bad_phi # inverse of 'idx_good_phi'

# The fraction of area occupied by "good" phi in TightFiducialPhi (out of 2\pi)
global TightFiducialPhiAreaFraction;
global fraction_bad_phi; # same as above but inverse

# acceptance correction as a function of pion \phi angle
# integral "height" of the acceptance correction for pi+ and pi-, based on GEMC
global AccCorrecHeight, AccCorrecHeight_err;
# acceptance corrcetion in each \phi bin
global AccCorrec,       AccCorrec_err;
# acceptance correction in "good" \phi regions (=0 in "bad" \phi regions)
global AccCorrecTightFiducial, AccCorrecTightFiducial_err;





e_e_pi, e_e_pi_n                                   = dict(),dict()
e_e_pi_pass_cuts, e_e_pi_n_pass_cuts               = dict(),dict()
e_e_pi_GEMC                                        = dict()
e_e_pi_GEMC_pass_cuts                              = dict()
h, h_err                                           = dict(),dict()
AccCorrecHeight, AccCorrecHeight_err               = dict(),dict() 
AccCorrec, AccCorrec_err                           = dict(),dict()
AccCorrecTightFiducial, AccCorrecTightFiducial_err = dict(),dict()
TightFiducialPhi, TightFiducialPhiAreaFraction     = dict(),dict()
idx_good_phi,idx_bad_phi                           = dict(),dict()
fraction_bad_phi                                   = dict()









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
        if fdebug>1: print('Subdividing into %d subsets up to %d events'%(Nsubsets,NMaxPerSubset))
        
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
        e_e_pi_pass_cuts         = e_e_pi_after_Mx_cut;

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
        #}
        print(' ')
    #}
    if (e_e_pi_n=={}) is False:#{
        print("(e,e'πn)")
        
        # tagged (e,e'\pi) SIDIS data
        if doApply_minPn_cut:
            e_e_pi_n_after_minPn_cut   = apply_minPn_cut( e_e_pi_n )
        else:
            e_e_pi_n_after_p_theta_cut = dict()
            for pi_ch,pi_print in zip(pi_charge_names,pi_prints):#{
                df = e_e_pi_n[pi_ch]
                e_e_pi_n_after_p_theta_cut[pi_ch] = df[0:NeventsMax]
            #}

            
        if doAcceptanceMatchingCut:
            e_e_pi_n_after_p_theta_cut = apply_p_theta_acceptance_cut( e_e_pi_n_after_minPn_cut,
                                                                      NeventsMax=NeventsMax,
                                                                      NMaxPerSubset=NMaxPerSubset,
                                                                      fdebug=fdebug )
        else:
            e_e_pi_n_after_p_theta_cut = e_e_pi_n_after_minPn_cut
            
        e_e_pi_n_pass_cuts         = e_e_pi_n_after_p_theta_cut;

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
def load_SIDIS_ratio_DataFrame(z_bins   = np.arange(0.3,0.8,0.1),
                               z_widths = 0.01*np.ones(5),
                               prefix   = 'Untagged_SIDIS_ratio_',
                               suffix   = '',
                               doPlotResults=False,
                               data_path= '/Users/erezcohen/Desktop/data/BAND/Results/'):
    '''
    Load SIDIS ratio results
    last update July-12, 2022
    
    
    '''
    
    SIDIS_results = dict()
    for z_bin,z_width in zip(z_bins,z_widths):
        z_min,z_max = z_bin-z_width,z_bin+z_width
        filelabel = 'z_%.2f-%.2f'%(z_bin-z_width,z_bin+z_width)
        filename  =  data_path + prefix + filelabel + suffix  + '.csv'
        df = pd.read_csv(filename)
        SIDIS_results[prefix + filelabel + suffix] = df
        
        
    if doPlotResults:#{
        fig = plt.figure(figsize=(9,6))
        ax  = fig.add_subplot(1,1,1)
        for z_bin,z_width in zip(z_bins,z_widths):

            z_min,z_max = z_bin-z_width,z_bin+z_width
            filelabel = 'z_%.2f-%.2f'%(z_bin-z_width,z_bin+z_width)
            filename  =  prefix + filelabel + suffix

            df = SIDIS_results[filename]
            y    = df['$R$']
            y_err= (df['$\Delta R_{+}$'],df['$\Delta R_{-}$'])
            # plot
            l=ax.errorbar(x=x, xerr=x_err,  y=y, yerr=y_err,
                        marker='o',markeredgecolor='k',
                        label='$z=%.2f\pm%.2f$'%(z_bin,z_width))

        set_axes(ax,xlabel,"$N(e,e'\pi^+)/N(e,e'\pi^-)$",
                 title="$\pi^+/\pi^-$ ratio as a function of $x_B$ without a tagged neutron",
                 do_add_grid=True, do_add_legend=True, fontsize=18,
                );
        plt.legend(bbox_to_anchor=(1,1.05),loc='best',fontsize=18)
        
    return SIDIS_results
# ------------------------------------------------------------------------------------------------ #





# ------------------------------------------------------------------------------------------------ #
def save_SIDIS_ratio_DataFrame(df_dict  = None,
                               x_var    = 'xB' ,
                               x_bins   = np.linspace(0.2,0.6,11),
                               z_bins   = np.arange(0.3,0.8,0.1),
                               z_widths = 0.01*np.ones(5),
                               data_path= '/Users/erezcohen/Desktop/data/BAND/Results/',
                               fdebug   = 0,
                               prefix   = 'Untagged_SIDIS_ratio_',
                               suffix   = ''):
    '''
    Save SIDIS ratio results
    last update July-12, 2022
    
    
    '''
    
    x        = (x_bins[1:] + x_bins[:-1])/2
    x_err    = (x_bins[1:] - x_bins[:-1])/2
    results_data_path = data_path + '/' + 'Results' + '/'
    for z_bin,z_width in zip(z_bins,z_widths):
        z_min,z_max = z_bin-z_width,z_bin+z_width
        (R,R_err_up,R_err_dw,
         N_pips,N_pims) = compute_ratio_pips_to_pims(df_dict= df_dict ,
                                                     var    = x_var,
                                                     bins   = x_bins,
                                                     z_min  = z_min,
                                                     z_max  = z_max)

        df_to_save = pd.DataFrame({"$x_B$":x,
                                   "$\Delta x_B$":x_err,
                                   '$N(\pi_{+})$':N_pips,
                                   '$N(\pi_{-})$':N_pims,
                                   '$R$':R,
                                   '$\Delta R_{+}$':R_err_up,
                                   '$\Delta R_{-}$':R_err_dw})
        filelabel = 'z_%.2f-%.2f'%(z_bin-z_width,z_bin+z_width)
        filename  =  data_path + prefix + filelabel + suffix  + '.csv'
        df_to_save.to_csv(filename)
        print('saved',filename)
        if fdebug>1:
            print('$z=%.2f\pm%.2f$'%(z_bin,z_width))
            print(filename)
            display(df_to_save)
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
def compute_ratio_pips_to_pims(df_dict,
                               var='xB', bins=np.linspace(0,1,10),
                               weight_option = '',
                               z_min=0, z_max=1):#{
    '''
    last edit July-8, 2022
    
    weight_option: None, 'Acc. correction as f(phi)'
    '''
    # z_min,z_max are z limits on the pion outgoing momentum
    df_pips = df_dict['piplus']
    df_pims = df_dict['piminus']
    
    # cut on z
    df_pips = df_pips[ (z_min < df_pips.Zpi) & (df_pips.Zpi < z_max) ]
    df_pims = df_pims[ (z_min < df_pims.Zpi) & (df_pims.Zpi < z_max) ]

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
            
            pips_in_bin  = pips[ (x_min < pips) & (pips < x_max) ]
            Npips_in_bin = len(pips_in_bin)
            pims_in_bin  = pims[ (x_min < pims) & (pims < x_max) ]
            Npims_in_bin = len(pims_in_bin)
            
        #}

        R     = Npips_in_bin / np.max([Npims_in_bin,1])
        R_err = R * np.sqrt( 1./np.max([1,Npips_in_bin]) + 1./np.max([1,Npims_in_bin]) )

        N_pips            .append(Npips_in_bin)
        N_pims            .append(Npims_in_bin)
        R_pips_to_pims    .append(R)
        R_pips_to_pims_err.append(R_err)
    #}
    R_pips_to_pims_errup,R_pips_to_pims_errdw = get_err_up_dw(R_pips_to_pims, R_pips_to_pims_err)
    
    return np.array(R_pips_to_pims),np.array(R_pips_to_pims_errup),np.array(R_pips_to_pims_errdw),np.array(N_pips), np.array(N_pims)
#}
# ------------------------------------------------------------------------------------------------ #












# ------------------------------------------------------------------------------------------------ #
def compute_acceptance_correction_as_function_of_phi_from_MC( ):#{
    '''    
    compute_acceptance_correction_as_function_of_phi_from_MC()
    last update May-26, 2022
    
    
    Comments:
    
    We omit regions where the acceptance-correction is greater than the "plateau value"
    as we want to avoid large acceptance correction uncertainty.
    This effectively means tightening the fiducial cut
    We do this by,
    first finding this plateau by computing the meidan of the acceptance correction
    and defining the "good" regions as the ones with a correction
    smaller than the median + std/3
    and assigning an acceptance correction of 0 to these regions
  
    '''
  


    
    global h, h_err, AccCorrec, AccCorrec_err
    global TightFiducialPhi, AccCorrecTightFiducial, AccCorrecTightFiducial_err
    global e_e_pi_GEMC, e_e_pi_GEMC_pass_cuts

    var_gen      = 'pi_Phi_g'
    label        = '$\phi$'
    units        = '[deg.]'
    scale_factor = r2d

    # (1) compute histograms of generated, reconstructed and accepted events asna function of $\phi$
    fig = plt.figure(figsize=(14,6));
    for pi_ch,pi_charge_label,pi_idx in zip(pi_charge_names,pi_labels,range(2)):
        df_gen = e_e_pi_GEMC[pi_ch];
        df_rec = e_e_pi_GEMC[pi_ch][e_e_pi_GEMC[pi_ch].pi_reconstructed==1];
        df_acc = e_e_pi_GEMC_pass_cuts[pi_ch][e_e_pi_GEMC_pass_cuts[pi_ch].pi_passed_cuts==1];    
        x_gen = df_gen[var_gen]
        x_rec = df_rec[var_gen]    
        x_acc = df_acc[var_gen]    


        ax = fig.add_subplot(1,2,pi_idx+1)
        x,h[pi_ch+'gen'],x_err,h_err[pi_ch+'gen'] = plot_step_hist( x_gen*scale_factor, phi_bins, color='k',           label='generated')
        x,h[pi_ch+'rec'],x_err,h_err[pi_ch+'rec'] = plot_step_hist( x_rec*scale_factor, phi_bins, color='royalblue',   label='reconsutrcted')
        x,h[pi_ch+'acc'],x_err,h_err[pi_ch+'acc'] = plot_step_hist( x_acc*scale_factor, phi_bins, color='forestgreen', label='accepted')            

        # ax.set_yscale('log')
        set_axes(ax,label + ' ' + units,'counts' if pi_idx==0 else '',
                    do_add_legend=True if pi_idx==1 else False,
                     title='$'+pi_charge_label+'$ acceptance as a function of '+label, fontsize=24, 
                 do_add_grid=True,xlim=phi_xlim,xticks=phi_xticks)

    plt.tight_layout();
    
    # (2) compute acceptance correction
    for pi_ch,pi_idx in zip(pi_charge_names,range(2)):#{
        h[pi_ch+'eff'] = h[pi_ch+'acc']/h[pi_ch+'gen']
        h_err[pi_ch+'eff'] = h[pi_ch+'eff']*np.sqrt( np.square(h_err[pi_ch+'acc']/h[pi_ch+'acc']) + np.square(h_err[pi_ch+'gen']/h[pi_ch+'gen'])  )
        AccCorrec[pi_ch]     = 1./h[pi_ch+'eff']
        AccCorrec_err[pi_ch] = h_err[pi_ch+'eff']/np.square( h[pi_ch+'eff'] )
    #}
    
    # (3) omit regions where the acceptance-correction is greater than the "plateau value"
    # (3.1) first find this plateau
    median_correction_phi_Plateau,std_correction_phi_Plateau = dict(), dict()
    for pi_ch,pi_idx in zip(pi_charge_names,range(2)):#{
        median_correction_phi_Plateau[pi_ch] = np.median( AccCorrec[pi_ch] )
        std_correction_phi_Plateau[pi_ch]  = np.std( AccCorrec[pi_ch] )
    #}
    # (3.2) now compute all correction
    for pi_ch,pi_idx in zip(pi_charge_names,range(2)):#{        
        AccCorrecTightFiducial[pi_ch],AccCorrecTightFiducial_err[pi_ch] = np.zeros(Nphi_pts),np.zeros(Nphi_pts)
        TightFiducialPhi[pi_ch] = np.zeros(Nphi_pts)
        NgoodPhi = 0.0
        for phi_idx in range(Nphi_pts):#{
            if AccCorrec[pi_ch][phi_idx] > (median_correction_phi_Plateau[pi_ch] + std_correction_phi_Plateau[pi_ch]/3): #{
                AccCorrecTightFiducial_err[pi_ch][phi_idx] = 1e-5
                TightFiducialPhi[pi_ch][phi_idx]           = 0
                AccCorrecTightFiducial[pi_ch][phi_idx]     = 0            
            #}
            else :#{
                AccCorrecTightFiducial_err[pi_ch][phi_idx] = AccCorrec_err[pi_ch][phi_idx]
                AccCorrecTightFiducial[pi_ch][phi_idx]     = AccCorrec[pi_ch][phi_idx]
                TightFiducialPhi[pi_ch][phi_idx]           = 1
                NgoodPhi = NgoodPhi + 1
            #}
        #}
        TightFiducialPhiAreaFraction[pi_ch] = float(NgoodPhi)/Nphi_pts
    #}

    # (4) plot acceptance correction weight
    color,capsize,capthick,marker,linewidth = 'k', 2, 2, 'o', 2    
    fig = plt.figure(figsize=(14,12));
    for pi_ch,pi_charge_label,pi_idx in zip(pi_charge_names,pi_labels,range(2)):
        ax = fig.add_subplot(2,2,pi_idx+1)
        plt.step ( x, h[pi_ch+'eff']*100., color=color, where='mid', label=None )
        plt.errorbar ( x=x, xerr=x_err, y=h[pi_ch+'eff']*100., yerr=h_err[pi_ch+'eff']*100.,
                      color=color, marker='o', linestyle='None',label=label,
                      capsize=capsize, capthick=capthick, linewidth=linewidth )
        
        plt.fill_between ( x, np.zeros(len(x)), h[pi_ch+'eff']*100.*TightFiducialPhi[pi_ch], color='k', alpha=0.1)
        
        set_axes(ax,label + ' ' + units,'Efficiency [%]' if pi_idx==0 else '',
            do_add_legend=False,
            title='$'+pi_charge_label+'$ acceptance as a function of '+label, 
                 fontsize=24, do_add_grid=True,xlim=phi_xlim, xticks=phi_xticks)

    # for pi_ch,pi_charge_label,pi_idx in zip([eepips_GEMC,eepims_GEMC],pi_charge_names,pi_labels,range(2)):
        ax = fig.add_subplot(2,2,pi_idx+3)
        plt.step ( x, AccCorrec[pi_ch], '--', color=color, where='mid', label=None )
        plt.errorbar ( x=x, xerr=x_err, y=AccCorrecTightFiducial[pi_ch], yerr=AccCorrecTightFiducial_err[pi_ch],
                      color=color, marker='o', linestyle='None',label=label,
                      capsize=capsize, capthick=capthick, linewidth=linewidth )
        plt.fill_between ( x, np.zeros(len(x)), AccCorrecTightFiducial[pi_ch], color='k', alpha=0.1)
        set_axes(ax,label + ' ' + units,'Acceptance correction weight' if pi_idx==0 else '',
            do_add_legend=False,
            title='$'+pi_charge_label+'$ acceptance correction (%.1f'%(100.*TightFiducialPhiAreaFraction[pi_ch]) + '% covered)', 
                 fontsize=24, do_add_grid=True,xlim=phi_xlim, xticks=phi_xticks, ylim=(0.5,np.max(ax.get_ylim()))
                )
        # ymax = np.max(ax.get_ylim());
        # plt.step ( x, AccCorrecTightFiducial[pi_ch], color='royalblue', where='mid', label=None )
        # plt.fill_between ( x, np.zeros(len(x)), TightFiducialPhi[pi_ch]*ymax, color='royalblue', alpha=0.1)
        # plt.errorbar ( x=x, xerr=x_err, y=AccCorrec[pi_ch], yerr=AccCorrec_err[pi_ch],
        #               color=color, marker='o', linestyle='None',label=label,
        #               capsize=capsize, capthick=capthick, linewidth=linewidth )
    plt.tight_layout(); 
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
def ComputeAverageAcceptanceInPhiPlateau( do_plot_reconstructed=False, do_add_legend=False ):#{
    '''
    Comments:
    
    Compute the average acceptance of pi+ and pi- in phi-plateau 
    We compute this since CLAS12 MC is missing some local ingredients,
    like dead wires in certain sectors and local ineffiencies,
    but we believe that the integral pion detection effieciency is not completely off.
    In addition, we are studying the super-ratio of tagged/untagged pi+/pi- production ratio,
    and the integral height of the plateau, used for acceptance correction,
    will cancel out in this super-ratio 
  
    '''
    global h, h_err, AccCorrecHeight, AccCorrecHeight_err

    var_gen      = 'pi_Phi_g'
    # scale_factor = r2d
    label        = '$\phi$'
    units        = '[deg.]'    


    # (1) compute histograms of generated, reconstructed and accepted events asna function of $\phi$
    fig = plt.figure(figsize=(14,6));
    for pi_ch,pi_charge_label,pi_idx in zip(pi_charge_names,pi_labels,range(2)):
        df_gen = e_e_pi_GEMC[pi_ch];
        df_rec = e_e_pi_GEMC[pi_ch][e_e_pi_GEMC[pi_ch].pi_reconstructed==1];
        df_acc = e_e_pi_GEMC_pass_cuts[pi_ch][e_e_pi_GEMC_pass_cuts[pi_ch].pi_passed_cuts==1];    
        x_gen = df_gen[var_gen]
        x_rec = df_rec[var_gen]    
        x_acc = df_acc[var_gen]    

        N_gen = len(x_gen)/Nphi_pts;
        ax = fig.add_subplot(1,2,pi_idx+1)
        
        # compute (and plot) histograms of generated and accepted (after cuts)
        x,h[pi_ch+'gen'],x_err,h_err[pi_ch+'gen'] = plot_step_hist( x_gen*r2d, phi_bins,     color='k',           label='Generated',             ScaleFactor=100./N_gen)
        x,h[pi_ch+'acc'],x_err,h_err[pi_ch+'acc'] = plot_step_hist( x_acc*r2d, phi_bins,     color='forestgreen', label='Accepted passed cuts',  ScaleFactor=100./N_gen)
        if do_plot_reconstructed: 
            x,h[pi_ch+'rec'],x_err,h_err[pi_ch+'rec'] = plot_step_hist( x_rec*r2d, phi_bins, color='royalblue',   label='Reconsutrcted (no cuts)',  ScaleFactor=100./N_gen)


        # compute (and plot) average height in plateau of generated and accepted (after cuts)
        indices_in_plateau = np.where( h[pi_ch+'acc'] > np.median(h[pi_ch+'acc']) - 0.2*np.std(h[pi_ch+'acc']) )
        average_gen_in_plateau = np.mean( h[pi_ch+'gen'][indices_in_plateau] );
        average_acc_in_plateau = np.mean( h[pi_ch+'acc'][indices_in_plateau] );
        ax.plot( x , average_gen_in_plateau*np.ones(len(x)) , '-', color='k',           label=None , alpha=0.2, linewidth=10)
        ax.plot( x , average_acc_in_plateau*np.ones(len(x)) , '-', color='forestgreen', label=None , alpha=0.2, linewidth=10)

        AccCorrecHeight[pi_ch]     = average_gen_in_plateau/average_acc_in_plateau
        # The uncertainty needs to be rescaed by N_gen since we 'pealed' it off from the histograms
        AccCorrecHeight_err[pi_ch] = AccCorrecHeight[pi_ch] * np.sqrt(1./(N_gen*average_gen_in_plateau) + 1./(N_gen*average_acc_in_plateau))
            
        # plot a text with average acceptance
        plt.text(-120, 1.2*average_acc_in_plateau, 
                 'Average acceptance $%.1f\pm%.1f$'%(100.*1./AccCorrecHeight[pi_ch],
                                                     100.*AccCorrecHeight_err[pi_ch]/(AccCorrecHeight[pi_ch]*AccCorrecHeight[pi_ch]))+'%',
                 color='k',fontsize=20)
        
        # cosmetics
        set_axes(ax,label + ' ' + units,'Fraction of events [%]' if pi_idx==0 else '',
                    do_add_legend=do_add_legend if pi_idx==1 else False,
                     title='$'+pi_charge_label+'$ acceptance as a function of '+label, fontsize=24, 
                 do_add_grid=True,xlim=phi_xlim,xticks=phi_xticks)

    plt.tight_layout();

# ------------------------------------------------------------------------------------------------ #








# ------------------------------------------------------------------------------------------------ #
def load_MC_data_for_acceptance_correction(sim_configuration = 'gcard_default/1M_events'):#{
    global e_e_pi_GEMC
    
    sim_data_path = '/Users/erezcohen/Desktop/data/BAND/AcceptanceCorrection/GEMCimulationOuputFiles/'
    e_e_pi_GEMC['piplus']  = pd.read_csv( sim_data_path + sim_configuration + '/eepips_p_uniform_distribution_1M_events.csv')
    e_e_pi_GEMC['piminus'] = pd.read_csv( sim_data_path + sim_configuration + '/eepims_p_uniform_distribution_1M_events.csv')

    for pi_ch in pi_charge_names:#{
        df = e_e_pi_GEMC[pi_ch]
        print(pi_ch+': %d generated events'%len(df))


    # Omit events in which the electron is not reconstructed
    e_e_pi_GEMC['piplus']  = e_e_pi_GEMC['piplus'] [e_e_pi_GEMC['piplus'] .e_passed_cuts==1]
    e_e_pi_GEMC['piminus'] = e_e_pi_GEMC['piminus'][e_e_pi_GEMC['piminus'].e_passed_cuts==1]

    for pi_ch in pi_charge_names:#{
        print(pi_ch)
        df = e_e_pi_GEMC[pi_ch]
        print('%d events in which electron was reconsutrcted'%len(df[df.e_passed_cuts==1]))
        Nevents = len(df)
        N = len(df[df.pi_reconstructed==1]);
        print('%.1f'%(100.*N/Nevents),'% (','%d events) include pi reconstructed'%N)
        N = len(df[df.pi_passed_fiducial_cuts==1])
        print('%.1f'%(100.*N/Nevents),'% (','%d events) include pi passed fiducial cuts'%N)
        N = len(df[df.pi_passed_PID_cuts==1])
        print('%.1f'%(100.*N/Nevents),'% (','%d events) include pi passed PID cuts'%N)
        N = len(df[df.pi_passed_cuts==1])
        print('%.1f'%(100.*N/Nevents),'% (','%d events) include pi passed all cuts'%N)
    #}
#}
# ------------------------------------------------------------------------------------------------ #





# ------------------------------------------------------------------------------------------------ #
def Compute_acceptance_correction_weight( pi_charge_name, phi ):#{
    '''
    last update May-4, 2022
    '''
    phi_bin = Find_phi_bin( phi )
    acceptance_correction_weight     = AccCorrec[pi_charge_name][phi_bin]
    acceptance_correction_weight_err = AccCorrec_err[pi_charge_name][phi_bin]
    # areaCoverageCorrection           = 1./TightFiducialPhiAreaFraction[pi_charge_name]
    # acceptance_correction_weight     = AccCorrecTightFiducial[pi_charge_name][phi_bin]    
    # acceptance_correction_weight_err = AccCorrecTightFiducial_err[pi_charge_name][phi_bin]
    return acceptance_correction_weight
#}
# ------------------------------------------------------------------------------------------------ #





# ------------------------------------------------------------------------------------------------ #
def Find_phi_bin( phi ):#{
    '''
    input 
    -------
    phi in [deg.]
    
    return
    -------
    phi bin index (/indices)
    '''        
    #
    # non-working version:
    # idx_arr = np.digitize(phi,phi_centers) 
    #
    # working version:
    idx_arr = np.digitize(phi,phi_bins) 
    if np.isscalar(phi):
        if idx_arr > Nphi_pts-1:
            idx_arr = Nphi_pts-1
    else:
        idx_arr[ idx_arr>Nphi_pts-1 ] = Nphi_pts-1
    return idx_arr-1
#}
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



