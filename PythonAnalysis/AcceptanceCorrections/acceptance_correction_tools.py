import numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib as mpl
import sys; sys.path.insert(0, '/Users/erezcohen/Desktop/Software/mySoftware/Python/'); 
from my_tools               import *; 
from plot_tools             import *;
from my_data_analysis_tools import *;
# global eepips_GEMC, eepims_GEMC;



r2d = 180./np.pi

pi_charge_names  = ['piplus'   ,'piminus'  ]
pi_labels        = ['\pi^{+}'  ,'\pi^{-}'  ]
pi_colors        = ['royalblue','salmon'   ]

NphiBins    = 180
phi_min     = -180 # deg.
phi_max     = 180  # deg.
phi_bins    = np.linspace(phi_min, phi_max,NphiBins+1)
phi_centers = (phi_bins[1:]+phi_bins[:-1])/2
Nphi_pts    = len(phi_centers)
phi_xticks  = [-180,-60,60,180]
phi_xlim    = [-180,180] 


e_e_pi, e_e_pi_n                                   = dict(),dict()
e_e_pi_GEMC                                        = dict()
h, h_err                                           = dict(),dict()
AccCorrec, AccCorrec_err                           = dict(),dict()
AccCorrecTightFiducial, AccCorrecTightFiducial_err = dict(),dict()
TightFiducialPhi, TightFiducialPhiAreaFraction     = dict(),dict()
# TightFiducialPhi = the are in phi by "good" phi, which we want to keep (0 or 1 for each of phi_centers)
# TightFiducialPhiAreaFraction = the fraction of are occupied by "good" phi in TightFiducialPhi (out of 2\pi)




# ------------------------------------------------------------------------------------------------ #
def load_SIDIS_data(runs_filename = "good_runs_10-2.txt", 
                    main_data_path = '/Users/erezcohen/Desktop/data/BAND/',
                    Nruns = 1, 
                    fdebug = 2):#{
    '''
    Fill e_e_pi and e_e_pi_n with data
    
    e_e_pi & e_e_pi_n are dict() 
    each has 'piplus' and 'piminus' keys,
    and e.g. e_e_pi['piplus'] is a pandas.DataFrame object containing the (e,e'p) events data
    
    '''
    global e_e_pi, e_e_pi_n;

    e_e_pi_data_path   = main_data_path + 'SIDIS_skimming/'
    e_e_pi_n_data_path = main_data_path + 'merged_SIDIS_and_BAND_skimming/'
    runfile_path = "/Users/erezcohen/Desktop/Software/CLAS12/BAND/SIDIS_at_BAND/macros/runlists/";
    
    # Using readlines()
    runs_file     = open( runfile_path + runs_filename, 'r')
    run_fileLines = runs_file.readlines()
    runs = []
    for line in run_fileLines[0:Nruns]:#{
        run = int(line.strip())
        runs.append(run)
    #}
    runs = np.array(runs)        
    
    
    for runnum,runIdx in zip(runs,range(len(runs))):
        print('Run number ',runnum,'(%d/%d runs)'%(runIdx+1,len(runs)))
        for pi_charge_name in pi_charge_names:
            eepi   = pd.read_csv(e_e_pi_data_path                                                
                                 +'skimmed_SIDIS_inc_00%d_e_%s_selected_eepi_kinematics.csv'%
                                 (runnum,pi_charge_name))
            
            eepin = pd.read_csv(e_e_pi_n_data_path                                            
                                +'skimmed_SIDIS_and_BAND_inc_00%d_e_%s_n.csv'%
                                (runnum,pi_charge_name))        


            if fdebug>1: print('Loaded',len(eepi),'(e,e',pi_charge_name,
                               ') events and ',                           
                               len(eepin),'(e,e',pi_charge_name,'n), events')
            
            # large momentum neutrons ( Pn > 275 MeV/c)
            eepin = eepin[eepin['n_P']>0.275]
            if fdebug>1: print('retained',len(eepin),'(e,e',pi_charge_name,'n), events with Pn > 275 MeV/c')            
            
            if runIdx==0:             
                e_e_pi[pi_charge_name]   = eepi
                e_e_pi_n[pi_charge_name] = eepin                
            else:
                e_e_pi[pi_charge_name]   = pd.concat([e_e_pi[pi_charge_name],  eepi ])
                e_e_pi_n[pi_charge_name] = pd.concat([e_e_pi_n[pi_charge_name],eepin])                
    print('Done loading files.')    
    #}
#}
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def load_MC_data_for_acceptance_correction():#{
    global e_e_pi_GEMC
    
    sim_data_path = '/Users/erezcohen/Desktop/data/BAND/AcceptanceCorrection/GEMCimulationOuputFiles/'
    e_e_pi_GEMC['piplus']  = pd.read_csv( sim_data_path + 'eepips_p_uniform_distribution_1M_events.csv')
    e_e_pi_GEMC['piminus'] = pd.read_csv( sim_data_path + 'eepims_p_uniform_distribution_1M_events.csv')

    # Omit events in which the electron is not reconstructed
    print('%d generated events'%len(e_e_pi_GEMC))
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
def ComputeAcceptanceCorrectionAsFunctionOfPhi( max_AccCorrec = 1.5  # omit regions where the acceptance-correction is greater than 50%
                                              ):#{
    global h, h_err, AccCorrec, AccCorrec_err
    global TightFiducialPhi, AccCorrecTightFiducial, AccCorrecTightFiducial_err


    var_gen      = 'pi_Phi_g'
    label        = '$\phi$'
    units        = '[deg.]'
    scale_factor = r2d

    # (1) compute histograms of generated, reconstructed and accepted events asna function of $\phi$
    fig = plt.figure(figsize=(14,6));
    for eepi_GEMC,pi_charge_name,pi_charge_label,pi_idx in zip([eepips_GEMC,eepims_GEMC],pi_charge_names,pi_labels,range(2)):
        pi_ch = pi_charge_name
        df_gen = eepi_GEMC;
        df_rec = eepi_GEMC[eepi_GEMC.pi_reconstructed==1];
        df_acc = eepi_GEMC[eepi_GEMC.pi_passed_cuts  ==1];    
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
    for eepi_GEMC,pi_charge_name,pi_charge_label,pi_idx in zip([eepips_GEMC,eepims_GEMC],pi_charge_names,pi_labels,range(2)):#{
        pi_ch = pi_charge_name
        h[pi_ch+'eff'] = h[pi_ch+'acc']/h[pi_ch+'gen']
        h_err[pi_ch+'eff'] = h[pi_ch+'eff']*np.sqrt( np.square(h_err[pi_ch+'acc']/h[pi_ch+'acc']) + np.square(h_err[pi_ch+'gen']/h[pi_ch+'gen'])  )

        AccCorrec[pi_ch]     = 1./h[pi_ch+'eff']
        AccCorrec_err[pi_ch] = h_err[pi_ch+'eff']/np.square( h[pi_ch+'eff'] )
    #}
    
    # (3) omit regions where the acceptance-correction is greater than 50%
    # as we want to avoid large acceptance correction uncertainty.
    # This effectively means tightening the fiducial cut
    # We do this by assigning an acceptance correction of 0 to these regions
    for eepi_GEMC,pi_ch,pi_charge_label,pi_idx in zip([eepips_GEMC,eepims_GEMC],pi_charge_names,pi_labels,range(2)):#{        
        AccCorrecTightFiducial[pi_ch],AccCorrecTightFiducial_err[pi_ch] = np.zeros(Nphi_pts),np.zeros(Nphi_pts)
        TightFiducialPhi[pi_ch] = np.zeros(Nphi_pts)
        NgoodPhi = 0.0
        for phi_idx in range(Nphi_pts):#{
            if AccCorrec[pi_ch][phi_idx] > max_AccCorrec: #{
                AccCorrecTightFiducial_err[pi_ch][phi_idx] = AccCorrec_err[pi_ch][phi_idx]
                TightFiducialPhi[pi_ch][phi_idx]           = 0
                AccCorrecTightFiducial[pi_ch][phi_idx]     = 0            
            #}
            else :#{
                AccCorrecTightFiducial_err[pi_ch][phi_idx] = 1e-5
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
    for eepi_GEMC,pi_charge_name,pi_charge_label,pi_idx in zip([eepips_GEMC,eepims_GEMC],pi_charge_names,pi_labels,range(2)):
        pi_ch = pi_charge_name
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

    for eepi_GEMC,pi_ch,pi_charge_label,pi_idx in zip([eepips_GEMC,eepims_GEMC],pi_charge_names,pi_labels,range(2)):
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
def Compute_acceptance_correction_weight( pi_charge_name, phi ):#{
    phi_bin = Find_phi_bin( phi )
    # acceptance_correction_weight     = AccCorrec[pi_charge_name][phi_bin]
    # acceptance_correction_weight_err = AccCorrec_err[pi_charge_name][phi_bin]
    # areaCoverageCorrection           = 1./TightFiducialPhiAreaFraction[pi_charge_name]
    acceptance_correction_weight     = AccCorrecTightFiducial[pi_charge_name][phi_bin]    
    acceptance_correction_weight_err = AccCorrecTightFiducial_err[pi_charge_name][phi_bin]
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
    df_dict_after_cut = dict()
    for pi_charge_name in pi_charge_names:#{
        df = df_dict[pi_charge_name]
        df = df[ (Mx_min < df.M_X) & (df.M_X < Mx_max)]    
        df_dict_after_cut[pi_charge_name] = df
    #}
    return df_dict_after_cut
#}
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def apply_p_theta_acceptance_cut( df_dict=None ):    
    '''
        df_dict_after_cut = apply_p_theta_acceptance_cut(df_dict)
        
        Apply a pi+/pi- acceptance matching cut on the in p-\theta plane
        
    '''    
    import numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib as mpl
    
    df_dict_after_cut = dict()
    for pi_charge_name in pi_charge_names:
        good_indices = np.array([])
        df = df_dict[pi_charge_name]
        
        for sector in range(1,7):#{
            df_in_sector = df[df.pi_DC_sector == sector]
            pi_P      = np.array(df_in_sector.pi_P)
            pi_Theta  = np.array(df_in_sector.pi_Theta)*r2d
            theta_min_pips = pi_min_theta_cut( pi_charge = 'pips', sector=sector, p=pi_P )
            theta_min_pims = pi_min_theta_cut( pi_charge = 'pims', sector=sector, p=pi_P )
            
            # print(len(df),sector,len(df_in_sector))
            good_indices_in_sector = []
            for evt_idx,step_idx in zip(df_in_sector.index,range(len(df_in_sector))):#{
                # print( evt_idx, pi_P[step_idx], pi_Theta[step_idx] , theta_min[step_idx] )
                if (pi_Theta[step_idx] > theta_min_pips[step_idx]) and (pi_Theta[step_idx] > theta_min_pims[step_idx]):#{
                    good_indices_in_sector.append(evt_idx)
                #}
            #}
            # print (pi_charge_name,sector,good_indices_in_sector)
            good_indices = np.concatenate([good_indices,good_indices_in_sector])
        #}
        # print (pi_charge_name,good_indices)
        df = df.loc[good_indices]
        # print ('defined df of good indices in ', pi_charge_name)
        df_dict_after_cut[pi_charge_name] = df
        
    return df_dict_after_cut
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def pi_min_theta_cut( pi_charge = 'pips', sector=1, p=2 ):
    '''
    p-theta cut optimized by Alex K., Feb-2022
    
    input:
    --------
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

    if pi_charge=='pips' or pi_charge=='piplus':
        parameters = pips_parameters
            
    elif pi_charge=='pims'or pi_charge=='piminus':
        parameters = pims_parameters

    theta_min = parameters[sector-1][0] + parameters[sector-1][1]/p
    return theta_min
# ------------------------------------------------------------------------------------------------ #