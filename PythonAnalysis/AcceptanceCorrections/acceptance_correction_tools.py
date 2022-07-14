#import numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib as mpl
#import sys; sys.path.insert(0, '/Users/erezcohen/Desktop/Software/mySoftware/Python/');
#import sys; sys.path.insert(0, '/Users/erezcohen/Desktop/Software/CLAS12/BAND/SIDIS_at_BAND/PythonAnalysis/python_auxiliary/');
#from my_tools               import *;
#from plot_tools             import *;
#from my_data_analysis_tools import *;
#from event_selection_tools  import *;

#
#r2d = 180./np.pi
#
#pi_charge_names  = ['piplus'   ,'piminus'  ]
#pi_labels        = ['\pi^{+}'  ,'\pi^{-}'  ]
#pi_prints        = ['π+'       ,'π-'       ]
#pi_colors        = ['royalblue','salmon'   ]
#
#NphiBins    = 180
#phi_min     = -180 # deg.
#phi_max     = 180  # deg.
#phi_bins    = np.linspace(phi_min, phi_max,NphiBins+1)
#phi_centers = (phi_bins[1:]+phi_bins[:-1])/2
#Nphi_pts    = len(phi_centers)
#phi_xticks  = [-180,-60,60,180]
#phi_xlim    = [-180,180] 
#
## (e,e'π) and (e,e'πn) events 
#global e_e_pi, e_e_pi_n; 
## (e,e'π) and (e,e'πn) events after all selection cuts
## (p-theta acceptance mathcing, Mx)
#global e_e_pi_pass_cuts, e_e_pi_n_pass_cuts;

## (e,e'π) simulated events, generated from uniform electron direction and uniform pion direction and momoentum
#global e_e_pi_GEMC;
#global e_e_pi_GEMC_pass_cuts; 
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





#e_e_pi, e_e_pi_n                                   = dict(),dict()
#e_e_pi_pass_cuts, e_e_pi_n_pass_cuts               = dict(),dict()
#e_e_pi_GEMC                                        = dict()
#e_e_pi_GEMC_pass_cuts                              = dict()
h, h_err                                           = dict(),dict()
AccCorrecHeight, AccCorrecHeight_err               = dict(),dict() 
AccCorrec, AccCorrec_err                           = dict(),dict()
AccCorrecTightFiducial, AccCorrecTightFiducial_err = dict(),dict()
TightFiducialPhi, TightFiducialPhiAreaFraction     = dict(),dict()
idx_good_phi,idx_bad_phi                           = dict(),dict()
fraction_bad_phi                                   = dict()




















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


