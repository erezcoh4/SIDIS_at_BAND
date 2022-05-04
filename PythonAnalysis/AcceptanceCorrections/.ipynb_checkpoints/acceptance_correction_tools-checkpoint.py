import numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib as mpl
import sys; sys.path.insert(0, '/Users/erezcohen/Desktop/Software/mySoftware/Python/'); 
from my_tools               import *; 
from plot_tools             import *;
from my_data_analysis_tools import *;


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
def ComputeAcceptanceCorrectionFromBANDDataRun( run=6421 ):#{
    '''
    last update May-4, 2022
    
    Comments:
    
    Compute the acceptance correction for pi+ and pi- using data from a BAND data run 
  
    input: data run
    
    '''
    global e_e_pi_pass_cuts
    global idx_good_phi, idx_bad_phi
    global fraction_bad_phi
    med_AccCorr_phi, std_AccCorr_phi = dict(),dict()
     
    fig = plt.figure(figsize=(16,6));

    
    for pi_ch,pi_label,pi_color,pi_idx in zip(pi_charge_names,pi_labels,pi_colors,[1,2]):#{
        ax = fig.add_subplot(1,2,pi_idx)

        df  = e_e_pi_pass_cuts[pi_ch];
        df = df[df.runnum==run]
        phi = np.array(df.pi_Phi)*r2d
        x,histo,x_err,histo_err = plot_step_hist( x_arr=phi,  bins=phi_bins , label='Data',  
                                                 color='royalblue', 
                                                 do_plot_errorbar=False, density=False, linewidth=3)
        h['data-%d'%run]     = histo
        h_err['data-%d'%run] = histo_err
        
        # compute average "height" of the data in the plateau region
        # here for the "plateau" we take the median and not median - 0.2std, 
        # since the data is much "sharper" than the simulation
        indices_in_plateau = np.where( histo > np.median(histo) - 0.*np.std(histo) )
        mean_in_plateau = np.mean( histo[indices_in_plateau] );
        ax.plot( x , mean_in_plateau*np.ones(len(x)) , '-', color='royalblue',
                label=None , alpha=0.2, linewidth=10)

        
        # move average "height" of the data in the plateau region to average
        mean_in_plateau_corrected = mean_in_plateau * AccCorrecHeight[pi_ch]
        ax.plot( x , mean_in_plateau_corrected*np.ones(len(x)) , '-', color='k',
                label=None , alpha=0.2, linewidth=10)

        
        # compute acceptance ccorrection weight in each \phi bin
        AccCorrec[pi_ch]     = mean_in_plateau_corrected * (1./ h['data-%d'%run])
        AccCorrec_err[pi_ch] = mean_in_plateau_corrected * (h_err['data-%d'%run] / (h['data-%d'%run]*h['data-%d'%run]))
        
        # print(AccCorrec[pi_ch])
        # print(AccCorrec_err[pi_ch])
        
        # omit regions of the sharp drop in efficiency 
        med_AccCorr_phi[pi_ch] = np.median( AccCorrec[pi_ch] )
        std_AccCorr_phi[pi_ch] = np.std   ( AccCorrec[pi_ch] )
        idx_bad_phi[pi_ch]     = np.where ( AccCorrec[pi_ch] > med_AccCorr_phi[pi_ch] + 0.2*std_AccCorr_phi[pi_ch] )
        idx_bad_phi[pi_ch]     = np.array (idx_bad_phi[pi_ch]).flatten()
        
        for i in idx_bad_phi[pi_ch]: 
            AccCorrec[pi_ch][i]     = 0
            AccCorrec_err[pi_ch][i] = 0
            
        fraction_bad_phi[pi_ch] = float(len(idx_bad_phi[pi_ch])) / Nphi_pts
        
        # and account for this fraction as another scalar factor in acceptance correction
        AccCorrec[pi_ch]     = AccCorrec[pi_ch]     / (1 - fraction_bad_phi[pi_ch])
        AccCorrec_err[pi_ch] = AccCorrec_err[pi_ch] / (1 - fraction_bad_phi[pi_ch])
        
        # plot acceptance correction in right-axis
        axR = ax.twinx()
        # axR.errorbar( x , AccCorrec[pi_ch], AccCorrec_err[pi_ch], marker='.', color='salmon', linecolor=None )
        axR.step( x , AccCorrec[pi_ch], color='salmon')
        
        # indicate "bad" (omitted) regions by shadow
        ymax = 1.*np.max(AccCorrec[pi_ch])        
        for i in idx_bad_phi[pi_ch]:
            phi_min = phi_bins[i]
            phi_max = phi_bins[i+1]
            axR.fill_between( [phi_min,phi_max] , [ymax,ymax] , color='salmon' , alpha=0.2, edgecolor=None )

                    
        # cosmetics
        set_axes(ax, '$\phi$ [deg.]', 'Counts [a.u.]' if pi_idx==1 else '',
                 title="$(e,e'"+pi_label+')$',
                 do_add_grid=True,
                 remove_ticks_y=False,
                 do_add_legend=False if pi_idx==2 else False,
                 xlim=phi_xlim,
                 xticks=phi_xticks)   

        set_axes(axR, '', 'Acceptance correction' if pi_idx==2 else '', remove_ticks_y=True)   

#}
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def compute_ratio_pips_to_pims(df_dict, 
                               var='xB', bins=np.linspace(0,1,10), 
                               weight_option = 'Acc. correction as f(phi)',
                               z_min=0, z_max=1):#{
    '''
    last edit Apr-28, 2022
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
            # no weight, no acceptance corrcetion            
            
            pips_in_bin  = pips[ (x_min < pips) & (pips < x_max) ]
            Npips_in_bin = len(pips_in_bin)
            pims_in_bin  = pims[ (x_min < pims) & (pims < x_max) ]
            Npims_in_bin = len(pims_in_bin)    
            
        #}

        R     = Npips_in_bin / np.max([Npims_in_bin,1])
        R_err = R * np.sqrt( 1./np.max([1,Npips_in_bin]) + 1./np.max([1,Npims_in_bin]) )


        R_pips_to_pims    .append(R)
        R_pips_to_pims_err.append(R_err)
    #}    
    R_pips_to_pims_errup,R_pips_to_pims_errdw = get_err_up_dw(R_pips_to_pims, R_pips_to_pims_err)
    
    return np.array(R_pips_to_pims),np.array(R_pips_to_pims_errup),np.array(R_pips_to_pims_errdw)
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
def apply_further_selection_cuts_to_data(fdebug=2):#{
    '''
    e_e_pi_pass_cuts, e_e_pi_n_pass_cuts, e_e_pi_GEMC_pass_cuts = apply_further_selection_cuts_to_data(fdebug=2)
    
    Apply selection cuts not previously imposed
    
    1. pi+/pi- acceptance matching cut in p-theta plane 
    2. Missing mass cut on (e,e'\pi) events
    
    '''
    global e_e_pi, e_e_pi_n, e_e_pi_GEMC
    global e_e_pi_pass_cuts, e_e_pi_n_pass_cuts, e_e_pi_GEMC_pass_cuts
    
    
    # (e,e'\pi n) SIDIS data complete this -  need to add sector ID in the (e,e'\pi n) data 
    
    # print number of events retained on every cut
    if fdebug < 1: return
    Nevents      = dict()
    frac_Nevents = dict()
    
    # (1) Data
    if (e_e_pi=={}) is False:#{
        
        # (e,e'\pi) SIDIS data
        e_e_pi_after_p_theta_cut = apply_p_theta_acceptance_cut( e_e_pi )
        e_e_pi_after_Mx_cut      = apply_Mx_cut( e_e_pi_after_p_theta_cut )
        e_e_pi_pass_cuts         = e_e_pi_after_Mx_cut;

        for pi_ch in pi_charge_names:#{
            print('(e,e',pi_ch,')')

            Nevents[pi_ch + ' original'] = len(e_e_pi[pi_ch])
            frac_Nevents[pi_ch + ' original'] = 1        
            print(Nevents[pi_ch + ' original'],'events before cut')    

            Nevents[pi_ch +' p-theta cut'] = len(e_e_pi_after_p_theta_cut[pi_ch])
            frac_Nevents[pi_ch + ' p-theta cut'] = float(Nevents[pi_ch +' p-theta cut'])/ Nevents[pi_ch + ' original']
            print(Nevents[pi_ch +' p-theta cut'],'events after p-theta cut (%.1f'%(100.*frac_Nevents[pi_ch + ' p-theta cut']),'%)')    


            Nevents[pi_ch +' Mx cut'] = len(e_e_pi_after_Mx_cut[pi_ch])
            frac_Nevents[pi_ch + ' Mx cut'] = float(Nevents[pi_ch +' Mx cut'])/Nevents[pi_ch + ' original']
            print(Nevents[pi_ch +' Mx cut'],'events after M_X cut (%.1f'%(100.*frac_Nevents[pi_ch + ' Mx cut']),'%)')
        #}
    #}
    
    # (2) MC
    if (e_e_pi_GEMC=={}) is False:#{
        
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
        for pi_ch in pi_charge_names:#{
            print('(e,e',pi_ch,') in uniform GEMC simulation')

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
    #}
    return e_e_pi_pass_cuts, e_e_pi_n_pass_cuts, e_e_pi_GEMC_pass_cuts
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
    global e_e_pi
    global e_e_pi_n;

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
def ComputeAcceptanceCorrectionAsFunctionOfPhi( ):#{
    '''
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
    import numpy as np
    
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