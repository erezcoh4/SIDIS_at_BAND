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

# d(e,e'π), d(e,e'πn) and  p(e,e'π) events before all selection cuts
global e_e_pi          , e_e_pi_n,          e_e_pi_FreeP;

# d(e,e'π), d(e,e'πn) and  p(e,e'π) events after all selection cuts
global e_e_pi_pass_cuts, e_e_pi_n_pass_cuts,e_e_pi_FreeP_pass_cuts;

# (e,e'π) simulated events, generated from uniform electron direction and uniform pion direction and momoentum
global e_e_pi_GEMC     , e_e_pi_GEMC_pass_cuts;


e_e_pi, e_e_pi_n, e_e_pi_FreeP                               = dict(),dict(),dict()
e_e_pi_pass_cuts, e_e_pi_n_pass_cuts, e_e_pi_FreeP_pass_cuts = dict(),dict(),dict()
e_e_pi_GEMC                                                  = dict()
e_e_pi_GEMC_pass_cuts                                        = dict()



# Compute weight to each run of 1/beam-charge
# such that events are counted in 1/nC
beam_charge_all_runs = pd.read_csv('/Users/erezcohen/Desktop/data/BAND/metaData/beam_charge_all_runs.csv');
weight_per_run = np.zeros(np.max(beam_charge_all_runs.runnum)+1)
for run in beam_charge_all_runs.runnum:
    weight_per_run[run] = 1./float(beam_charge_all_runs[beam_charge_all_runs.runnum==run].beam_charge);



# ----------------------- #
def extract_SIDIS_ratio(df_dict  = None,
                        specific_run_number=None,
                        x_var    = 'xB' ,
                        x_bins   = np.linspace(0.2,0.6,11),
                        zvar     = "Zpi",
                        z_bins   = np.arange(0.3,0.8,0.1),
                        z_widths = 0.01*np.ones(5),
                        data_path= '/Users/erezcohen/Desktop/data/BAND/Results/',
                        fdebug   = 0,
                        prefix   = 'Untagged_SIDIS_ratio_',
                        suffix   = '',
                        M_x_min  = 0,
                        M_x_max  = np.inf,
                        W_min    = 0,
                        W_max    = np.inf,
                        Q2_min   = 0,
                        Q2_max   = np.inf,
                        Mx_d_min = 0,
                        weight_option = 'beam-charge',):
    '''
    Extract SIDIS results,
    the number of d(e,e'π+) and d(e,e'π-) events,
    and save them to a CSV file
    
    last update Nov-26, 2022
    
    
    input
    ---------
    zvar            "Zpi" / "zeta_pi"
    
    '''
    
    x        = (x_bins[1:] + x_bins[:-1])/2
    x_err    = (x_bins[1:] - x_bins[:-1])/2
    for z_bin,z_width in zip(z_bins,z_widths):
        z_min,z_max = z_bin-z_width, z_bin+z_width
        
        (R,
         R_err_up,
         R_err_dw,
         N_pips,
         N_pims,
         Zavg_pips,
         Zavg_pims,
         N_pips_err,
         N_pims_err) = compute_ratio_pips_to_pims(df_dict = df_dict ,
                                                 specific_run_number=specific_run_number,
                                                 var     = x_var,
                                                 bins    = x_bins,
                                                  zvar   = zvar,
                                                 z_min   = z_min,
                                                 z_max   = z_max,
                                                 M_x_min = M_x_min,
                                                 M_x_max = M_x_max,
                                                 W_min   = W_min,
                                                 W_max   = W_max,
                                                 Mx_d_min= Mx_d_min,
                                                  Q2_min = Q2_min,
                                                  Q2_max = Q2_max,
                                                  weight_option=weight_option,
                                                 fdebug  = fdebug)

        df_to_save = pd.DataFrame({"$x_B$":x,
                                   "$\Delta x_B$":x_err,
                                   '$N(\pi_{+})$':N_pips,
                                   '$N(\pi_{-})$':N_pims,
                                   '$\Delta N(\pi_{+})$':N_pips_err,
                                   '$\Delta N(\pi_{-})$':N_pims_err,
                                   '$R$':R,
                                   '$\Delta R_{+}$':R_err_up,
                                   '$\Delta R_{-}$':R_err_dw})
        
        filelabel = '%s_min%.3f_%s_mean_pips%.3f_pims%.3f_%s_max%.3f'%(zvar,z_min,zvar,Zavg_pips,Zavg_pims,zvar,z_max)
        filename  =  data_path + prefix + filelabel + suffix  + '.csv'
        df_to_save.to_csv(filename)
        if fdebug:
            print('saved',filename)
            if fdebug>1:
                print('$%s=%.3f\pm%.3f$'%(zvar,z_bin,z_width))
                print(filename)
                if fdebug>2: display(df_to_save)
# ----------------------- #





# ----------------------- #
def compute_ratio_pips_to_pims(df_dict,
                               specific_run_number=None,
                               var='xB',
                               bins=np.linspace(0,1,10),
                               weight_option = 'beam-charge',
                               zvar="Zpi",
                               z_min=0,    z_max=1,
                               M_x_min=0,  M_x_max=np.inf,
                               W_min=0,    W_max=np.inf,
                               Q2_min = 0, Q2_max= np.inf,
                               Mx_d_min=0, fdebug=0,
                               cutoff = 1.e-8 ):#{
    '''
    last edit Oct-19, 2022
    
    [R_pips_to_pims, R_pips_to_pims_errup, R_pips_to_pims_errdw,
     N_pips, N_pims,
     Zavg_pips, Zavg_pims,
     N_pips_err, N_pims_err] = compute_ratio_pips_to_pims(df_dict,
                               specific_run_number=None,
                               var='xB',
                               bins=np.linspace(0,1,10),
                               weight_option = 'beam-charge',
                               z_min=0,   z_max=1,
                               M_x_min=0, M_x_max=np.inf,
                               W_min=0,   W_max=np.inf,
                               Q2_min = 0, Q2_max= np.inf,
                               Mx_d_min=0, fdebug=0 )
    
    
    input:
    -------
    M_x_min               float         minimal M_x
    M_x_max               float         maximal M_x
    
    
    return:
    -------
    R_pips_to_pims                np.array()   number of π+ events in each x-bin / number of π-
    R_pips_to_pims_errup          np.array()   err-up in number of π+ events in each x-bin / number of π-
    R_pips_to_pims_errdw          np.array()   err-dw in number of π+ events in each x-bin / number of π-

    R_normed_pips_to_pims         np.array()   R_pips_to_pims normalized by beam-charge
    R_normed_pips_to_pims_errup   np.array()   R_pips_to_pims_errup normalized by beam-charge
    R_normed_pips_to_pims_errdw   np.array()   R_pips_to_pims_errdw normalized by beam-charge

    N_pips                        np.array()   number of π+ events in each x-bin
    N_pims                        np.array()   number of π- events in each x-bin
    Zavg_pips                     float        mean z-value in the range z_min < z < z_max for π+
    Zavg_pims                     float        mean z-value in the range z_min < z < z_max for π-
    N_pips_err                    np.array()   uncertainty number of π+ events in each x-bin
    N_pims_err                    np.array()   uncertainty in the number of π- events in each x-bin

        
    comments:
    -------
    weight_option: None, 'beam-charge', 'Acc. correction as f(phi)'
    if weight is 'beam-charge', the results are given in number of events per nC
    
    '''
    # z_min,z_max are z limits on the pion outgoing momentum
    df_pips = df_dict['piplus']
    df_pims = df_dict['piminus']
    
    # cut on z and other variables for pi+
    df_pips = df_pips[  (z_min   < df_pips[zvar]) & (df_pips[zvar] < z_max  )
                      & (W_min   < df_pips.W  )   & (df_pips.W   < W_max  )   ]
    
    if 0 < M_x_min or M_x_max < np.inf:
        df_pips = df_pips[ (M_x_min < df_pips.M_x) & (df_pips.M_x < M_x_max) ]
        
    if 0 < Mx_d_min:
        df_pips = df_pips[ Mx_d_min < df_pips.M_x_d ]
        
    if 0 < Q2_min or Q2_max < np.inf:
        df_pips = df_pips[ (Q2_min < df_pips.Q2) & (df_pips.Q2 < Q2_max) ]

        
    if specific_run_number is not None:
        if fdebug>1:  print('df_pips: before run %d filter: %d events'%(specific_run_number,len(df_pips)))
        df_pips = df_pips[df_pips.runnum == specific_run_number]
        if fdebug>1:  print('after run %d filter: %d events'%(specific_run_number,len(df_pips)))
                           
    Zavg_pips = np.mean( np.array(df_pips[zvar])  )
    
    # cut on z and other variables for pi-
    df_pims = df_pims[  (z_min   < df_pims[zvar]) & (df_pims[zvar] < z_max  )
                      & (W_min   < df_pims.W  )   & (df_pims.W   < W_max  )   ]
    
    if 0 < M_x_min or M_x_max < np.inf:
        df_pims = df_pims[ (M_x_min < df_pims.M_x) & (df_pims.M_x < M_x_max) ]
    
    if 0 < Mx_d_min:
        df_pims = df_pims[Mx_d_min < df_pims.M_x_d]
        
    if 0 < Q2_min or Q2_max < np.inf:
        df_pims = df_pims[ (Q2_min < df_pims.Q2) & (df_pims.Q2 < Q2_max) ]

        
    if specific_run_number is not None:
        df_pims = df_pims[df_pims.runnum == specific_run_number]

    Zavg_pims = np.mean( np.array(df_pims[zvar])  )


    pips = df_pips[var]
    pims = df_pims[var]
    if weight_option == 'Acc. correction as f(phi)':#{
        phi_pips = np.array( df_pips.pi_Phi )*r2d
        phi_pims = np.array( df_pims.pi_Phi )*r2d
        
    elif weight_option == 'beam-charge':#{
        w_pips = np.array( df_pips.weight )
        w_pims = np.array( df_pims.weight )
    #}
        
    R_pips_to_pims, R_pips_to_pims_err = [],[]
    N_pips, N_pims                     = [],[]
    N_pips_err, N_pims_err             = [],[]
    if fdebug>1: print('%.2f<%s<%.2f: '%(z_min,zvar,z_max), len(pips),'π+ and',len(pims),'π-')
        

    for x_min,x_max in zip(bins[:-1],bins[1:]):#{
        
        if fdebug>1: print('\t%.2f<%s<%.2f:'%(x_min,var,x_max),
                           len(pips[ (x_min < pips) & (pips < x_max) ]),'π+ and',
                           len(pims[ (x_min < pims) & (pims < x_max) ]),'π-')
        
        if weight_option == 'Acc. correction as f(phi)':#{
            # each event is weighted by the acceptance correction weight
            phi_pips_in_bin = phi_pips[ (x_min < pips) & (pips < x_max) ]
            W_pips_in_bin   = [ Compute_acceptance_correction_weight( 'piplus' , phi ) for phi in phi_pips_in_bin ]
            Npips_in_bin    = np.sum( W_pips_in_bin )
            
            phi_pims_in_bin = phi_pims[ (x_min < pims) & (pims < x_max) ]
            W_pims_in_bin   = [ Compute_acceptance_correction_weight( 'piminus', phi ) for phi in phi_pims_in_bin ]
            Npims_in_bin    = np.sum( W_pims_in_bin )
            
        elif weight_option == 'beam-charge':#{
            
            # weighted by 1/beam-charge
            W_pips_in_bin   = w_pips[ (x_min < pips) & (pips < x_max) ]
            Npips_in_bin    = np.sum( W_pips_in_bin )
            Npips_in_bin_err= np.sqrt(np.sum( np.square(W_pips_in_bin )))
            
            W_pims_in_bin   = w_pims[ (x_min < pims) & (pims < x_max) ]
            Npims_in_bin    = np.sum( W_pims_in_bin )
            Npims_in_bin_err= np.sqrt(np.sum( np.square(W_pims_in_bin )))
            
            # cutoff = np.min(weight_per_run)
            # print(cutoff,np.max([Npims_in_bin,cutoff]))
            if Npims_in_bin < cutoff:
                R = 0
            else:
                R = Npips_in_bin / np.max([Npims_in_bin,cutoff])
            R_err = R * np.sqrt( np.square(Npips_in_bin_err/np.max([Npips_in_bin,cutoff]))
                                + np.square(Npims_in_bin_err/np.max([Npims_in_bin,cutoff]) ) )

           
        else:
            # no weight, no acceptance correction
            pips_in_bin      = pips[ (x_min < pips) & (pips < x_max) ]
            Npips_in_bin     = len(pips_in_bin)
            Npips_in_bin_err = sqrt(len(pips_in_bin))
            
            pims_in_bin      = pims[ (x_min < pims) & (pims < x_max) ]
            Npims_in_bin     = len(pims_in_bin)
            Npims_in_bin_err = sqrt(len(pims_in_bin))

            
            R     = Npips_in_bin / np.max([Npims_in_bin,1])
            R_err = R * np.sqrt( 1./np.max([1,Npips_in_bin]) + 1./np.max([1,Npims_in_bin]) )
        #}


        N_pips            .append(Npips_in_bin)
        N_pips_err        .append(Npips_in_bin_err)

        N_pims            .append(Npims_in_bin)
        N_pims_err        .append(Npims_in_bin_err)

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
            Zavg_pims,
            np.array(N_pips_err),
            np.array(N_pims_err)]
#}
# ----------------------- #









# ----------------------- #
def load_SIDIS_ratio(xlabel   = "Bjorken $x$",
                     fdebug   = 0,
                     prefix   = 'Untagged_SIDIS_ratio_',
                     suffix   = '',
                     doPlotResults=False,
                     axPlot   = [],
                     Nzbins2Plot = 5,
                     titlePlot="$\pi^+/\pi^-$ ratio as a function of $x_B$ without a tagged neutron",
                     zvar     = "Zpi",
                     data_path= '/Users/erezcohen/Desktop/data/BAND/Results/'):
    '''
    Load SIDIS ratio results
    last update Jan-2, 2023
    
    input
    -------
    zvar      "Zpi","zeta_pi"
    
    '''
    
    SIDIS_results = dict()
    
    z_min_arr, z_max_arr, Zavg_pips_arr, Zavg_pims_arr = [],[],[],[]
    if fdebug: print('Reading files from ' + data_path)
    filelist = os.listdir(data_path)
    for filename in filelist:
        if prefix in filename and suffix in filename:
            filenameparts = filename.split('_')
            if fdebug>2: print(filename)
            if fdebug>3: print(filenameparts)
            if zvar=="Zpi":
                if "Zpi" not in filenameparts: continue
                z_min     = float(filenameparts[4][4:8])
                # Zavg_pips = float(filenameparts[7][5:9])
                Zavg_pips_str = filenameparts[7][5:9]
                if 'an' in Zavg_pips_str: 
                    Zavg_pips = 0.000;
                else:                     Zavg_pips = float(Zavg_pips_str)
                # Zavg_pims = float(filenameparts[8][5:9])
                Zavg_pims_str = filenameparts[8][5:9]
                if 'an' in Zavg_pims_str: 
                    Zavg_pims = 0.000;
                else:                     Zavg_pims = float(Zavg_pims_str) 
                z_max     = float(filenameparts[10][4:8])

            elif zvar=="zeta_pi":
                if "zeta" not in filenameparts: continue
                z_min     = float(filenameparts[5][4:8])
                Zavg_pips = float(filenameparts[9][5:9])
                Zavg_pims = float(filenameparts[10][5:9])
                z_max     = float(filenameparts[13][4:8])

            filelabel = '%s_min%.3f_%s_mean_pips%.3f_pims%.3f_%s_max%.3f'%(zvar,z_min,zvar,Zavg_pips,Zavg_pims,zvar,z_max)

            df = pd.read_csv( data_path + '/' + filename )
            SIDIS_results[filelabel] = df
            z_min_arr.append(z_min)
            z_max_arr.append(z_max)
            Zavg_pips_arr.append(Zavg_pips)
            Zavg_pims_arr.append(Zavg_pims)

        
    z_min_arr = np.sort(z_min_arr)
    z_max_arr = np.sort(z_max_arr)
    Zavg_pips_arr = np.sort(Zavg_pips_arr)
    Zavg_pims_arr = np.sort(Zavg_pims_arr)
    
    if doPlotResults:#{
        x     = np.array(df["$x_B$"])
        x_err = np.array(df["$\Delta x_B$"])

        if axPlot==[]:
            fig = plt.figure(figsize=(9,6))
            axPlot  = fig.add_subplot(1,1,1)
            
        for z_min,z_max,Zavg_pips,Zavg_pims in zip( z_min_arr, z_max_arr, Zavg_pips_arr, Zavg_pims_arr[0:Nzbins2Plot] ):
            Zavg = (Zavg_pips+Zavg_pims)/2.
            filelabel = '%s_min%.3f_%s_mean_pips%.3f_pims%.3f_%s_max%.3f'%(zvar,z_min,zvar,Zavg_pips,Zavg_pims,zvar,z_max)
            
            df = SIDIS_results[filelabel]
            y    = df['$R$']
            y_err= (df['$\Delta R_{+}$'],df['$\Delta R_{-}$'])
            # plot
            l=axPlot.errorbar(x=x, xerr=x_err,  y=y, yerr=y_err,
                        marker='o',markeredgecolor='k',
                        label='$%.3f<z<%.3f, \\bar{z}=%.3f$'%(z_min,z_max,Zavg))

        set_axes(axPlot,xlabel,"$N(e,e'\pi^+)/N(e,e'\pi^-)$",
                 title=titlePlot,
                 do_add_grid=True, do_add_legend=True,fontsize=18);
        # plt.legend(bbox_to_anchor=(1,1.05),loc='best',fontsize=18)
        
    if fdebug: print('Done.')
    return SIDIS_results
# ----------------------- #




# NEED TO CHANGE THE "NAN" to something like 0.000 with 4 digits





# ----------------------- #
def apply_p_theta_acceptance_cut_single_set( df_dict=None,
                                 NeventsMax=-1,
                                 fdebug=1):
    '''
        df_dict_after_cut = apply_p_theta_acceptance_cut(df_dict)
        
        Apply a π+/π- acceptance matching cut on the in p-\theta plane
        last update Sep-22, 2022
        
    '''
    import numpy as np
    print("Apply a π+/π- acceptance matching cut on the in p-theta plane")
    df_dict_after_cut = dict()
        
    for pi_ch in pi_charge_names:
        if fdebug: print('Applying p-theta on cut for '+pi_ch+' which includes %d events'%len(df_dict[pi_ch]))
        if NeventsMax > 0: NeventsMax = np.min( [NeventsMax, len(df_dict[pi_ch])] )
        else:              NeventsMax = len(df_dict[pi_ch])
        df = df_dict[pi_ch][0:NeventsMax]
        if fdebug: print('Applying p-theta on cut for '+pi_ch+' on %d events'%NeventsMax)
        # good_indices = np.array([])
        df_after_cut = pd.DataFrame();
        for sector in range(1,7):#{
            df_in_sector   = df[df.pi_DC_sector == sector]
            if fdebug: print(len(df_in_sector),'in sector',sector)
            theta_min_pi   = pi_min_theta_cut( pi_charge = 'any', sector=sector, p=np.array(df_in_sector.pi_P) )
            df_in_sector_pass_cut = df_in_sector[ df_in_sector.pi_Theta*r2d > theta_min_pi ]
            df_after_cut = df_after_cut.append(df_in_sector_pass_cut);
        #}
        df_dict_after_cut[pi_ch] = df_after_cut
    return df_dict_after_cut
# ----------------------- #

#
#
## ----------------------- #
#def apply_p_theta_acceptance_cut_single_set( df_dict=None,
#                                 NeventsMax=-1,
#                                 NMaxPerSubset = 500000,
#                                 fdebug=1):
#    '''
#        df_dict_after_cut = apply_p_theta_acceptance_cut(df_dict)
#
#        Apply a π+/π- acceptance matching cut on the in p-\theta plane
#        last update Sep-22, 2022
#
#    '''
#    import numpy as np
#    print("Apply a π+/π- acceptance matching cut on the in p-theta plane")
#    df_dict_after_cut = dict()
#
#    for pi_ch in pi_charge_names:
#
#        if NeventsMax > 0: NeventsMax = np.min( [NeventsMax, len(df_dict[pi_ch])] )
#        else:              NeventsMax = len(df_dict[pi_ch])
#        df = df_dict[pi_ch][0:NeventsMax]
#        if fdebug: print('Applying p-theta on cut for '+pi_ch+' on %d events'%NeventsMax)
#        good_indices = np.array([])
#        for sector in range(1,7):#{
#            df_in_sector   = df[df.pi_DC_sector == sector]
#            if fdebug: print(len(df_in_sector),'in sector',sector)
#            theta_min_pi   = pi_min_theta_cut( pi_charge = 'any', sector=sector, p=np.array(df_in_sector.pi_P) )
#            good_indices_in_sector = []
#            good_indices_in_sector = df_in_sector[ df_in_sector.pi_Theta*r2d > theta_min_pi ].index
#            good_indices = np.concatenate([good_indices,good_indices_in_sector])
#        #}
#        good_indices = np.unique(good_indices)
#        df_after_cut = df.loc[good_indices]
#
#        df_dict_after_cut[pi_ch] = df_after_cut
#    return df_dict_after_cut
## ----------------------- #


# ----------------------- #
def load_SIDIS_data(runs_filename   = "good_runs_10-2-final.txt",
                    main_data_path  = '/Users/erezcohen/Desktop/data/BAND/',
                    Nruns           = 1,
                    do_e_e_pi       = True,
                    do_e_e_pi_n     = True,
                    do_e_e_pi_FreeP = True,
                    do_all_vars     = False,
                    fdebug          = 2,
                    prefix          = "sidisdvcs",
                    subdirname      = "",
                    taggedsubdirname= "",
                    FreeP_prefix    = "ntupleNew"):#{
    '''
    e_e_pi, e_e_pi_n, e_e_pi_FreeP = load_SIDIS_data()
    Load SIDIS data, and fill e_e_pi and e_e_pi_n with data
    last update Nov-25, 2022
    
    input:
    -------------
    do_e_e_pi       flag to read d(e,e'π) data  from RGB - takes much time for a large number of runs
    do_e_e_pi_n     flag to read d(e,e'πn) data from RGB - takes less time
    do_e_e_pi_FreeP flag to read p(e,e'π) data from RGA  - takes much time
    prefix          "sidisdvcs" / "inc"      - inclusive skimming train
    subdirname      "With_W0.5cut" / "With_W2.5cut"
    
    
    Comments:
    -------------
    e_e_pi, e_e_pi_n, e_e_pi_FreeP       dict(['piplus','piminus'])
    e.g. :
    e_e_pi['piplus'] = pandas.DataFrame( (e,e'π) events data )
    
    '''
    global e_e_pi, e_e_pi_n, e_e_pi_FreeP;

    if taggedsubdirname=="": taggedsubdirname = subdirname;
    e_e_pi_data_path       = main_data_path + 'SIDIS_skimming/' + prefix + '/' + subdirname + '/'
    e_e_pi_n_data_path     = main_data_path + 'merged_SIDIS_and_BAND_skimming/' + prefix + '/' + taggedsubdirname + '/'
    e_e_pi_FreeP_data_path = main_data_path + 'RGA_Free_proton/'

    runs = read_run_nunmbers(runs_filename=runs_filename,Nruns=Nruns)
    e_e_pi, e_e_pi_n, e_e_pi_FreeP = dict(),dict(),dict()
    
    for runnum,runIdx in zip(runs,range(len(runs))):#{
        if fdebug>1: print('Run number ',runnum,'(%d/%d runs)'%(runIdx+1,len(runs)))
        for pi_charge_name,pi_print in zip(pi_charge_names,pi_prints):
            if do_e_e_pi:#{
                # e_e_pi[pi_charge_name] = []
                if do_all_vars:
                    eepi   = pd.read_csv(e_e_pi_data_path
                                     +'skimmed_SIDIS_'
                                     +prefix + '_'
                                     +'00%d_e_%s_selected_eepi_kinematics.csv'%(runnum,pi_charge_name))

                else: # more economic
                    eepi   = pd.read_csv(e_e_pi_data_path
                                     +'skimmed_SIDIS_'
                                     +prefix + '_'
                                     +'00%d_e_%s_selected_eepi_kinematics.csv'%(runnum,pi_charge_name),
                                     usecols=['runnum','evnum',
                                              'e_P','e_Theta','e_Phi',
                                              'pi_P', 'pi_Theta', 'pi_Phi',
                                              'Q2', 'W',
                                              'xB', 'Zpi',
                                              'M_x', 'e_DC_sector',
                                              'pi_DC_sector','pi_qFrame_pT','pi_qFrame_pL'],
                                     dtype={'runnum':int,'evnum': int,
                                            'e_DC_sector':int, 'pi_DC_sector':int,
                                            'e_P':np.half,'e_Theta':np.half,'e_Phi':np.half,
                                            'pi_P':np.half,'pi_Theta':np.half, 'pi_Phi':np.half,
                                            'Q2':np.half,  'W':np.half,
                                            'xB':np.half, 'Zpi':np.half,
                                            'M_x':np.half,
                                            'pi_qFrame_pT':np.half,'pi_qFrame_pL':np.half})
                
                if runIdx==0: e_e_pi[pi_charge_name] = eepi
                else:         e_e_pi[pi_charge_name] = pd.concat([e_e_pi[pi_charge_name],eepi])

                
                if fdebug>1: print('Loaded',len(eepi)," d(e,e'"+pi_print+") events")

            #}
            if do_e_e_pi_n:#{
                eepin = pd.read_csv(e_e_pi_n_data_path
                                    + 'skimmed_SIDIS_and_BAND_'
                                    + prefix + '_'
                                    + '00%d_e_%s_n.csv'%(runnum,pi_charge_name))
                
                if fdebug>1: print('Loaded',len(eepin)," d(e,e'"+pi_print+"n) events")

                if runIdx==0: e_e_pi_n[pi_charge_name] = eepin
                else:         e_e_pi_n[pi_charge_name] = pd.concat([e_e_pi_n[pi_charge_name],eepin])
            #}
        #}
    #}
    if do_e_e_pi_FreeP:#{
        for pi_charge_name,pi_print in zip(pi_charge_names,pi_prints):
            if do_all_vars:
                eepi   = pd.read_csv(e_e_pi_FreeP_data_path
                                 +FreeP_prefix
                                 +'_e_e_%s_selected_eepi_kinematics.csv'%(pi_charge_name))

            else: # more economic
                eepi   = pd.read_csv(e_e_pi_FreeP_data_path
                                 +FreeP_prefix
                                 +'_e_e_%s_selected_eepi_kinematics.csv'%(pi_charge_name),
                                 usecols=['runnum','evnum',
                                          'e_P','e_Theta','e_Phi',
                                          'pi_P', 'pi_Theta', 'pi_Phi',
                                          'Q2', 'W',
                                          'xB', 'Zpi',
                                          'M_x', 'e_DC_sector', 'pi_DC_sector','pi_qFrame_pT','pi_qFrame_pL'],
                                 dtype={'runnum':int,'evnum': int,
                                        'e_DC_sector':int, 'pi_DC_sector':int,
                                        'e_P':np.half,'e_Theta':np.half,'e_Phi':np.half,
                                        'pi_P':np.half,'pi_Theta':np.half, 'pi_Phi':np.half,
                                        'Q2':np.half,  'W':np.half,
                                        'xB':np.half, 'Zpi':np.half,
                                        'M_x':np.half,
                                        'pi_qFrame_pT':np.half,'pi_qFrame_pL':np.half})
            
            # Aug-2022: in e_e_pi_FreeP we only have 1 data-file
            e_e_pi_FreeP[pi_charge_name] = eepi
            #            if runIdx==0: e_e_pi_FreeP[pi_charge_name] = eepi
            #            else:         e_e_pi_FreeP[pi_charge_name] = pd.concat([e_e_pi_FreeP[pi_charge_name],eepi])
            if fdebug>1: print('Loaded',len(eepi)," p(e,e'"+pi_print+") events")
        #}
    #}
       
    print('Done loading files.')
    
    if fdebug>0:
        print('')
        print('Total statistics:')
        for pi_charge_name,pi_print in zip(pi_charge_names,pi_prints):
            if do_e_e_pi:       print(len(e_e_pi[pi_charge_name])      ," d(e,e'"+pi_print+")  events")
            if do_e_e_pi_n:     print(len(e_e_pi_n[pi_charge_name])    ," d(e,e'"+pi_print+"n) events")
            if do_e_e_pi_FreeP: print(len(e_e_pi_FreeP[pi_charge_name])," p(e,e'"+pi_print+")  events")
    #}
    return e_e_pi, e_e_pi_n, e_e_pi_FreeP
#}
# ----------------------- #











# ----------------------- #
def runnum_weight( runnumbers ):
    return weight_per_run[runnumbers]
# ----------------------- #


# ----------------------- #
def apply_cuts_to_e_e_pi_n(fdebug=2,
                             NeventsMax=-1,
                             NMaxPerSubset = 500000,
                             doAcceptanceMatchingCut = True,
                             doApply_minPn_cut       = True,
                             doApply_Mx_cut          = True):#{
    '''
    e_e_pi_n_pass_cuts = apply_cuts_to_e_e_pi_n(fdebug,
                                         NeventsMax,
                                         NMaxPerSubset,
                                         doAcceptanceMatchingCut,
                                         doApply_minPn_cut,
                                         doApply_Mx_cut)
                                         
    Aug-29, 2022
    
    Tagged d(e,e'\pi) SIDIS data
    '''
    global e_e_pi_n, e_e_pi_n_pass_cuts


    if doApply_minPn_cut:#{
        e_e_pi_n_after_minPn_cut   = apply_minPn_cut( e_e_pi_n )
    else:
        e_e_pi_n_after_minPn_cut = dict()
        for pi_ch,pi_print in zip(pi_charge_names,pi_prints):#{
            df = e_e_pi_n[pi_ch]
            e_e_pi_n_after_minPn_cut[pi_ch] = df[0:NeventsMax]
        #}
    #}


    if doAcceptanceMatchingCut:#{
        e_e_pi_n_after_p_theta_cut = apply_p_theta_acceptance_cut_single_set( e_e_pi_n_after_minPn_cut,
                                                                  NeventsMax=NeventsMax,
#                                                                  NMaxPerSubset=NMaxPerSubset,
                                                                  fdebug=fdebug )
    else:
        e_e_pi_n_after_p_theta_cut = e_e_pi_n_after_minPn_cut
    #}

    if doApply_Mx_cut:  #{
        e_e_pi_n_after_Mx_cut      = apply_Mx_cut( e_e_pi_n_after_p_theta_cut )
    else:
        e_e_pi_n_after_Mx_cut      = e_e_pi_n_after_p_theta_cut
    #}

    e_e_pi_n_after_Kinematical_cuts = apply_Kinematical_cuts( e_e_pi_n_after_Mx_cut )
    e_e_pi_n_pass_cuts  = e_e_pi_n_after_Kinematical_cuts;
    # add beam-charge weight
    for pi_ch,pi_print in zip(pi_charge_names,pi_prints):#{
        runnumbers = np.array(e_e_pi_n_pass_cuts[pi_ch].runnum).astype(int);
        e_e_pi_n_pass_cuts[pi_ch]['weight'] = runnum_weight( runnumbers )
    #}


    Nevents      = dict()
    frac_Nevents = dict()
    for pi_ch,pi_print in zip(pi_charge_names,pi_prints):#{
        print('(e,e',pi_print,')')
        Nevents,frac_Nevents = dict(),dict()
        if NeventsMax < 0:
            Nevents,frac_Nevents = get_Nevents(pi_ch, 'original',  e_e_pi_n, Nevents, frac_Nevents);
        else:
            Nevents[pi_ch + 'original cut'] = NeventsMax
            frac_Nevents[pi_ch + 'original cut'] = 1

        Nevents,frac_Nevents = get_Nevents(pi_ch, 'p-theta',      e_e_pi_n_after_p_theta_cut,      Nevents, frac_Nevents);
        Nevents,frac_Nevents = get_Nevents(pi_ch, 'Kinematical',  e_e_pi_n_after_Kinematical_cuts, Nevents, frac_Nevents);
    #}
    print(' ')
    return e_e_pi_n_pass_cuts
#}


# ----------------------- #
def apply_cuts_to_e_e_pi(fdebug=0,
                         NeventsMax=-1,
                         NMaxPerSubset = 500000,
                         doAcceptanceMatchingCut = True,
                         doApply_minPn_cut       = True,
                         doApply_Mx_cut          = True,
                         W_min                   = None):#{
    '''
    e_e_pi_pass_cuts = apply_cuts_to_e_e_pi(fdebug,
                                         NeventsMax,
                                         NMaxPerSubset,
                                         doAcceptanceMatchingCut,
                                         doApply_minPn_cut,
                                         doApply_Mx_cut)
                                         
    Sep-15, 2022
    '''
    # d(e,e'\pi) SIDIS data
    global e_e_pi, e_e_pi_pass_cuts

    
    if doAcceptanceMatchingCut:#{
        # e_e_pi_after_p_theta_cut = apply_p_theta_acceptance_cut( e_e_pi,
#                                                                NeventsMax=NeventsMax,
#                                                                NMaxPerSubset=NMaxPerSubset,
#                                                                fdebug=fdebug )
        
        e_e_pi_after_p_theta_cut = apply_p_theta_acceptance_cut_single_set( e_e_pi,
                                                                  NeventsMax=NeventsMax,
                                                                  fdebug=fdebug )

    #}
    else:#{
        e_e_pi_after_p_theta_cut = dict()
        for pi_ch,pi_print in zip(pi_charge_names,pi_prints):#{
            df = e_e_pi[pi_ch]
            e_e_pi_after_p_theta_cut[pi_ch] = df[0:NeventsMax]
        #}
    #}

    if doApply_Mx_cut:  #{
        e_e_pi_after_Mx_cut      = apply_Mx_cut( e_e_pi_after_p_theta_cut )
    #}
    else: #{
        e_e_pi_after_Mx_cut      = e_e_pi_after_p_theta_cut
    #}

    e_e_pi_after_Kinematical_cuts = apply_Kinematical_cuts( e_e_pi_after_Mx_cut, W_min=W_min )
    e_e_pi_pass_cuts         = e_e_pi_after_Kinematical_cuts;

    Nevents      = dict()
    frac_Nevents = dict()
    for pi_ch,pi_print in zip(pi_charge_names,pi_prints):#{
        print('(e,e',pi_print,')')
        Nevents,frac_Nevents = dict(),dict()
        if NeventsMax < 0:
            Nevents,frac_Nevents = get_Nevents(pi_ch, 'original',  e_e_pi, Nevents, frac_Nevents);
        else:
            Nevents[pi_ch + ' original cut'] = NeventsMax
            frac_Nevents[pi_ch + ' original cut'] = 1

        Nevents,frac_Nevents = get_Nevents(pi_ch, 'p-theta',      e_e_pi_after_p_theta_cut,      Nevents, frac_Nevents);
        Nevents,frac_Nevents = get_Nevents(pi_ch, 'Mx',           e_e_pi_after_Mx_cut,           Nevents, frac_Nevents);
        Nevents,frac_Nevents = get_Nevents(pi_ch, 'Kinematical',  e_e_pi_after_Kinematical_cuts, Nevents, frac_Nevents);
        
        # add beam-charge weight
        runnumbers = np.array(e_e_pi_pass_cuts[pi_ch].runnum).astype(int);
        e_e_pi_pass_cuts[pi_ch]['weight'] = runnum_weight( runnumbers )
    #}
    print(' ')
    return e_e_pi_pass_cuts
#}



# ----------------------- #
def apply_cuts_to_e_e_pi_FreeP(fdebug=2,
                             NeventsMax=-1,
                             NMaxPerSubset = 500000,
                             doAcceptanceMatchingCut = True,
                             doApply_minPn_cut       = True,
                             doApply_Mx_cut          = True):#{
    '''
    e_e_pi_FreeP_pass_cuts = apply_cuts_to_e_e_pi_FreeP(fdebug,
                                         NeventsMax,
                                         NMaxPerSubset,
                                         doAcceptanceMatchingCut,
                                         doApply_minPn_cut,
                                         doApply_Mx_cut)
                                         
    Aug-29, 2022
    
    p(e,e'\pi) SIDIS data
    '''
    
    global e_e_pi_FreeP, e_e_pi_FreeP_pass_cuts

    if doAcceptanceMatchingCut:#{
        e_e_pi_FreeP_after_p_theta_cut = apply_p_theta_acceptance_cut_single_set( e_e_pi_FreeP,
                                                                NeventsMax=NeventsMax,
#                                                                NMaxPerSubset=NMaxPerSubset,
                                                                fdebug=fdebug )
    #}
    else:#{
        e_e_pi_FreeP_after_p_theta_cut = dict()
        for pi_ch,pi_print in zip(pi_charge_names,pi_prints):#{
            df = e_e_pi_FreeP[pi_ch]
            e_e_pi_FreeP_after_p_theta_cut[pi_ch] = df[0:NeventsMax]
        #}
    #}

    if doApply_Mx_cut:  #{
        e_e_pi_FreeP_after_Mx_cut      = apply_Mx_cut( e_e_pi_FreeP_after_p_theta_cut )
    #}
    else: #{
        e_e_pi_FreeP_after_Mx_cut      = e_e_pi_FreeP_after_p_theta_cut
    #}
    e_e_pi_FreeP_after_Kinematical_cuts = apply_Kinematical_cuts( e_e_pi_FreeP_after_Mx_cut )
    e_e_pi_FreeP_pass_cuts         = e_e_pi_FreeP_after_Kinematical_cuts;
    
    
    # add beam-charge weight
    for pi_ch,pi_print in zip(pi_charge_names,pi_prints):#{
    # For free-proton (RGA) data we do not have run-numbers or beam charge info
    # (Sep-2022)
    #        runnumbers = np.array(e_e_pi_FreeP_pass_cuts[pi_ch].runnum).astype(int);
    #        e_e_pi_FreeP_pass_cuts[pi_ch]['weight'] = runnum_weight( runnumbers )
        e_e_pi_FreeP_pass_cuts[pi_ch]['weight'] = 1;
    #}



    Nevents      = dict()
    frac_Nevents = dict()
    for pi_ch,pi_print in zip(pi_charge_names,pi_prints):#{
        print('p(e,e',pi_print,')')
        Nevents,frac_Nevents = dict(),dict()
        if NeventsMax < 0:
            Nevents,frac_Nevents = get_Nevents(pi_ch, 'original',  e_e_pi_FreeP, Nevents, frac_Nevents);
        else:
            Nevents[pi_ch + ' original cut'] = NeventsMax
            frac_Nevents[pi_ch + ' original cut'] = 1
            
            
        Nevents,frac_Nevents = get_Nevents(pi_ch, 'p-theta',      e_e_pi_FreeP_after_p_theta_cut,      Nevents, frac_Nevents);
        Nevents,frac_Nevents = get_Nevents(pi_ch, 'Mx',           e_e_pi_FreeP_after_Mx_cut,           Nevents, frac_Nevents);
        Nevents,frac_Nevents = get_Nevents(pi_ch, 'Kinematical',  e_e_pi_FreeP_after_Kinematical_cuts, Nevents, frac_Nevents);
    #}
    print(' ')
    
    return e_e_pi_FreeP_pass_cuts
#}



# ----------------------- #
def apply_cuts_to_e_e_pi_GEMC(fdebug=2,
                             NeventsMax=-1,
                             NMaxPerSubset = 500000,
                             doAcceptanceMatchingCut = True,
                             doApply_minPn_cut       = True,
                             doApply_Mx_cut          = True):#{
    '''
    e_e_pi_GEMC_pass_cuts = apply_cuts_to_e_e_pi_GEMC(fdebug,
                                         NeventsMax,
                                         NMaxPerSubset,
                                         doAcceptanceMatchingCut,
                                         doApply_minPn_cut,
                                         doApply_Mx_cut)
                                         
    June-29, 2022
    
    (e,e'\pi) - (uniform) MC for acceptance correction (uniform in e and \pi)
    '''
    
    global e_e_pi_GEMC, e_e_pi_GEMC_pass_cuts

    e_e_pi_GEMC_after_eepi_cuts       = dict()

    # Apply (e,e'pi) SIDIS kinematical cuts while asking if pion was accepted,
    # externally (here, and not in the CLAS12ROOT script) since we
    # want to retain and record also the events that did not pass these cuts, in the simulation
    # whereas in data we just omit events that did not pass these cuts
    for pi_ch in pi_charge_names:#{
        e_e_pi_GEMC_after_eepi_cuts[pi_ch] = e_e_pi_GEMC[pi_ch][(e_e_pi_GEMC[pi_ch].pi_passed_cuts==1) & (e_e_pi_GEMC[pi_ch].eepiPastKinematicalCuts==1)];
    #}
    e_e_pi_GEMC_after_p_theta_cut = apply_p_theta_acceptance_cut_single_set( e_e_pi_GEMC_after_eepi_cuts )
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
    return e_e_pi_GEMC_pass_cuts
#}




# ----------------------- #
def apply_further_selection_cuts_to_data(fdebug=0,
                                         NeventsMax=-1,
                                         NMaxPerSubset = 500000,
                                         doAcceptanceMatchingCut = True,
                                         doApply_minPn_cut       = True,
                                         doApply_Mx_cut          = True,
                                         W_min                   = None):#{
    '''
    e_e_pi_pass_cuts, e_e_pi_n_pass_cuts, e_e_pi_FreeP_pass_cuts, e_e_pi_GEMC_pass_cuts = apply_further_selection_cuts_to_data(fdebug=2)
    last edit Jan-7, 2023
    
    Apply selection cuts not previously imposed
    
    The cuts applied for d(e,e'π), d(e,e'πn) and p(e,e'π) events:
    1. pi+/pi- acceptance matching cut in p-theta plane
    2. Missing mass cut
    3. ...
    
    
    input:
    --------
    doApply_*_cut      flag to apply the cut or not
    W_min              default is read-off automatically from macros/cuts/BANDcutValues.csv
        
    return:
    --------
    e_e_pi_pass_cuts, e_e_pi_n_pass_cuts, e_e_pi_FreeP_pass_cuts, e_e_pi_GEMC_pass_cuts
    
    '''
    global e_e_pi, e_e_pi_n, e_e_pi_FreeP, e_e_pi_GEMC
    global e_e_pi_pass_cuts, e_e_pi_n_pass_cuts, e_e_pi_FreeP_pass_cuts, e_e_pi_GEMC_pass_cuts
    
    
    # (e,e'\pi n) SIDIS data complete this -  need to add sector ID in the (e,e'\pi n) data
    print('Applying selection cuts not previously imposed')
    
    # print number of events retained on every cut
    # if fdebug < 1: return
    Nevents      = dict()
    frac_Nevents = dict()
    
    # (1) d(e,e'π) Data
    if (e_e_pi=={}) is False:#{
        if fdebug: print("(1) Applying cuts to d(e,e'π) data")
        e_e_pi_pass_cuts = apply_cuts_to_e_e_pi(fdebug, NeventsMax, NMaxPerSubset,
                                                doAcceptanceMatchingCut,
                                                doApply_minPn_cut,
                                                doApply_Mx_cut,
                                                W_min = W_min)
    # (2) d(e,e'πn) data
    if (e_e_pi_n=={}) is False:#{
        if fdebug: print("(2) Applying cuts to d(e,e'πn) data")
        e_e_pi_n_pass_cuts = apply_cuts_to_e_e_pi_n(fdebug, NeventsMax, NMaxPerSubset,
                                                    doAcceptanceMatchingCut,
                                                    doApply_minPn_cut,
                                                    doApply_Mx_cut)
    # (3) p(e,e'π) Data
    if (e_e_pi_FreeP=={}) is False:#{
        if fdebug: print("(3) Applying cuts to p(e,e'π) data")
        e_e_pi_FreeP_pass_cuts = apply_cuts_to_e_e_pi_FreeP(fdebug, NeventsMax,NMaxPerSubset,
                                                    doAcceptanceMatchingCut,
                                                    doApply_minPn_cut,
                                                    doApply_Mx_cut)
    # (4) MC
    if (e_e_pi_GEMC=={}) is False:#{
        if fdebug: print('(4) MC')
        e_e_pi_GEMC_pass_cuts = apply_cuts_to_e_e_pi_GEMC(fdebug, NeventsMax, NMaxPerSubset,
                                         doAcceptanceMatchingCut,
                                         doApply_minPn_cut,
                                         doApply_Mx_cut)
    print('Done applying event-selection cuts')
    return e_e_pi_pass_cuts, e_e_pi_n_pass_cuts, e_e_pi_FreeP_pass_cuts, e_e_pi_GEMC_pass_cuts
#}
# ----------------------- #


# ----------------------- #
def read_run_nunmbers(runfile_path = "/Users/erezcohen/Desktop/Software/CLAS12/BAND/SIDIS_at_BAND/macros/runlists/",
                      runs_filename  = "good_runs_10-2-final.txt",
                      Nruns=-1):#{
    '''
    runs = read_run_nunmbers(runs_filename  = "good_runs_10-2-final.txt")
    
    last edit Aug-26, 2022
    Read run numbers using readlines()
    '''
    runs_file     = open( runfile_path + runs_filename, 'r')
    run_fileLines = runs_file.readlines()
    if Nruns==-1: Nruns = len(run_fileLines)
    runs = []
    for line in run_fileLines[0:Nruns]:#{
        run = int(line.strip())
        runs.append(run)
    #}
    runs = np.array(runs)
    return runs;
#}










# ----------------------- #
def apply_Kinematical_cuts( df_dict,
                           Q2_min=None,     Q2_max =None,
                           Pe_min=None,     Pe_max =None,
                           Ppi_min=None,    Ppi_max=None ,
                           W_min=None,      W_max  =None,
                           cuts_filepath = "/Users/erezcohen/Desktop/Software/CLAS12/BAND/SIDIS_at_BAND/macros/cuts/",
                           cuts_filename = "BANDcutValues.csv"):#{
    '''
        df_dict_after_cut = apply_Kinematical_cuts(df_dict)
        
        Apply kinematical cuts to match BAND neutron skimming
        last update Sep-17, 2022
        
        cuts are read-off from cuts_filename,
        unless specifically stated in function input paraemters
        
    '''
    
    cuts = pd.read_csv( cuts_filepath + cuts_filename );
    if W_min is None: W_min = float(cuts[cuts.parameter=="W_min"].value)
    if W_max is None: W_max = float(cuts[cuts.parameter=="W_max"].value)
    if Q2_min is None: Q2_min = float(cuts[cuts.parameter=="Q2_min"].value)
    if Q2_max is None: Q2_max = float(cuts[cuts.parameter=="Q2_max"].value)

    if Pe_min is None: Pe_min = float(cuts[cuts.parameter=="Pe_min"].value)
    if Pe_max is None: Pe_max = float(cuts[cuts.parameter=="Pe_max"].value)

    if Ppi_min is None: Ppi_min = float(cuts[cuts.parameter=="Ppi_min"].value)
    if Ppi_max is None: Ppi_max = float(cuts[cuts.parameter=="Ppi_max"].value)


    
    
    df_dict_after_cut = dict()
    for pi_charge_name in pi_charge_names:#{
        df = df_dict[pi_charge_name]
        df = df[ (W_min < df.W)                 & (df.W < W_max)     ]
        df = df[ (Q2_min < df.Q2)               & (df.Q2 < Q2_max)   ]
        df = df[ (Pe_min < df.e_P)              & (df.e_P < Pe_max)  ]
        df = df[ (Ppi_min < df.pi_P)            & (df.pi_P < Ppi_max)]
        df_dict_after_cut[pi_charge_name] = df
    #}
    return df_dict_after_cut
#}
# ----------------------- #











# ----------------------- #
def get_Nevents(pi_ch='piplus',
                cut_name='p-theta',
                df_dict = None,
                Nevents = dict(),
                frac_Nevents = dict()):
    '''
    last update July-26, 2022
    
    input:
    --------
    cut name
    
    return:
    --------
    Nevents,frac_Nevents
    
    '''
    df_after_cut = df_dict[pi_ch]
    if pi_ch+' original cut' not in Nevents.keys(): #{
        Nevents[pi_ch + ' original cut']      = len(df_after_cut)
        frac_Nevents[pi_ch + ' original cut'] = 1
    #}

    label = pi_ch +' ' + cut_name + ' cut'
    Nevents[label] = len(df_after_cut)
    frac_Nevents[label] = float(Nevents[label])/ Nevents[pi_ch + ' original cut']
    print(Nevents[label],'events after '+cut_name+' cut (%.1f'%(100.*frac_Nevents[label]),'%)')
    return Nevents,frac_Nevents
# ----------------------- #





# ----------------------- #
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
# ----------------------- #























#
#
#
## ----------------------- #
#def apply_p_theta_acceptance_cut( df_dict=None,
#                                 NeventsMax=-1,
#                                 NMaxPerSubset = 500000,
#                                 fdebug=1):
#    '''
#        df_dict_after_cut = apply_p_theta_acceptance_cut(df_dict)
#
#        Apply a π+/π- acceptance matching cut on the in p-\theta plane
#        last update July-12, 2022
#
#    '''
#    import numpy as np
#    print("Apply a π+/π- acceptance matching cut on the in p-theta plane")
#    df_dict_after_cut = dict()
#
#
#    for pi_ch in pi_charge_names:
#
#        if NeventsMax > 0: NeventsMax = np.min( [NeventsMax, len(df_dict[pi_ch])] )
#        else:              NeventsMax = len(df_dict[pi_ch])
#        df = df_dict[pi_ch][0:NeventsMax]
#        if fdebug: print('Applying p-theta on cut for '+pi_ch+' on %d events'%NeventsMax)
#
#        # subdivide samples into sub-samples of no more than 1M events each
#        Nsubsets  = np.max([1,(int)(NeventsMax/NMaxPerSubset)]);
#        if fdebug: print('Subdividing into %d subsets up to %d events'%(Nsubsets,NMaxPerSubset))
#
#        for subset_idx in range(Nsubsets):#{
#            i_min_subset = NMaxPerSubset*subset_idx
#            i_max_subset = np.min([NMaxPerSubset*(subset_idx+1), NeventsMax])
#            subset_df = df[i_min_subset:i_max_subset]
#            good_indices = np.array([])
#            if fdebug>1:
#                print('subset %d of index %d-%d'%(subset_idx,i_min_subset,i_max_subset-1))
#
#            # focus on each sector and select the events that pass the p-theta cut
#            for sector in range(1,7):#{
#                df_in_sector   = subset_df[subset_df.pi_DC_sector == sector]
#                pi_P           = np.array(df_in_sector.pi_P)
#                theta_min_pi   = pi_min_theta_cut( pi_charge = 'any', sector=sector, p=pi_P )
#                good_indices_in_sector = []
#                good_indices_in_sector = df_in_sector[ df_in_sector.pi_Theta*r2d > theta_min_pi ].index
#                good_indices = np.concatenate([good_indices,good_indices_in_sector])
#            #}
#            good_indices = np.unique(good_indices)
#            subset_df_in_cut = subset_df.loc[good_indices]
#            if subset_idx==0: df_after_cut = subset_df_in_cut
#            else:             df_after_cut = pd.concat([df_after_cut , subset_df_in_cut])
#        #}
#        df_dict_after_cut[pi_ch] = df_after_cut
#    return df_dict_after_cut
## ----------------------- #







# ----------------------- #
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
# ----------------------- #





# ----------------------- #
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
# ----------------------- #


# ----------------------- #
def apply_Mx_cut( df_dict=None, Mx_min = 1.3, Mx_max = 5 ): #{
    '''
        df_dict_after_cut = apply_Mx_cut(df_dict)
        
        Apply a cut on the missing mass of a (e,e'\pi) reaction
        
    '''
    print("Apply a cut on the missing mass of a (e,e'π) reaction: %.1f<Mx<%.1f GeV/c2"%(Mx_min,Mx_max))
    df_dict_after_cut = dict()
    for pi_charge_name in pi_charge_names:#{
        df = df_dict[pi_charge_name]
        df = df[ (Mx_min < df.M_x) & (df.M_x < Mx_max)]
        df_dict_after_cut[pi_charge_name] = df
    #}
    return df_dict_after_cut
#}
# ----------------------- #





# ------------------------------------------
#
#      Bad functions that were corrected
#
# ------------------------------------------



# ----------------------- #
def apply_further_selection_cuts_to_data_before_correction(fdebug=2,
                                         NeventsMax=-1,
                                         NMaxPerSubset = 500000,
                                         doAcceptanceMatchingCut = True,
                                         doApply_minPn_cut       = True,
                                         doApply_Mx_cut          = True,
                                         W_min                   = None):#{
    '''
    e_e_pi_pass_cuts, e_e_pi_n_pass_cuts, e_e_pi_GEMC_pass_cuts = apply_further_selection_cuts_to_data(fdebug=2)
    last edit Aug-26, 2022
    
    Apply selection cuts not previously imposed
    
    The cuts applied for d(e,e'π), d(e,e'πn) and p(e,e'π) events:
    1. pi+/pi- acceptance matching cut in p-theta plane
    2. Missing mass cut
    3. ...
    
    
    input:
    --------
    doApply_*_cut      flag to apply the cut or not
    W_min              default is read-off automatically from macros/cuts/BANDcutValues.csv
        
    return:
    --------
    e_e_pi_pass_cuts, e_e_pi_n_pass_cuts, e_e_pi_FreeP_pass_cuts, e_e_pi_GEMC_pass_cuts
    
    '''
    global e_e_pi, e_e_pi_n, e_e_pi_FreeP, e_e_pi_GEMC
    global e_e_pi_pass_cuts, e_e_pi_n_pass_cuts, e_e_pi_FreeP_pass_cuts, e_e_pi_GEMC_pass_cuts
    
    
    # (e,e'\pi n) SIDIS data complete this -  need to add sector ID in the (e,e'\pi n) data
    print('Applying selection cuts not previously imposed')
    
    # print number of events retained on every cut
    # if fdebug < 1: return
    Nevents      = dict()
    frac_Nevents = dict()
    
    # (1) d(e,e'π) Data
    if (e_e_pi=={}) is False:#{
        print("(1) Applying cuts to d(e,e'π) data")
        e_e_pi_pass_cuts = apply_cuts_to_e_e_pi_before_correction(fdebug, NeventsMax, NMaxPerSubset,
                                                doAcceptanceMatchingCut,
                                                doApply_minPn_cut,
                                                doApply_Mx_cut,
                                                W_min = W_min)
        print(e_e_pi_pass_cuts.keys())
        
    # (2) d(e,e'πn) data
    if (e_e_pi_n=={}) is False:#{
        print("(2) Applying cuts to d(e,e'πn) data")
        e_e_pi_n_pass_cuts = apply_cuts_to_e_e_pi_n(fdebug, NeventsMax, NMaxPerSubset,
                                                    doAcceptanceMatchingCut,
                                                    doApply_minPn_cut,
                                                    doApply_Mx_cut)
    # (3) p(e,e'π) Data
    if (e_e_pi_FreeP=={}) is False:#{
        print("(3) Applying cuts to p(e,e'π) data")
        e_e_pi_FreeP_pass_cuts = apply_cuts_to_e_e_pi_FreeP(fdebug, NeventsMax,NMaxPerSubset,
                                                    doAcceptanceMatchingCut,
                                                    doApply_minPn_cut,
                                                    doApply_Mx_cut)
    # (4) MC
    if (e_e_pi_GEMC=={}) is False:#{
        print('(4) MC')
        e_e_pi_GEMC_pass_cuts = apply_cuts_to_e_e_pi_GEMC(fdebug, NeventsMax, NMaxPerSubset,
                                         doAcceptanceMatchingCut,
                                         doApply_minPn_cut,
                                         doApply_Mx_cut)
    print('Done applying event-selection cuts')
    print('e_e_pi_pass_cuts.keys():',e_e_pi_pass_cuts.keys())
    # print(len(e_e_pi_pass_cuts['piplus']))
    return e_e_pi_pass_cuts, e_e_pi_n_pass_cuts, e_e_pi_FreeP_pass_cuts, e_e_pi_GEMC_pass_cuts
#}
# ----------------------- #



# ----------------------- #
def apply_p_theta_acceptance_cut_before_correction( df_dict=None,
                                NeventsMax=-1,
                                NMaxPerSubset = 500000,
                                fdebug=1):
    '''
       df_dict_after_cut = apply_p_theta_acceptance_cut(df_dict)

       Apply a π+/π- acceptance matching cut on the in p-\theta plane
       last update Sep-22, 2022

    '''
    import numpy as np
    print("Apply a π+/π- acceptance matching cut on the in p-theta plane")
    df_dict_after_cut = dict()
    for pi_ch in pi_charge_names:

        if NeventsMax > 0: NeventsMax = np.min( [NeventsMax, len(df_dict[pi_ch])] )
        else:              NeventsMax = len(df_dict[pi_ch])
        df = df_dict[pi_ch][0:NeventsMax]
        if fdebug: print('Applying p-theta on cut for '+pi_ch+' on %d events'%NeventsMax)
        good_indices = np.array([])
        for sector in range(1,7):#{
            df_in_sector   = df[df.pi_DC_sector == sector]
            if fdebug: print(len(df_in_sector),'in sector',sector)
            theta_min_pi   = pi_min_theta_cut( pi_charge = 'any', sector=sector, p=np.array(df_in_sector.pi_P) )
            good_indices_in_sector = []
            good_indices_in_sector = df_in_sector[ df_in_sector.pi_Theta*r2d > theta_min_pi ].index
            good_indices = np.concatenate([good_indices,good_indices_in_sector])
        #}
        good_indices = np.unique(good_indices)
        df_after_cut = df.loc[good_indices]

        df_dict_after_cut[pi_ch] = df_after_cut
    return df_dict_after_cut
# ----------------------- #



# ----------------------- #
def apply_cuts_to_e_e_pi_before_correction(fdebug=2,
                         NeventsMax=-1,
                         NMaxPerSubset = 500000,
                         doAcceptanceMatchingCut = True,
                         doApply_minPn_cut       = True,
                         doApply_Mx_cut          = True,
                         W_min                   = None):#{
    '''
    e_e_pi_pass_cuts = apply_cuts_to_e_e_pi(fdebug,
                                         NeventsMax,
                                         NMaxPerSubset,
                                         doAcceptanceMatchingCut,
                                         doApply_minPn_cut,
                                         doApply_Mx_cut)
                                         
    Sep-22, 2022
    '''
    # d(e,e'\pi) SIDIS data
    global e_e_pi, e_e_pi_pass_cuts

    
    if doAcceptanceMatchingCut:#{
        e_e_pi_after_p_theta_cut = apply_p_theta_acceptance_cut_before_correction( e_e_pi,
                                                                  NeventsMax=NeventsMax,
                                                                  fdebug=fdebug )
        if fdebug: print('e_e_pi_after_p_theta_cut.keys():',e_e_pi_after_p_theta_cut.keys())

    #}
    else:#{
        e_e_pi_after_p_theta_cut = dict()
        for pi_ch,pi_print in zip(pi_charge_names,pi_prints):#{
            df = e_e_pi[pi_ch]
            e_e_pi_after_p_theta_cut[pi_ch] = df[0:NeventsMax]
        #}
    #}

    if doApply_Mx_cut:  #{
        e_e_pi_after_Mx_cut      = apply_Mx_cut( e_e_pi_after_p_theta_cut )
    #}
    else: #{
        e_e_pi_after_Mx_cut      = e_e_pi_after_p_theta_cut
    #}

    e_e_pi_after_Kinematical_cuts = apply_Kinematical_cuts( e_e_pi_after_Mx_cut, W_min=W_min )
    e_e_pi_pass_cuts         = e_e_pi_after_Kinematical_cuts;
    if fdebug: print('e_e_pi_pass_cuts.keys():',e_e_pi_pass_cuts.keys())

    Nevents      = dict()
    frac_Nevents = dict()
    for pi_ch,pi_print in zip(pi_charge_names,pi_prints):#{
        if fdebug: print('(e,e',pi_print,')')
        Nevents,frac_Nevents = dict(),dict()
        if NeventsMax < 0:
            Nevents,frac_Nevents = get_Nevents(pi_ch, 'original',  e_e_pi, Nevents, frac_Nevents);
        else:
            Nevents[pi_ch + ' original cut'] = NeventsMax
            frac_Nevents[pi_ch + ' original cut'] = 1

        Nevents,frac_Nevents = get_Nevents(pi_ch, 'p-theta',      e_e_pi_after_p_theta_cut,      Nevents, frac_Nevents);
        Nevents,frac_Nevents = get_Nevents(pi_ch, 'Mx',           e_e_pi_after_Mx_cut,           Nevents, frac_Nevents);
        Nevents,frac_Nevents = get_Nevents(pi_ch, 'Kinematical',  e_e_pi_after_Kinematical_cuts, Nevents, frac_Nevents);
        
        # add beam-charge weight
        runnumbers = np.array(e_e_pi_pass_cuts[pi_ch].runnum).astype(int);
        e_e_pi_pass_cuts[pi_ch]['weight'] = runnum_weight( runnumbers )
    #}
    if fdebug: print(' ')
    return e_e_pi_pass_cuts
#}



