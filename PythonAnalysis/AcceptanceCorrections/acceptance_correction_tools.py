import numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib as mpl
import sys; sys.path.insert(0, '/Users/erezcohen/Desktop/Software/mySoftware/Python/'); 
from my_tools               import *; 
from plot_tools             import *;
from my_data_analysis_tools import *;
r2d = 180./np.pi

pi_charge_names  = ['piplus'   ,'piminus'  ]
pi_labels        = ['\pi^{+}'  ,'\pi^{-}'  ]
pi_colors        = ['royalblue','black'   ]

def apply_p_theta_acceptance_cut( df_dict=None ):    
    # df_dict_after_cut = apply_p_theta_acceptance_cut(df_dict)
    
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
