# -*- coding: utf-8 -*-
"""
Produce a LUND file of (e,e'\pi) events to simulate with GEMC.
                        We are interested in pion accetpance,
                        so we simulate the same electron multiple times
                        
                        
 LUND output:
2   1   1    0.0   0.0 PDG_beam   Ebeam   1       1      event_weight
1  -1.  1     11   0    0  -0.9830   0.0981  9.6502  9.7007  0.0005  0.0000 0.0000  -0.8072
2   1.  1   2212   0    0   0.7333   0.1126  0.6391  1.3560  0.9380  0.0000 0.0000  -0.8072
                       

"""
#%% Imports and definitions
# =============================================================================
# Imports and definitions
# =============================================================================
import numpy as np, matplotlib.pyplot as plt, pandas as pd
import sys; sys.path.insert(0, '/Users/erezcohen/Desktop/Software/mySoftware/Python/');
from my_tools               import *;
from plot_tools             import *;
from my_data_analysis_tools import *;
#%config InlineBackend.figure_format = 'retina'
    
Nevents = 10000;

particles = ['e-','pi'];
r2d       = 180./np.pi;
d2r       = np.pi/180.;
p_min     = 1.25 # GeV/c
p_max     = 5.00 # GeV/c
theta_min = 0    # deg.
theta_max = 180  # deg.
phi_min   = -180 # deg.
phi_max   = 180  # deg.
vz_min    = -10; # cm
vz_max    = 10;  # cm
fdebug    = 3
m_pi      = 0.13957 # GeV/c2
m_e       = 0.511e-3 # GeV/c2

PDG_e       = 11;
PDG_beam    = PDG_e;
PDG_pips    = 211;
PDG_pims    = -211;
Ebeam       = 10.2;
event_weight= 1;
Nparticles  = 2;



#%% Initialize
# =============================================================================
# Initialize
# =============================================================================

M           = np.zeros(Nparticles)
PDG         = np.zeros(Nparticles)
Px,Py,Pz,E  = np.zeros(Nparticles),np.zeros(Nparticles),np.zeros(Nparticles),np.zeros(Nparticles)
vx,vy,vz    = np.zeros(Nparticles),np.zeros(Nparticles),np.zeros(Nparticles)
p_arr       = np.zeros(Nevents)
theta_arr   = np.zeros(Nevents)
phi_arr     = np.zeros(Nevents)

M[0]  = m_e
PDG[0]= PDG_e
# electron momentum from example DVCS event [https://gemc.jlab.org/gemc/html/documentation/generator/lund.html]
Px[0] = -0.9830
Py[0] = 0.0981
Pz[0] = 9.6502
Pe    = np.sqrt(np.square(Px[0])+np.square(Py[0])+np.square(Pz[0]))
E[0]  = np.sqrt( np.square(Pe) + np.square(m_e) )
# electron vertex - at origin
vx[0] = 0
vy[0] = 0
vz[0] = 0;
 



M           = np.zeros(Nparticles)
PDG         = np.zeros(Nparticles)
Px,Py,Pz,E  = np.zeros(Nparticles),np.zeros(Nparticles),np.zeros(Nparticles),np.zeros(Nparticles)
vx,vy,vz    = np.zeros(Nparticles),np.zeros(Nparticles),np.zeros(Nparticles)

M[0]  = m_e
PDG[0]= PDG_e
# electron momentum from example DVCS event [https://gemc.jlab.org/gemc/html/documentation/generator/lund.html]
Px[0] = -0.9830
Py[0] = 0.0981
Pz[0] = 9.6502
Pe    = np.sqrt(np.square(Px[0])+np.square(Py[0])+np.square(Pz[0]))
E[0]  = np.sqrt( np.square(Pe) + np.square(m_e) )
# electron vertex - at origin
vx[0] = 0
vy[0] = 0
vz[0] = 0;
 


#%% Sample pion momentum and print to file
# do the same process of positive and negative pions
fdebug=1
for pi_charge,pi_label,pi_PDG in zip(['pips','pims'],['\pi^+','\pi^-'],[PDG_pips,PDG_pims]):
    
    p_arr       = np.zeros(Nevents)
    theta_arr   = np.zeros(Nevents)
    phi_arr     = np.zeros(Nevents)

    PDG[1]= pi_PDG
    M[1]  = m_pi

    print("Generating %d (e,e'%s) events."%(Nevents,pi_label))    
    outputfile = open("/Users/erezcohen/Desktop/data/BAND/AcceptanceCorrection/InputFiles/"
                      +"/ee%s_p_uniform_distribution.dat"%(pi_charge), "w")
    for n in range(Nevents):
        # sample pion momentum uniformly in p, cos(theta), and phi
        p           = np.random.uniform(p_min    ,p_max    )
        cos_theta   = np.random.uniform( np.cos(theta_min*d2r) ,np.cos(theta_max*d2r) )
        theta       = np.arccos(cos_theta)
        phi         = np.random.uniform(phi_min*d2r  ,phi_max*d2r  )
        
        # record for monitoring
        p_arr[n]    = p
        theta_arr[n]= theta
        phi_arr[n]  = phi
        
        # move to Cartesian coordinates
        Px[1] = p*np.sin(theta)*np.cos(phi)
        Py[1] = p*np.sin(theta)*np.sin(phi)
        Pz[1] = p*np.cos(theta)
        E[1]  = np.sqrt( np.square(p) + np.square(m_pi) )
        
        # sample pion vertex 
        vx[1] = 0
        vy[1] = 0
        vz[1] = np.random.uniform( vz_min, vz_max );
        
        if fdebug:  
            if n%(Nevents/10)==0: print( '%.1f'%(100.*n/Nevents)+'%')
            if fdebug>1: print('event',n,', p =',p,'GeV/c, \\theta =',theta, ', \phi =',phi );
        
        
        event_header_str = ('%d \t 1 \t 1 \t 0.0 \t 0.0 %d \t %.3f \t 1 \t 1 \t %d\n'%
                      (Nparticles, PDG_beam, Ebeam, event_weight))
        outputfile.write( event_header_str )
        if fdebug>2: print( event_header_str )
        
        for particle, j in zip(particles,range(len(particles))):
            particle_str = ('%d \t %.3f \t %d \t %d \t %d \t %d \t %.3f \t %.3f \t %.3f \t %.3f \t %.3f \t %.3f \t %.3f \t %.3f \n'%
                            (j+1, 0.0, 1, PDG[j], 0, 0, Px[j], Py[j], Pz[j], E[j], M[j], vx[j], vy[j], vz[j])) 
            outputfile.write( particle_str )
            if fdebug>2: print(particle_str)
    
        
    print("done generating %d (e,e'%s) events."%(Nevents,pi_label))
    outputfile.close();
    print("Saved output file \n%s"%("/Users/erezcohen/Desktop/data/BAND/AcceptanceCorrection/InputFiles/"
                                    +"/ee%s_p_uniform_distribution.dat"%(pi_charge)))    
print('done.')

#%% Monitoring plots
# =============================================================================
# Monitoring plots
# =============================================================================

fig = plt.figure(figsize=(16,6))
for x,label,units,bins,scale_factor,subplot_idx in zip([p_arr,theta_arr,phi_arr],
                                    ['$p$','$\\theta$', '$\phi$'], 
                                    ['[GeV/c]','[deg.]', '[deg.]'], 
                                    [np.linspace(p_min, p_max,30),np.linspace(theta_min, theta_max,30),np.linspace(phi_min, phi_max,30)],
                                    [1,r2d,r2d],
                                    range(3)):
    ax = fig.add_subplot(1,3,subplot_idx+1)
    ax.hist( x*scale_factor, bins , edgecolor='k')
    set_axes(ax,label + ' ' + units,'counts' if subplot_idx==0 else '',title=label, fontsize=18)

plt.tight_layout();



#%% arxiv
# =============================================================================
# archive
# =============================================================================












