# -*- coding: utf-8 -*-
"""
Produce a LUND file of (e,e'\pi) events to simulate with GEMC.
                        We are interested in pion accetpance,
                        so we simulate the same electron multiple times
                        
   last update Jan-26, 2022
                        
 LUND output:
2   1   1    0.0   0.0 PDG_beam   Ebeam   1       1      event_weight
1  -1.  1     11   0    0  -0.9830   0.0981  9.6502  9.7007  0.0005  0.0000 0.0000  -0.8072
2   1.  1   2212   0    0   0.7333   0.1126  0.6391  1.3560  0.9380  0.0000 0.0000  -0.8072
      
First attempt:                 
The electron is taken from an accepted and selected (e,e'pi-) event from run 6420
/Users/erezcohen/Desktop/data/BAND/SIDIS_skimming/skimmed_SIDIS_inc_006420_e_piminus_selected_eepi_kinematics.csv

runnum,evnum,e_P,e_Theta,e_Phi,e_Vz
6420,993,6.33416,0.175525,1.18165,-3.39518

Second attempt: 
Electron is taken from a simulated event that worked....
MC particle PDG code 11, p: 1.345 GeV/c, , theta: 30.9 deg, , phi: -163.4 deg, , V(z): 7.758 cm


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

#%% Electron information 

M[0]  = m_e
PDG[0]= PDG_e
# electron momentum from example DVCS event [https://gemc.jlab.org/gemc/html/documentation/generator/lund.html]
# Px[0] = -0.9830 #Py[0] = 0.0981 #Pz[0] = 9.6502
# Pe    = np.sqrt(np.square(Px[0])+np.square(Py[0])+np.square(Pz[0]))
# E[0]  = np.sqrt( np.square(Pe) + np.square(m_e) )

# electron taken from an accepted and selected (e,e'pi-) event from run 6420 event 993
# e_P    = 6.33416
# e_theta = 0.175525 # = 10.1 deg.
# e_phi = 1.18165    # = 67.4 deg.
# e_Vz  = -3.39518

# electron taken from a simulated event that passed reconstruction stage with default clas12 yaml card
# MC particle PDG code 11, p: 1.345 GeV/c, , theta: 30.9 deg, , phi: -163.4 deg, , V(z): 7.758 cm
e_P    = 1.345
e_theta= 30.9 * r2d # = 30.9 deg.
e_phi = -163.4 * r2d   # = 67.4 deg.
e_Vz  = 7.758 

Px[0] = e_P*np.sin(e_theta)*np.cos(e_phi)
Py[0] = e_P*np.sin(e_theta)*np.sin(e_phi)
Pz[0] = e_P*np.cos(e_theta)
E[0]  = np.sqrt( np.square(e_P) + np.square(m_e) )

# electron vertex - at origin
vx[0] = 0
vy[0] = 0
vz[0] = e_Vz;
 



#%% Sample pion momentum and print to file

# do the same process of positive and negative pions
fdebug=1
for pi_charge,pi_label,pi_PDG in zip(['pips','pims'],['\pi^+','\pi^-'],[PDG_pips,PDG_pims]):
    
    p_arr       = np.zeros(Nevents)
    theta_arr   = np.zeros(Nevents)
    phi_arr     = np.zeros(Nevents)

    px_arr       = np.zeros(Nevents)
    py_arr       = np.zeros(Nevents)
    pz_arr       = np.zeros(Nevents)
    
    PDG[1]= pi_PDG
    M[1]  = m_pi

    print("Generating %d (e,e'%s) events."%(Nevents,pi_label))    
    outputfile = open("/Users/erezcohen/Desktop/data/BAND/AcceptanceCorrection/InputFiles/"
                      +"/ee%s_p_uniform_distribution.dat"%(pi_charge), "w")
    for n in range(Nevents):
        
        # # Sample electron momentum to find a good one that is accepted
        # p           = np.random.uniform(p_min    ,p_max    )
        # cos_theta   = np.random.uniform( np.cos(theta_min*d2r) ,np.cos(theta_max*d2r) )
        # theta       = np.arccos(cos_theta)
        # phi         = np.random.uniform(phi_min*d2r  ,phi_max*d2r  )        
        # Px[0] = p*np.sin(theta)*np.cos(phi)
        # Py[0] = p*np.sin(theta)*np.sin(phi)
        # Pz[0] = p*np.cos(theta)
        # E[0]  = np.sqrt( np.square(p) + np.square(m_e) )
        # vz[0] = np.random.uniform( vz_min, vz_max );

        
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
        
        px_arr[n]    = Px[1]
        py_arr[n]    = Py[1]
        pz_arr[n]    = Pz[1]
        
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
            particle_str = ('%d \t %.3f \t %d \t %d \t %d \t %d \t %.5f \t %.5f \t %.5f \t %.5f \t %.5f \t %.5f \t %.5f \t %.3f \n'%
                            (j+1, 0.0, 1, PDG[j], 0, 0, Px[j], Py[j], Pz[j], E[j], M[j], vx[j], vy[j], vz[j])) 
            outputfile.write( particle_str )
            if fdebug>2: print(particle_str)
    
        
    print("done generating %d (e,e'%s) events."%(Nevents,pi_label))
    outputfile.close();
    print("Saved output file \n%s"%("/Users/erezcohen/Desktop/data/BAND/AcceptanceCorrection/InputFiles/"
                                    +"/ee%s_p_uniform_distribution.dat"%(pi_charge)))    
    print('scp -r /Users/erezcohen/Desktop/data/BAND/AcceptanceCorrection/InputFiles/ee%s_p_uniform_distribution.dat cohen@ftp.jlab.org:/volatile/clas12/users/ecohen/GEMC/LUND/10.2/AcceptanceCorrection/%s/'%(pi_charge,pi_charge))
print('done.')
print(' ')

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

fig = plt.figure(figsize=(16,6))
for x,label,units,bins,scale_factor,subplot_idx in zip([px_arr,py_arr,pz_arr],
                                    ['$p_x$','$p_y$', '$p_z$'], 
                                    ['[GeV/c]','[GeV/c]', '[GeV/c]'], 
                                    [np.linspace(-p_max, p_max,30),np.linspace(-p_max, p_max,30),np.linspace(-p_max, p_max,30)],
                                    [1,1,1],
                                    range(3)):
    ax = fig.add_subplot(1,3,subplot_idx+1)
    ax.hist( x*scale_factor, bins , edgecolor='k')
    set_axes(ax,label + ' ' + units,'counts' if subplot_idx==0 else '',title=label, fontsize=18)

plt.tight_layout();



#%% arxiv
# =============================================================================
# archive
# =============================================================================













