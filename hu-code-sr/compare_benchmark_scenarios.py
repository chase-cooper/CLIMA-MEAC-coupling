"""
Purpose of this code is to reproduce Figure 4 and Table 7 of Hu+2012
"""
########################
###Import useful libraries
########################
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import pdb

########################
###Define useful constants, all in CGS (via http://www.astro.wisc.edu/~dolan/constants.html)
########################

#Unit conversions
km2m=1.e3 #1 km in m
km2cm=1.e5 #1 km in cm
cm2km=1.e-5 #1 cm in km
amu2g=1.66054e-24 #1 amu in g
bar2atm=0.9869 #1 bar in atm
Pa2bar=1.e-5 #1 Pascal in bar
bar2Pa=1.e5 #1 bar in Pascal
Pa2Ba=10.#1 Pa in Barye
deg2rad=np.pi/180.
bar2barye=1.e+6 #1 Bar in Barye (the cgs unit of pressure)
barye2bar=1.e-6 #1 Barye in Bar
micron2m=1.e-6 #1 micron in m
micron2cm=1.e-4 #1 micron in cm
metricton2kg=1000. #1 metric ton in kg

#Fundamental constants
c=2.997924e10 #speed of light, cm/s
h=6.6260755e-27 #planck constant, erg/s
k=1.380658e-16 #boltzmann constant, erg/K
sigma=5.67051e-5 #Stefan-Boltzmann constant, erg/(cm^2 K^4 s)
R_earth=6371.*km2m#radius of earth in m
R_sun=69.63e9 #radius of sun in cm
AU=1.496e13#1AU in cm


########################
###Establish key
########################

#Corrected for Hu 1-indexing vs Python 0-indexing
ind_o=1-1 #O
ind_h=3-1 #H
ind_oh=4-1 #OH
ind_ph3 = 112-1 #PH3
ind_ph3d = 113-1 #PH3D

ind_so2=43-1
ind_so=42-1
ind_ch4=21-1
ind_h2s=45-1
ind_h2=53-1
ind_h2o=7-1

ind_no=12-1
ind_n2o=11-1

ind_co2=52-1
ind_n2=55-1


ind_s=40-1
ind_s2=41-1
ind_s8=79-1
ind_s8a=111-1
ind_ch4o=24-1
ind_c2h2=27-1

ind_ocs=49-1

ind_o1d=56-1
ind_co=20-1
ind_o2=54-1
ind_c2h6=31-1
ind_h2o2=6-1
ind_h2so4=73-1
ind_h2so4a=78-1
ind_ch2o=22-1
ind_o3=2-1
ind_ho2=5-1

def plot_comparison(h2_file, n2_file, co2_file, name):
    """
    """
    ########################
    ###Read in H2-dominated data
    ########################
    h2_data=np.genfromtxt(h2_file, skip_header=2, skip_footer=0, unpack=False) #Import mapping between numerical ID in code and species name.
    
    z_center_h2=h2_data[:,0] # Center of altitude bins, km 
    T_z_h2=h2_data[:,3] # Temperature(z), in K
    P_z_h2=h2_data[:,4]# Pressure(z), in Pa
    n_z_s_h2=h2_data[:,5:] #Number concentrations of the 111 chemical species, in cm**-3, as a function of (altitude, species)

    ###Get molar concentrations
    ###NOTE: May (probably) need to exclude condensed-phase species for molar concentration calculation...probably doesn't matter most of the time, but formally required and mioght matter in some weird edge cases.
    n_z_h2=P_z_h2*Pa2Ba/(k*T_z_h2)

    mc_z_s_h2=np.zeros(np.shape(n_z_s_h2))
    num_s=np.shape(n_z_s_h2)[1]

    for ind2 in range(0, num_s):
        mc_z_s_h2[:,ind2]=n_z_s_h2[:,ind2]/n_z_h2#molar concentration of each species.

    ########################
    ###Read in N2-dominated data
    ########################
    n2_data=np.genfromtxt(n2_file, skip_header=2, skip_footer=0, unpack=False) #Import mapping between numerical ID in code and species name.
    
    z_center_n2=n2_data[:,0] # Center of altitude bins, km 
    T_z_n2=n2_data[:,3] # Temperature(z), in K
    P_z_n2=n2_data[:,4] # Pressure(z), in Pa 
    n_z_s_n2=n2_data[:,5:] #Number concentrations of the 111 chemical species, in cm**-3, as a function of (altitude, species)

    ###Get molar concentrations
    ###NOTE: May (probably) need to exclude condensed-phase species for molar concentration calculation...probably doesn't matter most of the time, but formally required and mioght matter in some weird edge cases.
    n_z_n2=P_z_n2*Pa2Ba/(k*T_z_n2)

    mc_z_s_n2=np.zeros(np.shape(n_z_s_n2))
    num_s=np.shape(n_z_s_n2)[1]

    for ind2 in range(0, num_s):
        mc_z_s_n2[:,ind2]=n_z_s_n2[:,ind2]/n_z_n2#molar concentration of each species.
        
        
    ########################
    ###Read in CO2-dominated data
    ########################
    co2_data=np.genfromtxt(co2_file, skip_header=2, skip_footer=0, unpack=False) #Import mapping between numerical ID in code and species name.
    
    z_center_co2=co2_data[:,0] # Center of altitude bins, km 
    T_z_co2=co2_data[:,3] # Temperature(z), in K
    P_z_co2=co2_data[:,4] # Pressure(z), in Pa 
    n_z_s_co2=co2_data[:,5:] #Number concentrations of the 111 chemical species, in cm**-3, as a function of (altitude, species)

    ###Get molar concentrations
    ###NOTE: May (probably) need to exclude condensed-phase species for molar concentration calculation...probably doesn't matter most of the time, but formally required and mioght matter in some weird edge cases.
    n_z_co2=P_z_co2*Pa2Ba/(k*T_z_co2)

    mc_z_s_co2=np.zeros(np.shape(n_z_s_co2))
    num_s=np.shape(n_z_s_co2)[1]

    for ind2 in range(0, num_s):
        mc_z_s_co2[:,ind2]=n_z_s_co2[:,ind2]/n_z_co2#molar concentration of each species.
#    ########################
#    ###Print key parameters
#    ########################

    print('Column-Averaged Mixing ratios): (Hu+2012 Table 6)')
    print('H2')
    print('Emitted Gases')
    print('CO2. Hu+2012: 8.9E-5, this code: {0:1.1e}'.format((np.sum(mc_z_s_h2[:,ind_co2]*n_z_h2)/np.sum(n_z_h2))))
    print('SO2. Hu+2012: 9.9E-12, this code: {0:1.1e}'.format((np.sum(mc_z_s_h2[:,ind_so2]*n_z_h2)/np.sum(n_z_h2))))
    print('CH4. Hu+2012: 5.9E-6, this code: {0:1.1e}'.format((np.sum(mc_z_s_h2[:,ind_ch4]*n_z_h2)/np.sum(n_z_h2))))
    print('H2S. Hu+2012: 9.1E-10, this code: {0:1.1e}'.format((np.sum(mc_z_s_h2[:,ind_h2s]*n_z_h2)/np.sum(n_z_h2))))
    print('Photochemical Products')
    print('CO. Hu+2012: 8.0E-6, this code: {0:1.1e}'.format((np.sum(mc_z_s_h2[:,ind_co]*n_z_h2)/np.sum(n_z_h2))))
    print('C2H6. Hu+2012: 4.7E-10, this code: {0:1.1e}'.format((np.sum(mc_z_s_h2[:,ind_c2h6]*n_z_h2)/np.sum(n_z_h2))))
    print('S8. Hu+2012: 3.5E-10, this code: {0:1.1e}'.format((np.sum(mc_z_s_h2[:,ind_s8]*n_z_h2)/np.sum(n_z_h2) + np.sum(mc_z_s_h2[:,ind_s8a]*n_z_h2)/np.sum(n_z_h2))))
    print('CH2O. Hu+2012: 2.9E-10, this code: {0:1.1e}'.format((np.sum(mc_z_s_h2[:,ind_ch2o]*n_z_h2)/np.sum(n_z_h2))))
    print('CH4O. Hu+2012: 5.6E-11, this code: {0:1.1e}'.format((np.sum(mc_z_s_h2[:,ind_ch4o]*n_z_h2)/np.sum(n_z_h2))))
    print('Reactive Agents')
    print('H. Hu+2012: 2.1E-9, this code: {0:1.1e}'.format((np.sum(mc_z_s_h2[:,ind_h]*n_z_h2)/np.sum(n_z_h2))))
    print('OH. Hu+2012: 2.0E-14, this code: {0:1.1e}'.format((np.sum(mc_z_s_h2[:,ind_oh]*n_z_h2)/np.sum(n_z_h2))))
    print('O. Hu+2012: 1.2E-11, this code: {0:1.1e}'.format((np.sum(mc_z_s_h2[:,ind_o]*n_z_h2)/np.sum(n_z_h2))))
    print('O(1D). Hu+2012: 2.2E-21, this code: {0:1.1e}'.format((np.sum(mc_z_s_h2[:,ind_o1d]*n_z_h2)/np.sum(n_z_h2))))
    print('')
    print('N2')
    print('Emitted Gases')
    print('CO2. Hu+2012: 1.3E-4, this code: {0:1.1e}'.format((np.sum(mc_z_s_n2[:,ind_co2]*n_z_n2)/np.sum(n_z_n2))))
    print('H2. Hu+2012: 4.5E-4, this code: {0:1.1e}'.format((np.sum(mc_z_s_n2[:,ind_h2]*n_z_n2)/np.sum(n_z_n2))))
    print('SO2. Hu+2012: 8.9E-12, this code: {0:1.1e}'.format((np.sum(mc_z_s_n2[:,ind_so2]*n_z_n2)/np.sum(n_z_n2))))
    print('CH4. Hu+2012: 3.1E-5, this code: {0:1.1e}'.format((np.sum(mc_z_s_n2[:,ind_ch4]*n_z_n2)/np.sum(n_z_n2))))
    print('H2S. Hu+2012: 1.1E-14, this code: {0:1.1e}'.format((np.sum(mc_z_s_n2[:,ind_h2s]*n_z_n2)/np.sum(n_z_n2))))
    print('Photochemical Products')
    print('CO. Hu+2012: 1.7E-7, this code: {0:1.1e}'.format((np.sum(mc_z_s_n2[:,ind_co]*n_z_n2)/np.sum(n_z_n2))))
    print('C2H6. Hu+2012: 9.0E-9, this code: {0:1.1e}'.format((np.sum(mc_z_s_n2[:,ind_c2h6]*n_z_n2)/np.sum(n_z_n2))))
    print('CH4O. Hu+2012: 1.4E-9, this code: {0:1.1e}'.format((np.sum(mc_z_s_n2[:,ind_ch4o]*n_z_n2)/np.sum(n_z_n2))))
    print('O2. Hu+2012: 3.4E-10, this code: {0:1.1e}'.format((np.sum(mc_z_s_n2[:,ind_o2]*n_z_n2)/np.sum(n_z_n2))))
    print('S8. Hu+2012: 3.0E-10, this code: {0:1.1e}'.format((np.sum(mc_z_s_n2[:,ind_s8]*n_z_n2)/np.sum(n_z_n2) + np.sum(mc_z_s_n2[:,ind_s8a]*n_z_n2)/np.sum(n_z_n2))))
    print('CH2O. Hu+2012: 4.0E-11, this code: {0:1.1e}'.format((np.sum(mc_z_s_n2[:,ind_ch2o]*n_z_n2)/np.sum(n_z_n2))))
    print('C2H2. Hu+2012: 1.5E-11, this code: {0:1.1e}'.format((np.sum(mc_z_s_n2[:,ind_c2h2]*n_z_n2)/np.sum(n_z_n2))))
    print('Reactive Agents')
    print('H. Hu+2012: 1.2E-9, this code: {0:1.1e}'.format((np.sum(mc_z_s_n2[:,ind_h]*n_z_n2)/np.sum(n_z_n2))))
    print('OH. Hu+2012: 9.3E-14, this code: {0:1.1e}'.format((np.sum(mc_z_s_n2[:,ind_oh]*n_z_n2)/np.sum(n_z_n2))))
    print('O. Hu+2012: 6.5E-10, this code: {0:1.1e}'.format((np.sum(mc_z_s_n2[:,ind_o]*n_z_n2)/np.sum(n_z_n2))))
    print('O(1D). Hu+2012: 1.8E-20, this code: {0:1.1e}'.format((np.sum(mc_z_s_n2[:,ind_o1d]*n_z_n2)/np.sum(n_z_n2))))
    print('')
    print('CO2')
    print('Emitted Gases')
    print('H2. Hu+2012: 1.0E-3, this code: {0:1.1e}'.format((np.sum(mc_z_s_co2[:,ind_h2]*n_z_co2)/np.sum(n_z_co2))))
    print('SO2. Hu+2012: 1.6E-10, this code: {0:1.1e}'.format((np.sum(mc_z_s_co2[:,ind_so2]*n_z_co2)/np.sum(n_z_co2))))
    print('CH4. Hu+2012: 3.7E-5, this code: {0:1.1e}'.format((np.sum(mc_z_s_co2[:,ind_ch4]*n_z_co2)/np.sum(n_z_co2))))
    print('H2S. Hu+2012: 1.4E-10, this code: {0:1.1e}'.format((np.sum(mc_z_s_co2[:,ind_h2s]*n_z_co2)/np.sum(n_z_co2))))
    print('Photochemical Products')
    print('CO. Hu+2012: 7.7E-3, this code: {0:1.1e}'.format((np.sum(mc_z_s_co2[:,ind_co]*n_z_co2)/np.sum(n_z_co2))))
    print('O2. Hu+2012: 6.4E-7, this code: {0:1.1e}'.format((np.sum(mc_z_s_co2[:,ind_o2]*n_z_co2)/np.sum(n_z_co2))))
    print('C2H6. Hu+2012: 6.1E-10, this code: {0:1.1e}'.format((np.sum(mc_z_s_co2[:,ind_c2h6]*n_z_co2)/np.sum(n_z_co2))))
    print('H2O2. Hu+2012: 3.7E-10, this code: {0:1.1e}'.format((np.sum(mc_z_s_co2[:,ind_h2o2]*n_z_co2)/np.sum(n_z_co2))))
    print('H2SO4. Hu+2012: 5.0E-11, this code: {0:1.1e}'.format((np.sum(mc_z_s_co2[:,ind_h2so4]*n_z_co2)/np.sum(n_z_co2) + np.sum(mc_z_s_co2[:,ind_h2so4a]*n_z_co2)/np.sum(n_z_co2))))
    print('CH2O. Hu+2012: 2.5E-12, this code: {0:1.1e}'.format((np.sum(mc_z_s_co2[:,ind_ch2o]*n_z_co2)/np.sum(n_z_co2))))
    print('Reactive Agents')
    print('H. Hu+2012: 6.0E-11, this code: {0:1.1e}'.format((np.sum(mc_z_s_co2[:,ind_h]*n_z_co2)/np.sum(n_z_co2))))
    print('OH. Hu+2012: 7.8E-15, this code: {0:1.1e}'.format((np.sum(mc_z_s_co2[:,ind_oh]*n_z_co2)/np.sum(n_z_co2))))
    print('O. Hu+2012: 2.0E-8, this code: {0:1.1e}'.format((np.sum(mc_z_s_co2[:,ind_o]*n_z_co2)/np.sum(n_z_co2))))
    print('O(1D). Hu+2012: 3.0E-18, this code: {0:1.1e}'.format((np.sum(mc_z_s_co2[:,ind_o1d]*n_z_co2)/np.sum(n_z_co2))))

    ########################
    ###Plot
    ########################
    
    ###Initialize plot
    fig2, ax=plt.subplots(3,2, figsize=(7.5, 10.), sharex=True)
    markersizeval=5.
    
    ###Top plot:
    ax[0,0].plot(mc_z_s_h2[:,ind_h2], P_z_h2, linewidth=2, linestyle='-', color='red', label='H2')
    ax[0,0].plot(mc_z_s_h2[:,ind_h2o], P_z_h2, linewidth=2, linestyle='-', color='blue', label='H2O')
    ax[0,0].plot(mc_z_s_h2[:,ind_h], P_z_h2, linewidth=2, linestyle='--', color='red', label='H')
    ax[0,0].plot(mc_z_s_h2[:,ind_o2], P_z_h2, linewidth=2, linestyle='-', color='yellowgreen', label='O2')
    ax[0,0].plot(mc_z_s_h2[:,ind_o], P_z_h2, linewidth=2, linestyle='--', color='yellowgreen', label='O')
    ax[0,0].plot(mc_z_s_h2[:,ind_o3], P_z_h2, linewidth=2, linestyle=':', color='yellowgreen', label='O3')
    ax[0,0].plot(mc_z_s_h2[:,ind_oh], P_z_h2, linewidth=2, linestyle=':', color='blue', label='OH')
    #
    ax[0,0].set_title('H2')
    ax[0,0].set_yscale('log')
    ax[0,0].set_ylim([1.e+5, 1.e-1])
    ax[0,0].invert_yaxis()
    ax[0,0].set_ylabel('Pressure (Pa)')
    ax[0,0].set_xlabel('Mixing Ratio')  
    ax[0,0].legend(ncol=1, borderaxespad=0., fontsize=10, bbox_to_anchor=(2.6,1))    
    ax[0,0].tick_params(top=True, right=True)

    ax[0,1].plot(mc_z_s_h2[:,ind_n2], P_z_h2, linewidth=2, linestyle='-', color='skyblue', label='N2')
    ax[0,1].plot(mc_z_s_h2[:,ind_co2], P_z_h2, linewidth=2, linestyle='-', color='black', label='CO2')
    ax[0,1].plot(mc_z_s_h2[:,ind_co], P_z_h2, linewidth=2, linestyle='--', color='black', label='CO')
    ax[0,1].plot(mc_z_s_h2[:,ind_ch4], P_z_h2, linewidth=2, linestyle=':', color='black', label='CH4')
    ax[0,1].plot(mc_z_s_h2[:,ind_h2s], P_z_h2, linewidth=2, linestyle='-', color='brown', label='H2S')
    ax[0,1].plot(mc_z_s_h2[:,ind_so2], P_z_h2, linewidth=2, linestyle='--', color='brown', label='SO2')
    ax[0,1].plot(mc_z_s_h2[:,ind_h2so4a], P_z_h2, linewidth=2, linestyle=':', color='brown', label='H2SO4A')
    ax[0,1].plot(mc_z_s_h2[:,ind_s], P_z_h2, linewidth=2, linestyle='-', color='purple', label='S')
    ax[0,1].plot(mc_z_s_h2[:,ind_s2], P_z_h2, linewidth=2, linestyle='--', color='purple', label='S2')
    ax[0,1].plot(mc_z_s_h2[:,ind_s8a], P_z_h2, linewidth=2, linestyle=':', color='purple', label='S8A')
    #
    ax[0,1].set_title('H2')
    ax[0,1].set_yscale('log')
    ax[0,1].set_ylim([1.e+5, 1.e-1])
    ax[0,1].invert_yaxis()
    ax[0,1].set_ylabel('Pressure (Pa)')
    ax[0,1].set_xlabel('Mixing Ratio')  
#    ax[0,1].legend(loc='best', ncol=1, borderaxespad=0., fontsize=12)    
    ax[0,1].tick_params(top=True, right=True)
    
    ####Middle Plot:
    ax[1,0].plot(mc_z_s_n2[:,ind_h2], P_z_n2, linewidth=2, linestyle='-', color='red', label='H2')
    ax[1,0].plot(mc_z_s_n2[:,ind_h2o], P_z_n2, linewidth=2, linestyle='-', color='blue', label='H2O')
    ax[1,0].plot(mc_z_s_n2[:,ind_h], P_z_n2, linewidth=2, linestyle='--', color='red', label='H')
    ax[1,0].plot(mc_z_s_n2[:,ind_o2], P_z_n2, linewidth=2, linestyle='-', color='yellowgreen', label='O2')
    ax[1,0].plot(mc_z_s_n2[:,ind_o], P_z_n2, linewidth=2, linestyle='--', color='yellowgreen', label='O')
    ax[1,0].plot(mc_z_s_n2[:,ind_o3], P_z_n2, linewidth=2, linestyle=':', color='yellowgreen', label='O3')
    ax[1,0].plot(mc_z_s_n2[:,ind_oh], P_z_n2, linewidth=2, linestyle=':', color='blue', label='OH')
    #
    ax[1,0].set_title('N2')
    ax[1,0].set_yscale('log')
    ax[1,0].set_ylim([1.e+5, 1.e-1])
    ax[1,0].invert_yaxis()
    ax[1,0].set_ylabel('Pressure (Pa)')
    ax[1,0].set_xlabel('Mixing Ratio')  
#    ax[1,0].legend(loc='best', ncol=1, borderaxespad=0., fontsize=12)    
    ax[1,0].tick_params(top=True, right=True)

    ax[1,1].plot(mc_z_s_n2[:,ind_n2], P_z_n2, linewidth=2, linestyle='-', color='skyblue', label='N2')
    ax[1,1].plot(mc_z_s_n2[:,ind_co2], P_z_n2, linewidth=2, linestyle='-', color='black', label='CO2')
    ax[1,1].plot(mc_z_s_n2[:,ind_co], P_z_n2, linewidth=2, linestyle='--', color='black', label='CO')
    ax[1,1].plot(mc_z_s_n2[:,ind_ch4], P_z_n2, linewidth=2, linestyle=':', color='black', label='CH4')
    ax[1,1].plot(mc_z_s_n2[:,ind_h2s], P_z_n2, linewidth=2, linestyle='-', color='brown', label='H2S')
    ax[1,1].plot(mc_z_s_n2[:,ind_so2], P_z_n2, linewidth=2, linestyle='--', color='brown', label='SO2')
    ax[1,1].plot(mc_z_s_n2[:,ind_h2so4a], P_z_n2, linewidth=2, linestyle=':', color='brown', label='H2SO4A')
    ax[1,1].plot(mc_z_s_n2[:,ind_s], P_z_n2, linewidth=2, linestyle='-', color='purple', label='S')
    ax[1,1].plot(mc_z_s_n2[:,ind_s2], P_z_n2, linewidth=2, linestyle='--', color='purple', label='S2')
    ax[1,1].plot(mc_z_s_n2[:,ind_s8a], P_z_n2, linewidth=2, linestyle=':', color='purple', label='S8A')
    #
    ax[1,1].set_title('N2')
    ax[1,1].set_yscale('log')
    ax[1,1].set_ylim([1.e+5, 1.e-1])
    ax[1,1].invert_yaxis()
    ax[1,1].set_ylabel('Pressure (Pa)')
    ax[1,1].set_xlabel('Mixing Ratio')  
    ax[1,1].legend(ncol=1, borderaxespad=0., fontsize=10, bbox_to_anchor=(1.001,1))   
    ax[1,1].tick_params(top=True, right=True)

    ####Bottom Plot:
    ax[2,0].plot(mc_z_s_co2[:,ind_h2], P_z_co2, linewidth=2, linestyle='-', color='red', label='H2')
    ax[2,0].plot(mc_z_s_co2[:,ind_h2o], P_z_co2, linewidth=2, linestyle='-', color='blue', label='H2O')
    ax[2,0].plot(mc_z_s_co2[:,ind_h], P_z_co2, linewidth=2, linestyle='--', color='red', label='H')
    ax[2,0].plot(mc_z_s_co2[:,ind_o2], P_z_co2, linewidth=2, linestyle='-', color='yellowgreen', label='O2')
    ax[2,0].plot(mc_z_s_co2[:,ind_o], P_z_co2, linewidth=2, linestyle='--', color='yellowgreen', label='O')
    ax[2,0].plot(mc_z_s_co2[:,ind_o3], P_z_co2, linewidth=2, linestyle=':', color='yellowgreen', label='O3')
    ax[2,0].plot(mc_z_s_co2[:,ind_oh], P_z_co2, linewidth=2, linestyle=':', color='blue', label='OH')
    #
    ax[2,0].set_title('CO2')
    ax[2,0].set_yscale('log')
    ax[2,0].set_ylim([1.e+5, 1.e-1])
    ax[2,0].invert_yaxis()
    ax[2,0].set_ylabel('Pressure (Pa)')
    ax[2,0].set_xlabel('Mixing Ratio')  
#    ax[2,0].legend(loc='best', ncol=1, borderaxespad=0., fontsize=12)    
    ax[2,0].tick_params(top=True, right=True)

    ax[2,1].plot(mc_z_s_co2[:,ind_n2], P_z_co2, linewidth=2, linestyle='-', color='skyblue', label='N2')
    ax[2,1].plot(mc_z_s_co2[:,ind_co2], P_z_co2, linewidth=2, linestyle='-', color='black', label='CO2')
    ax[2,1].plot(mc_z_s_co2[:,ind_co], P_z_co2, linewidth=2, linestyle='--', color='black', label='CO')
    ax[2,1].plot(mc_z_s_co2[:,ind_ch4], P_z_co2, linewidth=2, linestyle=':', color='black', label='CH4')
    ax[2,1].plot(mc_z_s_co2[:,ind_h2s], P_z_co2, linewidth=2, linestyle='-', color='brown', label='H2S')
    ax[2,1].plot(mc_z_s_co2[:,ind_so2], P_z_co2, linewidth=2, linestyle='--', color='brown', label='SO2')
    ax[2,1].plot(mc_z_s_co2[:,ind_h2so4a], P_z_co2, linewidth=2, linestyle=':', color='brown', label='H2SO4A')
    ax[2,1].plot(mc_z_s_co2[:,ind_s], P_z_co2, linewidth=2, linestyle='-', color='purple', label='S')
    ax[2,1].plot(mc_z_s_co2[:,ind_s2], P_z_co2, linewidth=2, linestyle='--', color='purple', label='S2')
    ax[2,1].plot(mc_z_s_co2[:,ind_s8a], P_z_co2, linewidth=2, linestyle=':', color='purple', label='S8A')
    #
    ax[2,1].set_title('CO2')
    ax[2,1].set_yscale('log')
    ax[2,1].set_ylim([1.e+5, 1.e-1])
    ax[2,1].invert_yaxis()
    ax[2,1].set_ylabel('Pressure (Pa)')
    ax[2,1].set_xlabel('Mixing Ratio')  
#    ax[2,1].legend(loc='best', ncol=1, borderaxespad=0., fontsize=12)   
    ax[2,1].tick_params(top=True, right=True)
    
    ax[2,1].set_xscale('log')
    ax[2,1].set_xlim([1.0E-14, 1.0E0])
    
    ax[2,1].tick_params(top=True, right=True)
    plt.subplots_adjust(hspace=0.2, wspace=0.3)
    plt.savefig('./Plots/plot'+name+'.pdf', orientation='portrait', format='pdf') #,papertype='letter'
    plt.show()
    
plot_comparison('./scenario_library/Sun/H2-Full/ConcentrationSTD.dat', './scenario_library/Sun/N2-Full/ConcentrationSTD.dat', './scenario_library/Sun/CO2-Full/ConcentrationSTD.dat', 'Hu_Benchmark_Scenarios')

