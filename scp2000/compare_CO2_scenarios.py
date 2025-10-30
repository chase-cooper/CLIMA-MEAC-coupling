"""
Purpose of this code is to reproduce Figure 6 and part of Table 7 of Hu+2012
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
Pa2Ba=10. #1 Pa in Ba
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

def plot_comparison(fullemission_file, reducedemission_file, noemission_file, name):
    """
    """
    ########################
    ###Read in full emission (f) data
    ########################
    f_data=np.genfromtxt(fullemission_file, skip_header=2, skip_footer=0, unpack=False) #Import mapping between numerical ID in code and species name.
    
    z_center_f=f_data[:,0] # Center of altitude bins, km 
    T_z_f=f_data[:,3] # Temperature(z), in K
    P_z_f=f_data[:,4]# Pressure(z), in Pa
    n_z_s_f=f_data[:,5:] #Number concentrations of the 111 chemical species, in cm**-3, as a function of (altitude, species)

    ###Get molar concentrations
    ###NOTE: May (probably) need to exclude condensed-phase species for molar concentration calculation...probably doesn't matter most of the time, but formally required and mioght matter in some weird edge cases.
    n_z_f=P_z_f*Pa2Ba/(k*T_z_f)

    mc_z_s_f=np.zeros(np.shape(n_z_s_f))
    num_s=np.shape(n_z_s_f)[1]

    for ind2 in range(0, num_s):
        mc_z_s_f[:,ind2]=n_z_s_f[:,ind2]/n_z_f#molar concentration of each species.

    ########################
    ###Read in reduced emission (r) data
    ########################
    r_data=np.genfromtxt(reducedemission_file, skip_header=2, skip_footer=0, unpack=False) #Import mapping between numerical ID in code and species name.
    
    z_center_r=r_data[:,0] # Center of altitude bins, km 
    T_z_r=r_data[:,3] # Temperature(z), in K
    P_z_r=r_data[:,4] # Pressure(z), in Pa 
    n_z_s_r=r_data[:,5:] #Number concentrations of the 111 chemical species, in cm**-3, as a function of (altitude, species)

    ###Get molar concentrations
    ###NOTE: May (probably) need to exclude condensed-phase species for molar concentration calculation...probably doesn't matter most of the time, but formally required and mioght matter in some weird edge cases.
    n_z_r=P_z_r*Pa2Ba/(k*T_z_r)

    mc_z_s_r=np.zeros(np.shape(n_z_s_r))
    num_s=np.shape(n_z_s_r)[1]

    for ind2 in range(0, num_s):
        mc_z_s_r[:,ind2]=n_z_s_r[:,ind2]/n_z_r#molar concentration of each species.
        
        
    ########################
    ###Read in no emission (n) data
    ########################
    n_data=np.genfromtxt(noemission_file, skip_header=2, skip_footer=0, unpack=False) #Import mapping between numerical ID in code and species name.
    
    z_center_n=n_data[:,0] # Center of altitude bins, km 
    T_z_n=n_data[:,3] # Temperature(z), in K
    P_z_n=n_data[:,4] # Pressure(z), in Pa 
    n_z_s_n=n_data[:,5:] #Number concentrations of the 111 chemical species, in cm**-3, as a function of (altitude, species)

    ###Get molar concentrations
    ###NOTE: May (probably) need to exclude condensed-phase species for molar concentration calculation...probably doesn't matter most of the time, but formally required and mioght matter in some weird edge cases.
    n_z_n=P_z_n*Pa2Ba/(k*T_z_n)

    mc_z_s_n=np.zeros(np.shape(n_z_s_n))
    num_s=np.shape(n_z_s_n)[1]

    for ind2 in range(0, num_s):
        mc_z_s_n[:,ind2]=n_z_s_n[:,ind2]/n_z_n#molar concentration of each species.
#    ########################
#    ###Print key parameters
#    ########################
#
    print('Column-Averaged Mixing Ratios: O2, O3 (Hu+2012 Table 7)')
    print('O2 (full emission). Hu+2012: 6.4E-7, this code: {0:1.1e}'.format((np.sum(mc_z_s_f[:,ind_o2]*n_z_f)/np.sum(n_z_f))))
    print('O3 (full emission). Hu+2012: 7.0E-11, this code: {0:1.1e}'.format((np.sum(mc_z_s_f[:,ind_o3]*n_z_f)/np.sum(n_z_f))))
    print('O2 (reduced emission). Hu+2012: 3.8E-6, this code: {0:1.1e}'.format((np.sum(mc_z_s_r[:,ind_o2]*n_z_r)/np.sum(n_z_r))))
    print('O3 (reduced emission). Hu+2012: 3.7E-10, this code: {0:1.1e}'.format((np.sum(mc_z_s_r[:,ind_o3]*n_z_r)/np.sum(n_z_r))))
    print('O2 (no H2, CH4 emission). Hu+2012: 1.3E-3, this code: {0:1.1e}'.format((np.sum(mc_z_s_n[:,ind_o2]*n_z_n)/np.sum(n_z_n))))
    print('O3 (no H2, CH4 emission). Hu+2012: 1.3E-7, this code: {0:1.1e}'.format((np.sum(mc_z_s_n[:,ind_o3]*n_z_n)/np.sum(n_z_n))))

    ########################
    ###Plot
    ########################
    
    ###Initialize plot
    fig2, ax=plt.subplots(2, figsize=(8., 10.), sharex=True)
    markersizeval=5.
    
    ###Top plot:
    ax[0].plot(mc_z_s_f[:,ind_h2o], P_z_f, linewidth=2, linestyle='-', color='blue', label='H2O')
    ax[0].plot(mc_z_s_r[:,ind_h2o], P_z_r, linewidth=2, linestyle='--', color='blue', label='H2O')
    ax[0].plot(mc_z_s_n[:,ind_h2o], P_z_n, linewidth=2, linestyle=':', color='blue', label='H2O')
    #
    ax[0].plot(mc_z_s_f[:,ind_co], P_z_f, linewidth=2, linestyle='-', color='brown', label='CO')
    ax[0].plot(mc_z_s_r[:,ind_co], P_z_r, linewidth=2, linestyle='--', color='brown', label='CO')
    ax[0].plot(mc_z_s_n[:,ind_co], P_z_n, linewidth=2, linestyle=':', color='brown', label='CO')
    #
    ax[0].plot(mc_z_s_f[:,ind_h2], P_z_f, linewidth=2, linestyle='-', color='skyblue', label='H2')
    ax[0].plot(mc_z_s_r[:,ind_h2], P_z_r, linewidth=2, linestyle='--', color='skyblue', label='H2')
    ax[0].plot(mc_z_s_n[:,ind_h2], P_z_n, linewidth=2, linestyle=':', color='skyblue', label='H2')
    #
    ax[0].plot(mc_z_s_f[:,ind_ch4], P_z_f, linewidth=2, linestyle='-', color='darkorchid', label='CH4')
    ax[0].plot(mc_z_s_r[:,ind_ch4], P_z_r, linewidth=2, linestyle='--', color='darkorchid', label='CH4')
    ax[0].plot(mc_z_s_n[:,ind_ch4], P_z_n, linewidth=2, linestyle=':', color='darkorchid', label='CH4')   
    
    ax[0].set_title('Effect of Variable Outgassing on CO2-dominated scenario')
    ax[0].set_yscale('log')
    ax[0].set_ylim([1.e+5, 1.e-1])
    ax[0].invert_yaxis()
    ax[0].set_ylabel('Pressure (Pa)')
    ax[0].set_xscale('log')
    ax[0].set_xlabel('Mixing Ratio')  
    ax[0].set_xlim([1.e-14, 1.e0])
    ax[0].legend(loc='best', ncol=1, borderaxespad=0., fontsize=12)    
   
    ####Bottom Plot:
    ax[1].plot(mc_z_s_f[:,ind_o2], P_z_f, linewidth=2, linestyle='-', color='black', label='O2')
    ax[1].plot(mc_z_s_r[:,ind_o2], P_z_r, linewidth=2, linestyle='--', color='black', label='O2')
    ax[1].plot(mc_z_s_n[:,ind_o2], P_z_n, linewidth=2, linestyle=':', color='black', label='O2')
    #
    ax[1].plot(mc_z_s_f[:,ind_o], P_z_f, linewidth=2, linestyle='-', color='yellowgreen', label='O')
    ax[1].plot(mc_z_s_r[:,ind_o], P_z_r, linewidth=2, linestyle='--', color='yellowgreen', label='O')
    ax[1].plot(mc_z_s_n[:,ind_o], P_z_n, linewidth=2, linestyle=':', color='yellowgreen', label='O')    
    #
    ax[1].plot(mc_z_s_f[:,ind_o3], P_z_f, linewidth=2, linestyle='-', color='red', label='O3')
    ax[1].plot(mc_z_s_r[:,ind_o3], P_z_r, linewidth=2, linestyle='--', color='red', label='O3')
    ax[1].plot(mc_z_s_n[:,ind_o3], P_z_n, linewidth=2, linestyle=':', color='red', label='O3')      
    
    ax[1].legend(loc=2, ncol=1, borderaxespad=0., fontsize=12)    
    ax[1].set_yscale('log')
    ax[1].set_ylim([1.e+5, 1.e-1])
    ax[1].invert_yaxis()
    ax[1].set_ylabel('Pressure (Pa)')
    ax[1].set_xscale('log')
    ax[1].set_xlabel('Mixing Ratio')  
    ax[1].set_xlim([1.e-14, 1.e-0])
   

    
    plt.savefig('./Plots/plot'+name+'.pdf', orientation='portrait',papertype='letter', format='pdf')
    plt.show()
    
plot_comparison('./scenario_library/scenario_CO2r_emi/ConcentrationSTD.dat', './scenario_library/scenario_CO2r/ConcentrationSTD.dat', './scenario_library/scenario_CO2r_NoE/ConcentrationSTD.dat', 'Hu_CO2_Outgassing')

