"""
Purpose of this code is to see if we have reproduced Hu+2012 Figure 3 and Table 4 (i.e. if we correctly reproduce their modern Mars calculations)
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
deg2rad=np.pi/180.
bar2barye=1.e6 #1 Bar in Barye (the cgs unit of pressure)
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

###For Plotting
ind_co=20-1
ind_o2=54-1
ind_o3=2-1
ind_o=1-1
ind_h2o=7-1
ind_oh=4-1 #OH
ind_h2o2=6-1
ind_ho2=5-1
ind_h2=53-1
ind_h= 3-1




def plot_modMars_comparison(base_file, new_file, name):
    """
    #Base file. Eventually: original Hu+2012 result, for exact comparison.
    #Title of plot and name of file. 
    """
    
    ###Initialize plot
    fig2, ax=plt.subplots(1, figsize=(8., 8.))
    markersizeval=5.

    ########################
    ###Read in base data
    ########################
    base_data=np.genfromtxt(base_file, skip_header=2, skip_footer=0, unpack=False) #Import mapping between numerical ID in code and species name.
    
    z_center_base=base_data[:,0] # Center of altitude bins, km 
    deltazs=base_data[:,2]-base_data[:,1] #width of altitude bins, km
    T_z_base=base_data[:,3] # Temperature(z), in K
    P_z_base=base_data[:,4]*Pa2bar*bar2barye # Pressure(z), in Pa converted to Barye
    n_z_s_base=base_data[:,5:] #Number concentrations of the 111 chemical species, in cm**-3, as a function of (altitude, species)

    ###Get molar concentrations
    ###NOTE: May (probably) need to exclude condensed-phase species for molar concentration calculation...probably doesn't matter most of the time, but formally required and mioght matter in some weird edge cases.
    n_z_base=P_z_base/(k*T_z_base)

    mc_z_s_base=np.zeros(np.shape(n_z_s_base))
    num_s=np.shape(n_z_s_base)[1]

    for ind2 in range(0, num_s):
        mc_z_s_base[:,ind2]=n_z_s_base[:,ind2]/n_z_base#molar concentration of each species.

    ########################
    ###Read in new data
    ########################    
    new_data=np.genfromtxt(new_file, skip_header=2, skip_footer=0, unpack=False) #Import mapping between numerical ID in code and species name.

    z_center_new=new_data[:,0] # Center of altitude bins, km
    T_z_new=new_data[:,3] # Temperature(z), in K
    P_z_new=new_data[:,4]*Pa2bar*bar2barye # Pressure(z), in Pa converted to Barye
    n_z_s_new=new_data[:,5:] #Number concentrations of the 111 chemical species, in cm**-3, as a function of (altitude, species)
    
    
    ###Get molar concentrations
    n_z_new=P_z_new/(k*T_z_new)

    mc_z_s_new=np.zeros(np.shape(n_z_s_new))
    num_s=np.shape(n_z_s_new)[1]

    for ind2 in range(0, num_s):
        mc_z_s_new[:,ind2]=n_z_s_new[:,ind2]/n_z_new#molar concentration of each species.

    ########################
    ###Plot
    ########################
    ax.plot(mc_z_s_base[:,ind_co], z_center_base, linewidth=2, linestyle='-', color='black', label='CO')
    ax.plot(mc_z_s_base[:,ind_o2], z_center_base, linewidth=2, linestyle='-', color='lightgreen', label='O2')
    ax.plot(mc_z_s_base[:,ind_o3], z_center_base, linewidth=2, linestyle='-', color='green', label='O3')
    ax.plot(mc_z_s_base[:,ind_o], z_center_base, linewidth=2, linestyle='-', color='red', label='O')
    ax.plot(mc_z_s_base[:,ind_h2o], z_center_base, linewidth=2, linestyle='-', color='pink', label='H2O')
    ax.plot(mc_z_s_base[:,ind_oh], z_center_base, linewidth=2, linestyle='-', color='blue', label='OH')
    ax.plot(mc_z_s_base[:,ind_h2o2], z_center_base, linewidth=2, linestyle='-', color='magenta', label='H2O2')
    ax.plot(mc_z_s_base[:,ind_ho2], z_center_base, linewidth=2, linestyle='-', color='brown', label='HO2')
    ax.plot(mc_z_s_base[:,ind_h2], z_center_base, linewidth=2, linestyle='-', color='cyan', label='H2')
    ax.plot(mc_z_s_base[:,ind_h], z_center_base, linewidth=2, linestyle='-', color='orange', label='H')
    
    ax.plot(mc_z_s_new[:,ind_co], z_center_new, linewidth=2, linestyle='--', color='black')
    ax.plot(mc_z_s_new[:,ind_o2], z_center_new, linewidth=2, linestyle='--', color='lightgreen')
    ax.plot(mc_z_s_new[:,ind_o3], z_center_new, linewidth=2, linestyle='--', color='green')
    ax.plot(mc_z_s_new[:,ind_o], z_center_new, linewidth=2, linestyle='--', color='red')
    ax.plot(mc_z_s_new[:,ind_h2o], z_center_new, linewidth=2, linestyle='--', color='pink')
    ax.plot(mc_z_s_new[:,ind_oh], z_center_new, linewidth=2, linestyle='--', color='blue')
    ax.plot(mc_z_s_new[:,ind_h2o2], z_center_new, linewidth=2, linestyle='--', color='magenta')
    ax.plot(mc_z_s_new[:,ind_ho2], z_center_new, linewidth=2, linestyle='--', color='brown')
    ax.plot(mc_z_s_new[:,ind_h2], z_center_new, linewidth=2, linestyle='--', color='cyan')
    ax.plot(mc_z_s_new[:,ind_h], z_center_new, linewidth=2, linestyle='--', color='orange')
    
    ax.legend(ncol=1, loc='upper right')
    ax.set_xscale('log')
    ax.set_xlabel('Mixing Ratio')
    ax.set_xlim([1.e-15, 1.])
    ax.set_yscale('linear')
    ax.set_ylabel('Altitude (km)')
    ax.set_ylim([0., 120.])
    
    
    plt.savefig('./Plots/plot_'+name+'.pdf', orientation='portrait', format='pdf') # ,papertype='letter'
    plt.show()
    
    ###Get column-averaged mixing ratios for comparison
    print('Column-Averaged Mixing Ratios:')
    print('O2 (ppm). Measured: 1200-2000, Hu+2012: 1545, this code: {0:4.0f}'.format(1.e+6*np.sum(mc_z_s_new[:,ind_o2]*n_z_new)/np.sum(n_z_new)))
    print('CO (ppm). Measured: 800, Hu+2012: 572, this code: {0:4.0f}'.format(1.e+6*np.sum(mc_z_s_new[:,ind_co]*n_z_new)/np.sum(n_z_new)))
    print('H2 (ppm). Measured: 17, Hu+2012: 23, this code: {0:4.0f}'.format(1.e+6*np.sum(mc_z_s_new[:,ind_h2]*n_z_new)/np.sum(n_z_new)))
    print('H2O2 (ppb). Measured: 0-40, Hu+2012: 18, this code: {0:4.0f}'.format(1.e+9*np.sum(mc_z_s_new[:,ind_h2o2]*n_z_new)/np.sum(n_z_new)))
    print('O3 (ppb). Measured: 0-120, Hu+2012: 18, this code: {0:4.0f}'.format(1.e+9*np.sum(mc_z_s_new[:,ind_o3]*n_z_new)/np.sum(n_z_new)))


plot_modMars_comparison('./scenario_library/Mars/ConcentrationSTD_base.dat', './scenario_library/Mars/ConcentrationSTD.dat', 'Modern_Mars') #
