"""
Purpose of this code is to see if we have reproduced Hu+2012 Figure 2 and Table 3 (i.e. if we correctly reproduce their modern Earth calculations)
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
ind_oh=4-1 #OH
ind_o3=2-1
ind_n2o=11-1
ind_h2o=7-1
ind_ho2=5-1
ind_no=12-1
ind_no2=13-1
ind_hno3=18-1
ind_ch4=21-1

###For surface mixing ratios
ind_co= 20-1
ind_nh3= 9-1
ind_so2= 43-1
ind_ocs= 49-1
ind_h2s= 45-1
ind_h2so4= 73-1
ind_h2so4a= 78-1




def plot_modEarth_comparison(base_file, new_file, name):
    """
    #Base file. Eventually: original Hu+2012 result, for exact comparison. For now: duplicate new file. 
    #New file. Our calculation
    #Title of plot and name of file. 
    """
    
    ###Initialize plot
    fig2, ax=plt.subplots(3,3, figsize=(8., 8.))
    markersizeval=5.

    ########################
    ###Read in base data
    ########################
    base_data=np.genfromtxt(base_file, skip_header=2, skip_footer=0, unpack=False) #Import mapping between numerical ID in code and species name.
    
    z_center_base=base_data[:,0] # Center of altitude bins, km 
    T_z_base=base_data[:,3] # Temperature(z), in K
    P_z_base=base_data[:,4]*Pa2bar*bar2barye # Pressure(z), in Pa converted to Barye
    n_z_s_base=base_data[:,5:] #Number concentrations of the 111 chemical species, in cm**-3, as a function of (altitude, species)

    ###Get molar concentrations
    ###NOTE: May (probably) need to exclude condensed-phase species for molar concentration calculation...probably doesn't matter most of the time, but formally required and mioght matter in some weird edge cases.
    n_z_base=P_z_base/(k*T_z_base) #sum number densities across species. This is a profile for the whole atmosphere.

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
    n_z_new=P_z_new/(k*T_z_new) #sum number densities across species. This is a profile for the whole atmosphere.

    mc_z_s_new=np.zeros(np.shape(n_z_s_new))
    num_s=np.shape(n_z_s_new)[1]

    for ind2 in range(0, num_s):
        mc_z_s_new[:,ind2]=n_z_s_new[:,ind2]/n_z_new#molar concentration of each species.

    ########################
    ###Plot
    ########################
    linestyles=np.array(['-',':'])
    ###Plot O3
    ax[0,0].plot(mc_z_s_base[:,ind_o3], z_center_base, linewidth=2, linestyle=linestyles[0], color='black')
    ax[0,0].plot(mc_z_s_new[:,ind_o3], z_center_new, linewidth=2, linestyle=linestyles[1], color='red')
    ax[0,0].set_title('O3')
    ax[0,0].set_xscale('log')
    ax[0,0].set_xlabel('Mixing Ratio')
    ax[0,0].set_xlim([1.e-8, 1.e-4])
    ax[0,0].set_yscale('linear')
    ax[0,0].set_ylabel('Altitude (km)')
    ax[0,0].set_ylim([0., 80.])
    
    ###Plot N2O
    ax[0,1].plot(mc_z_s_base[:,ind_n2o], z_center_base, linewidth=2, linestyle=linestyles[0], color='black')
    ax[0,1].plot(mc_z_s_new[:,ind_n2o], z_center_new, linewidth=2, linestyle=linestyles[1], color='red')
    ax[0,1].set_title('N2O')
    ax[0,1].set_xscale('log')
    ax[0,1].set_xlabel('Mixing Ratio')
    ax[0,1].set_xlim([1.e-10, 1.e-6])
    ax[0,1].set_yscale('linear')
    ax[0,1].set_ylabel('Altitude (km)')
    ax[0,1].set_ylim([0., 80.]) 
    
    ###Plot CH4
    ax[0,2].plot(mc_z_s_base[:,ind_ch4], z_center_base, linewidth=2, linestyle=linestyles[0], color='black')
    ax[0,2].plot(mc_z_s_new[:,ind_ch4], z_center_new, linewidth=2, linestyle=linestyles[1], color='red')
    ax[0,2].set_title('CH4')
    ax[0,2].set_xscale('log')
    ax[0,2].set_xlabel('Mixing Ratio')
    ax[0,2].set_xlim([1.e-8, 1.e-5])
    ax[0,2].set_yscale('linear')
    ax[0,2].set_ylabel('Altitude (km)')
    ax[0,2].set_ylim([0., 80.])     
    
    ###Plot H2O
    ax[1,0].plot(n_z_s_base[:,ind_h2o], z_center_base, linewidth=2, linestyle=linestyles[0], color='black')
    ax[1,0].plot(n_z_s_new[:,ind_h2o], z_center_new, linewidth=2, linestyle=linestyles[1], color='red')
    ax[1,0].set_title('H2O')
    ax[1,0].set_xscale('log')
    ax[1,0].set_xlabel('Number Density (cm**-3)')
    ax[1,0].set_xlim([1.e+8, 1.e+18])
    ax[1,0].set_yscale('linear')
    ax[1,0].set_ylabel('Altitude (km)')
    ax[1,0].set_ylim([0., 80.])

     ###Plot OH
    ax[1,1].plot(n_z_s_base[:,ind_oh], z_center_base, linewidth=2, linestyle=linestyles[0], color='black')
    ax[1,1].plot(n_z_s_new[:,ind_oh], z_center_new, linewidth=2, linestyle=linestyles[1], color='red')
    ax[1,1].set_title('OH')
    ax[1,1].set_xscale('log')
    ax[1,1].set_xlabel('Number Density (cm**-3)')
    ax[1,1].set_xlim([1.e+5, 1.e+8])
    ax[1,1].set_yscale('linear')
    ax[1,1].set_ylabel('Altitude (km)')
    ax[1,1].set_ylim([0., 80.])

     ###Plot HO2
    ax[1,2].plot(n_z_s_base[:,ind_ho2], z_center_base, linewidth=2, linestyle=linestyles[0], color='black')
    ax[1,2].plot(n_z_s_new[:,ind_ho2], z_center_new, linewidth=2, linestyle=linestyles[1], color='red')
    ax[1,2].set_title('HO2')
    ax[1,2].set_xscale('log')
    ax[1,2].set_xlabel('Number Density (cm**-3)')
    ax[1,2].set_xlim([1.e+5, 1.e+9])
    ax[1,2].set_yscale('linear')
    ax[1,2].set_ylabel('Altitude (km)')
    ax[1,2].set_ylim([0., 80.])    
    
    ###Plot NO
    ax[2,0].plot(mc_z_s_base[:,ind_no], z_center_base, linewidth=2, linestyle=linestyles[0], color='black')
    ax[2,0].plot(mc_z_s_new[:,ind_no], z_center_new, linewidth=2, linestyle=linestyles[1], color='red')
    ax[2,0].set_title('NO')
    ax[2,0].set_xscale('log')
    ax[2,0].set_xlabel('Mixing Ratio')
    ax[2,0].set_xlim([1.e-12, 1.e-7])
    ax[2,0].set_yscale('linear')
    ax[2,0].set_ylabel('Altitude (km)')
    ax[2,0].set_ylim([0., 80.])
    
    ###Plot NO2
    ax[2,1].plot(mc_z_s_base[:,ind_no2], z_center_base, linewidth=2, linestyle=linestyles[0], color='black')
    ax[2,1].plot(mc_z_s_new[:,ind_no2], z_center_new, linewidth=2, linestyle=linestyles[1], color='red')
    ax[2,1].set_title('NO2')
    ax[2,1].set_xscale('log')
    ax[2,1].set_xlabel('Mixing Ratio')
    ax[2,1].set_xlim([1.e-14, 1.e-7])
    ax[2,1].set_yscale('linear')
    ax[2,1].set_ylabel('Altitude (km)')
    ax[2,1].set_ylim([0., 80.])
    
    ###Plot HNO3
    ax[2,2].plot(mc_z_s_base[:,ind_hno3], z_center_base, linewidth=2, linestyle=linestyles[0], color='black')
    ax[2,2].plot(mc_z_s_new[:,ind_hno3], z_center_new, linewidth=2, linestyle=linestyles[1], color='red')
    ax[2,2].set_title('HNO3')
    ax[2,2].set_xscale('log')
    ax[2,2].set_xlabel('Mixing Ratio')
    ax[2,2].set_xlim([1.e-18, 1.e-6])
    ax[2,2].set_yscale('linear')
    ax[2,2].set_ylabel('Altitude (km)')
    ax[2,2].set_ylim([0., 80.])
    
    plt.subplots_adjust(wspace=0.4, hspace=0.4)
    
    plt.savefig('./Plots/plot_'+name+'.pdf', orientation='portrait', format='pdf') #,papertype='letter'
    plt.show()
    
    print('Surface Mixing Ratios:')
    print('CO (ppb). Measured: 40-200, Hu+2012: 101-113, Current code: {0:3.0f}'.format(1.e+9*mc_z_s_new[0,ind_co]))
    print('CH4 (ppb). Measured: 700-1745, Hu+2012: 1235-1939, Current code: {0:3.0f}'.format(1.e+9*mc_z_s_new[0,ind_ch4]))
    print('NH3 (ppb). Measured: 0.1-10, Hu+2012: .24, Current code: {0:3.2f}'.format(1.e+9*mc_z_s_new[0,ind_nh3]))
    print('N2O (ppb). Measured: 276-315, Hu+2012: 290-302, Current code: {0:3.0f}'.format(1.e+9*mc_z_s_new[0,ind_n2o]))
    print('NO (ppb). Measured: 0.02-10, Hu+2012: 0.024-0.025, Current code: {0:3.3f}'.format(1.e+9*mc_z_s_new[0,ind_no]))
    print('SO2 (ppt). Measured: 30-260, Hu+2012: 237-239, Current code: {0:3.0f}'.format(1.e+12*mc_z_s_new[0,ind_so2]))
    print('OCS (ppt). Measured: 510, Hu+2012: 185-188, Current code: {0:3.0f}'.format(1.e+12*mc_z_s_new[0,ind_ocs]))
    print('H2S (ppt). Measured: 1-13, Hu+2012: 3.62-3.92, Current code: {0:3.2f}'.format(1.e+12*mc_z_s_new[0,ind_h2s]))
    print('H2SO4 (ppt). Measured*: 26-170, Hu+2012: 127, Current code: {0:3.1f}'.format(1.e+12*(mc_z_s_new[0,ind_h2so4]+mc_z_s_new[0,ind_h2so4a]))) #Tracing back the 5-70 ppt, I think there is an error in the source Hu cited -- the underlying source appears to be Warneck+1988, which gives 20-130 ng S m^-3 (as PARTICULATE SO4--) in the free troposphere. This converts to about 26-170 ppt PARTICULATE SO4--. That is within a factor of 2. 

plot_modEarth_comparison('./scenario_library/Earth/ConcentrationSTD_base_Earth.dat', './scenario_library/Earth/ConcentrationSTD.dat', 'Modern_Earth') #
