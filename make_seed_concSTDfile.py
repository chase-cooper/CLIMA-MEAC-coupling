"""
Purpose of this script is to generate new seed files for MEAC to enable climate-photochemistry coupling. Specifically, the goal is to take an old ConcentrationSTD file, a new TP profile, and generate a new concentrationSTD file with the old file's molar concentrations and the new TP profile. 
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
bar2barye=1.e+6 #1 Bar in Barye (the cgs unit of pressure)
barye2bar=1.e-6 #1 Barye in Bar
micron2m=1.e-6 #1 micron in m
micron2cm=1.e-4 #1 micron in cm
metricton2kg=1000. #1 metric ton in kg

#Fundamental constants
c=2.997924e10 #speed of light, cm/s
h=6.6260755e-27 #planck constant, erg/s
k=1.380658e-16 #boltzmann constant, erg/K #Exact same as in MEAC. 
sigma=5.67051e-5 #Stefan-Boltzmann constant, erg/(cm^2 K^4 s)
R_earth=6371.*km2m#radius of earth in m
R_sun=69.63e9 #radius of sun in cm
AU=1.496e13#1AU in cm

#Mean molecular masses, in amu #From the CRC handbook
m_h2o=18.015
m_co2=44.010
m_o3=47.998
m_ch4=16.043
m_o2=31.998
m_no2=46.006
m_n2=28.014
m_c2h6=30.069
m_h2=2.016


########################
###Establish key
########################

#Corrected for Hu 1-indexing vs Python 0-indexing
ind_o=1-1 #O
ind_h=3-1 #H
ind_oh=4-1 #OH

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
ind_cho=61-1


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
ind_n2o5=15-1
ind_hno4=70-1

def MEAC_interp(x_new, x_old, y_old):
    """
    Linear interpolation as nearly exactly as the MEAC call in Interpolation.c does it as I can easily make it.
    The x_new must be a subset (closed) of x_old.
    x_old, x_new are assumed to be strictly nondecreasing. 
    """
    y_new=np.zeros(np.shape(x_new))
    for ind in range(0, len(x_new)):
        x_new_val=x_new[ind] #this is the value we are interpolating to
        if x_new_val==x_old[0]: #if the interpolation value is at the lower limit of x_old
            y_new[ind]=y_old[0]
        elif x_new_val==x_old[-1]: #if the interpolation value is at the upper limit of x_old
            y_new[ind]=y_old[-1]
        else: #In this case, the interpolation value is somewhere in the middle. There is a value of x_old that is strictly below AND one that is strictly above this value. 
            j=np.squeeze(np.where(x_old<x_new_val)) #this is the index at which the old data is just under the interpolation value. so for [0, 1, 2, 3] and 0.5, it would return 0. 
            if j.ndim ==0: j=j.item()      # The line above sometimes returns a 0d array, this extracts the value. CC2025
            else: j=j[-1]
            # pdb.set_trace()
            y_new[ind]=y_old[j] + (x_new_val-x_old[j])*((y_old[j+1]-y_old[j])/(x_old[j+1]-x_old[j]))
    return y_new

def generate_MEACzTP_from_cloudyClima(cloudyClimafile, newfilename, zmin, zmax, zbin, g_planet):
    """
    cloudyClimafile: CloudyClima .atm file (including address) giving atmospheric zTP and mixing ratios.
    newfilename: new MEAC zTP file (including address) which will be generated. 
    zmin: lower altitude limit in km. Should match planet_scenario file of target MEAC simulation. 
    zmax: upper altitude limit in km. Should match planet_scenario file of target MEAC simulation. 
    zbin: number of altitude bins. Should match planet_scenario file of target MEAC simulation. 
    g_planet: gravity on planet in cm s^-2 (cgs).
    """
    ###Import CLIMA file
    data=np.genfromtxt(cloudyClimafile, skip_header=2, skip_footer=0, unpack=False) #Pa, K, km, v/v
    ##NOTE this is in REVERSE altitude order. 
    reversed_data=np.flip(data, axis=0)#Get it increasing altitude
    
    clima_P_z=reversed_data[:,0] #Pa
    clima_z=reversed_data[:,1] #km
    clima_T_z=reversed_data[:,2] #K
    clima_mr_h2o_z=reversed_data[:,3] #v/v
    clima_mr_co2_z=reversed_data[:,4] #v/v
    clima_mr_o3_z=reversed_data[:,5] #v/v
    clima_mr_ch4_z=reversed_data[:,6] #v/v
    clima_mr_o2_z=reversed_data[:,7] #v/v
    clima_mr_no2_z=reversed_data[:,8] #v/v
    clima_mr_n2_z=reversed_data[:,9] #v/v
    clima_mr_c2h6_z=reversed_data[:,10] #v/v
    clima_mr_h2_z=reversed_data[:,11] #v/v
    
    ##Remember that MEAC works in units of DRY pressure. So need to correct out the effects of H2O.
    #P_tot = P_dry + P_H2O = P_dry(1+PH2O/P_dry) = P_dry(1+mr_H2O) ---> P_dry=P_tot/(1+mr_H2O)
    clima_Pdry_z=clima_P_z/(1.0+clima_mr_h2o_z)
    #WARNING this does not pass order 0 by-hand test...something is weird in terms of how MR is being defined in CLIMA. Need to understand cloudy-Clima source code better. 
    
    ##Useful parameters in case of extrapolation
    clima_z_lowest=clima_z[0]
    clima_T_lowest=clima_T_z[0]# T at lowest altitude in grid, temperature extrapolation assumed isothermal from here.
    clima_Pdry_lowest=clima_Pdry_z[0]# Dry P in Pa at lowest altitude in grid
    clima_mu_lowest=(clima_mr_h2o_z[0]*m_h2o + clima_mr_co2_z[0]*m_co2 + clima_mr_o3_z[0]*m_o3 + clima_mr_ch4_z[0]*m_ch4 + clima_mr_o2_z[0]*m_o2 + clima_mr_no2_z[0]*m_no2 + clima_mr_n2_z[0]*m_n2 + clima_mr_c2h6_z[0]*m_c2h6 + clima_mr_h2_z[0]*m_h2)*amu2g#mmm in g of lowest altitude in clima grid
    H_lowest_km=(k*clima_T_lowest/(clima_mu_lowest*g_planet))*cm2km #scale height at bottom of atmosphere, converted to km. 
    
    clima_z_highest=clima_z[-1]
    clima_T_highest=clima_T_z[-1]# T at highest altitude in grid, temperature extrapolation assumed isothermal from here.
    clima_Pdry_highest=clima_Pdry_z[-1]# Dry P in Pa at highest altitude in grid
    clima_mu_highest=(clima_mr_h2o_z[-1]*m_h2o + clima_mr_co2_z[-1]*m_co2 + clima_mr_o3_z[-1]*m_o3 + clima_mr_ch4_z[-1]*m_ch4 + clima_mr_o2_z[-1]*m_o2 + clima_mr_no2_z[-1]*m_no2 + clima_mr_n2_z[-1]*m_n2 + clima_mr_c2h6_z[-1]*m_c2h6 + clima_mr_h2_z[-1]*m_h2)*amu2g#mmm in g of highest altitude in clima grid
    H_highest_km=(k*clima_T_highest/(clima_mu_highest*g_planet))*cm2km #scale height at bottom of atmosphere, converted to km. 

    
    ###Generate z_axis, initialize array to hold outputs
    z_grid_new=np.linspace(zmin, zmax, num=zbin+1, endpoint=True)
    towrite_tp=np.zeros([len(z_grid_new),3])
    
    ###Populate temperatures
    T_grid_new=np.interp(z_grid_new, clima_z, clima_T_z, left=clima_T_lowest, right=clima_T_highest) #assume isothermal at last value when extrapolating.
    
    ###Populate pressures
    ##This one is tricky because interpolation etc needs to be performed in log(P) space. Fortulately MEAC zTP is in log(P in Pa) anyway.
    clima_logPdry_z=np.log10(clima_Pdry_z) #log10(Pa)
    logPdry_grid_new=np.interp(z_grid_new, clima_z, clima_logPdry_z, left=0, right=0) #do not extrapolate. Need to do this carefully.
    
    ##It's possible we need to extrapolate if the new z grid is below or above the limits of the clima file. 
    #First, let's consider the case that the z_grid extends deeper than the clima file.

    belowinds=np.where(z_grid_new<clima_z[0]) #inds where z_grid_new is below the lowest altitude in clima. Should never be trigged, but let's be careful. 
    aboveinds=np.where(z_grid_new>clima_z[-1]) #inds where z_grid_new is above the highest altitude in clima. May well be trigged. 
    
    logPdry_grid_new[belowinds]=np.log10(clima_Pdry_lowest*np.exp(-(z_grid_new[belowinds]-clima_z_lowest)/H_lowest_km))
    logPdry_grid_new[aboveinds]=np.log10(clima_Pdry_highest*np.exp(-(z_grid_new[aboveinds]-clima_z_highest)/H_highest_km))
    
    #Populate final outputs
    towrite_tp[:,0]=z_grid_new
    towrite_tp[:,1]=logPdry_grid_new
    towrite_tp[:,2]=T_grid_new
    np.savetxt(newfilename, towrite_tp, delimiter=' ', fmt='%1.6f %1.6f %3.6f', newline='\n')


# g_Earth=(6.674E-11*5.9376E+24/(6371000.0)**2.0)*100.0  #m s^-2 --> cm s^-2. Using same constants as from MEAC to be sure of match. 
# generate_MEACzTP_from_cloudyClima('./Clima/001bar_CO2.atm', './Data/TP_clima_001barCO2.dat', 0.0, 120.0, 600, g_Earth)
# generate_MEACzTP_from_cloudyClima('./Clima/01bar_CO2.atm', './Data/TP_clima_01barCO2.dat', 0.0, 120.0, 600, g_Earth)
# generate_MEACzTP_from_cloudyClima('./Clima/1bar_CO2.atm', './Data/TP_clima_1barCO2.dat', 0.0, 120.0, 600, g_Earth)

def plot_climatoMEACTP(file_list, label_list, color_list):
    
    z_center={} # Center of altitude bins, km
    logP_z={} #log10(Pa(z) in Pa)
    T_z={} # Temperature(z), in K
    fig, ax=plt.subplots(2, figsize=(8., 10))
    markersizeval=5.
    
    for ind in range(0, len(label_list)):
        file=file_list[ind]
        label=label_list[ind]
        color=color_list[ind]
        
        z_center[label], logP_z[label], T_z[label]=np.genfromtxt(file, skip_header=0, skip_footer=0, unpack=True)
        
        ax[0].plot(T_z[label], z_center[label], label=label, linestyle='-', color=color)
        ax[1].plot(T_z[label], 10**(logP_z[label])*Pa2bar, label=label, linestyle='-', color=color)
    ax[0].legend(loc='best', fontsize=12)
    ax[0].set_xlabel('Temperature (K)')
    ax[0].set_ylabel('Altitude (km)')

    ax[1].set_ylabel('Dry Pressure (bar)')
    ax[1].set_xlabel('Temperature (K)')
    ax[1].set_yscale('log')
    ax[1].invert_yaxis()
    
    plt.show()
        
        

# plot_climatoMEACTP(['./Data/TPStd175288CO2N2_extended.dat', './Data/TPStd200288N2.dat', './Data/TP_clima_001barCO2.dat', './Data/TP_clima_01barCO2.dat', './Data/TP_clima_1barCO2.dat'],['MEAC_CO2N2','MEAC_N2','Clima_0.01barCO2','Clima_0.1barCO2','Clima_1barCO2'], ['gold', 'green', 'blue', 'purple', 'red'])

    
def generate_new_concentrationSTD(old_concSTD, new_zTP, zmin, zmax, zbin,new_concSTD, fixed_species_stdn=[], fixed_species_mrs=[]):
    """
    old_concSTD: old MEAC ConcentrationSTD.dat file (including address), from which mixing ratios will be taken.
    new_zTP: new MEAC zTP profile (including address), from which the new T-P profile will be taken (remember: dry pressure for MEAC!)
    zmin: lower altitude limit in km. Should match planet_scenario file of target MEAC simulation. 
    zmax: upper altitude limit in km. Should match planet_scenario file of target MEAC simulation. 
    zbin: number of altitude bins. Should match planet_scenario file of target MEAC simulation. 
    new_concSTD: new MEAC ConcentrationSTD.dat file (including address) to seed new simulation. 
    fixed_species_stdn: tuple containing standard numbers of chemical species to be assigned constant mixing ratio (e.g. 0.99 N2 throughout atm-->[55])
    fixed_species_mrs: tuple containing mixing ratios of chemical species to be assigned constant mixing ratio e.g. 0.99 N2 throughout atm-->[0.99] )
    """

    ########################
    ###Read in old ConcentrationSTD file
    ########################
    base_data=np.genfromtxt(old_concSTD, skip_header=2, skip_footer=0, unpack=False) #Import mapping between numerical ID in code and species name.
    
    z_center_base=base_data[:,0] # Center of altitude bins, km 
    T_z_base=base_data[:,3] # Temperature(z), in K
    P_z_base=base_data[:,4] # Pressure(z), in Pa 
    n_z_s_base=base_data[:,5:] #Number concentrations of the 111 chemical species, in cm**-3, as a function of (altitude, species)

    ###Get molar concentrations
    n_z_bulkatm_base=(P_z_base*Pa2bar*bar2barye)/(k*T_z_base) #cm^-3

    mr_z_s_base=np.zeros(np.shape(n_z_s_base))

    num_s=np.shape(n_z_s_base)[1]

    for ind2 in range(0, num_s):
        mr_z_s_base[:,ind2]=n_z_s_base[:,ind2]/n_z_bulkatm_base#mixing ratio of each species relative to bulk atmosphere.
    
    #If we are manually adjusting any of the mixing ratios in the new file, e.g. because we are moving to a new bulk atmospheric compostion.
    if len(fixed_species_stdn)>0:
        for ind in range(0, len(fixed_species_stdn)):
            stdn=fixed_species_stdn[ind]
            mr=fixed_species_mrs[ind]
            mr_z_s_base[:,stdn-1]*=0#zero the mixing ratios
            mr_z_s_base[:,stdn-1]+=mr#adjust it to the new mixing ratio.         
            
    ########################
    ###Read in new zTP file
    ########################    
    new_zTP_data=np.genfromtxt(new_zTP, skip_header=0, skip_footer=0, unpack=False) #Import new zTP profile

    z_zTP=new_zTP_data[:,0] # altitude, km
    logP_zTP=new_zTP_data[:,1] # Pressure(z), in log10(Pa)
    T_zTP=new_zTP_data[:,2] # Temperature(z), in K

    ########################
    ###Define new altitude grid, interpolated in EXACTLY the same way MEAC does. (lines 149-186)
    ########################
    #First, generate the new altitude grid.
    z_grid_new=np.linspace(zmin, zmax, num=zbin+1, endpoint=True)
    
    #Second, interpolate the new zTP to the new altitude grid.
    # T_new=np.interp(z_grid_new, z_zTP, T_zTP) 
    # logP_new=np.interp(z_grid_new, z_zTP, logP_zTP) #interpolation for pressure is performed in log(p) space.
    T_new=MEAC_interp(z_grid_new, z_zTP, T_zTP) 
    logP_new=MEAC_interp(z_grid_new, z_zTP, logP_zTP) #interpolation for pressure is performed in log(p) space.   
    
    P_Pa_new=10.0**logP_new #now, convert to pressure in Pa. 
    
    #Third, get layer-center values.
    z_left_new=np.zeros(len(z_grid_new)-1)
    z_center_new=np.zeros(len(z_grid_new)-1)
    z_right_new=np.zeros(len(z_grid_new)-1)
    P_z_new=np.zeros(len(z_grid_new)-1)
    T_z_new=np.zeros(len(z_grid_new)-1)
    n_z_bulkatm_new=np.zeros(len(z_grid_new)-1)
    
    for ind in range(0, len(z_grid_new)-1):
        z_left_new[ind]=z_grid_new[ind]
        z_right_new[ind]=z_grid_new[ind+1]
        z_center_new[ind]=0.5*(z_grid_new[ind]+z_grid_new[ind+1])
        T_z_new[ind]=0.5*(T_new[ind]+T_new[ind+1])
        P_z_new[ind]=np.sqrt(P_Pa_new[ind]*P_Pa_new[ind+1])
        
    # T_z_new=T_z_base
    # P_z_new=P_z_base
    n_z_bulkatm_new=(P_z_new*Pa2bar*bar2barye)/(k*T_z_new) #cm^-3
    
    #Now, input into concSTD file. 
    concSTD_new=np.zeros((len(z_center_new), 5+num_s))
    concSTD_new[:,0]=z_center_new
    concSTD_new[:,1]=z_left_new
    concSTD_new[:,2]=z_right_new
    concSTD_new[:,3]=T_z_new
    concSTD_new[:,4]=P_z_new

    for ind2 in range(0, num_s):
        concSTD_new[:,5+ind2]=n_z_bulkatm_new*np.interp(z_center_new, z_center_base, mr_z_s_base[:,ind2],left=mr_z_s_base[0,ind2], right=mr_z_s_base[-1,ind2])
    
    #Format and save
    header='z\t\t z0\t\t z1\t\t T\t\t P\t\t'
    fmt='%5.6f %5.6f %5.6f %5.6f %5.6e'
    for i in range(1, 111+1):
        header += (str(i)+'\t\t')
        fmt += ' %5.6e'
    header += '\n km\t\t km\t\t km\t\t K\t\t Pa\t\t'
    np.savetxt(new_concSTD, concSTD_new, delimiter='\t', newline='\n', header=header, fmt=fmt) #Print checkfile for reaction rates.

# generate_new_concentrationSTD('./scenario_library/Sun/N2-Full/ConcentrationSTD.dat', './Data/TPStd200288N2.dat', 0.0, 100.0, 50, './scenario_library/Sun/N2-Full-changePT/ConcentrationSTD.dat')
# generate_new_concentrationSTD('./scenario_library/Sun/N2-Full/ConcentrationSTD.dat', './Data/TPStd200288N2.dat', 0.0, 100.0, 50, './scenario_library/Sun/N2_CO2-Full-seeded/ConcentrationSTD.dat', [55, 52], [0.99, 0.01])

# generate_new_concentrationSTD('./scenario_library/Sun/N2_CO2-Full-seeded/ConcentrationSTD.dat', './Data/TP_clima_001barCO2.dat', 0.0, 100.0, 50, './scenario_library/Sun/N2_CO2-Full-climaTP-seeded/ConcentrationSTD.dat', [55, 52], [0.99, 0.01])

# generate_new_concentrationSTD('./scenario_library/Sun/N2_CO2-Full-climaTP-seeded/ConcentrationSTD.dat', './Data/TP_clima_01barCO2.dat', 0.0, 100.0, 50, './scenario_library/Sun/N2_CO2-Full-climaTP-seeded-01CO2/ConcentrationSTD.dat', [55, 52], [0.9, 0.1])

# generate_new_concentrationSTD('./scenario_library/Sun/N2_CO2-Full-climaTP-seeded-01CO2/ConcentrationSTD.dat', './Data/TP_clima_1barCO2.dat', 0.0, 100.0, 50, './scenario_library/Sun/N2_CO2-Full-climaTP-seeded-1CO2/ConcentrationSTD.dat', [55, 52], [0.0, 1.0])
