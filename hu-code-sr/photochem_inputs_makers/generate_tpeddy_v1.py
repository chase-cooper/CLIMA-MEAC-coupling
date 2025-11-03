"""
Purpose of this script is to generate new T-P-Kzz profiles for use with Hu code.

This will, at least initially, need to follow the basic profiles laid down by the Hu code. 
"""

import numpy as np
import matplotlib.pyplot as plt
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
Pascal2bar=1.e-5 #1 Pascal in bar
Pa2bar=1.e-5 #1 Pascal in bar
bar2Pa=1.e5 #1 bar in Pascal
deg2rad=np.pi/180.
bar2barye=1.e6 #1 Bar in Barye (the cgs unit of pressure)
barye2bar=1.e-6 #1 Barye in Bar
micron2m=1.e-6 #1 micron in m
micron2cm=1.e-4 #1 micron in cm
metricton2kg=1000. #1 metric ton in kg
AU=1.496e13#1AU in cm
m2cm=1.0E2 #1 m in cm
kg2g=1.0E3 #1 kg in g

#Fundamental constants
c=2.997924e10 #speed of light, cm/s
h=6.6260755e-27 #planck constant, erg/s
k_boltzmann=1.380658e-16 #boltzmann constant, erg/K
sigma=5.67051e-5 #Stefan-Boltzmann constant, erg/(cm^2 K^4 s)
G_const=6.6726e-8 #cm^3 g^-1 s^-2

#Planet parameters
g_earth=981. #surface gravity of Earth, cm s**-2
R_earth=6.3781e6*m2cm #cm; Pinned post from Ryan MacDonald #6371.*km2cm#radius of earth in cm
R_sun=69.63e9 #radius of sun in cm
M_earth=5.972e24*kg2g #5.96e27 #mass of Earth, in g. From https://www.prl.res.in/~snaik/const.html


#Mean molecular masses
m_co2=44.01*amu2g #co2, in g
m_h2=2.02*amu2g #h2, in g
m_n2=28.01*amu2g #n2, in g
m_h2o=18.02*amu2g #h2o, in g
m_o2=16.00*amu2g #o2, in g
m_ar=39.95*amu2g #ar, in g
m_air=28.964*amu2g # air, from CRC pg 14-19

#c_p, specific heat capacity at constant pressure. All at 0C/1 bar (PPC pg. 92, Table 2.1), converted from J Kg**-1 K**-1 to erg g**-1 K**-1.
c_p_co2= 820.*1.e4 #
c_p_h2= 14230.*1.e4 #
c_p_n2= 1037.*1.e4 #
c_p_o2=916.*1.e4 #
c_p_air=1.006*1.e7 #From CRC pg. 6-18; 280K, 1 bar, converged from kJ kg**-1 K**-1

########################
###Define core functions
########################
def T_z_dry_adiabat(z, T_0, g, c_p):
    """
    Returns temperature as a function of altitude for the dry adiabat. 
    
    Inputs:
    ---z: altitude in cm
    ---T_0: temperature at z=0, in K
    ---g: acceleration due to gravity, in cm s**-2
    ---c_p: specific heat capacity of atmosphere at constant pressure, in units of erg g**-1 K**-1
    
    Outputs: T(z), in K
    
    Uses: Equation 16.13 of Seinfeld & Pandis 2016
    """
    
    return T_0 - (g/c_p)*z

def P_T_dry_adiabat(T, T_0, P_0, c_p, mmm):
    """
    Returns pressure as a function of temperature for the dry adiabat.
    Inputs:
    ---T(z), in K
    ---T(0), in K
    ---P(0), in barye
    ---c_p: specific heat capacity of atmosphere at constant pressure in units of erg g**-1 K**-1.
    ---mmm: mean molecular mass of atmosphere, in g
    
    Uses: Equation 16.9 of Seinfeld & Pandis 2016
    """
    R_d=k_boltzmann/mmm #specific gas constant, in units of erg g**-1 K**-
    return P_0*(T/T_0)**(c_p/R_d)

def ntp_dryadiabat_isothermalstrat(z_list_km, mr_co2, mr_n2, mr_h2,mr_o2, P_0_bar, T_0, T_strat, grav):
    """
    z_list_km: list of altitudes, in km
    mr_co2: CO2 mixing ratio
    mr_n2: N2 mixing ratio
    mr_h2: H2 mixing ratio
    mr_o2: O2 mixing ratio.
    P_0_bar: surface pressure in bar
    T_0: surface temperature in K
    T_strat: stratospheric temperature in K
    grav: gravitational acceleration, in cm s**-2.
    """
    #Basic prep work.
    mmm = mr_co2*m_co2 + mr_n2*m_n2 + mr_h2*m_h2 + mr_o2*m_o2#Mean molecular mass in grams.
    c_p = mr_co2*c_p_co2 + mr_n2*c_p_n2 + mr_h2*c_p_h2 + mr_o2*c_p_o2 #specific heat capacity at constant pressure, ergs g**-1 K**-1.
    p_0=P_0_bar*bar2barye #surface pressure, converted to barye
    z_list=z_list_km*km2cm
    
    #Calculate stratosphere
    z_strat_0=(c_p/grav)*(T_0-T_strat) #At this altitude, transition to isothermal stratosphere (cm)
    p_strat_0=P_T_dry_adiabat(T_strat, T_0, p_0, c_p, mmm) #This is the pressure at which the transition to isothermal stratosphere begins (barye) 
    H_strat=k_boltzmann*T_strat/(mmm*grav) #scale height in cm

    #Create all-adiabat atmosphere
    T_z_atm=T_z_dry_adiabat(z_list, T_0, grav, c_p)
    T_z_atm[T_z_atm<0]=0 #Temp cannot go below absolute 0. 
    p_z_atm=P_T_dry_adiabat(T_z_atm, T_0, p_0, c_p, mmm)
    
    #Insert stratosphere
    strat_inds=np.where(T_z_atm<=T_strat)
    
    T_z_atm[strat_inds]=T_strat
    p_z_atm[strat_inds]=p_strat_0*np.exp(-(z_list[strat_inds]-z_strat_0)/H_strat)
    n_z_atm=p_z_atm/(k_boltzmann*T_z_atm)
        
    return n_z_atm, T_z_atm, p_z_atm

def make_earth_eddy(mr_co2, mr_n2, mr_h2, mr_o2, new_atm_z, new_atm_logp):
    """
    NOTE: THIS CODE ADAPTED FROM CODE WRITTEN BY SANGITA MANDAL UNDER GUIDANCE OF SUKRIT RANJAN DURING PROJECT INITIATED AT NISER IN FALL 2021
    
    new_atm_z, new_atm_logp: Because of how np.interp works, new_atm_logp needs to be strictly increasing. So flipud as needed.
    """
    ###Get mean molecular mass of specified atmospheric composition
    m_atm=mr_co2*m_co2 + mr_n2*m_n2 + mr_h2*m_h2 + mr_o2*m_o2
    

    ###Import TP profiles for the modern earth
    TP_modern_earth_z,TP_modern_earth_logp,TP_modern_earth_T=np.genfromtxt("./../Data/TP1986.dat", skip_header=0, skip_footer=0,usecols=(0,1,2), unpack=True) #units: km, log(Pa),K
    eddy_modern_earth_z,eddy_modern_earth_K=np.genfromtxt("./../Data/EddyEarth.dat", skip_header=0, skip_footer=0,usecols=(0,1), unpack=True)# km,cm^2sec^-1
    

    ###Scale the Kz(P), convert it to Kz(z) 
    eddy_modern_earth_logP=np.interp(eddy_modern_earth_z,TP_modern_earth_z,TP_modern_earth_logp,left=None, right=None, period=None) #logp at each gridpoint at which Kz(z) is known for modern Earth
    new_atm_eddy=(m_air/m_atm)*np.interp(new_atm_logp, np.flipud(eddy_modern_earth_logP), np.flipud(eddy_modern_earth_K), left=eddy_modern_earth_K[-1], right=eddy_modern_earth_K[0])
    # eddy_atm_K=(m_air/m_atm)*eddy_modern_earth_K #cm^2sec^-1 #Eddy diffusion scaled as a function of pressure
    # eddy_atm_z=np.interp(eddy_modern_earth_logP,new_atm_logp,new_atm_z,left=None, right=None, period=None) #km,km,km,logP ###Pressures at which eddy diffusion is reported, mapped to altitude.
    # #Add in a high-altitude point as well to ensure we have coverage
    # eddy_N2_CO2_z_padded=np.append(eddy_N2_CO2_z_padded,np.array([100.0])) #km
    # eddy_N2_CO2_K_padded=np.append(eddy_N2_CO2_K_padded,np.array(eddy_N2_CO2_K_padded[-1])) #km
    # return eddy_atm_z, eddy_atm_K
    return new_atm_eddy

    
def make_new_ztpkzz(z_list_km, mr_co2, mr_n2, mr_h2,mr_o2, P_0_bar, T_0, T_strat, grav, filename):

    ######
    ###ZTP profile
    ######
    ###Make the ZTP profile
    n_z_atm, T_z_atm, p_z_atm=ntp_dryadiabat_isothermalstrat(z_list_km, mr_co2, mr_n2, mr_h2,mr_o2, P_0_bar, T_0, T_strat, grav)
        
    ###Output it in MEAC form. 
    towrite_tp=np.zeros([len(z_list_km),3])
    #Add BOA manually
    towrite_tp[:,0]=z_list_km
    towrite_tp[:,1]=np.log10(p_z_atm*barye2bar*bar2Pa) 
    towrite_tp[:,2]=T_z_atm
    np.savetxt('./photochem_outputs/TP_'+filename+'.dat', towrite_tp, delimiter=' ', fmt='%1.6f %1.6f %3.6f', newline='\n')
    
    ######
    ###Eddy diffusion profile
    ######
    ###Make the scaled profile
    # eddy_atm_z, eddy_atm_K=make_earth_eddy(mr_co2, mr_n2, mr_h2, mr_o2, np.flipud(z_list_km), np.flipud(towrite_tp[:,1]))
    Kz_z_atm=make_earth_eddy(mr_co2, mr_n2, mr_h2, mr_o2, z_list_km, np.log10(p_z_atm*barye2bar*bar2Pa) )

    ###Output the scaled profile
    towrite_eddy=np.zeros([len(z_list_km),2])
    towrite_eddy[:,0]=z_list_km #altitude in km
    towrite_eddy[:,1]= Kz_z_atm# #cm^2sec^-1
    np.savetxt('./photochem_outputs/Eddy_'+filename+'.dat', towrite_eddy, delimiter=' ', fmt='%1.6f %1.6f', newline='\n')
    
    ########################
    ###Import reference data (MEAC/Hu), for comparison/scaling.
    ########################
    ###Temperature-pressure
    tp_CO2_N2_z, tp_CO2_N2_logp, tp_CO2_N2_T=np.genfromtxt('./../Data/TPStd175288CO2N2.dat', skip_header=0, skip_footer=0,usecols=(0,1,2), unpack=True) #units: km, log(Pa), K
    tp_N2_z, tp_N2_logp, tp_N2_T=np.genfromtxt('./../Data/TPStd200288N2.dat', skip_header=0, skip_footer=0,usecols=(0,1,2), unpack=True) #units: km, log(Pa), K
    tp_H2_N2_z, tp_H2_N2_logp, tp_H2_N2_T=np.genfromtxt('./../Data/TPStd160288H2N2.dat', skip_header=0, skip_footer=0,usecols=(0,1,2), unpack=True) #units: km, log(Pa), K
    tp_Earth_z, tp_Earth_logp, tp_Earth_T=np.genfromtxt('./../Data/TP1986.dat', skip_header=0, skip_footer=0,usecols=(0,1,2), unpack=True) #units: km, log(Pa), K
    
    
    ###Eddy diffusion
    eddy_CO2_N2_z,eddy_CO2_N2_K=np.genfromtxt("./../Data/EddyCO2N2.dat", skip_header=0, skip_footer=0,usecols=(0,1), unpack=True)# km,cm^2sec^-1
    eddy_N2_z,eddy_N2_K=np.genfromtxt("./../Data/EddyN2.dat", skip_header=0, skip_footer=0,usecols=(0,1), unpack=True)# km,cm^2sec^-1
    eddy_H2_N2_z,eddy_H2_N2_K=np.genfromtxt("./../Data/EddyH2N2.dat", skip_header=0, skip_footer=0,usecols=(0,1), unpack=True)# km,cm^2sec^-1
    eddy_Earth_z, eddy_Earth_K=np.genfromtxt("./../Data/EddyEarth.dat", skip_header=0, skip_footer=0,usecols=(0,1), unpack=True)# km,cm^2sec^-1
    
    
    ######
    ###Plot, to make sure we are sensible/not batshit crazy.
    ######
    fig, (ax1, ax2, ax3)=plt.subplots(3, figsize=(8,10))

    ### Plot 1
    ax1.plot(T_z_atm, p_z_atm*barye2bar, color='black',linestyle='-', label='This Atm.')
    ax1.plot(tp_CO2_N2_T, 10.0**(tp_CO2_N2_logp)*Pa2bar, color='red',linestyle='-', label='Hu+2012 CO2-N2 dominated')
    ax1.plot(tp_Earth_T, 10.0**(tp_Earth_logp)*Pa2bar, color='green',linestyle='-', label='Hu+2012 Modern Earth')
    ax1.plot(tp_N2_T, 10.0**(tp_N2_logp)*Pa2bar, color='blue',linestyle='--', label='Hu+2012 N2 dominated')
    # ax1.plot(tp_H2_N2_T, 10.0**(tp_H2_N2_logp)*Pa2bar, color='pink',linestyle='-', label='Hu+2012 H2-N2 dominated')
    ax1.legend(loc='upper right', fontsize=12)
    ax1.set_ylabel('Dry Pressure (bar)')
    ax1.set_xlabel('Temperature (K)')
    ax1.set_yscale('log')
    ax1.invert_yaxis()
    ax1.set_ylim([1.0, 1.0E-9])
 
    ### Plot 2
    ax2.set_title('T(z)')
    ax2.plot(T_z_atm, z_list_km, color='black',linestyle='-', label='This Atm.')
    ax2.plot(tp_CO2_N2_T, tp_CO2_N2_z, color='red',linestyle='-', label='Hu+2012 CO2-N2 dominated')
    ax2.plot(tp_Earth_T, tp_Earth_z, color='green',linestyle='-', label='Hu+2012 Modern Earth')
    ax2.plot(tp_N2_T, tp_N2_z, color='blue',linestyle='--', label='Hu+2012 N2 dominated')
    # ax2.plot(tp_H2_N2_T, tp_H2_N2_z, color='pink',linestyle='-', label='Hu+2012 H2-N2 dominated')
    ax2.set_ylim([0.0, 100.0])

    ax2.set_ylabel('Altitude (km)')
    ax2.set_xlabel('Temperature (K)')
    
    ax3.plot(Kz_z_atm,z_list_km ,color='black',linestyle='-', label='This Atm')
    ax3.plot(eddy_CO2_N2_K,eddy_CO2_N2_z ,color='red',linestyle='-', label='Hu+2012 CO2-N2 dominated')
    ax3.plot(eddy_Earth_K,eddy_Earth_z ,color='green',linestyle='-', label='Hu+2012 Modern Earth')
    ax3.plot(eddy_N2_K,eddy_N2_z,color='blue',linestyle='--', label='Hu+2012 N2 dominated')
    # ax3.plot(eddy_H2_N2_K,eddy_H2_N2_z ,color='hotpink',linestyle='-', label='Hu+2012 H2-N2 dominated')
     
    ax3.set_ylim([0.0, 100.0])
    ax3.legend(loc='upper left', fontsize=12)
    ax3.set_ylabel('altitude(km)')
    ax3.set_xlabel('Eddy diffusion coefficient(cm^2sec^-1)')
    #ax1.invert_yaxis()
    
    plt.subplots_adjust(wspace=0., hspace=0.4)
    
    plt.savefig('./Plots/newzTPKzz_'+filename+'.pdf', orientation='portrait', format='pdf')

### make_new_ztpkzz(z_list_km_0, 0.0, 1.0, 0.0,0.0, 1.0, 288.0, 200.0, grav_0, 'test')
z_list_km_0=np.linspace(0.0, 150.0, num=151)
grav_0=G_const*(0.772*M_earth)/(0.910*R_earth)**2.0

pCO2_list=np.array([0.01, 1.0, 0.3])
pN2_list=np.array([0.99, 0.1, 1.0])
pH2_list=np.array([0.0,0.0,0.0])
pO2_list=np.array([0.0, 1.0, 0.0])

for ind in range(0, len(pCO2_list)):
    pCO2=pCO2_list[ind]
    pN2=pN2_list[ind]
    pH2=pH2_list[ind]
    pO2=pO2_list[ind]
    
    ptot=pCO2+pN2+pH2+pO2
    make_new_ztpkzz(z_list_km_0, pCO2/ptot, pN2/ptot, pH2/ptot,pO2/ptot, ptot, 288.0, 180.0, grav_0, 'T1e_pCO2={0}_pN2={1}_pH2={2}_pO2={3}_Tstrat=180K'.format(pCO2, pN2, pH2, pO2))
