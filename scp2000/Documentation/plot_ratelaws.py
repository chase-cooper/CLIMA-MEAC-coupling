#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 16:53:48 2024

@author: sukrit

Purpose of this code is to plot the updated reactions for MEAC and make sure they are not insane.
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
R = 8.314  # Universal gas constant in J/(mol·K)

########################
###Define useful constants, all in CGS (via http://www.astro.wisc.edu/~dolan/constants.html)
########################
def MM(P,T):
    """
    Returns number density at a given pressure, temperature
    Pressure in bar
    Temperature in K
    """
    P_barye=P*bar2barye #get into CGS units
    return P_barye/(k*T)

P_list=np.logspace(-10, 6, base=10, num=100) #pressure in bars


### function to compute reaction rate of the reverse reaction when the forward reaction is known
### product/reactant refer to the forward reaction
def ReverseRate(mu, n_prod, n_react, a1, a2, a3, a4, a5, a6, a7, T, k_f):
    n_total = n_prod + n_react
    delta_n = n_prod - n_react

    delta_a1 = 0.0
    delta_a2 = 0.0
    delta_a3 = 0.0
    delta_a4 = 0.0
    delta_a5 = 0.0
    delta_a6 = 0.0
    delta_a7 = 0.0

    # Loop over species
    for j in range(n_total):
        delta_a1 += mu[j] * a1[j]
        delta_a2 += mu[j] * a2[j]
        delta_a3 += mu[j] * a3[j]
        delta_a4 += mu[j] * a4[j]
        delta_a5 += mu[j] * a5[j]
        delta_a6 += mu[j] * a6[j]
        delta_a7 += mu[j] * a7[j]

    # Calculate equilibrium constant K_c
    K_c = (k*1.0e-6 * T) ** (-delta_n) * np.exp(
        delta_a1 * (np.log(T) - 1) +
        delta_a2 * T / 2 +
        delta_a3 * T**2 / 6 +
        delta_a4 * T**3 / 12 +
        delta_a5 * T**4 / 20 -
        delta_a6 / T +
        delta_a7
    )

    # Reverse rate constant
    k_r = k_f / K_c
    return k_r


# ########################
# ###Wogan updates: H + H2CO + M => CH3O + M
# #This WAS R147 H + H2CO --> CH3O + M
# #It is NOW M88 H + CH2O + M --> CH3O + M
# #Let's plot the effective bimolecular rate over a range of temperatures and pressures.
# ########################
# def r147(T):
#     """
#     T in K
#     """
#     k=4.0E-11*np.exp(-2068.7/T)
#     return k

# def m88(T,MM):
#     """
#     T in K
#     MM in cm^-3
#     Based on fit performed by Wogan (N. Wogan, personal communication, 12/18/2023) to Xu et al. 2015
#     """
#     LH=1.22E-23*pow(T,-3.0)*np.exp(-2900.0/T)
#     RH=6.56E3*pow(T,-5.0)*np.exp(-4000.0/T)
#     kkM=MM*LH/(1.0+(MM*LH/RH))
#     return kkM

# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(r147(150.0)*np.ones(np.shape(P_list))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m88(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='blue')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(r147(300.0)*np.ones(np.shape(P_list))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m88(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='blue')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(r147(600.0)*np.ones(np.shape(P_list))), linewidth=2, linestyle='-', color='red')
# ax[2].plot(P_list, np.log10(m88(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='blue')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].set_ylim([-21, -12])

########################
###Ranjan+2020 updates: CO + OH ratelaws.
##Earlier, this was represented by:
#R511 OH + CO --> CO2 + H
#M39  OH + CO + M --> CHO2 + M
##We will now represent it by the JPL 2020 expressions. 
########################
# def r511(T):
#     """
#     T in K
#     """
#     k=5.4E-14*pow(T/298.0,1.5)*np.exp(250.0/T)
#     return k

# def m39(T,MM):
#     """
#     T in K
#     MM in cm^-3
#     Based on JPL2015. 
#     """
#     RH=1.1E-12*pow(T/300.0, 1.3);
#     LH=5.9E-33*pow(T/300.0, -1.0); 
#     ind=1.0/(1.0+pow(np.log10(LH*MM/RH),2.0))
#     kkM=LH*MM/(1+LH*MM/RH)*pow(0.6,ind)
    
#     return kkM

# def oh_co_co2_h_jpl2015(T, MM):
#     """
#     Full JPL2015 expression for the CO + OH +M --> CO2 + H + M
#     """   
#     LH=1.5E-13
#     RH=2.1E9*pow(T/300.0, 6.1)
#     ind=1.0/(1.0+pow(np.log10(LH/(RH/MM)),2.0))
    
#     kkM=LH/(1.0+LH/(RH/MM))*pow(0.6,ind)
#     return kkM
    
    

# def oh_co_co2_h_jpl2020(T, MM):
#     """
#     JPL 2020 expression for CO + OH + M --> CO2 + H +M
#     """
#     K0=6.9E-33*pow((298.0/T), 2.1)
#     Kinf=1.1E-12*pow((298.0/T), -1.3)
#     ind=1.0/(1.0+pow(np.log10(K0*MM/Kinf),2.0))
#     Kf=(K0*MM/(1.0+K0*MM/Kinf))*pow(0.6, ind)

#     kint=1.85E-13*np.exp(-65.0/T)
    
#     kkM=kint*(1.0-Kf/Kinf)
#     return kkM
    
    
# def oh_co_ch2o_jpl2020(T, MM):
#     """
#     JPL 2020 expression for CO + OH + M --> CH2O + M
#     """
    
#     K0=6.9E-33*pow((298.0/T), 2.1)
#     Kinf=1.1E-12*pow((298.0/T), -1.3)
#     ind=1.0/(1.0+pow(np.log10(K0*MM/Kinf),2.0))
    
#     kkM=(K0*MM/(1.0+K0*MM/Kinf))*pow(0.6, ind)
#     return kkM

# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(r511(150.0)*np.ones(np.shape(P_list))), linewidth=2, linestyle='-', color='red', label='OH+CO-->CO2+H (R511, Hu2012)')
# ax[0].plot(P_list, np.log10(oh_co_co2_h_jpl2015(150.0, MM(P_list,150.0))), linewidth=2, linestyle='--', color='red', label='OH+CO+M-->CO2+H+M (JPL2015)')
# ax[0].plot(P_list, np.log10(oh_co_co2_h_jpl2020(150.0, MM(P_list,150.0))), linewidth=2, linestyle=':', color='red', label='OH+CO+M-->CO2+H+M (JPL2020)')
# ax[0].plot(P_list, np.log10(m39(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='blue', label='OH+CO+M-->CHO2+M (M39, JPL15)')
# ax[0].plot(P_list, np.log10(oh_co_ch2o_jpl2020(150.0, MM(P_list,150.0))), linewidth=2, linestyle=':', color='blue', label='OH+CO+M-->CO2+H+M (JPL2020)')

# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(r511(300.0)*np.ones(np.shape(P_list))), linewidth=2, linestyle='-', color='red', label='OH+CO-->CO2+H (R511, Hu2012)')
# ax[1].plot(P_list, np.log10(oh_co_co2_h_jpl2015(300.0, MM(P_list,300.0))), linewidth=2, linestyle='--', color='red', label='OH+CO+M-->CO2+H+M (JPL2015)')
# ax[1].plot(P_list, np.log10(oh_co_co2_h_jpl2020(300.0, MM(P_list,300.0))), linewidth=2, linestyle=':', color='red', label='OH+CO+M-->CO2+H+M (JPL2020)')
# ax[1].plot(P_list, np.log10(m39(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='blue', label='OH+CO+M-->CHO2+M (M39, JPL15))')
# ax[1].plot(P_list, np.log10(oh_co_ch2o_jpl2020(300.0, MM(P_list,300.0))), linewidth=2, linestyle=':', color='blue', label='OH+CO+M-->CO2+H+M (JPL2020)')

# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(r511(600.0)*np.ones(np.shape(P_list))), linewidth=2, linestyle='-', color='red', label='OH+CO-->CO2+H (R511, Hu2012)')
# ax[2].plot(P_list, np.log10(oh_co_co2_h_jpl2015(600.0, MM(P_list,600.0))), linewidth=2, linestyle='--', color='red', label='OH+CO+M-->CO2+H+M (JPL2015)')
# ax[2].plot(P_list, np.log10(oh_co_co2_h_jpl2020(600.0, MM(P_list,600.0))), linewidth=2, linestyle=':', color='red', label='OH+CO+M-->CO2+H+M (JPL2020)')
# ax[2].plot(P_list, np.log10(m39(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='blue', label='OH+CO+M-->CHO2+M (M39, JPL15))')
# ax[2].plot(P_list, np.log10(oh_co_ch2o_jpl2020(300.0, MM(P_list,300.0))), linewidth=2, linestyle=':', color='blue', label='OH+CO+M-->CO2+H+M (JPL2020)')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].set_ylim([-21, -12])
# ax[2].legend()

########################
###Update 3-body reactions in code. 
########################
# def m3_old(T,MM):
#     """
#     T in K
#     MM in cm^-3
#     Based on JPL2015. 
#     """
#     kkM=1.68E-24*pow(T/298.0,-7.0)*np.exp(-1390.0/T)*MM
    
#     return kkM

# def m3_lindemann(T,MM):
#     """
#     T in K
#     MM in cm^-3
#     Based on JPL2015. 
#     """
#     k_0_M=1.68E-24*pow(T/298.0,-7.0)*np.exp(-1390.0/T)*MM
#     k_inf=6.0E-11
#     kkM=k_0_M*k_inf/(k_0_M+k_inf)
    
#     return kkM

# def m3_baulch1994full(T,MM):
#     """
#     T in K
#     MM in cm^-3
#     Based on JPL2015. 
#     """
#     k_0_M=1.68E-24*pow(T/298.0,-7.0)*np.exp(-1390.0/T)*MM
#     k_inf=6.0E-11
#     F_c=0.38*np.exp(-T/73.0)+0.62*np.exp(-T/1180.0)
#     M_M_c=k_0_M/k_inf
#     N=0.75-1.27*np.log10(F_c)
#     F=pow(F_c, 1.0/(1.0+pow((np.log10(M_M_c)/N),2.0)))
#     kkM=(k_0_M*k_inf/(k_0_M+k_inf))*F
    
    
#     return kkM

# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m3_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red', label='CH3+CH3_M-->C2H6+M (M3, Hu2012)')
# ax[0].plot(P_list, np.log10(m3_lindemann(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple', label='CH3+CH3_M-->C2H6+M (M3, Lindemann)')
# ax[0].plot(P_list, np.log10(m3_baulch1994full(150.0, MM(P_list,150.0))), linewidth=2, linestyle='--', color='blue', label='CH3+CH3_M-->C2H6+M (M3, Baulch94)')


# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m3_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red', label='CH3+CH3_M-->C2H6+M (M3, Hu2012)')
# ax[1].plot(P_list, np.log10(m3_lindemann(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple', label='CH3+CH3_M-->C2H6+M (M3, Lindemann)')
# ax[1].plot(P_list, np.log10(m3_baulch1994full(300.0, MM(P_list,300.0))), linewidth=2, linestyle='--', color='blue', label='CH3+CH3_M-->C2H6+M (M3, Baulch94)')

# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m3_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='CH3+CH3_M-->C2H6+M (M3, Hu2012)')
# ax[2].plot(P_list, np.log10(m3_lindemann(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='CH3+CH3_M-->C2H6+M (M3, Lindemann)')
# ax[2].plot(P_list, np.log10(m3_baulch1994full(600.0, MM(P_list,600.0))), linewidth=2, linestyle='--', color='blue', label='CH3+CH3_M-->C2H6+M (M3, Baulch94)')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# # ax[2].set_ylim([-21, -12])
# ax[2].legend()


# def m5_old(T,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kkM=4.1E-30*pow(T/298.0,-2.1)*MM
    
#     return kkM

# def m5_new(T,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M=4.1E-30*pow(T/298.0,-2.1)*MM
#     k_inf=3.00E-11*pow(T/298.0, -0.90)
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
    
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m5_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red', label='M5, Hu2012')
# ax[0].plot(P_list, np.log10(m5_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple', label='M5, NIST+Lindemann')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m5_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red', label='M5, Hu2012')
# ax[1].plot(P_list, np.log10(m5_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple', label='M5, NIST+Lindemann')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m5_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='M5, Hu2012')
# ax[2].plot(P_list, np.log10(m5_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='M5, NIST+Lindemann')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# # ax[2].set_ylim([-21, -12])
# ax[2].legend()

# def m6_old(T,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kkM=7.83E-29*pow(T/298.0,-7.56)*np.exp(-5490.7/T)*MM
    
#     return kkM

# def m6_new(T,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M=7.83E-29*pow(T/298.0,-7.56)*np.exp(-5490.7/T)*MM
#     k_inf=1.9E-13*pow(T/298.0, 2.25)*np.exp(-3033.4/T)
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
    
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m6_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red', label='M6, Hu2012')
# ax[0].plot(P_list, np.log10(m6_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple', label='M6, NIST+Lindemann')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m6_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red', label='M6, Hu2012')
# ax[1].plot(P_list, np.log10(m6_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple', label='M6, NIST+Lindemann')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m6_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='M6, Hu2012')
# ax[2].plot(P_list, np.log10(m6_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='M6, NIST+Lindemann')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()


# def m7_old(T,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kkM=3.31E-30*np.exp(-740.0/T)*MM
    
#     return kkM

# def m7_new(T,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M=6.14E-28*pow(T/298.0, -4.6)*np.exp(-2237.2/T)*MM
#     k_inf=9.01E-12*np.exp(-1190.8/T)
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
    
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m7_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m7_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m7_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m7_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m7_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='M7, Hu2012')
# ax[2].plot(P_list, np.log10(m7_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='M7, NIST+Lindemann')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

# def m8_old(T,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kkM=1.38E-30*MM

#     return kkM

# def m8_new(T,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M=1.38E-30*MM
#     k_inf=3.68E-12*pow(T/298.0, 1.61)*np.exp(-1321.8/T)
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
    
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m8_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m8_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m8_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m8_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m8_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='M8, Hu2012')
# ax[2].plot(P_list, np.log10(m8_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='M8, NIST+Lindemann')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()


# def m9_old(T,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kkM=7.69E-30*np.exp(-380.0/T)*MM

#     return kkM

# def m9_hu2021(T,MM):
#     """
#     T in K
#     MM in cm^-3
#     """

#     k_0_M = np.max([1.3e-29*np.exp(-380/T), 3.7e-30])*MM
#     k_inf = 6.6e-15*np.pow(T,1.28)*np.exp(-650/T)*np.ones(len(MM))
#     Fc = 0.24*np.exp(-T/40) + 0.76*np.exp(-T/1025)
#     kkM = k_0_M / ((1.0 + k_0_M / k_inf)*Fc)

#     return kkM


# def m9_new(T,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M=7.09E-23*pow(T/298.0, -6.20)*np.exp(-2394.6/T)*MM
#     k_inf=8.51E-12*pow(T/298.0, 1.87)*np.exp(-586.4/T)
#     # k_inf=8.51E-12*pow(T/298.0, 1.87)*np.exp(-947.7/T) #NIST expression, which does not match source paper.

#     kkM=k_0_M/(1.0+k_0_M/k_inf)
    
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m9_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m9_hu2021(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='blue')
# ax[0].plot(P_list, np.log10(m9_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m9_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m9_hu2021(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='blue')
# ax[1].plot(P_list, np.log10(m9_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m9_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='M9, Hu2012')
# ax[2].plot(P_list, np.log10(m9_hu2021(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='blue', label='M9, Hu2021')
# ax[2].plot(P_list, np.log10(m9_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='M9, NIST+Lindemann')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()


# def m10_old(T,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kkM=3.21E-30*pow(T/298.0,-2.57)*np.exp(-215.0/T)*MM

#     return kkM

# def m10_new_NIST(T,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M=3.21E-30*pow(T/298.0,-2.57)*np.exp(-215.0/T)*MM
#     k_inf=7.77E-14*np.exp(2280./T)
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
    
#     return kkM

# def m10_new_RMG(T,MM):
#     """
#     T in K
#     MM in cm^-3
#     Just took Result #1 under assumption it is somehow the "best"
#     """
    
#     ###The below are from Result #1 - HFO-1234yf-seed/34 
#     alpha=0.7824
#     T3=271.0
#     T1=2755.0
#     T2=6570.0
    
#     k_inf=1.09E6*pow(T,0.48)*np.exp(1087.84/(8.314*T))*1E6/6.02E23 #m^3/mol to cm^3/molecule
#     k0=1.35E12*pow(T, -2.57)*np.exp(-5962.20/(8.314*T))*1E12/(6.02E23)**2.0 #m^6/mol^2 to cm^6/molecule^2
#     P_r=k0*MM/k_inf

#     F_cent=(1-alpha)*np.exp(-T/T3) + alpha*np.exp(-T/T1) + np.exp(-T2/T) + np.exp(-T2/T)
    
#     c=-0.4-0.67*np.log(F_cent)
#     n=0.75-1.27*np.log(F_cent)
#     d=0.14
#     F=np.exp(np.log(F_cent)*pow(1.0+pow((np.log(P_r)+c)/(n-d*(np.log(P_r)+c)),2),-1))
    
#     kkM=k_inf*(P_r/(1.0+P_r))*F
    
#     # ###The ones after are from Result #6 - Klippenstein_Glarborg2016/36
#     # k_inf=7.52621E6*pow(T,0.41)*np.exp(3963.11/(8.314*T))*1E6/6.02E23 #convert from m^3/mol to cm^3/molecule
#     # k_0=1.22463E21*pow(T, -5.09)*np.exp(-22002.80/(8.314*T))*1E12/(6.02E23)**2.0 #convert from m^6/mol^2 to cm^6/molecule^2
#     # P_r=k_0*MM/k_inf
#     # kkM=k_inf*(P_r/(1.0+P_r))
    
#     return kkM

# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m10_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m10_new_NIST(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[0].plot(P_list, np.log10(m10_new_RMG(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='blue')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m10_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m10_new_NIST(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].plot(P_list, np.log10(m10_new_RMG(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='blue')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m10_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='M10, Hu2012')
# ax[2].plot(P_list, np.log10(m10_new_NIST(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='M10, NIST+Lindemann')
# ax[2].plot(P_list, np.log10(m10_new_RMG(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='blue', label='M10, RMG')

# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

# def m11_old(T,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kkM=9.35E-30*pow(T/298.0,-2.0)*np.exp(-521.0/T)*MM

#     return kkM

# def m11_new(T,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M=8.63E-30*pow(T/298.0, -2.20)*np.exp(-567/T)*MM
#     k_inf=1.73E-10*pow(T/298,-0.50)
#     M_M_c=k_0_M/k_inf
#     F_c=0.95-(1.0E-4)*T
#     F=pow(F_c, 1.0/(1.0+pow(np.log10(M_M_c),2.0)))
    
#     kkM=k_0_M/(1.0+k_0_M/k_inf)*F
    
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m11_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m11_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m11_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m11_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m11_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='M11, Hu2012')
# ax[2].plot(P_list, np.log10(m11_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='M11, NIST+Troe')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

# def m12_old(T,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kkM=5.29E-34*np.exp(-370.0/T)*MM

#     return kkM

# def m12_new(T,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M=5.29E-34*np.exp(-370.0/T)*MM
#     k_inf=1.96E-13*np.exp(-1370.0/T)

    
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
    
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m12_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m12_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m12_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m12_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m12_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='M12, Hu2012')
# ax[2].plot(P_list, np.log10(m12_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='M12, NIST+Lindemann')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

# def m13_old(T,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kkM=6.04E-33*pow(T/298.0, -1.0)*MM

#     return kkM

# def m13_hu2021(T,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M = np.min([8.85e-33*(T/298)**(-0.6), 1.0e-32])*MM
#     k_inf = 1.0e-11
    
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
    
#     return kkM

# def m13_new(T,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M=6.04E-33*pow(T/298.0, -1.0)*MM
#     k_inf=1.81E-13*np.exp(-754.0/T)

    
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
    
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m13_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m13_hu2021(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='blue')
# ax[0].plot(P_list, np.log10(m13_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m13_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m13_hu2021(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='blue')
# ax[1].plot(P_list, np.log10(m13_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m13_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='M13, Hu2012')
# ax[2].plot(P_list, np.log10(m13_hu2021(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='blue', label='M13, Hu2021')
# ax[2].plot(P_list, np.log10(m13_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='M13, NIST+Lindemann')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()


# def m14_old(T,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kkM=3.0E-30*MM

#     return kkM

# def m14_new(T,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M=4.40E-30*pow(T/298,-1.76)*MM
#     k_inf=2.6E-10
    
#     #F_cent=0.5
#     #N=0.75-1.27*np.log(F_cent)
    
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
    
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m14_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m14_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m14_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m14_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m14_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='M13, Hu2012')
# ax[2].plot(P_list, np.log10(m14_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='M13, Altinay')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

# def m15_old(T,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     RH=2.44E-10*pow(T/298.0, -0.41)
#     LH=1.34E-31*pow(T/298.0, -1.32)*np.exp(-370.5/T)
#     ind=1.0/(1.0+pow(np.log10(LH*MM/RH),2))
#     kkM=LH*MM/(1.0+LH*MM/RH)*pow(0.6,ind)

#     return kkM

# def m15_new(T,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_inf=2.44E-10*pow(T/298.0, -0.41)
#     k_0_M=1.34E-31*pow(T/298.0, -1.32)*np.exp(-370.5/T)*MM
#     F_c=0.82
#     F=pow(F_c,1.0/(1.0+pow(np.log(k_0_M/k_inf),2)))
    
#     #F_cent=0.5
#     #N=0.75-1.27*np.log(F_cent)
    
#     kkM=k_0_M/(1.0+k_0_M/k_inf)*F
    
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m15_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m15_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m15_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m15_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m15_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='M15, Hu2012')
# ax[2].plot(P_list, np.log10(m15_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='M15, Altinay')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

# def m16_old(T,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
    
#     kkM=4.36E-32*pow(T/298.0,-1.0)*MM

#     return kkM

# def m16_new(T,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M=4.36E-32*pow(T/298.0,-1.0)*MM
#     k_inf=1.0E-11
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m16_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m16_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m16_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m16_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m16_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='m16, Hu2012')
# ax[2].plot(P_list, np.log10(m16_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='M17, Altinay')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()


# def m17_old(T,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     RH=7.5E-11*pow(T/300.0, 0.2);
#     LH=4.4E-32*pow(T/300.0, -1.3);
#     ind=1.0/(1.0+pow(np.log10(LH*MM/RH),2.0));
#     kkM=LH*MM/(1.0+LH*MM/RH)*pow(0.6,ind)

#     return kkM

# def m17_new(T,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_inf=9.5E-11*pow(T/298.0, 0.4)
#     k_0_M=5.3E-32*pow(T/298.0, -1.8)*MM
#     ind=1.0/(1.0+pow(np.log10(k_0_M/k_inf),2.0))    
#     kkM=(k_0_M/(1+k_0_M/k_inf))*pow(0.6, ind)
    
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m17_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m17_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m17_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m17_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m17_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='M17, Hu2012')
# ax[2].plot(P_list, np.log10(m17_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='M17, Altinay')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()


# def m18_old(T,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
    
#     kkM=6.87E-31*pow(T/298.0, -2.0)*MM

#     return kkM

# def m18_new(T,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M=6.87E-31*pow(T/298.0, -2.0)*MM
#     k_inf=4.17E-11*pow(T,0.234)*np.exp(57.5/T) #from Sellevag et al. 2008. Note expression little bit different than NIST extraction.
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m18_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m18_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m18_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m18_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m18_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='m18, Hu2012')
# ax[2].plot(P_list, np.log10(m18_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='M17, Altinay')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()


# def m19_old(T,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
    
#     kkM=1.7E-33*np.exp(1000.0/T)*MM

#     return kkM

# def m19_new(T,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M=2.1E-33*np.exp(920.0/T)*MM
#     k_inf=7.0E-12*np.exp(620.0/T) #This is almost certainly wrong -- it is in the Zahnle network but based on chatting with Nick it is probably an error (Two of them, (1) typo in pulling from Atkinson+2005 and (2) this is probably an entirely separate reaction. DEFINITELY needs to be fixed.)
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m19_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m19_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m19_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m19_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m19_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='m19, Hu2012')
# ax[2].plot(P_list, np.log10(m19_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='M17, Altinay')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

# def m20_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
    
#     RH=2.9E-12*pow(tl/300.0, -1.1)
#     LH=2.0E-31*pow(tl/300.0, -3.4)
#     ind=1.0/(1.0+pow(np.log10(LH*MM/RH),2.0))
#     kkM=LH*MM/(1.0+LH*MM/RH)*pow(0.6,ind)

#     return kkM

# def m20_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     RH=4.0E-12*pow(tl/298.0, -0.3)
#     LH=1.9E-31*pow(tl/298.0, -3.4)
#     ind=1.0/(1.0+pow(np.log10(LH*MM/RH),2.0))
#     kkM=LH*MM/(1.0+LH*MM/RH)*pow(0.6,ind)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m20_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m20_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m20_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m20_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m20_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='m20, Hu2012')
# ax[2].plot(P_list, np.log10(m20_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='M20, New')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

# def m21_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kkM=9.18E-34*pow(tl/298.0, -1.69)*MM

#     return kkM

# def m21_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M=9.18E-34*pow(tl/298.0, -1.69)*MM
#     k_inf=2.0E-10*pow(tl/298.0,0.31)
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m21_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m21_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m21_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m21_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m21_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='m21, Hu2012')
# ax[2].plot(P_list, np.log10(m21_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='m21, New')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

# def m22_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kkM=9.4E-33*MM

#     return kkM

# def m22_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M=9.4E-33*MM
#     # k_inf=8.0E-19 #This is from Antipov+2009, which is NOT the high-pressure limit. 
#     k_inf=1.0E-10 #This is PURELY AD-HOC VALUE TO PREVENT INSANITY IN HIGH-PRESSURE LIMIT
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m22_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m22_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m22_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m22_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m22_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='m22, Hu2012')
# ax[2].plot(P_list, np.log10(m22_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='m22, New')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

# def m23_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kkM=5.0E-32*MM

#     return kkM

# def m23_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M=5.0E-32*MM
#     k_inf=1.0E-10 #This is PURELY AD-HOC VALUE TO PREVENT INSANITY IN HIGH-PRESSURE LIMIT
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m23_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m23_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m23_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m23_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m23_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='m23, Hu2012')
# ax[2].plot(P_list, np.log10(m23_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='m23, New')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()


# def m24_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kkM=1.0E-36*MM

#     return kkM

# def m24_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M=8.3E-38*MM
#     k_inf=1.94E-20
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m24_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m24_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m24_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m24_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m24_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='m24, Hu2012')
# ax[2].plot(P_list, np.log10(m24_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='m24, New')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

# def m25_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kkM=1.38E-33*np.exp(502.7/tl)*MM

#     return kkM

# def m25_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M=1.38E-33*np.exp(502.7/tl)*MM
#     k_inf=5.0E-16
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m25_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m25_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m25_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m25_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m25_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='m25, Hu2012')
# ax[2].plot(P_list, np.log10(m25_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='m25, New')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()


# def m26_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kkM=5.46E-33*pow(tl/298.0,-1.29)*MM

#     return kkM

# def m26_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M=5.46E-33*np.exp(155.2/tl)*MM
#     k_inf=1.0E-10
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m26_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m26_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m26_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m26_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m26_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='m26, Hu2012')
# ax[2].plot(P_list, np.log10(m26_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='m26, New')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

# def m29_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     RH=1.4E-12*pow(tl/300.0, -0.7)
#     LH=2.0E-30*pow(tl/300.0, -4.4)
#     ind=1/(1.0+pow(np.log10(LH*MM/RH),2.0));
#     kkM=LH*MM/(1.0+LH*MM/RH)*pow(0.6,ind);
#     return kkM

# def m29_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     RH=1.6E-12*pow(tl/298.0, 0.1)
#     LH=2.4E-30*pow(tl/298.0, -3.0)
#     ind=1/(1.0+pow(np.log10(LH*MM/RH),2.0));
#     kkM=LH*MM/(1.0+LH*MM/RH)*pow(0.6,ind);
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m29_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m29_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m29_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m29_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m29_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='m29, Hu2012')
# ax[2].plot(P_list, np.log10(m29_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='m29, New')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()


# def m30_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kkM=2.0E-34*MM
#     return kkM

# def m30_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M=2.0E-34*MM
#     k_inf=4.82E-15*pow((tl/300.0),-1)
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m30_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m30_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m30_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m30_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m30_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='m30, Hu2012')
# ax[2].plot(P_list, np.log10(m30_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='m30, New')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()


# def m31_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kkM=1.7E-33*np.exp(-1509.0/tl)*MM
#     return kkM

# def m31_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M=1.7E-33*np.exp(-1509.0/tl)*MM
#     k_inf=1.0E-14*np.exp(-1630.0/tl)
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m31_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m31_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m31_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m31_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m31_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='m31, Hu2012')
# ax[2].plot(P_list, np.log10(m31_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='m31, New')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

# ########################
# ###M33: Association reaction O + NO2 (Burkholder+2020)
# ##Earlier, this was represented by:
# #R442 O + NO2 --> NO + O2
# #M33  O + NO2 + M --> NO3 + M
# ##We will now represent it by the JPL 2020 expressions:
# #M33 O + NO2 + M--> NO3 + M
# #M90 O + NO2 + M--> NO + O2 + M
# #R442: zeroed, excised. 
# ########################
# def r442(tl):
#     """
#     tl in K
#     """
#     k=5.1E-12*np.exp(210.0/tl)
#     return k

# def m33_old(tl,MM):
#     """
#     tl in K
#     MM in cm^-3
#     """    
#     RH=2.2E-11*pow(tl/300.0, -0.7)
#     LH=2.5E-31*pow(tl/300.0, -1.8)
#     ind=1.0/(1.0+pow(np.log10(LH*MM/RH),2.0))
#     kkM=LH*MM/(1.0+LH*MM/RH)*pow(0.6,ind) 
#     return kkM

# def m33_new(tl, MM):
#     """
#     JPL 2020 expression for O + NO2 + M --> NO3 + M
#     """
    
#     K0=3.4E-31*pow((298.0/tl), 1.6)
#     Kinf=2.3E-11*pow((298.0/tl), 0.2)
#     ind=1.0/(1.0+pow(np.log10(K0*MM/Kinf),2))
    
#     kkM=(K0*MM/(1.0+K0*MM/Kinf))*pow(0.6, ind)
#     return kkM

# def m90_new(tl, MM):
#     """
#     JPL 2020 expression for O + NO2 + M --> O2 + NO +M
#     """
#     K0=3.4E-31*pow((298.0/tl), 1.6)
#     Kinf=2.3E-11*pow((298.0/tl), 0.2)
#     ind=1.0/(1.0+pow(np.log10(K0*MM/Kinf),2))
#     Kf=(K0*MM/(1.0+K0*MM/Kinf))*pow(0.6, ind)

#     kint=5.3E-12*np.exp(200.0/tl)
    
#     kkM=kint*(1.0-Kf/Kinf)
#     return kkM
    
    


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(r442(150.0)*np.ones(np.shape(P_list))), linewidth=2, linestyle='-', color='red', label='O+NO2-->O2+NO (R442, Hu2012)')
# ax[0].plot(P_list, np.log10(m90_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle=':', color='red', label='O+NO2+M-->O2+NO+M (JPL2020, M90_new)')
# ax[0].plot(P_list, np.log10(m33_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='blue', label='O+NO2+M-->NO3+M (M33, Hu2012)')
# ax[0].plot(P_list, np.log10(m33_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle=':', color='blue', label='O+NO2+M-->NO3+M (JPL2020, M33_new)')

# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(r442(300.0)*np.ones(np.shape(P_list))), linewidth=2, linestyle='-', color='red', label='O+NO2-->O2+NO (R442, Hu2012)')
# ax[1].plot(P_list, np.log10(m90_new(300.0, MM(P_list,150.0))), linewidth=2, linestyle=':', color='red', label='O+NO2+M-->O2+NO+M (JPL2020, M90_new)')
# ax[1].plot(P_list, np.log10(m33_old(300.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='blue', label='O+NO2+M-->NO3+M (M33, Hu2012)')
# ax[1].plot(P_list, np.log10(m33_new(300.0, MM(P_list,150.0))), linewidth=2, linestyle=':', color='blue', label='O+NO2+M-->NO3+M (JPL2020, M33_new)')

# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(r442(600.0)*np.ones(np.shape(P_list))), linewidth=2, linestyle='-', color='red', label='O+NO2-->O2+NO (R442, Hu2012)')
# ax[2].plot(P_list, np.log10(m90_new(600.0, MM(P_list,150.0))), linewidth=2, linestyle=':', color='red', label='O+NO2+M-->O2+NO+M (JPL2020, M90_new)')
# ax[2].plot(P_list, np.log10(m33_old(600.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='blue', label='O+NO2+M-->NO3+M (M33, Hu2012)')
# ax[2].plot(P_list, np.log10(m33_new(600.0, MM(P_list,150.0))), linewidth=2, linestyle=':', color='blue', label='O+NO2+M-->NO3+M (JPL2020, M33_new)')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].set_ylim([-21, -10])
# ax[2].legend()
# plt.show()


# def m34_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kkM=5.21E-35*np.exp(900.0/tl)*MM
#     return kkM

# def m34_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M=5.21E-35*np.exp(900.0/tl)*MM
#     k_inf=1.21E-11*pow(tl/300,-2.0)#1.0E-12
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m34_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m34_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m34_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m34_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m34_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='m34, Hu2012')
# ax[2].plot(P_list, np.log10(m34_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='m34, New')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()


# plt.show()


# def m35_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kkM=6.0E-34*pow(tl/300.0, -2.4)*MM
#     return kkM

# def m35_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M=6.0E-34*pow(tl/298.0, -2.4)*MM
#     k_inf=2.8E-11
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m35_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m35_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m35_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m35_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m35_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='m35, Hu2012')
# ax[2].plot(P_list, np.log10(m35_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='m35, New')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

# def m36_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kkM=2.8E-36*pow(tl/300.0, -0.9)*MM
#     return kkM

# def m36_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M=2.8E-36*pow(tl/298.0, -0.9)*MM
#     k_inf=3.4e-16
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m36_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m36_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m36_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m36_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m36_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='m36, Hu2012')
# ax[2].plot(P_list, np.log10(m36_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='m36, New')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()


# def m38_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     RH=7.5E-12*pow(tl/300.0, -0.85)
#     LH=1.0E-28*pow(tl/300.0, -4.5)
#     ind=1.0/(1.0+pow(np.log10(LH*MM/RH),2))
#     kkM=LH*MM/(1.0+LH*MM/RH)*pow(0.6,ind)
#     return kkM

# def m38_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     RH=8.5E-12*pow(tl/298.0, -1.75)
#     LH=1.1E-28*pow(tl/298.0, -3.5)
#     ind=1.0/(1.0+pow(np.log10(LH*MM/RH),2))
#     kkM=LH*MM/(1.0+LH*MM/RH)*pow(0.6,ind)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m38_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m38_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m38_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m38_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m38_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='m38, Hu2012')
# ax[2].plot(P_list, np.log10(m38_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='m38, New')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

# def m40_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kkM=1.8E-27*pow(tl/298.0,-3.85)*MM
#     return kkM

# def m40_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M=2.0E-23*pow(tl/1.0, -1.3)*np.exp(-362.0/tl)*MM
#     k_inf=1.50E-10
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m40_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m40_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m40_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m40_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m40_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='m40, Hu2012')
# ax[2].plot(P_list, np.log10(m40_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='m40, New')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()


# def m42_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     RH=3.6E-11*pow(tl/300.0, -0.1)
#     LH=7.0E-31*pow(tl/300.0, -2.6)
#     ind=1.0/(1.0+pow(np.log10(LH*MM/RH),2.0))
#     kkM=LH*MM/(1.0+LH*MM/RH)*pow(0.6,ind)
#     return kkM

# def m42_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     RH=3.6E-11*pow(tl/298.0, -0.1)
#     LH=7.1E-31*pow(tl/298.0, -2.6)
#     ind=1.0/(1.0+pow(np.log10(LH*MM/RH),2.0))
#     kkM=LH*MM/(1.0+LH*MM/RH)*pow(0.6,ind)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m42_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m42_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m42_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m42_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m42_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='m42, Hu2012')
# ax[2].plot(P_list, np.log10(m42_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='m42, New')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

# def m43_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     RH=2.8E-11
#     LH=1.8E-30*pow(tl/300.0, -3.0)
#     ind=1.0/(1.0+pow(np.log10(LH*MM/RH),2.0))
#     kkM=LH*MM/(1.0+LH*MM/RH)*pow(0.6,ind)
#     return kkM

# def m43_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     RH=2.8E-11
#     LH=1.8E-30*pow(tl/298.0, -3.0)
#     ind=1.0/(1.0+pow(np.log10(LH*MM/RH),2.0))
#     kkM=LH*MM/(1.0+LH*MM/RH)*pow(0.6,ind)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m43_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m43_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m43_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m43_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m43_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='m43, Hu2012')
# ax[2].plot(P_list, np.log10(m43_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='m43, New')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()


# def m45_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kkM=5.0E-32*np.exp(900.0/tl)*MM
#     return kkM

# def m45_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M=1.98E-33*np.exp(206.0/tl)*MM
#     k_inf=2.26E-14*np.exp(415.0/tl)
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m45_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m45_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m45_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m45_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m45_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='m45, Hu2012')
# ax[2].plot(P_list, np.log10(m45_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='m45, New')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

# def m46_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kkM=4.82E-31*pow(tl/298.0,-2.17)*MM
#     return kkM

# def m46_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M=4.82E-31*pow(tl/298.0,-2.17)*MM
#     k_inf=5.30E-11
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m46_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m46_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m46_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m46_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m46_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='m46, Hu2012')
# ax[2].plot(P_list, np.log10(m46_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='m46, New')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()


# def m47_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kkM=6.45E-29*pow(tl/298.0,-3.48)*np.exp(-490.0/tl)*MM
#     return kkM

# def m47_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M=6.45E-29*pow(tl/298.0,-3.48)*np.exp(-490.0/tl)*MM
#     k_inf=4.47E-11*pow(tl/298.0,0.50)*np.exp(199.7/tl)
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m47_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m47_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m47_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m47_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m47_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='m47, Hu2012')
# ax[2].plot(P_list, np.log10(m47_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='m47, New')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()


# def m48_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kkM=5.61E-30*pow(tl/298.0,-5.19)*np.exp(-2271.0/tl)*MM
#     return kkM

# def m48_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M=5.61E-30*pow(tl/298.0,-5.19)*np.exp(-2271.0/tl)*MM
#     k_inf=7.58E-12*pow(tl/298.0,1.59)*np.exp(-1243.6/tl)
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM

# def m48_alt(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M=5.74E-31*pow(tl/298.0,-3.69)*np.exp(-2410.3/tl)*MM
#     k_inf=2.31E-11*pow(tl/298.0,0.62)*np.exp(-1819.7/tl)
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM



# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m48_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m48_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[0].plot(P_list, np.log10(m48_alt(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='blue')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m48_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m48_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].plot(P_list, np.log10(m48_alt(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='blue')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m48_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='m48, Hu2012')
# ax[2].plot(P_list, np.log10(m48_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='m48, New')
# ax[2].plot(P_list, np.log10(m48_alt(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='blue', label='m48, alt')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

# def m49_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     RH=4.2E-14*pow(tl/300.0, 1.8)
#     LH=1.8E-33*pow(tl/300.0, 2.0)
#     ind=1.0/(1.0+pow(np.log10(LH*MM/RH),2.0))
#     kkM=LH*MM/(1.0+LH*MM/RH)*pow(0.6,ind)
#     return kkM

# def m49_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     RH=4.1E-14*pow(tl/298.0, 1.8)
#     LH=1.8E-33*pow(tl/298.0, 2.0)
#     ind=1.0/(1.0+pow(np.log10(LH*MM/RH),2.0))
#     kkM=LH*MM/(1.0+LH*MM/RH)*pow(0.6,ind)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m49_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m49_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m49_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m49_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m49_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='m49, Hu2012')
# ax[2].plot(P_list, np.log10(m49_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='m49, New')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

# def m50_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     RH=1.6E-12
#     LH=3.3E-31*pow(tl/300.0, -4.3)
#     ind=1.0/(1.0+pow(np.log10(LH*MM/RH),2.0))
#     kkM=LH*MM/(1.0+LH*MM/RH)*pow(0.6,ind)
#     return kkM

# def m50_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     RH=1.7E-12*pow(tl/298.0,0.2)
#     LH=2.9E-31*pow(tl/298.0, -4.1)
#     ind=1.0/(1.0+pow(np.log10(LH*MM/RH),2.0))
#     kkM=LH*MM/(1.0+LH*MM/RH)*pow(0.6,ind)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m50_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m50_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m50_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m50_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m50_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='m50, Hu2012')
# ax[2].plot(P_list, np.log10(m50_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='m50, New')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()


# def m51_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     RH=7.5E-11
#     LH=5.7E-32*pow(tl/300.0, -1.6)
#     ind=1.0/(1.0+pow(np.log10(LH*MM/RH),2.0))
#     kkM=LH*MM/(1.0+LH*MM/RH)*pow(0.6,ind)
#     return kkM

# # def m51_new(tl,MM): ##NOT ADOPTED
# #     """
# #     T in K
# #     MM in cm^-3
# #     """
# #     #This is from Wogan+2023
# #     k_0_M=7.431425e-34*pow(tl,0.28)*MM
# #     k_inf=1.0E-11
# #     kkM=k_0_M/(1.0+k_0_M/k_inf)
# #     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m51_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m51_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m51_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m51_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m51_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='m51, Hu2012')
# ax[2].plot(P_list, np.log10(m51_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='m51, New')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()


# def m52_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kkM=1.1E-30*pow(tl/300.0, -2.0)*MM
#     return kkM

# def m52_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M=1.1E-30*pow(tl/300.0, -2.0)*MM
#     k_inf=8.0E-11
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m52_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m52_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m52_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m52_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m52_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='m52, Hu2012')
# ax[2].plot(P_list, np.log10(m52_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='m52, New')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

# def m53_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kkM=1.1E-30*pow(tl/300.0, -2.0)*MM
#     return kkM

# def m53_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M=1.1E-30*pow(tl/300.0, -2.0)*MM
#     k_inf=8.0E-11
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m53_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m53_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m53_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m53_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m53_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='m53, Hu2012')
# ax[2].plot(P_list, np.log10(m53_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='m53, New')
# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

# def m54_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kkM=4.0E-31*np.exp(900.0/tl)*MM
#     return kkM

# def m54_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M=4.0E-31*np.exp(900.0/tl)*MM
#     k_inf=3.0E-11
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM

# def m54_alt(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     #This is the Wogan+2023 prescription
#     k_0_M=9.0E-26*pow(tl, -2.0)*MM
#     k_inf=3.0E-11
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m54_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m54_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[0].plot(P_list, np.log10(m54_alt(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='blue')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m54_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m54_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].plot(P_list, np.log10(m54_alt(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='blue')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m54_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='m54, Hu2012')
# ax[2].plot(P_list, np.log10(m54_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='m54, New')
# ax[2].plot(P_list, np.log10(m54_alt(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='blue', label='m54, Wogan+2023')

# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()


# def m55_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kkM=4.0E-31*np.exp(900.0/tl)*MM
#     return kkM

# def m55_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M=4.0E-31*np.exp(900.0/tl)*MM
#     k_inf=1.0E-11
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM

# def m55_alt(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     #This is the Wogan+2023 prescription
#     k_0_M=9.0E-26*pow(tl, -2.0)*MM
#     k_inf=1.0E-11
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m55_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m55_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[0].plot(P_list, np.log10(m55_alt(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='blue')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m55_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m55_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].plot(P_list, np.log10(m55_alt(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='blue')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m55_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='m55, Hu2012')
# ax[2].plot(P_list, np.log10(m55_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='m55, New')
# ax[2].plot(P_list, np.log10(m55_alt(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='blue', label='m55, Wogan+2023')

# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

# def m56_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     RH=2.37E-12*np.exp(523.0/tl)
#     LH=5.8E-30*np.exp(355.0/tl)
#     ind=1.0/(1.0+pow(np.log10(LH*MM/RH),2.0))
#     kkM=LH*MM/(1.0+LH*MM/RH)*pow(0.6,ind)
#     return kkM

# def m56_alt(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     # Wogan+2023
#     k_0_M=1.8E-24*pow(tl,-2.0)*MM
#     k_inf=1.8E-5*pow(tl, -2.0)
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m56_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m56_alt(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m56_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m56_alt(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m56_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='m56, Hu2012')
# ax[2].plot(P_list, np.log10(m56_alt(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='m56, alt')

# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()


# def m57_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kkM=6.2E-29*pow(tl/298.0,-1.8)*MM
#     return kkM

# def m57_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M=6.2E-29*pow(tl/298.0,-1.8)*MM
#     k_inf=3.5E-10
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM

# def m57_alt(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     #This is the Wogan+2023 prescription
#     k_0_M=1.55214E-23*pow(tl, -2.17)*MM
#     k_inf=3.5E-10
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM

# def m57_hu2021(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M = 6.0e-29*np.max([tl/298, 1.0])**-1.8*MM
#     k_inf = 1.92e-8*np.max([tl,110])**-0.5*np.exp(-400/np.max([tl,110]))
#     Fc = 0.3 + 0.58*np.exp(-tl/800)
#     kkM = k_0_M / ((1.0 + k_0_M / k_inf)*Fc)

#     return kkM

# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m57_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m57_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[0].plot(P_list, np.log10(m57_alt(150.0, MM(P_list,150.0))), linewidth=2, linestyle='--', color='blue')
# ax[0].plot(P_list, np.log10(m57_hu2021(150.0, MM(P_list,150.0))), linewidth=2, linestyle='--', color='green')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m57_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m57_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].plot(P_list, np.log10(m57_alt(300.0, MM(P_list,300.0))), linewidth=2, linestyle='--', color='blue')
# ax[1].plot(P_list, np.log10(m57_hu2021(300.0, MM(P_list,300.0))), linewidth=2, linestyle='--', color='green')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m57_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='m57, Hu2012')
# ax[2].plot(P_list, np.log10(m57_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='m57, New')
# ax[2].plot(P_list, np.log10(m57_alt(600.0, MM(P_list,600.0))), linewidth=2, linestyle='--', color='blue', label='m57, Wogan+2023')
# ax[2].plot(P_list, np.log10(m57_hu2021(600.0, MM(P_list,600.0))), linewidth=2, linestyle='--', color='green', label='m57, Hu+2021')

# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

# def m58_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     RH=1.5E-10
#     LH=5.5E-23*pow(tl,-2.0)*np.exp(-1040/tl)
#     ind=1.0/(1.0+pow(np.log10(LH*MM/RH),2.0))
#     kkM=LH*MM/(1.0+LH*MM/RH)*pow(0.6,ind)
#     return kkM
      
# def m58_alt(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     # Rimmer+2019
#     #H       C2H5            C2H6                                 0   2.00e-28-1.50e+00 0.00e-00 B  XX         *XXXX*
#     #H       C2H5            C2H6                                 0   1.70e-10 0.00e-00 0.00e-00 B  XX         *XXXX*
#     #Key from Rimmer+2016 equation 9-10
#     k_0_M=2.00E-28*pow(tl,-1.5)*MM
#     k_inf=1.7E-10
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m58_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m58_alt(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m58_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m58_alt(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m58_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='m58, Hu2012')
# ax[2].plot(P_list, np.log10(m58_alt(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='m58, alt')

# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()


# def m73_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     RH=2.37E-12*np.exp(523.0/tl)
#     LH=5.8E-30*np.exp(355.0/tl)
#     ind=1.0/(1.0+pow(np.log10(LH*MM/RH),2.0))
#     kkM=LH*MM/(1.0+LH*MM/RH)*pow(0.6,ind)
#     return kkM

# def m73_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     #This is from Moses+2011
#     k_0_M=9.0E-31*np.exp(550.0/tl)*MM
#     k_inf=8.55E-11*pow(tl,0.15)
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM

# def m73_alt(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     #This is the Wogan+2024 prescription
#     k_0_M=3.387566E-25*pow(tl, -2.3)*MM
#     k_inf=8.500839E-11*pow(tl, 0.15)
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m73_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m73_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[0].plot(P_list, np.log10(m73_alt(150.0, MM(P_list,150.0))), linewidth=2, linestyle='--', color='blue')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m73_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m73_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].plot(P_list, np.log10(m73_alt(300.0, MM(P_list,300.0))), linewidth=2, linestyle='--', color='blue')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m73_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='m73, Hu2012')
# ax[2].plot(P_list, np.log10(m73_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='m73, New')
# ax[2].plot(P_list, np.log10(m73_alt(600.0, MM(P_list,600.0))), linewidth=2, linestyle='--', color='blue', label='m73, Wogan+2023')

# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

# def m85_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kkM=1.96E-29*pow(tl/298.0,-3.9)*MM
#     return kkM

# def m85_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     #This is from Moses+2011
#     k_0_M=4.48E-14*pow(tl,-5.49)*np.exp(-1000.0/tl)*MM
#     k_inf=9.33E-10*pow(tl,-0.414)*np.exp(-33/tl)
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM

# def m85_alt(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     #This is the Wogan+2024 prescription
#     k_0_M=9.158047E-20*pow(tl, -3.9)*MM
#     k_inf=2.167928E-11*pow(tl, 0.3)
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM


# def m85_dai(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     #This is from Dai+2020
#     k_0_M=2.33E-19*pow(tl/298.0, -14.40)*np.exp(-6914.4/tl)*MM
#     k_inf=1.04E-8*pow(tl/298.0,-4.17)*np.exp(-1555.1/tl)
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM

# def m85_fit(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     #This is from NIST fit to Dai+2020 and Klippenstein+2009
#     k_0_M=1.34E-23*pow(tl/298.0, -9.8)*np.exp(-4025.5/tl)*MM
#     k_inf=2.17E-10*pow(tl/298.0,-1.4)*np.exp(-300.7/tl)
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM

# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m85_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m85_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[0].plot(P_list, np.log10(m85_alt(150.0, MM(P_list,150.0))), linewidth=2, linestyle='--', color='blue')
# ax[0].plot(P_list, np.log10(m85_dai(150.0, MM(P_list,150.0))), linewidth=2, linestyle='--', color='green')
# ax[0].plot(P_list, np.log10(m85_fit(150.0, MM(P_list,150.0))), linewidth=2, linestyle='--', color='hotpink')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m85_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m85_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].plot(P_list, np.log10(m85_alt(300.0, MM(P_list,300.0))), linewidth=2, linestyle='--', color='blue')
# ax[1].plot(P_list, np.log10(m85_dai(300.0, MM(P_list,300.0))), linewidth=2, linestyle='--', color='green')
# ax[1].plot(P_list, np.log10(m85_fit(300.0, MM(P_list,300.0))), linewidth=2, linestyle='--', color='hotpink')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m85_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='m85, Hu2012')
# ax[2].plot(P_list, np.log10(m85_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='m85, New')
# ax[2].plot(P_list, np.log10(m85_alt(600.0, MM(P_list,600.0))), linewidth=2, linestyle='--', color='blue', label='m85, Wogan+2023')
# ax[2].plot(P_list, np.log10(m85_dai(600.0, MM(P_list,600.0))), linewidth=2, linestyle='--', color='green', label='Dai+2020')
# ax[2].plot(P_list, np.log10(m85_fit(600.0, MM(P_list,600.0))), linewidth=2, linestyle='--', color='hotpink', label='NIST Fit')

# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()


# def m86_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kkM=1.8E-27*pow(tl/298.0,-3.85)*MM
#     return kkM

# def m86_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     #Underlying reference is Jodkowski+1995.
#     k_0_M=1.8E-27*pow(tl/298.0,-3.85)*MM
#     k_inf=1.30E-10*pow(tl/298.0,0.42)
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM

# def m86_alt(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     #This is from Rimmer+2019
#     # CH3     H2N             CH5N                                 0   5.25e-11 0.00e-00 1.77e+04 G  XX         *XXXX*
#     # CH3     H2N             CH5N                                 0   6.92e+10 0.00e-00 2.42e+04 G  XX         *XXXX*

#     k_0_M=5.25E-11*np.exp(-1.77E4/tl)*MM
#     k_inf=6.92E10*np.exp(-2.42E4/tl)
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m86_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m86_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[0].plot(P_list, np.log10(m86_alt(150.0, MM(P_list,150.0))), linewidth=2, linestyle='--', color='blue')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m86_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m86_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].plot(P_list, np.log10(m86_alt(300.0, MM(P_list,300.0))), linewidth=2, linestyle='--', color='blue')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m86_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='m86, Hu2012')
# ax[2].plot(P_list, np.log10(m86_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='m86, New')
# ax[2].plot(P_list, np.log10(m86_alt(600.0, MM(P_list,600.0))), linewidth=2, linestyle='--', color='blue', label='m86, Rimmer+2021')

# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()


# def m87_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kkM=2.2E-33*np.exp(-1780.0/tl)*MM
#     return kkM

# def m87_alt(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     # Zahnle+2016
#     k_0_M=3.6E-34*pow(tl/298.0,-0.57)*MM
#     k_inf=3.0E-14
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(m87_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(m87_alt(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(m87_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(m87_alt(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(m87_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='m87, Hu2012')
# ax[2].plot(P_list, np.log10(m87_alt(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='m87, alt')

# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()


########################
###Update 2-body reactions in code. 
########################

# def r415_old(tl):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kk=np.zeros(np.shape(tl))+1.4E-10
#     return kk

# def r415_alt(tl):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kk=np.zeros(np.shape(tl))+9.0E-11
#     return kk


# fig, ax=plt.subplots(1, figsize=(8.,6.), sharex=True, sharey=True)
# markersizeval=5.
# T_list=np.linspace(200, 1000, num=1000)

# ax.set_title('R415')
# ax.plot(T_list, np.log10(r415_old(T_list)), linewidth=2, linestyle='-', color='red', label='R415, Hu2012')
# ax.plot(T_list, np.log10(r415_alt(T_list)), linewidth=2, linestyle='-', color='blue', label='R415, Wogan+2024')
        
# ax.set_ylabel('log10(k)')
# ax.set_xscale('linear')
# ax.set_xlabel('Temperature (K)')
# ax.set_yscale('linear')
# ax.legend()

# def r148_old(tl):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kk=2.14E-12*pow(tl/298.0, 1.62)*np.exp(-1090.0/tl)
#     return kk

# def r148_alt(tl):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kk=2.28E-19*pow(tl, 2.65)*np.exp(-766.5/tl)
#     return kk


# fig, ax=plt.subplots(1, figsize=(8.,6.), sharex=True, sharey=True)
# markersizeval=5.
# T_list=np.linspace(200, 1000, num=1000)

# ax.set_title('r148')
# ax.plot(T_list, np.log10(r148_old(T_list)), linewidth=2, linestyle='-', color='red', label='r148, Hu2012')
# ax.plot(T_list, np.log10(r148_alt(T_list)), linewidth=2, linestyle='-', color='blue', label='r148, Wogan+2024')
        
# ax.set_ylabel('log10(k)')
# ax.set_xscale('linear')
# ax.set_xlabel('Temperature (K)')
# ax.set_yscale('linear')
# ax.legend()


########################
###Update unimolecular reactions in code.
########################

# def t1_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kkM=7.16E-10*np.exp(-11200.0/tl)*MM
#     return kkM

# def t1_alt(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M=7.16E-10*np.exp(-11200.0/tl)*MM #Heimerl & Coffee; valid from 300-3000 K. 
#     k_inf=7.60E12*np.exp(-12268.0/tl) #Nominally Popovich et al. 1985, but that is in Russian so really just trusting NIST here. 
    
#     # # Peukert+2013
#     # k_0_M=5.82E-7*pow(tl/298.0,-4.37)*np.exp(-13711.0/tl)*MM
#     # k_inf=3.01E13*pow(tl/298.0, -0.67)*np.exp(-13110.0/tl)
    
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(t1_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(t1_alt(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(t1_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(t1_alt(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(t1_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='t1, Hu2012')
# ax[2].plot(P_list, np.log10(t1_alt(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='t1, alt')

# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()


# def t2_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kkM=2.41E-8*pow(tl/298.0, -1.18)*np.exp(-24415.0/tl)*MM
#     return kkM

# def t2_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M=2.41E-8*pow(tl/298.0, -1.18)*np.exp(-24415.0/tl)*MM #From Tsang+1987.
#     k_inf=5.82E11*pow(tl/300.0, -2.18)*np.exp(-24400.0/tl) #Adapted from Rimmer+2016. No justification for this value provided. Rimmer+2016 use a barrier of 24400, not 24415; insignificant, but adopted for consistency. 
    
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(t2_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(t2_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(t2_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(t2_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(t2_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='t2, Hu2012')
# ax[2].plot(P_list, np.log10(t2_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='t2, new')

# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()


# def t3_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kkM=2.01E-7*np.exp(-22852.0/tl)*MM
#     return kkM

# def t3_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M=8.43E-6*pow(tl/298.0, -2.30)*np.exp(-24536.9/tl)*MM #Yang+2021 via NIST,thorough review. 
#     k_inf=3.37E14*pow(tl/298.0, 0.90)*np.exp(-24536.9/tl) #Yang+2021 via NIST,thorough review. 
    
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(t3_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(t3_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(t3_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(t3_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(t3_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='t3, Hu2012')
# ax[2].plot(P_list, np.log10(t3_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='t3, new')

# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

# def t4_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kkM=7.07E-9*np.exp(-29467.0/tl)*MM
#     return kkM

# def t4_moses(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     Prescription from Moses+2011
#     """
#     k_0_M=5.88E-10*np.exp(-28265/tl)*MM #Breshears (1995) via Moses+2011
#     k_inf=1.3E11*np.exp(-30000.0/tl) #Tsang and Herron (1991) via Moses+2011
    
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM


# def t4_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     This is a tough one for selection. 
#     """
#     k_0_M=1.55E-9*np.exp(-30190.0/tl)*MM #Pham+2020, 700-5000
#     # k_0_M=4.73E-5*pow(tl/298.0,-3.96)*np.exp(-34760.6/tl)*MM #Kovac+2020, 700-3600 ###This is a tricky one. It leads to much lower rates at low temperatures. To avoid over-extrapolating (it is being used outside validity range, as usual), I disfavor it. Also, they only fit data (no theory) and did not compare against past work. So all in all, a little less good. 
        
#     k_inf=7.94E10*np.exp(-30911.7/tl) #Pham+2020, 925-2500
#     # k_inf=1.58E13*pow(tl/298.0,0.36)*np.exp(-34399.8/tl)*MM #Yu+2022, 300-4000 ###I Disfavor this because there doesn't seem to be a comparison to actual data  AFAIK. 
#     # k_inf=1.3E11*np.exp(-30000.0/tl) #Tsang and Herron (1991) via Moses+2011, 700-2500 ###Disfavored b/c old. Underlying data addressed by Pham.

#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(t4_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(t4_moses(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='blue')
# ax[0].plot(P_list, np.log10(t4_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(t4_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(t4_moses(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='blue')
# ax[1].plot(P_list, np.log10(t4_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(t4_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='t4, Hu2012')
# ax[2].plot(P_list, np.log10(t4_moses(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='blue', label='t4, Moses+2011')
# ax[2].plot(P_list, np.log10(t4_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='t4, new')

# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()


# def t5_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kkM=1.88E-4*pow(tl/298.0, -3.37)*np.exp(-37645.0/tl)*MM
#     return kkM

# def t5_rimmer(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     Prescription from Rimmer+2019 STAND2019, assuming format described in Rimmer+2016
#     Also in STAND2024 on GitHub.
    
#     NOTE: This is a thermochemically reversed reaction in STAND. 
#     """
    
#     #      # Thermo coefficients for [CO, H, CHO]
#     #      a1 = [0.35795335E+01, 0.25000000E+01, 4.23754610E+00]
#     #      a2 = [-0.61035369E-03, 0.00000000E+00, -3.32075257E-03]
#     #      a3 = [0.10168143E-05, 0.00000000E+00, 1.40030264E-05]
#     #      a4 = [0.90700586E-09, 0.00000000E+00, -1.34239995E-08]
#     #      a5 = [-0.90442449E-12, 0.00000000E+00, 4.37416208E-12]
#     #      a6 = [-0.14344086E+05, 0.25473660E+05, 3.87241185E+03]
#     #      a7 = [0.35084093E+01, -0.44668285E+00, 3.30834869E+00]

         
#     #      kkM = ReverseRate(n_total, mu, n_prod, n_react, a1, a2, a3, a4, a5, a6, a7, tl, k_f)
    
#     #forward reaction
#     k_0_M=9.00E-32*(tl/300.0)**-1.50*MM
#     k_inf=3.00E-11
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
    
#     n_prod=1
#     n_reac=2
#     n_tot=3
#     mu=[-1, -1, 1] #NO + O --> NO2 is the FORWARD reaction we are reversing
#     #thermo coefficents for [NO, O, NO2]. For 200-1000K. Via https://respecth.elte.hu/. /WARNING there are two sets of entries for NO, with (very slightly) different coefficients! We have chosen to use the later, recalculated version
#     a1=[4.0851799, 3.1682671, 3.9440312]
#     a2=[-0.00364693188, -0.00327931884, -0.001585429]
#     a3=[0.00000849607612, 0.00000664306396, 0.000016657812]
#     a4=[-0.00000000662405734, -0.00000000612806624, -0.000000020475426]
#     a5=[0.00000000000177647373, 0.00000000000211265971, 0.0000000000078350564]
#     a6=[9840.61267, 29122.2592, 2896.618]
#     a7=[2.83578236, 2.05193346, 6.3119919 ]
    
    
#     k_r=ReverseRate(mu, n_prod, n_reac, a1, a2, a3, a4, a5, a6, a7, tl, kkM)
#     return k_r


# def t5_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     This is a tough one for selection. The Tsang+1991 is a good review, but it is only valid at very high temperatures and it is being wildly extrapolated at lower temperatures. Moses+2011, Zahnle+2008, and Zahnle+2016 do not include this reaction. Rimmer+2019 do, and it is much faster...but it lacks citation support. 
#     """
#     k_0_M=1.88E-4*pow(tl/298.0, -3.37)*np.exp(-37645.0/tl)*MM #From Tsang+1991, but only valid 1400-2500. 
#     k_inf=5.48E15*pow(tl/298.0,-1.27)*np.exp(-36925.6/tl) #From Tsang+1991, but only valid 1400-2500. 

#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(t5_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(t5_rimmer(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='blue')
# ax[0].plot(P_list, np.log10(t5_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(t5_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(t5_rimmer(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='blue')
# ax[1].plot(P_list, np.log10(t5_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(t5_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='t5, Hu2012')
# ax[2].plot(P_list, np.log10(t5_rimmer(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='blue', label='t5, Rimmer+2019')
# ax[2].plot(P_list, np.log10(t5_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='t5, new')

# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

# def t6_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
    
#     Note this ONLY has the high-pressure limit. 
#     """
#     kkM=2.5E6*np.exp(-6100.0/tl)*np.ones(np.shape(MM))
#     return kkM

# def t6_rimmer_stand2024(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     Prescription from STAND2024 on GitHub
#     Note: I have no idea where the forward reaction came from. I don't see it in JPL or in NIST. 
#     """
    
#     #Forward reaction
#     k_0_M=2.81E-41*MM
#     k_inf=1.36E-21*pow(tl/300, -1.0)
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
    
#     n_prod=1
#     n_reac=2
#     n_tot=3
#     mu=[-1, -1, 1] #NO + O2 --> NO3 is the FORWARD reaction we are reversing
#     #thermo coefficents for [NO, O2, NO3]. For 200-1000K. Via https://respecth.elte.hu/
#     a1=[4.0851799, 3.78245636, 2.1735933]
#     a2=[-0.00364693188, -0.00299673416, 0.0104902685]
#     a3=[0.00000849607612, 0.00000984730201, 0.0000110472669]
#     a4=[-0.00000000662405734, -0.00000000968129509, -0.0000000281561867]
#     a5=[0.00000000000177647373, 0.00000000000324372837, 0.000000000013658396]
#     a6=[9840.61267, -1063.94356, 7812.90905]
#     a7=[2.83578236, 3.65767573, 14.602209]
    
    
#     k_r=ReverseRate(mu, n_prod, n_reac, a1, a2, a3, a4, a5, a6, a7, tl, kkM)   
    
#     return k_r

# def t6_rimmer_stand2015(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     Ref given is 24=Graham & Johnston 1978a, The photochemistry of NO3 and the kinetics of the N2O5-O3 system. HOWEVER, this reference absolutely does not have the high-pressure limit, which is weird. 
#     """
#     k_0_M=2.51E-14*np.exp(-1230/tl)*MM
#     k_inf=6.06E5*pow(tl/300, -1.0)*np.exp(1230/tl)
    
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM

# def t6_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     This is another discrepant one. We have adopted the values from the Johnston group via NIST, because they are from the same laboratory and good citation support. However, they are only valid for cool atmospheres. 
    
#     Therefore, it is noteworthy that they WILDLY diverge from STAND, being much higher at all pressures for T=150K, and higher for high pressures at higher temps. This is noteworthy because to my reading, the high-pressure limit is more trustworthy than the low-pressure limit.
    
#     Another potential pain point, worth checking in more detail.     
#     """
#     k_0_M=2.51E-14*np.exp(-1230.0/tl)*MM #This is the value listed in NIST, but reading the underlying references (Graham & Johnston 1978, commentary in Johnston+1986) it seems it is for the reaction NO2 + NO3 --> NO2 + NO + O2, which based on Johnston+1986 may be an independent reaction mechanism. Nevertheless, it matches Rimmer+2016, 2024 pretty well, so we are going to go with  it for now. 
#     k_inf=2.5E6*np.exp(-6100.0/tl) #This is from Johnston+1986, and it seems pretty solid. Notably, it diverges strongly from STAND, but I still trust it because it reproduces e.g Schott and Davidson 1958 pretty well. 

#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(t6_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(t6_rimmer_stand2015(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='blue')
# ax[0].plot(P_list, np.log10(t6_rimmer_stand2024(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='cyan')
# ax[0].plot(P_list, np.log10(t6_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(t6_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(t6_rimmer_stand2015(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='blue')
# ax[1].plot(P_list, np.log10(t6_rimmer_stand2024(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='cyan')
# ax[1].plot(P_list, np.log10(t6_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(t6_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='t6, Hu2012')
# ax[2].plot(P_list, np.log10(t6_rimmer_stand2015(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='blue', label='t6, Rimmer+2016')
# ax[2].plot(P_list, np.log10(t6_rimmer_stand2024(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='cyan', label='t6, Rimmer+2024')
# ax[2].plot(P_list, np.log10(t6_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='t6, new')

# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

# def t7_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
    
#     """
#     kkM=1.0E-3*pow(tl/298.0, -3.5)*np.exp(-11000.0/tl)*MM
#     return kkM

# def t7_reverseHu(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     Reverse the Hu reaction (M29) using the Burcat polynomials. 
#     Note: I have no idea where the forward reaction came from. I don't see it in JPL or in NIST. 
#     """
    
#     #Forward reaction NO2 + NO3 --> N2O5
#     RH=1.6E-12*pow(tl/298.0, 0.1)
#     LH=2.4E-30*pow(tl/298.0, -3.0)
#     ind=1/(1.0+pow(np.log10(LH*MM/RH),2.0));
#     kkM=LH*MM/(1.0+LH*MM/RH)*pow(0.6,ind);
#     n_prod=1
#     n_reac=2
#     n_tot=3
#     mu=[-1, -1, 1] #NO2 + NO3 --> N2O5 is the FORWARD reaction we are reversing
#     #thermo coefficents for [NO2, NO3, N2O5]. For 200-1000K. Via https://respecth.elte.hu/
#     a1=[3.9440312, 2.1735933, 3.68767456]
#     a2=[ -0.001585429, 0.0104902685, 0.0392120802]
#     a3=[0.000016657812, 0.0000110472669, -0.0000553770082]
#     a4=[-0.000000020475426, -0.0000000281561867, 0.0000000420097925]
#     a5=[0.0000000000078350564,  0.000000000013658396, -0.0000000000131260758]
#     a6=[2896.618, 7812.90905, -573.270648]
#     a7=[6.3119919, 14.602209, 12.1967861]
    
    
#     k_r=ReverseRate(mu, n_prod, n_reac, a1, a2, a3, a4, a5, a6, a7, tl, kkM)   
    
#     return k_r

# def t7_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     This is from Atkinson+2004. They have both high and low-pressure limits. Range given is 200-400. 
#     """
#     k_0_M=1.3E-3*pow(tl/300.0, -3.5)*np.exp(-11000.0/tl)*MM
#     k_inf=9.7E14*pow(tl/300.0, 0.1)*np.exp(-11080/tl)

#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(t7_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(t7_reverseHu(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='blue')
# ax[0].plot(P_list, np.log10(t7_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(t7_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(t7_reverseHu(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='blue')
# ax[1].plot(P_list, np.log10(t7_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(t7_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='t7, Hu2012')
# ax[2].plot(P_list, np.log10(t7_reverseHu(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='blue', label='t7, MEAC Reversed')
# ax[2].plot(P_list, np.log10(t7_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='t7, new')

# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

# def t8_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
    
#     """
#     kkM=5.48E-7*pow(tl/298.0, -1.24)*np.exp(-25312.0/tl)*MM
#     return kkM

# def t8_reverseHu(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     Reverse the Hu reaction (M15) using the Burcat polynomials. 
#     """
    
#     #Forward reaction H+NO-->HNO
#     RH=2.44E-10*pow(tl/298.0, -0.41)
#     LH=1.34E-31*pow(tl/298.0, -1.32)*np.exp(-370.5/tl)
#     ind=1.0/(1.0+pow(np.log10(LH*MM/RH),2));
#     kkM=LH*MM/(1.0+LH*MM/RH)*pow(0.6,ind)

#     n_prod=1
#     n_reac=2
#     mu=[-1, -1, 1] #H + NO --> HNO is the FORWARD reaction we are reversing
#     #thermo coefficents for [H, NO, HNO]. For 200-1000K. Via https://respecth.elte.hu/
#     a1=[2.5, 4.0851799, 4.53525574]
#     a2=[0.0, -0.00364693188, -0.00568543377]
#     a3=[0.0, 0.00000849607612, 0.000018519854]
#     a4=[0.0, -0.00000000662405734, -0.0000000171881225]
#     a5=[0.0,  0.00000000000177647373, 0.00000000000555818157]
#     a6=[25473.66, 9840.61267, 11618.3003]
#     a7=[-0.44668285, 2.83578236, 1.74315886]
    
    
#     k_r=ReverseRate(mu, n_prod, n_reac, a1, a2, a3, a4, a5, a6, a7, tl, kkM)   
    
#     return k_r

# def t8_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     This is from Tsang+1991.
#     """
#     k_0_M=0.01*pow(tl, -1.61)*np.exp(-25585.0/tl)*MM #This is for N2.
#     # k_0_M=6.4E-4*pow(tl, -1.24)*np.exp(-25172.0/tl)*MM #This is for CO2. 
#     k_inf=1.2E16*pow(tl, -0.43)*np.exp(-24922.0/tl)

#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(t8_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(t8_reverseHu(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='blue')
# ax[0].plot(P_list, np.log10(t8_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(t8_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(t8_reverseHu(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='blue')
# ax[1].plot(P_list, np.log10(t8_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(t8_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='t8, Hu2012')
# ax[2].plot(P_list, np.log10(t8_reverseHu(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='blue', label='t8, MEAC Reversed')
# ax[2].plot(P_list, np.log10(t8_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='t8, new')

# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

# def t9_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
    
#     """
#     kkM=1.98E-3*pow(tl/298.0, -3.8)*np.exp(-25257.0/tl)*MM
#     return kkM

# def t9_reverseHu(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     Reverse the Hu reaction (M42) using the Burcat polynomials. 
#     Note: M42 is only valid over narrow range of temperatures (JPL expression)
#     """
    
#     #Forward reaction OH+NO-->HNO2
    
#     RH=3.6E-11*pow(tl/298.0, -0.1)
#     LH=7.1E-31*pow(tl/298.0, -2.6)
#     ind=1.0/(1.0+pow(np.log10(LH*MM/RH),2.0))
#     kkM=LH*MM/(1.0+LH*MM/RH)*pow(0.6,ind)
    
#     n_prod=1
#     n_reac=2
#     mu=[-1, -1, 1] #OH + NO --> HNO2 is the FORWARD reaction we are reversing
#     #thermo coefficents for [OH, NO, HNO2]. For 200-1000K. Via https://respecth.elte.hu/
#     a1=[3.99198424, 4.0851799, 3.21415915]
#     a2=[-0.00240106655, -0.00364693188, 0.00812778066]
#     a3=[0.00000461664033, 0.00000849607612, 0.00000165998916 ]
#     a4=[-0.00000000387916306, -0.00000000662405734, -0.00000000952814708]
#     a5=[0.00000000000136319502,  0.00000000000177647373, 0.00000000000487131424]
#     a6=[3368.89836, 9840.61267, -10783.0727]
#     a7=[-0.103998477, 2.83578236, 9.82200056]
    
    
#     k_r=ReverseRate(mu, n_prod, n_reac, a1, a2, a3, a4, a5, a6, a7, tl, kkM)   
    
#     return k_r

# def t9_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     This is from Tsang+1991.
#     """
#     k_0_M=5.0E6*pow(tl, -3.8)*np.exp(-25340/tl)*MM #This is for N2.
#     # k_0_M=1.7E6*pow(tl, -3.59)*np.exp(-25250/tl)*MM #This is for CO2. 
#     k_inf=1.2E19*pow(tl, -1.23)*np.exp(-25010.0/tl)

#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(t9_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(t9_reverseHu(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='blue')
# ax[0].plot(P_list, np.log10(t9_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(t9_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(t9_reverseHu(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='blue')
# ax[1].plot(P_list, np.log10(t9_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(t9_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='t9, Hu2012')
# ax[2].plot(P_list, np.log10(t9_reverseHu(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='blue', label='t9, MEAC Reversed')
# ax[2].plot(P_list, np.log10(t9_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='t9, new')

# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()


# def t10_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
    
#     """
#     kkM=1.15E-6*np.exp(-23092.0/tl)*MM
#     return kkM

# def t10_reverseHu(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     Reverse the Hu reaction (M43) using the Burcat polynomials. 
#     Note: M43 is only valid over narrow range of temperatures (JPL expression)
#     Note: This is just the reaction for the "main" channel at Earth atmosphere conditions (5-15% according to JPL2020), forming nitrous acid (HONO2). A minor channel forms pernitrous acid (HOONO). This reaction is not in MEAC (but perhaps should be?)
#     """
    
#     #Forward reaction OH+NO2-->HNO3
    
#     RH=2.8E-11
#     LH=1.8E-30*pow(tl/298.0, -3.0)
#     ind=1.0/(1.0+pow(np.log10(LH*MM/RH),2.0))
#     kkM=LH*MM/(1.0+LH*MM/RH)*pow(0.6,ind)
    
#     n_prod=1
#     n_reac=2
#     mu=[-1, -1, 1] #OH + NO2 --> HNO3 is the FORWARD reaction we are reversing
#     #thermo coefficents for [OH, NO2, HNO3]. For 200-1000K. Via https://respecth.elte.hu/
#     a1=[3.99198424, 3.9440312, 1.69329154]
#     a2=[-0.00240106655, -0.001585429, 0.0190167702]
#     a3=[0.00000461664033, 0.000016657812, -0.00000825176697]
#     a4=[-0.00000000387916306, -0.000000020475426, -0.00000000606113827]
#     a5=[0.00000000000136319502,  0.0000000000078350564, 0.00000000000465236978]
#     a6=[3368.89836, 2896.618, -17419.8909]
#     a7=[-0.103998477, 6.3119919, 17.1839838]
    
    
#     k_r=ReverseRate(mu, n_prod, n_reac, a1, a2, a3, a4, a5, a6, a7, tl, kkM)   
    
#     return k_r

# def t10_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M=1.15E-6*np.exp(-23092.0/tl)*MM# From Chakraborty+1998, stated good 300-2000K
#     k_inf=9.33E15*np.exp(-24657.0/tl) #From Smith+1978, stated good 220-1000K

#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(t10_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(t10_reverseHu(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='blue')
# ax[0].plot(P_list, np.log10(t10_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(t10_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(t10_reverseHu(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='blue')
# ax[1].plot(P_list, np.log10(t10_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(t10_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='t10, Hu2012')
# ax[2].plot(P_list, np.log10(t10_reverseHu(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='blue', label='t10, MEAC Reversed')
# ax[2].plot(P_list, np.log10(t10_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='t10, new')

# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()


# def t11_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
    
#     """
#     kkM=8.0E-2*pow(tl/298.0, -6.55)*np.exp(-26099.0/tl)*MM
#     return kkM

# def t11_STAND2015(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     Reverse the version of this reaction present in STAND2015 (Rimmer+2016)
#     This is super weird...the parameters given for the FORWARD reaction in STAND (HO2+NO-->HNO3) are the same as those given for the REVERSE reaction in NIST (HNO3-->HO2+NO).
#     Huh, the reference given for the reaction in STAND is 169=Sander+2011. However, Sander+2011 again do not make a recommendation for this reaction. So, I don't know what's up with it. 
#     """
    
#     #Forward reaction HO2+NO-->HNO3
    
#     k_0_M=8.00E-2*pow(tl/300.0, -6.55)*np.exp(-26110.0/tl)*MM
#     k_inf=5.56E17*pow(tl/300.0, -2.27)*np.exp(-26300.0/tl)
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
    
#     n_prod=1
#     n_reac=2
#     mu=[-1, -1, 1] #HO2 + NO --> HNO3 is the FORWARD reaction we are reversing
#     #thermo coefficents for [HO2, NO, HNO3]. For 200-1000K. Via https://respecth.elte.hu/
#     a1=[4.30179807, 4.21859896, 1.69329154]
#     a2=[-0.00474912097, -0.00463988124, 0.0190167702]
#     a3=[0.0000211582905, 0.0000110443049, -0.00000825176697]
#     a4=[-0.0000000242763914, -0.00000000934055507, -0.00000000606113827]
#     a5=[0.00000000000929225225, 0.00000000000280554874, 0.00000000000465236978]
#     a6=[264.018485, 2896.618, 9845.09964]
#     a7=[3.7166622, 6.3119919, 2.28061001]
    
    
#     k_r=ReverseRate(mu, n_prod, n_reac, a1, a2, a3, a4, a5, a6, a7, tl, kkM)   
    
#     return k_r


# def t11_STAND2024(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     Reverse the version of this reaction present in STAND2024.
#     """
    
#     #Forward reaction HO2+NO-->HNO3
    
#     k_0_M=6.57E-33*pow(tl/300.0, 0.18)*MM
#     k_inf=4.76E-11*pow(tl/300.0, -0.82)
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
    
#     n_prod=1
#     n_reac=2
#     mu=[-1, -1, 1] #HO2 + NO --> HNO3 is the FORWARD reaction we are reversing
#     #thermo coefficents for [HO2, NO, HNO3]. For 200-1000K. Via https://respecth.elte.hu/
#     a1=[4.30179807, 4.21859896, 1.69329154]
#     a2=[-0.00474912097, -0.00463988124, 0.0190167702]
#     a3=[0.0000211582905, 0.0000110443049, -0.00000825176697]
#     a4=[-0.0000000242763914, -0.00000000934055507, -0.00000000606113827]
#     a5=[0.00000000000929225225, 0.00000000000280554874, 0.00000000000465236978]
#     a6=[264.018485, 2896.618, 9845.09964]
#     a7=[3.7166622, 6.3119919, 2.28061001]
    
    
#     k_r=ReverseRate(mu, n_prod, n_reac, a1, a2, a3, a4, a5, a6, a7, tl, kkM)   
    
#     return k_r

# def t11_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     #These are from Zhu+2003, stated good 295-1200K. HOWEVER, the forward reaction is not in MEAC (that is, HO2 + NO --> HNO3) AND the JPL evaluation casts doubt on whether the forward reaction occurs at all. 
    
#     Actually, reading Zhu+2003, the rate constants given are actually for HNO3 --> OH + NO2! This appears to be an error by NIST, pure and simple. 
    
#     Therefore, zeroing this reaction for now. 
#     WARNING this really needs to be investigated. 
#     """
#     k_0_M=8.0E-2*pow(tl/298.0, -6.55)*np.exp(-26099.0/tl)*MM
#     k_inf=5.56E17*pow(tl/298.0,-2.27)*np.exp(-26341.0/tl) 

#     #kkM=k_0_M/(1.0+k_0_M/k_inf)
#     kkM=0*MM
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(t11_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(t11_STAND2015(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='blue')
# ax[0].plot(P_list, np.log10(t11_STAND2024(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='cyan')
# ax[0].plot(P_list, np.log10(t11_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(t11_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(t11_STAND2015(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='blue')
# ax[1].plot(P_list, np.log10(t11_STAND2024(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='cyan')
# ax[1].plot(P_list, np.log10(t11_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(t11_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='t11, Hu2012')
# ax[2].plot(P_list, np.log10(t11_STAND2015(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='blue', label='t11, STAND2015')
# ax[2].plot(P_list, np.log10(t11_STAND2024(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='cyan', label='t11, STAND2024')
# ax[2].plot(P_list, np.log10(t11_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='t11, new')

# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()


# def t12_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     From Simonaitas+1989
#     """
#     kkM=6.31E17*np.exp(-13110.0/tl)*MM
#     return kkM

# def t12_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     This reaction is probably heterogenous, and should be excised. 
#     """
#     kkM=0*MM
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(t12_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(t12_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(t12_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(t12_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(t12_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='t12, Hu2012')
# ax[2].plot(P_list, np.log10(t12_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='t12, new')

# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

# def t13_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kkM=4.10E-5*np.exp(-10600.0/tl)*MM
#     return kkM

# def t13_reverseHu(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     Reverse the Hu reaction (M20) using the Burcat polynomials. 
#     """
    
#     #Forward reaction HO2+NO2-->HNO4
    
#     RH=4.0E-12*pow(tl/298.0, -0.3)
#     LH=1.9E-31*pow(tl/298.0, -3.4)
#     ind=1.0/(1.0+pow(np.log10(LH*MM/RH),2.0))
#     kkM=LH*MM/(1.0+LH*MM/RH)*pow(0.6,ind)
    
#     n_prod=1
#     n_reac=2
#     mu=[-1, -1, 1] #HO2 + NO2 --> HNO4 is the FORWARD reaction we are reversing
#     #thermo coefficents for [HO2, NO2, HNO4=HOONO2]. For 200-1000K. Via https://respecth.elte.hu/
#     a1=[4.30179807, 3.9440312, 2.44847749]
#     a2=[-0.00474912097, -0.001585429, 0.0285012019]
#     a3=[0.0000211582905, 0.000016657812, -0.0000293784944]
#     a4=[-0.0000000242763914, -0.000000020475426, 0.0000000150460407]
#     a5=[0.00000000000929225225,  0.0000000000078350564, -0.00000000000295996331]
#     a6=[264.018485, 2896.618, -8611.84484]
#     a7=[3.7166622, 6.3119919, 14.421696]
    
    
#     k_r=ReverseRate(mu, n_prod, n_reac, a1, a2, a3, a4, a5, a6, a7, tl, kkM)   
    
#     return k_r



# def t13_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     WARNING this really needs to be investigated. 
#     """
#     k_0_M=4.10E-5*np.exp(-10650.0/tl)*MM #Atkinson+2004, 260-300K 
#     k_inf=4.8E15*np.exp(-11170.0/tl) #Atkinson+2004, 260-300K
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(t13_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(t13_reverseHu(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[0].plot(P_list, np.log10(t13_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='blue')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(t13_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(t13_reverseHu(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].plot(P_list, np.log10(t13_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='blue')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(t13_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='t13, Hu2012')
# ax[2].plot(P_list, np.log10(t13_reverseHu(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='t13, M20-reversed')
# ax[2].plot(P_list, np.log10(t13_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='blue', label='t13, new')

# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

# #T14 CH2O --> CO + H2
# def t14_old(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     kkM=7.7E-9*np.exp(-33075.0/tl)*MM
#     return kkM

# def t14_alt(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
    
#     k_0_M=7.7E-9*np.exp(-33075.0/tl)*MM #Troe+2005 for Kr
#     k_inf=1.51E14*np.exp(-50878.0/tl) #de Martins+1998
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
    
#     return kkM



# def t14_new(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     k_0_M=9.4E-9*np.exp(-33140.0/tl)*MM #Troe+2005 for Ar
#     k_inf=3.7E13*np.exp(-36220.0/tl) #Troe+2005
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
#     return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(t14_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(t14_alt(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[0].plot(P_list, np.log10(t14_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='blue')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(t14_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(t14_alt(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].plot(P_list, np.log10(t14_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='blue')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(t14_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='t14, Hu2012')
# ax[2].plot(P_list, np.log10(t14_alt(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='t14, Alternative')
# ax[2].plot(P_list, np.log10(t14_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='blue', label='t14, new')

# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()


# t21
# def t21_old(tl,MM):
#       """
#      T in K
#      MM in cm^-3
#      """
#       kkM = 2.12e13*(tl/298.0)*1.22**np.exp(-43539.0/tl)*MM
#       return kkM

# def t21_new(tl,MM):
#       """
#      T in K
#      MM in cm^-3
#      Stand2020
#      """
#       kkM = 2.03e10*(tl/298.0)*1.22**np.exp(-43539.0/tl)*MM
#       return kkM

# # def t21_new(tl,MM):
# #      """
# #      T in K
# #      MM in cm^-3
# #      """
     
# #      # low pressure limit from Sebbar+2019 (theory)
# #      k_0_M = 3.3E-9*np.power(tl/298, 0.65)*np.exp(-38607.7/tl)*MM    

# #      # high pressure limit from RMG 
# #      k_inf = 2.85e19*np.power(tl,-1.52)*np.exp(-43137.0/tl) 

# #      kkM=k_0_M/(1.0+k_0_M/k_inf)

# #      return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(t21_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(t21_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(t21_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(t21_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(t21_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='t41, Hu2012')
# ax[2].plot(P_list, np.log10(t21_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='t41, new')

# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()



# def t22_old(tl,MM):
#       """
#      T in K
#      MM in cm^-3
#      """
#       kkM = 9.51e15*(tl/298.0)**-1.02*np.exp(-46185.0/tl)*MM
#       return kkM


# def t22_new(tl,MM):
#     # STAND 2015 (Dombrowsky +1991)

#      """
#      T in K
#      MM in cm^-3
#      """
     
#      # low pressure limit from  STAND 2015(theory)
#      k_0_M = 1.16e-8*(tl/298)**-0*np.exp(-33400/tl)*MM    

#      #  # high pressure limit from STAND 2015 (theory)
#      k_inf = 2.8e11*(tl/298.0)**(-1)*np.exp(-33400.0/tl) * np.ones(len(MM))

#      kkM= k_0_M/(1.0+k_0_M/k_inf)

#      return kkM

# def t22_new2(tl,MM):
#     # NIST(Jasper+2007)

#      """
#      T in K
#      MM in cm^-3
#      """
     
#      # low pressure limit from  (theory)
#      k_0_M = (8.03e4*(tl/300)**-10.2*np.exp(-52454/tl) + 6.245*(tl/300)**-6.577*np.exp(-48007/tl) )*MM    

#      #  # high pressure limit from Jasper+2007 (theory)
#      k_inf = (9.443e15*(tl/300)**-1.0117*np.exp(-46156/tl))* np.ones(len(MM))

#      kkM= k_0_M/(1.0+k_0_M/k_inf)

#      return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(t22_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(t22_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[0].plot(P_list, np.log10(t22_new2(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='blue')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(t22_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(t22_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].plot(P_list, np.log10(t22_new2(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='blue')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(t22_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='t22, Hu2012')
# ax[2].plot(P_list, np.log10(t22_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='STAND2020, new')
# ax[2].plot(P_list, np.log10(t22_new2(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='blue', label='JASPER2007, new')

# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

# plt.show()

# #t23
# def t23_old(tl,MM):
#       """
#      T in K
#      MM in cm^-3
#      """
#       kkM = 6.0e14*np.exp(-21288.0/tl)*np.ones(len(MM))
#       return kkM


# def t23_STAND2015(tl,MM):
#      """
#      T in K
#      MM in cm^-3
#      """
#      # rate laws from STAND+2015
#      k_0_M = 2.48e-5*(tl/298.0)**1*np.exp(-21300/tl)*MM
#      k_inf = 6.0e14*np.exp(-21288.0/tl)
#      kkM = k_0_M / (1.0 + k_0_M / k_inf)

#      return kkM

# def t23_JASPER2009(tl,MM):
#      """
#      T in K
#      MM in cm^-3
#      """
#      # rate laws from JASPER+2009
#      k_0_M = 1.773*(tl/298.0)**-7.502*np.exp(-23531/tl)*MM
#      k_inf = 5.684e16*(tl/298.0)**-1.153*np.exp(-22270/tl)
#      kkM = k_0_M / (1.0 + k_0_M / k_inf)

#      return kkM

# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(t23_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(t23_STAND2015(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[0].plot(P_list, np.log10(t23_JASPER2009(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='pink')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(t23_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(t23_STAND2015(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].plot(P_list, np.log10(t23_JASPER2009(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='pink')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(t23_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='t47, Hu2012')
# ax[2].plot(P_list, np.log10(t23_STAND2015(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='t47, Stand+2019')
# ax[2].plot(P_list, np.log10(t23_JASPER2009(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='pink', label='t47, new')

# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

# plt.show()

# #t25
# def t25_old(tl,MM):
#       """
#      T in K
#      MM in cm^-3
#      """
#       kkM = 0.19*(tl/298.0)**-7.5*np.exp(-22852.0/tl) *MM
#       return kkM


# def t25_Baulch2005(tl,MM):
#      """
#      T in K
#      MM in cm^-3
#      """
#      # rate laws from Baulch+2005
#      k_0_M = 4.3e3*(tl)**-3.4*np.exp(-18020/tl)*MM
#      k_inf = 3.9e8*(tl)**-1.62*np.exp(-18650/tl)
#      #k_0_M = 6.6e3*(tl)**-3.5*np.exp(-18070/tl)*MM
#      kkM = (k_0_M / (1.0 + k_0_M / k_inf))*7.37e-4*tl**0.8

#      return kkM

# def t25_Baulch1992(tl,MM):
#      """
#      T in K
#      MM in cm^-3
#      """
#      # rate laws from Baulch+1992
#      k_0_M = 6.9e17*(tl)**-7.5*np.exp(-22900/tl)*MM
#      k_inf = 2e14**np.exp(-20000/tl)
#      kkM = (k_0_M / (1.0 + k_0_M / k_inf))*0.35

#      return kkM



# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(t25_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(t25_Baulch2005(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[0].plot(P_list, np.log10(t25_Baulch1992(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='pink')
# #ax[0].plot(P_list, np.log10(t25_STAND2015(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='blue')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(t25_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(t25_Baulch2005(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].plot(P_list, np.log10(t25_Baulch1992(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='pink')
# #ax[1].plot(P_list, np.log10(t25_STAND2015(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='blue')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(t25_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='t47, Hu2012')
# ax[2].plot(P_list, np.log10(t25_Baulch2005(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='t47, Baulch+2005')
# ax[2].plot(P_list, np.log10(t25_Baulch1992(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='pink', label='t47, Baulch+1992')
# #ax[2].plot(P_list, np.log10(t25_STAND2015(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='blue')

# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

# plt.show()


# #t26
# def t26_old(tl,MM):
#       """
#      T in K
#      MM in cm^-3
#      """
#       kkM = 5.8e-8*np.exp(-35961.0/tl)*MM
#       return kkM


# def t26_STAND2015(tl,MM):
#      """
#      T in K
#      MM in cm^-3
#      """
#      # rate laws from Baulch+1994 + Tsang+1986
#      k_0_M = 5.8e-8*np.exp(-35961.0/tl)*MM
#      k_inf = 10**(12.9)*(tl)**0.44*np.exp(-44700/tl)
#      kkM = (k_0_M / (1.0 + k_0_M / k_inf))

#      return kkM


# def t26_BAULCH2005(tl,MM):
#      """
#      T in K
#      MM in cm^-3
#      """
#      # rate laws from Baulch+2005 + Tsang+1986
#      k_0_M = 3.4e-7*np.exp(-39390/tl)*MM
#      k_inf = 10**(12.9)*(tl)**0.44*np.exp(-44700/tl)
#      kkM = (k_0_M / (1.0 + k_0_M / k_inf))

#      return kkM

# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(t26_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(t26_STAND2015(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[0].plot(P_list, np.log10(t26_BAULCH2005(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='blue')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(t26_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(t26_STAND2015(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].plot(P_list, np.log10(t26_BAULCH2005(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='blue')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(t26_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='t47, Hu2012')
# ax[2].plot(P_list, np.log10(t26_STAND2015(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='t47, STAND+2015')
# ax[2].plot(P_list, np.log10(t26_BAULCH2005(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='blue',label = 't47, BAULCH+2005')

# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

# plt.show()

# t41
# def t41_old(tl,MM):
#       """
#      T in K
#      MM in cm^-3
#      """
#       kkM = 9.47e-7*np.exp(-40051.0/tl)*MM
#       return kkM


# def t41_new(tl,MM):
#      """
#      T in K
#      MM in cm^-3
#      """
     
#      # low pressure limit from Sebbar+2019 (theory)
#      k_0_M = 3.3E-9*np.power(tl/298, 0.65)*np.exp(-38607.7/tl)*MM    

#      # high pressure limit from RMG 
#      k_inf = 2.85e19*np.power(tl,-1.52)*np.exp(-43137.0/tl) 

#      kkM=k_0_M/(1.0+k_0_M/k_inf)

#      return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(t41_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(t41_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(t41_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(t41_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(t41_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='t41, Hu2012')
# ax[2].plot(P_list, np.log10(t41_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='t41, new')

# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

#t42
# def t42_old(tl,MM):
#       """
#      T in K
#      MM in cm^-3
#      """
#       kkM = 1.4e-8*np.exp(-29467.0/tl)*MM
#       return kkM


# def t42_new(tl,MM):
#      """
#      T in K
#      MM in cm^-3
#      """
     
#      # low pressure limit from Tsuchiya+1994 (estimated)
#      k_0_M = 1.40E-8*np.exp(-29500.0/tl)*MM    

#      # high pressure limit from Rimmer+2021 (estimated)
#      k_inf = 3.38E+11*(tl/300)**(-1.0)*np.exp(-29500.0/tl)

#      kkM=k_0_M/(1.0+k_0_M/k_inf)

#      return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(t42_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(t42_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(t42_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(t42_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(t42_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='t42, Hu2012')
# ax[2].plot(P_list, np.log10(t42_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='t42, new')

# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()


#t43 
# def t43_old(tl,MM):
#       """
#      T in K
#      MM in cm^-3
#      """
#       kkM = 4.16e-7*(tl/298.0)**(-3.29)*np.exp(-9610.0/tl)*MM
#       return kkM


# def t43_new(tl,MM):
#      """
#      T in K
#      MM in cm^-3
#      """
     
#      # low pressure limit from Goumri+1999 (theory)
#      k_0_M = 4.16e-7*(tl/298.0)**(-3.29)*np.exp(-9610.0/tl)*MM    

#      # high pressure limit from Goumri+1999 (theory)
#      k_inf = 2.03e11*tl**0.9*np.exp(-9240/tl)

#      kkM=k_0_M/(1.0+k_0_M/k_inf)

#      return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(t43_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(t43_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(t43_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(t43_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(t43_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='t43, Hu2012')
# ax[2].plot(P_list, np.log10(t43_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='t43, new')

# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

# plt.show()



# t46 
# def t46_old(tl,MM):
#       """
#      T in K
#      MM in cm^-3
#      """
#       kkM = 6.0e-11*np.exp(-7721.0/tl)*MM
#       return kkM


# def t46_new(tl,MM):
#      """
#      T in K
#      MM in cm^-3
#      """
     
#      # Input values
#      n_prod = 1
#      n_react = 2
#      n_total = 3
#      mu = [-1, -1, 1]  # CO, H, CH
     
#      # Thermo coefficients for [CO, H, CHO]
#      a1 = [0.35795335E+01, 0.25000000E+01, 4.23754610E+00]
#      a2 = [-0.61035369E-03, 0.00000000E+00, -3.32075257E-03]
#      a3 = [0.10168143E-05, 0.00000000E+00, 1.40030264E-05]
#      a4 = [0.90700586E-09, 0.00000000E+00, -1.34239995E-08]
#      a5 = [-0.90442449E-12, 0.00000000E+00, 4.37416208E-12]
#      a6 = [-0.14344086E+05, 0.25473660E+05, 3.87241185E+03]
#      a7 = [0.35084093E+01, -0.44668285E+00, 3.30834869E+00]
     
#      k_f_0_M = 1.40e-34 * np.exp(-100 / tl) * MM
#      k_f_inf = 1.96e-13 * np.exp(-1370 / tl)
#      k_f = k_f_0_M / (1.0 + k_f_0_M / k_f_inf)
     
#      kkM = ReverseRate(mu, n_prod, n_react, a1, a2, a3, a4, a5, a6, a7, tl, k_f)

#      return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(t46_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(t46_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(t46_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(t46_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(t46_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='t46, Hu2012')
# ax[2].plot(P_list, np.log10(t46_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='t46, new')

# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

# plt.show()

# t47
# def t47_old(tl,MM):
#       """
#      T in K
#      MM in cm^-3
#      """
#       kkM = np.ones_like(MM)*1.69e14*(tl/298.0)**(-0.39)*np.exp(-13230.0/tl)
#       return kkM


# def t47_stand2019(tl,MM):
#      """
#      T in K
#      MM in cm^-3
#      """
#      # rate laws from Stand+2019
#      k_0_M = 6.59e-05*(tl/300)**(-2.7)*np.exp(-1.54e+04/tl) * MM
#      k_inf = 2.87e+13*(tl/300)**1.31*np.exp(-1.58e+04/tl)
#      kkM = k_0_M / (1.0 + k_0_M / k_inf)

#      return kkM

# def t47_new(tl,MM):
#      """
#      T in K
#      MM in cm^-3
#      """
#      # rate laws from Hippler+2001
#      k_0_M = 1.17e-06*(tl/298)**(-3.0)*np.exp(-12268/tl) * MM
#      k_inf = 6.8e+13*np.exp(-13230/tl)
#      F_c=0.97-tl/1950
#      M_M_c=k_0_M/k_inf
#      N=0.75-1.27*np.log10(F_c)
#      F=np.power(F_c, 1.0/(1.0+pow((np.log10(M_M_c)/N),2.0)))
#      kkM=k_0_M*k_inf/(k_0_M+k_inf)*F

#      return kkM

# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(t47_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(t47_stand2019(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[0].plot(P_list, np.log10(t47_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='pink')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(t47_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(t47_stand2019(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].plot(P_list, np.log10(t47_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='pink')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(t47_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='t47, Hu2012')
# ax[2].plot(P_list, np.log10(t47_stand2019(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='t47, Stand+2019')
# ax[2].plot(P_list, np.log10(t47_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='pink', label='t47, new')

# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

# plt.show()

#t48
# def t48_old(tl,MM):
#       """
#      T in K
#      MM in cm^-3
#      """
#       kkM = 7.67e-6*(tl/298.0)**-1.89*np.exp(-17800.0/tl)*MM
#       return kkM


# def t48_Larson1988(tl,MM):
#      """
#      T in K
#      MM in cm^-3
#      """
#      # rate laws from Larson+1988
#      k_0_M = 7.67e-6*(tl/298.0)**-1.89*np.exp(-17750.0/tl)*MM
#      k_inf = 1.21e14*(tl/298.0)**0.53*np.exp(-17100.0/tl)
#      kkM = k_0_M / (1.0 + k_0_M / k_inf)

#      return kkM

# def t48_Golden1998(tl,MM):
#      """
#      T in K
#      MM in cm^-3
#      """
#      # rate laws from Golden+1998
#      k_0_M = 4.92e-29*(tl/298.0)**-2.4*np.exp(-18862.0/tl)*MM
#      k_inf = 2.985e1*(tl/298.0)**0.13*np.exp(-18349.0/tl)
#      kkM = k_0_M / (1.0 + k_0_M / k_inf)

#      return kkM

# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(t48_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(t48_Larson1988(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[0].plot(P_list, np.log10(t48_Golden1998(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='pink')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(t48_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(t48_Larson1988(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].plot(P_list, np.log10(t48_Golden1998(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='pink')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(t48_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='t48, Hu2012')
# ax[2].plot(P_list, np.log10(t48_Larson1988(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='t48, Larson+1988')
# ax[2].plot(P_list, np.log10(t48_Golden1998(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='pink', label='t48, Golden+1998')

# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

# # t49
# def t49_old(tl,MM):
#       """
#      T in K
#      MM in cm^-3
#      """
#       kkM = 1.25e-5*(tl/298.0)**-3.02*np.exp(-17560.0/tl)*MM
#       return kkM


# def t49_Larson1988(tl,MM):
#      """
#      T in K
#      MM in cm^-3
#      """
#      # rate laws from Larson+1988
#      k_0_M = 1.25e-5*(tl/298.0)**-3.02*np.exp(-17560.0/tl)*MM
#      k_inf = 1.e13*(tl/298.0)**0.31*np.exp(-16570.0/tl)
#      kkM = k_0_M / (1.0 + k_0_M / k_inf)

#      return kkM

# def t49_Golden1998(tl,MM):
#      """
#      T in K
#      MM in cm^-3
#      """
#      # rate laws from Golden+1998
#      k_0_M = 7.21e-31*(tl/298.0)**-3.15*np.exp(-18629.0/tl)*MM
#      k_inf = 1.25e2*(tl/298.0)**0.41*np.exp(-17783.0/tl)
#      kkM = k_0_M / (1.0 + k_0_M / k_inf)

#      return kkM

# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(t49_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(t49_Larson1988(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[0].plot(P_list, np.log10(t49_Golden1998(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='pink')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(t49_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(t49_Larson1988(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].plot(P_list, np.log10(t49_Golden1998(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='pink')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(t49_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='t49, Hu2012')
# ax[2].plot(P_list, np.log10(t49_Larson1988(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='t49, Larson+1988')
# ax[2].plot(P_list, np.log10(t49_Golden1998(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='pink', label='t49, Golden+1998')

# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

# plt.show()

# #t27
# def t27_old(tl,MM):
#       """
#      T in K
#      MM in cm^-3
#      """
#       kkM = 1.69e-6*np.exp(-16838.0/tl)*MM
#       return kkM


# def t27_baulch2005(tl,MM):
#      """
#      T in K
#      MM in cm^-3
#      """
#      # rate laws from Baulch+2005
#      k_0_M = 1.7e-6*np.exp(-16800/tl)*MM
#      k_inf = 8.2e13*np.exp(-20070/tl)
#      kkM = k_0_M / (1.0 + k_0_M / k_inf)

#      return kkM

# def t27_Curran2006(tl,MM):
#      """
#      T in K
#      MM in cm^-3
#      """
#      # rate laws from Curran+2006
#      k_0_M = 1.7e-6*np.exp(-16800/tl)*MM
#      k_inf = 3.06e10*(tl**0.95)*np.exp(-18589/tl)
#      kkM = k_0_M / (1.0 + k_0_M / k_inf)

#      return kkM

# def t27_feng1993(tl,MM):
#      """
#      T in K
#      MM in cm^-3
#      """
#      # rate laws from Curran+2006
#      k_0_M = 6.63e9*(tl**-4.99)*np.exp(-20130/tl)*MM
#      k_inf = 1.11e10*(tl**1.037)*np.exp(-18504/tl)
#      kkM = k_0_M / (1.0 + k_0_M / k_inf)

#      return kkM


# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(t27_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(t27_baulch2005(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[0].plot(P_list, np.log10(t27_Curran2006(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='blue')
# ax[0].plot(P_list, np.log10(t27_feng1993(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='green')
# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(t27_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(t27_baulch2005(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].plot(P_list, np.log10(t27_Curran2006(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='blue')
# ax[1].plot(P_list, np.log10(t27_feng1993(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='green')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(t27_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='t27, Hu2012')
# ax[2].plot(P_list, np.log10(t27_baulch2005(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='t27, Baulch+2005')
# ax[2].plot(P_list, np.log10(t27_Curran2006(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='blue', label='t27, Curran+2005')
# ax[2].plot(P_list, np.log10(t27_feng1993(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='green', label='t27, Feng+1993')

# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()

# plt.show()

# t52
# use lower limit from Stand+2019 and upper limit from Zaslonko+1993 to compute the reverse rate (CH3O + NO2 --> CH3NO3), 
# and compare with that implemented in MEAC to determine if the rates are reasonable

# def t52_reverse_stand2019(tl,MM):
# #      """
# #      T in K
# #      MM in cm^-3
# #      """
#     # reverse reaction using Stand+2019 lower limit and Zaslonko+1993 upper limit
#     # the reaction we are reversing: CH3NO3 --> CH3O + NO2
#     n_prod = 2
#     n_react = 1
#     mu = [-1, 1, 1]

#     # Thermo coefficients for [CH3NO3, CH3O, NO2] from https://respecth.elte.hu/burcat.php

#     a1 = [3.91363583, 3.71180502, 3.9440312]
#     a2 = [0.0152137945, -0.00280463306, -0.001585429]
#     a3 = [0.0000173479131, 0.0000376550971, 0.000016657812]
#     a4 = [-0.0000000337074473, -0.0000000473072089, -0.000000020475426]
#     a5 = [0.0000000000144322204, 0.000000000018658842, 0.0000000000078350564]
#     a6 = [-16610.3232, 1295.6976, 2896.618]
#     a7 = [9.44208392, 6.57240864, 6.3119919]

#     # forward reaction rate
#     k_f_0_M = 4.14e-7 * (tl/300) * np.exp(-16800 / tl) * MM
#     k_f_inf = 1.0e13 * np.exp(-16800 / tl)
#     k_f = k_f_0_M / (1.0 + k_f_0_M / k_f_inf)

#     kkM = ReverseRate(mu, n_prod, n_react, a1, a2, a3, a4, a5, a6, a7, tl, k_f)

#     return kkM

# def t52_reverse_MEAC(tl, MM):

#     # reverse reaction rate (M28: CH3O + NO2 --> CH3NO3) implemented in MEAC
#     RH=1.9E-11*np.power(tl/298.0, -1.8)
#     LH=5.5E-29*np.power(tl/298.0, -4.4)
#     ind=1.0/(1.0+pow(np.log10(LH*MM/RH),2.0))
#     kkM=LH*MM/(1.0+LH*MM/RH)*np.power(0.6,ind)

#     return kkM

# def t52_reverse_M28(tl, MM):
#     """
#     T in K
#     MM in cm^-3
#     """
#     # compute the reaction rate for t52: CH3NO3 --> CH3O + NO2 by reversing the forward reaction M28: CH3O + NO2 --> CH3NO3
#     # the reaction we are reversing: CH3O + NO2 --> CH3NO3
#     n_prod = 1
#     n_react = 2
#     mu = [-1, -1, 1]

#     # Thermo coefficients for [CH3O, NO2, CH3NO3] from https://respecth.elte.hu/burcat.php

#     a1 = [3.71180502, 3.9440312, 3.91363583]
#     a2 = [-0.00280463306, -0.001585429, 0.0152137945]
#     a3 = [0.0000376550971, 0.000016657812, 0.0000173479131]
#     a4 = [-0.0000000473072089, -0.000000020475426, -0.0000000337074473]
#     a5 = [0.000000000018658842, 0.0000000000078350564, 0.0000000000144322204]
#     a6 = [1295.6976, 2896.618, -16610.3232]
#     a7 = [6.57240864, 6.3119919, 9.44208392]

#     # forward reaction rate
#     k_f_0_M = 5.5e-29 * (tl/298.)**(-4.4) * MM
#     k_f_inf = 1.9e-11 * (tl/298.)**(-1.8)
#     k_f = k_f_0_M / (1.0 + k_f_0_M / k_f_inf)

#     kkM = ReverseRate(mu, n_prod, n_react, a1, a2, a3, a4, a5, a6, a7, tl, k_f)

#     return kkM

# def t52_old(tl,MM):
#      #      """
#      #      T in K
#      #      MM in cm^-3
#      #      """

#      kkM = 1.0e13*np.exp(-16838.0/tl)*np.ones_like(MM)

#      return kkM

# def t52_new(tl,MM):
#      #      """
#      #      T in K
#      #      MM in cm^-3
#      #      """

#      k_0_M = 4.14e-7 * (tl/300) * np.exp(-16800 / tl) * MM
#      k_inf = 1.0e13 * np.exp(-16800 / tl)
#      kkM = k_0_M / (1.0 + k_0_M / k_inf)

#      return kkM

# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=False)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].set_ylim(-55,-30)
# ax[0].plot(P_list, np.log10(t52_reverse_M28(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='blue')
# ax[0].plot(P_list, np.log10(t52_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(t52_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')

# ax[1].set_title('T=300 K')
# ax[1].set_ylim(-30,-5)
# ax[1].plot(P_list, np.log10(t52_reverse_M28(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='blue')
# ax[1].plot(P_list, np.log10(t52_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(t52_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')

# ax[2].set_title('T=600 K')
# ax[2].set_ylim(-10,5)
# ax[2].plot(P_list, np.log10(t52_reverse_M28(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='blue', label='t52, reversed from m28 ')
# ax[2].plot(P_list, np.log10(t52_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='t52, Hu+2012')
# ax[2].plot(P_list, np.log10(t52_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='t52, new')

# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()
# plt.show()

# # t53

# def t53_old(tl, MM):
#     #      """
#     #      T in K
#     #      MM in cm^-3
#     #      """

#     kkM = 3.16e-4*(tl/298.0)**(-4.72)*np.exp(-13591.0/tl)*MM

#     return kkM

# def t53_new(tl, MM):
#     #      """
#     #      T in K
#     #      MM in cm^-3
#     #      """

#     # both limits from Somnitz+2004
#     k_0_M = 3.16e-4*(tl/298.0)**(-4.72)*np.exp(-13591.0/tl)*MM
#     k_inf = 4.34e14*np.exp(-12999/tl)
#     kkM = k_0_M / (1.0 + k_0_M / k_inf)

#     return kkM

# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=False)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].set_ylim(-35, -10)
# ax[0].plot(P_list, np.log10(t53_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(t53_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')

# ax[1].set_title('T=300 K')
# ax[1].set_ylim(-20, 10)
# ax[1].plot(P_list, np.log10(t53_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(t53_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')

# ax[2].set_title('T=600 K')
# ax[2].set_ylim(-20, 15)
# ax[2].plot(P_list, np.log10(t53_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='t53, Hu2012')
# ax[2].plot(P_list, np.log10(t53_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='t53, new')


# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()
# plt.show()

#t28

# def t28_old(tl, MM):
#     #      """
#     #      T in K
#     #      MM in cm^-3
#     #      """

#     # both limits from Somnitz+2004
#     kkM = 8.11e17*(tl/298.0)**(-1.23)*np.exp(-51356.0/tl)*np.ones(len(MM))

#     return kkM

# def t28_new(tl, MM):
#     #      """
#     #      T in K
#     #      MM in cm^-3
#     #      """

#     # both limits from Somnitz+2004
#     k_0_M = 10**42.838*(tl**-6.431)*np.exp(-53938/tl)*MM
#     k_inf = 10**20.947*(tl**-1.228)*np.exp(-51439/tl)*np.ones(len(MM))
#     Fc = 47.61*np.exp(-16182/tl) + np.exp(-tl/3371)
#     kkM = k_0_M / ((1.0 + k_0_M / k_inf)*Fc)

#     return kkM


# def t28_rimmer(tl,MM):
#     """
#     T in K
#     MM in cm^-3
#     Prescription from Rimmer+2019 STAND2019, assuming format described in Rimmer+2016
#     Also in STAND2024 on GitHub.
    
#     NOTE: This is a thermochemically reversed reaction in STAND. 
#     """
    
#     #      # Thermo coefficients for [CO, H, CHO]
#     #      a1 = [0.35795335E+01, 0.25000000E+01, 4.23754610E+00]
#     #      a2 = [-0.61035369E-03, 0.00000000E+00, -3.32075257E-03]
#     #      a3 = [0.10168143E-05, 0.00000000E+00, 1.40030264E-05]
#     #      a4 = [0.90700586E-09, 0.00000000E+00, -1.34239995E-08]
#     #      a5 = [-0.90442449E-12, 0.00000000E+00, 4.37416208E-12]
#     #      a6 = [-0.14344086E+05, 0.25473660E+05, 3.87241185E+03]
#     #      a7 = [0.35084093E+01, -0.44668285E+00, 3.30834869E+00]

         
#     #      kkM = ReverseRate(n_total, mu, n_prod, n_react, a1, a2, a3, a4, a5, a6, a7, tl, k_f)
    
#     #forward reaction
#     k_0_M=2e-28*(tl/300.0)**-1.50*MM
#     k_inf=1.7e-10
#     kkM=k_0_M/(1.0+k_0_M/k_inf)
    
#     n_prod=1
#     n_reac=2
#     n_tot=3
#     mu=[-1, -1, 1] #NO + O --> NO2 is the FORWARD reaction we are reversing
#     #thermo coefficents for [NO, O, NO2]. For 200-1000K. Via https://respecth.elte.hu/. /WARNING there are two sets of entries for NO, with (very slightly) different coefficients! We have chosen to use the later, recalculated version
#     a1 = [4.24185905, 2.5,         4.29142572]
#     a2 = [-0.00356905235, 0.0,     -0.00550154901]
#     a3 = [4.82667202e-05, 0.0,     5.99438458e-05]
#     a4 = [-5.85401009e-08, 0.0,    -7.08466469e-08]
#     a5 = [2.25804514e-11, 0.0,     2.68685836e-11]
#     a6 = [12969.0344, 25473.66,   -11522.2056]
#     a7 = [4.44703782, -0.44668285, 2.66678994]
    
    
#     k_r=ReverseRate(mu, n_prod, n_reac, a1, a2, a3, a4, a5, a6, a7, tl, kkM)
#     return k_r

# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(t28_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(t28_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
# ax[0].plot(P_list, np.log10(t28_rimmer(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='green')

# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(t28_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(t28_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
# ax[1].plot(P_list, np.log10(t28_rimmer(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='green')
# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(t28_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='t28, Hu2012')
# ax[2].plot(P_list, np.log10(t28_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='t28, new')
# ax[2].plot(P_list, np.log10(t28_rimmer(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='green', label='t28, rimmer')


# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()
# plt.show()

# t54
# def t54_old(tl, MM):
#          """
#          T in K
#          MM in cm^-3
#          """
         
#          kkM = 1.0e-12*MM+1.5E3
         
#          return kkM

# def t54_new(tl, MM):
#          """
#          T in K
#          MM in cm^-3
#          """
         
#          k_0_M = 1.0e-12*MM+1.5E3
#          k_inf = 1.0E10
#          kkM = k_0_M / (1.0 + k_0_M / k_inf)
         
#          return kkM

# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(t54_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(t54_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')

# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(t54_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(t54_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')

# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(t54_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='t54, Hu2012')
# ax[2].plot(P_list, np.log10(t54_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='t54, new')


# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()
# plt.show()

#t55
# def t55_old(tl, MM):
#          """
#          T in K
#          MM in cm^-3
#          """
         
#          kkM = 1.0e-11*MM+2.2E4
#          return kkM

# def t55_new(tl, MM):
#          """
#          T in K
#          MM in cm^-3
#          """
         
#          k_0_M = 1.0e-11*MM+2.2E4
#          k_inf = 1.0E10
#          kkM = k_0_M / (1.0 + k_0_M / k_inf)
         
#          return kkM

# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(t55_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(t55_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')

# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(t55_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(t55_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')

# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(t55_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='t55, Hu2012')
# ax[2].plot(P_list, np.log10(t55_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='t55, new')


# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()
# plt.show()

#t56
# def t56_old(tl, MM):
#          """
#          T in K
#          MM in cm^-3
#          """
         
#          kkM = 1.5e-13*MM+1.13E3
#          return kkM

# def t56_new(tl, MM):
#          """
#          T in K
#          MM in cm^-3
#          """
         
#          k_0_M = 1.5e-13*MM+1.13E3
#          k_inf = 1.0E10
#          kkM = k_0_M / (1.0 + k_0_M / k_inf)
         
#          return kkM

# fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
# markersizeval=5.

# ax[0].set_title('T=150 K')
# ax[0].plot(P_list, np.log10(t56_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
# ax[0].plot(P_list, np.log10(t56_new(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')

# ax[1].set_title('T=300 K')
# ax[1].plot(P_list, np.log10(t56_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
# ax[1].plot(P_list, np.log10(t56_new(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')

# ax[2].set_title('T=600 K')
# ax[2].plot(P_list, np.log10(t56_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='t56, Hu2012')
# ax[2].plot(P_list, np.log10(t56_new(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='t56, new')


# ax[0].set_ylabel('log10(k)')
# ax[1].set_ylabel('log10(k)')
# ax[2].set_ylabel('log10(k)')
# ax[2].set_xscale('log')
# ax[2].set_xlabel('Pressure (bar)')
# ax[2].set_yscale('linear')
# ax[2].legend()
# plt.show()

#t57
def t57_old(tl, MM):
         """
         T in K
         MM in cm^-3
         """
         
         kkM = 1.93e-4*(tl/298.0)**-2.44*np.exp(-62782.1/tl) * MM
         return kkM

def t57_NIST(tl, MM):
         """
         T in K
         MM in cm^-3
         """
         
         # reference from NIST (Tsang+1991)
         k_0_M = 1.93e-4*(tl/298.0)**-2.44*np.exp(-62782.1/tl) * MM
         k_inf = 4.15E15*(tl/298)**-0.93*np.exp(-62294/tl)
         Fc = 0.875 - 0.5e-4*tl
         kkM = k_0_M / ((1.0 + k_0_M / k_inf)*Fc)
         
         return kkM

def t57_STAND2019(tl, MM):
         """
         T in K
         MM in cm^-3
         """

         k_0_M = 6.14e-06*(tl/300)**-1.58*np.exp(-6.15e+04/tl)
         k_inf = 2.55e+12*np.ones_like(MM)
         kkM = k_0_M / (1.0 + k_0_M / k_inf)
         
         return kkM

fig, ax=plt.subplots(3, figsize=(8.,10.), sharex=True, sharey=True)
markersizeval=5.

ax[0].set_title('T=150 K')
ax[0].plot(P_list, np.log10(t57_old(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='red')
ax[0].plot(P_list, np.log10(t57_NIST(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='purple')
ax[0].plot(P_list, np.log10(t57_STAND2019(150.0, MM(P_list,150.0))), linewidth=2, linestyle='-', color='blue')

ax[1].set_title('T=300 K')
ax[1].plot(P_list, np.log10(t57_old(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='red')
ax[1].plot(P_list, np.log10(t57_NIST(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='purple')
ax[1].plot(P_list, np.log10(t57_STAND2019(300.0, MM(P_list,300.0))), linewidth=2, linestyle='-', color='blue')

ax[2].set_title('T=600 K')
ax[2].plot(P_list, np.log10(t57_old(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='red', label='t57, Hu2012')
ax[2].plot(P_list, np.log10(t57_NIST(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='purple', label='t57, NIST')
ax[2].plot(P_list, np.log10(t57_STAND2019(600.0, MM(P_list,600.0))), linewidth=2, linestyle='-', color='blue', label='t57, STAND2019')


ax[0].set_ylabel('log10(k)')
ax[1].set_ylabel('log10(k)')
ax[2].set_ylabel('log10(k)')
ax[2].set_xscale('log')
ax[2].set_xlabel('Pressure (bar)')
ax[2].set_yscale('linear')
ax[2].legend()
plt.show()