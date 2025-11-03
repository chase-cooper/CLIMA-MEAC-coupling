# -*- coding: iso-8859-1 -*-
"""
The purpose of this code is to process spectra of M-dwarf and K-dwarf exoplanets for use with the Hu+2012 photochemical code. 

We normalize all spectra to have the same total flux (integrated over wavelength) as the Sun at Earth, following Renyu's past conventions.
"""


########################
###Import useful libraries
########################
import numpy as np
import scipy.integrate
import scipy.optimize
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import pdb
from scipy import interpolate as interp
from astropy.io import fits

########################
###Define useful constants, all in CGS (via http://www.astro.wisc.edu/~dolan/constants.html)
########################

#Unit conversions
amu2g=1.66054e-24 #1 amu in g
bar2atm=0.9869 #1 bar in atm
Pa2bar=1.e-5 #1 Pascal in bar
bar2Pa=1.e5 #1 bar in Pascal
bar2barye=1.e6 #1 Bar in Barye (the cgs unit of pressure)
barye2bar=1.e-6 #1 Barye in Bar


#Fundamental constants
c=2.997924e10 #speed of light, cm/s
h=6.6260755e-27 #planck constant, erg/s
k=1.380658e-16 #boltzmann constant, erg/K

#Astronomical quantities
AU=1.496e13 #1 AU in cm
pc=3.086e18 #1 pc in cm
Lsun=3.9e33 #solar luminosity in erg/s

#masses
m_n2=28.01*amu2g #molecular mass of n2, converted to g
m_co2=44.01*amu2g #molecular mass of co2, converted to g


########################
###Define subfunction to integrate data to coarser bin
########################

def integrate_data(abscissa, data, leftedges, rightedges):
	"""
	Takes: abscissa and corresponding data to be integrated, left edges of new bins, right edges of new bins (abscissa, leftedges and rightedges must have same units)
	
	Returns: data integrated to new bins
	
	Method: functionalizes data using stepwise linear interpolation, then integrates using gaussian quadrature technique
	"""
	data_func=interp.interp1d(abscissa, data, kind='linear')
	
	num_bin=len(leftedges)
	
	data_integrated=np.zeros(num_bin)
	
	for ind in range(0, num_bin):
		data_integrated[ind]=scipy.integrate.quad(data_func, leftedges[ind], rightedges[ind], epsabs=0., epsrel=1.e-3, limit=100000)[0]/(rightedges[ind]-leftedges[ind])
	
	return data_integrated


########################################################################
########################################################################
########################################################################
###READ IN DATA
########################################################################
########################################################################
########################################################################

########################
###Initialize dicts to hold all data
#########################
#unbinned data0
wav={} #in nm
flux={} #in Watt m^-2 nm^-1

#binned data
wav_binned={}  #in nm
flux_binned={} #in Watt m^-2 nm^-1


########################
###Import quiescent Sun data
########################
###Modern Sun (Hu+2012) for comparison
hu_sun_wav, hu_sun_toa=np.genfromtxt('../Data/solar.txt', skip_header=0, skip_footer=0,usecols=(0,1), unpack=True) #units: nm, W m**-2 nm**-1
    
wav['Sun']=hu_sun_wav
flux['Sun']=hu_sun_toa
solar_constant=np.trapz(flux['Sun'], x=wav['Sun']) #solar constant. Integrated over nm, should be 1362 W m^-2...but isn't, quite, just a hair high. Could it be due to stellar activity?

print('TOA flux (Sun, 67 nm-10 um): {0:1.3e}'.format(solar_constant))

###Also get M-dwarf spectra from Seager+2013ab for comparison
hu_ma_wav, hu_ma_toa=np.genfromtxt('../Data/GJ1214New.txt', skip_header=0, skip_footer=0,usecols=(0,1), unpack=True) #units: nm, W m**-2 nm**-1
hu_mq_wav, hu_mq_toa=np.genfromtxt('../Data/T3000.txt', skip_header=0, skip_footer=0,usecols=(0,1), unpack=True) #units: nm, W m**-2 nm**-1

wav['Ma']=hu_ma_wav
flux['Ma']=hu_ma_toa
print('TOA flux (Ma, 67 nm-10 um): {0:1.3e}'.format(np.trapz(flux['Ma'], x=wav['Ma'])))

wav['Mq']=hu_mq_wav
flux['Mq']=hu_mq_toa
print('TOA flux (Mq, 67 nm-10 um): {0:1.3e}'.format(np.trapz(flux['Mq'], x=wav['Mq'])))

#########################
####Import VPL data (AD Leo)
#########################
##AD Leo (quiescence) [Segura et al 2005] NOTE: May lack LyA peak (see comments in file)
#adleo_wav_um, adleo_flux_dist_units=np.genfromtxt('./Raw_Stellar_Spectra/Steady-State/VPL/adleo_dat.txt', skip_header=175, skip_footer=1,usecols=(0,1), unpack=True) #um, Watt/cm2/um; fluxes are at Earth-star distance
#
#adleo_wav=adleo_wav_um*1.e3 #convert um to nm
#adleo_flux_dist = adleo_flux_dist_units*10. #Convert Watt/cm2/um to Watt/m2/nm; 
#total_flux=np.trapz(adleo_flux_dist, x=adleo_wav) #"stellar constant" at Earth-star separation
#
#adleo_flux=adleo_flux_dist*solar_constant/total_flux
#
#wav['adleo']= adleo_wav
#flux['adleo']=adleo_flux
#
####The AD Leo data feature many negative fluxes, which are artefacts of reduction; they must be binned. Further, the data are unevenly spaced; this is not intrinsically a problem, but certainly is untidy. Let's put them instead onto a grid of 1 nm spacing, approximately (but not exactly) matching Hu. 
#wav_left=np.arange(117.0, 10000.0, step=1.0) #nm
#wav_right=np.arange(118.0, 10001.0, step=1.0) #nm
#adleo_wav_rebinned=0.5*(wav_left+wav_right)
#
#adleo_flux_rebinned=integrate_data(adleo_wav, adleo_flux, wav_left, wav_right)
####This still leaves a couple of negative data points. Let's deal with them in the following (kludgy) way: the 2 negative data points are surrounded by positive data points. Add the negative data to the neighboring bins; set the negative data to zero. In this way total flux is conserved but negative data are avoided.
#
##This is really trivial but good practice. 
#neg_data_inds=np.where(adleo_flux_rebinned<0)[0]
#for ind in range(0, len(neg_data_inds)):
#    neg_data_ind=neg_data_inds[ind]
#    adleo_flux_rebinned[neg_data_ind+1] += adleo_flux_rebinned[neg_data_ind]
#    adleo_flux_rebinned[neg_data_ind] = 0
#
#
#wav_binned['adleo']= adleo_wav_rebinned
#flux_binned['adleo']=adleo_flux_rebinned
#
#print('TOA flux (Ad Leo Quiescent, 117-10 um): {0:1.3e}'.format(np.trapz(adleo_flux_rebinned, x=adleo_wav_rebinned))) #Looks to be a hair low, possibly due to omission of some of the longward flux?
#
#
#
#########################
####Import MUSCLES Data
#########################
#muscles_starnames=np.array(['gj1214','gj876','gj436','gj581','gj667c','gj176','gj832','hd85512','hd40307','hd97658','v-eps-eri', 'gj551']) #NOTE: Prox Cen=GJ 551. The HD stars and the epsilon eridani are K-dwarfs (not M-dwarfs)
#
#for muscles_starname in muscles_starnames:
#    filename='./Raw_Stellar_Spectra/Steady-State/MUSCLES/hlsp_muscles_multi_multi_'+muscles_starname+'_broadband_v22_adapt-const-res-sed.fits'
#    spec= fits.getdata(filename,1)
#    spec_wav=spec['WAVELENGTH'] #wavelength scale in A
#    spec_flux=spec['FLUX'] #fluxes in erg/s/cm2/A
#
#    header=fits.getheader(filename,0)
#    boloflux=header['BOLOFLUX'] #units of ergs/s/cm2 
#
#    boloflux_units=boloflux*0.001 #convert from units of erg/s/cm2 to units of Watts m^-2
#    
#    star_wav=spec_wav*0.1 #Convert from A to nm
#    spec_flux_units=spec_flux*0.01 #Convert from erg/s/cm2/A to Watts m^-2 nm^-1
#
##    checkboloflux=np.trapz(spec_flux_units, x=star_wav) #integrate the stellar spectrum to make sure the boloflux is what it says it is
##    print('For {0}, BOLOFLUX is {1:1.3e} whereas the integrated flux is {2:1.3e}'.format(muscles_starname, boloflux_units, checkboloflux))
##    ##checkboloflux will not match boloflux_units because the spectra terminate at 5.5 microns, while a Rayleigh-Jeans tail continues. This is included in BOLOFLUX and accounts for the difference. 
#   
#    star_flux=spec_flux_units*solar_constant/boloflux_units #Scale to match solar instellation
#    
#    wav[muscles_starname]=star_wav #nm
#    flux[muscles_starname]=star_flux #W m^-2 nm^-1
#    
####Integrate to 1 nm resolution for inclusion in Hu code. Can easily be something else since the MUSCLES data are sanitized and binned such that they lack negative fluxes.
#wav_left=np.arange(100.0, 5000.0, step=1.0) #nm
#wav_right=np.arange(101.0, 5001.0, step=1.0) #nm
#wav_rebinned=0.5*(wav_left+wav_right)
#
#for muscles_starname in muscles_starnames:
#    star_wav=wav[muscles_starname]
#    star_flux=flux[muscles_starname]
#    
#    wav_binned[muscles_starname]=wav_rebinned
#    flux_binned[muscles_starname]=integrate_data(star_wav, star_flux, wav_left, wav_right)
#    
#    print('TOA flux ({0}, 100 nm-5 um): {1:1.3e}'.format(muscles_starname, np.trapz(flux_binned[muscles_starname], x=wav_binned[muscles_starname]))) #Looks to be a hair low, possibly due to omission of some of the longward flux?

	

########################
###Import Peacock+2019 TRAPPIST-1 Models
########################
####Currently using the version from Eddie, which averages Peacock+2019 models 1, 2A, 2B. In future, probably best to use model 1 which matches HST reconstruction of LyA flux and seems to be a "best guess" model. Have emailed S. Peacock to request. 
####Use of this requires appropriate credits to Peacock, Lincowski, and Schwieterman.
#    
#trappist_wav_um, trappist_flux_units_unscaled=np.genfromtxt('./Raw_Stellar_Spectra/Steady-State/Peacock/TRAPPIST-1.dat', skip_header=5, skip_footer=0,usecols=(0,1), unpack=True) #um, Watt/m2/um; fluxes are at 1 AU, need to be scaled.
#
#trappist_wav=trappist_wav_um*1.0e3 #convert um to nm
#
#trappist_flux_unscaled=trappist_flux_units_unscaled*1.0e-3 #conver from W/m^2/um to W/m^2/nm
#
#trappist_flux_bolo=np.trapz(trappist_flux_unscaled, x=trappist_wav)#bolometric flux. Goes out to 95 microns so probably has everything!
#trappist_flux=trappist_flux_unscaled*solar_constant/trappist_flux_bolo
#
#wav['trappist-1']=trappist_wav
#flux['trappist-1']=trappist_flux
#
####Bin down to usable resolution
#wav_left=np.arange(100.0, 5000.0, step=1.0) #nm
#wav_right=np.arange(101.0, 5001.0, step=1.0) #nm
#wav_rebinned=0.5*(wav_left+wav_right)
#trappist_flux_rebinned=integrate_data(trappist_wav, trappist_flux, wav_left, wav_right)
#
#wav_binned['trappist-1']=wav_rebinned
#flux_binned['trappist-1']=trappist_flux_rebinned
#
#print('TOA flux (TRAPPIST-1, 100 nm-5 um): {0:1.3e}'.format(np.trapz(flux_binned[trappist-1], x=wav_binned[trappist-1]))) #Looks to be a hair low, possibly due to omission of some of the longward flux?
    
    
###Model 1 from Peacock+2019, synthetic spectrum of TRAPPIST-1 that is "best guess" for its spectrum (2A, 2B are bracketing cases)
#Use requires credit to Peacock+2019, and thanks to Peacock.
    
trappist_wav_A, trappist_flux_units_unscaled=np.genfromtxt('./Raw_Stellar_Spectra/Steady-State/Peacock/TRAPPIST1_1A_fullres.txt', skip_header=1, skip_footer=0,usecols=(0,1), unpack=True) #A, erg/s/cm2/A; fluxes are at stellar surface, need to be scaled.

trappist_wav=trappist_wav_A*0.1 #convert A to nm
trappist_flux_unscaled=trappist_flux_units_unscaled*0.01 #conver from erg/s/cm2/A to W/m^2/nm

trappist_flux_bolo=np.trapz(trappist_flux_unscaled, x=trappist_wav)#bolometric flux. Goes out to 95 microns so probably has everything!
trappist_flux=trappist_flux_unscaled*solar_constant/trappist_flux_bolo

wav['trappist-1']=trappist_wav
flux['trappist-1']=trappist_flux

###Bin down to usable resolution
wav_left=np.arange(67.0, 5000.0, step=1.0) #nm
wav_right=np.arange(68.0, 5001.0, step=1.0) #nm
wav_rebinned=0.5*(wav_left+wav_right)
trappist_flux_rebinned=integrate_data(trappist_wav, trappist_flux, wav_left, wav_right)

wav_binned['trappist-1']=wav_rebinned
flux_binned['trappist-1']=trappist_flux_rebinned

print('TOA flux (TRAPPIST-1, 67 nm-100 um): {0:1.3e}'.format(np.trapz(flux['trappist-1'], x=wav['trappist-1']))) 

print('TOA flux (TRAPPIST-1, 67 nm-5 um): {0:1.3e}'.format(np.trapz(flux_binned['trappist-1'], x=wav_binned['trappist-1']))) 


# ###################
# ###
# ###################

# ###New TRAPPIST-1 spectrum used by Eddie. To see if we can reproduce their findings. 
# trappist_new_wav, trappist_new_flux_units_unscaled=np.genfromtxt('./atmos-inputs/trappist-1_AVG_NEW.txt', skip_header=1, skip_footer=0,usecols=(0,1), unpack=True) #units: nm, mW m**-2 nm**-1. *should* be at 1-AU equivalent. 

# trappist_new_flux_unscaled=trappist_new_flux_units_unscaled*1.e-3 #convert from mW m**-2 nm**-1to W/m^2/nm

# # trappist_new_flux_bolo=np.trapz(trappist_new_flux_unscaled, x=trappist_new_wav)#bolometric flux. Goes out to 95 microns so probably has everything!

# # trappist_new_flux=trappist_new_flux_unscaled*solar_constant/trappist_flux_bolo
# trappist_new_flux=trappist_new_flux_unscaled
# wav['trappist-1_new']=trappist_new_wav
# flux['trappist-1_new']=trappist_new_flux

# ###Bin down to usable resolution
# wav_new_left=np.arange(67.0, 2000.0, step=1.0) #nm
# wav_new_right=np.arange(68.0, 2001.0, step=1.0) #nm
# wav_new_rebinned=0.5*(wav_new_left+wav_new_right)
# trappist_new_flux_rebinned=integrate_data(trappist_new_wav, trappist_new_flux, wav_new_left, wav_new_right)

# wav_binned['trappist-1_new']=wav_new_rebinned
# flux_binned['trappist-1_new']=trappist_new_flux_rebinned

# print('TOA flux (TRAPPIST-1, New, 1 nm-2.4 um): {0:1.3e}'.format(np.trapz(flux['trappist-1_new'], x=wav['trappist-1_new']))) 

# print('TOA flux (TRAPPIST-1, 67 nm-2 um): {0:1.3e}'.format(np.trapz(flux_binned['trappist-1_new'], x=wav_binned['trappist-1_new']))) 
#########################################################################
#########################################################################
#########################################################################
####Outputs: plots & input files.
#########################################################################
#########################################################################
#########################################################################
#dataset_list=np.array(['v-eps-eri', 'hd97658','hd40307','hd85512','gj832','gj176','adleo','gj667c','gj581','gj436','gj876','gj551', 'gj1214', 'trappist-1'])

# dataset_list=np.array(['trappist-1_new'])
dataset_list=np.array(['trappist-1'])


########################
###Plot binned vs. unbinned data, to make sure it worked right
########################    
for dataset in dataset_list:
    fig, ax=plt.subplots(1, figsize=(8, 6))
    ax.plot(wav[dataset], flux[dataset], linewidth=2, linestyle='-', color='black', label='Unbinned')
    ax.plot(wav_binned[dataset], flux_binned[dataset], linewidth=2, linestyle='--', color='red', label='Binned')
    ax.set_yscale('log')
    ax.set_ylim([1.0e-6, 1.0e+1])
    ax.set_xscale('log')
    ax.set_xlim([67., 1000.])
    ax.legend(loc='best', fontsize=12, borderaxespad=0)
    plt.tight_layout()
    plt.savefig('./Plots/'+dataset+'.pdf', orientation='portrait',papertype='letter', format='pdf')

########################
###Print binned data to Hu code compatible file
######################## 

for dataset in dataset_list:
    wav=wav_binned[dataset]
    flux=flux_binned[dataset]
    composite_toprint=np.zeros((len(wav), 2))
    composite_toprint[:,0]=wav #wavelengths in nm
    composite_toprint[:,1]=flux #TOA flux W m**-2 nm**-1
    
    np.savetxt('./photochem_outputs/'+dataset+'.txt', composite_toprint, fmt='%5.6f\t%1.6e', newline='\n')

