#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 16:36:13 2019

@author: sukrit

The purpose of this script is to calculate redox balance in the atmosphere and in the atmosphere-ocean system, following the treatment of Harman+2015.

The objectives of these calculations are twofold:
    --Redox balance in the atmosphere: this should be ~0, in order to conserve mass.
    --Redox balance in the atmosphere-ocean system ("global redox balance"): Under the assumption that in the abiotic case there is no burial of reductants or oxidants (oxidative weathering), then there should be redox balance in the atmosphere-ocean system. This needs to be checked as well.
    
How to check redox balance for atmospheres with fixed mixing ratio boundary conditions for molecules with non-zero redox state (e.g., O2, H2)
Such species don't have rainout or emission (or they shouldn't...) but they DO have escape and they DO influence atmospheric composition.
"""
###########################
###Import useful libraries
###########################
import matplotlib as plt
import numpy as np
import pandas as pd
import gspread

###########################
###Define useful constants
###########################


###########################
###Step 1: import mapping between species name, species standard number, and the "redox value" of each species.
###########################
df=pd.read_excel('SpeciesNameRedox.xlsx', sheet_name='Sheet1') #The nth entry in this sheet corresponds to species with std. no n+1 (zero-indexing). So, MUST BE IN ORDER

##Check the above statement:
#print(0==np.max(df['No.']-(np.arange(0, len(df['No.']))+1))) #should be 0, so should output "True". Can probably be commented out eventually once we trust it.

###########################
###Step 2: Import global balance file
###########################
globalbalance_table=np.genfromtxt('./Result/Sun/N2/GlobalBalance.dat', skip_header=1, skip_footer=0, names='Std_No, Outgassing, ChemProd, ChemLoss, DryDepo, WetDepo, Escape, Condensation, Net_Flux	, Global_Change', delimiter='\t') #genfromtxt does not like spaces or periods in column names. # ./GlobalBalance-NoE_fromRenyu.dat ##

type_1_species_stdnos=np.array([52, 55]) #N2, CO2

###########################
###Step 3: Calculate atmospheric redox balance
##Big picture: track whether electrons are being conserved well in the atmosphere, by converting to H equivalents. 
##Draw a box around the atmosphere and check whether the  net flux of electrons through the system is 0 (to reasonable numerical precision)
#This amounts to checking how well our model satisfies the equation 
#phi_outgassing(reductants) + phi_rainout(oxidants) = phi_esc(reductants)+phi_rainout(reductants) (Catling & Kasting 2017 eqn 8.34; Harman et al eqn 1). However, I'm not sure the signs work out here -- if one includes the sign on the redox coefficient, 
#Note that this formalism neglects outgassing and escape of oxidants. 
###########################
num_species=len(globalbalance_table['Std_No'])
redox_weighted_Outgassing=np.zeros(num_species)
redox_weighted_DryDepo=np.zeros(num_species)
redox_weighted_WetDepo=np.zeros(num_species)
redox_weighted_Escape=np.zeros(num_species)

#Loop over all lines in global balance table.
for ind in range(0, num_species):
    std_no=globalbalance_table['Std_No'][ind] #Standard number in Hu code identifying the species whose information is being read.
    
    redox_coefficient=df['Redox Coefficient (H, Hu)'][std_no-1] #redox coefficient corresponding to species. signed value: + for reductant, -for oxidant.
    
#    #Check to make sure we have gotten the right redox coefficient. Can probably be commented out eventually once we trust it.
#    print(0==(df['No.'][std_no-1] - std_no)) #If true, we have accessed the correct redox coefficient. If false, we have accessed the incorrect redox coefficient
    
    #Outgassing: positive in GlobalBalance file because SUPPLIED to atmosphere
    if (std_no in type_1_species_stdnos): #have to exclude boundary condition type 1 (fixed MR at base). Is this the right way to treat???
        redox_weighted_Outgassing[ind]=0
    else:
        redox_weighted_Outgassing[ind]=globalbalance_table['Outgassing'][ind] * redox_coefficient
        
    #Wet deposition, dry deposition, escape: negative because REMOVED from atmosphere
    redox_weighted_DryDepo[ind]=globalbalance_table['DryDepo'][ind] * redox_coefficient
    redox_weighted_WetDepo[ind]=globalbalance_table['WetDepo'][ind] * redox_coefficient
    redox_weighted_Escape[ind]=globalbalance_table['Escape'][ind] * redox_coefficient
    
total_redox_weighted_Outgassing=np.sum(redox_weighted_Outgassing)
total_redox_weighted_DryDepo=np.sum(redox_weighted_DryDepo)
total_redox_weighted_WetDepo=np.sum(redox_weighted_WetDepo)
total_redox_weighted_Escape=np.sum(redox_weighted_Escape)

total_redox_change=total_redox_weighted_Outgassing+total_redox_weighted_DryDepo+total_redox_weighted_WetDepo+total_redox_weighted_Escape #Enters atmosphere via outgassing; exist via wet, dry deposition, escape.

print('Atmospheric redox balance: {0:1.1e}'.format(total_redox_change))

###########################
###Step 4: Calculate atmosphere-ocean redox balance
###########################

#print('Outgassing and escape redox balance: {0:1.1e}'.format(total_redox_weighted_Outgassing+total_redox_weighted_Escape)) #this is the parameter adopted by Harman, etc. If nonzero, global redox balance requires a return flux. 

print('Global redox balance (AKA oceanic redox balance): {0:1.1e}'.format(total_redox_weighted_WetDepo+total_redox_weighted_DryDepo)) #this is the parameter adopted by Harman, etc. Should be 0.

###########################
###Step 4: Calculate atmosphere-ocean redox balance
###########################
