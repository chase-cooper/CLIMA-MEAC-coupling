"""
Generate int.rates.out.dat from colintrxnrates.dat

Version 4. Advances over v3 by taking a single directory as input, not the individual filenames.

Authors: Sukrit Ranjan, Zhuchang Zhan,  Sara Seager
"""

#Unit conversions
km2cm=1.e5 #1 km in cm
cm2km=1.e-5 #1 cm in km

import re
import sys
import numpy as np
import pandas as pd
from collections import Counter
import pdb

def return_colint_reaction_rates(concSTD, ChemReac, name):
    """
    Takes: concentrationSTD.dat and ChemicalRate.dat file from the Hu code
    From concentrationSTD.dat file extracts z scale. From ChemicalRate file takes reaction rate (cm**-3 s**-1).
    Prints: reaction ID, column integrated reaction rate for each reaction ID (cm**-2 s**-1)
    """
    
    ###Get z scale
    tp=np.genfromtxt(concSTD, skip_header=2, skip_footer=0, unpack=False) #Import simulation TPZ profile
    deltazs=(tp[:,2]-tp[:,1]) #width of altitude bins, km    

    ###Get Chemical Reaction Rates
    chemlabels=np.genfromtxt(ChemReac, skip_header=1, skip_footer=0, unpack=True, usecols=(0), dtype=None, encoding='UTF-8') #Labels of all filenames
    chemrates=np.genfromtxt(ChemReac, skip_header=1, skip_footer=0, unpack=False, dtype=float)[:,1:]
    
    ###Extract relevant chemical reaction rates.    
    colint_chemrates=np.zeros(len(chemlabels), dtype=float) #
    # chemlabels_decoded=np.zeros(len(chemlabels), dtype='U4')
    for ind in range(0, len(chemlabels)): #Step through each of the chemical reactions in the FILE
        # chemlabels_decoded[ind]=chemlabels[ind].decode('UTF-8') #decode the indth chemical reaction 
        colint_chemrates[ind]=np.sum(chemrates[ind,:]*deltazs*km2cm) #cm**-2 s**-1

    toprint=np.zeros(len(chemlabels), dtype=[('col1', 'U4'),('col2', float)])    
    toprint['col1']=chemlabels
    toprint['col2']=colint_chemrates
    np.savetxt(name+'colintrxnrates.dat', toprint, fmt='%s\t%3.6e', newline='\n', header='Reaction ID \t Column-Integrated Reaction Rate (cm**-2 s**-1)')

def read_rate_input(input_file):
    
    with open(input_file) as f:
        data = {}
        lines = [x.split("\t") for x in f.read().split("\n")[1:]]
        for line in lines:
            try:
                key,value = line
            except:
                continue
            data[key] = value

    return data

def read_reaction_r(reaction_rate_data,input_reaction_file,species_file):

    item = "R"
    mapping = pd.read_excel(species_file)
    reaction_data = np.genfromtxt(input_reaction_file)
    
    def maps(item):
        return "" if item == 0 else mapping[mapping["No."]==item]["Species Name"].values[0]
    
    molecules = []
    reactions = []
    for i,line in enumerate(reaction_data):
        r1m, r2m, ctm, p1m, p2m, p3m = [maps(x) for x in line]
        r1, r2, ct, p1, p2, p3 = line
        
        try:
            reaction_number = "%s%s"%(item,i+1)
            reaction_rate = reaction_rate_data[reaction_number]
            reactions.append([reaction_number,r1m, r2m, p1m, p2m, p3m,reaction_rate])
            molecules.extend([r1,r2,p1,p2,p3])
        except KeyError: # not all reaction is used in the input reaction_rate_data
            continue


    all_molecules = [maps(x) for x in list(set(molecules))]
    
    
    
    if "" in all_molecules:
        all_molecules.remove("")
        
    return all_molecules,reactions

def read_reaction_m(reaction_rate_data,input_reaction_file,species_file):

    item = "M"
    mapping = pd.read_excel(species_file)
    reaction_data = np.genfromtxt(input_reaction_file)
    
    def maps(item):
        return "" if item == 0 else mapping[mapping["No."]==item]["Species Name"].values[0]
    
    molecules = []
    reactions = []
    for i,line in enumerate(reaction_data):
        r1m, r2m, p1m, p2m= [maps(x) for x in line]
        r1, r2, p1, p2 = line
        
        try:
            reaction_number = "%s%s"%(item,i+1)
            reaction_rate = reaction_rate_data[reaction_number]
            reactions.append([reaction_number,r1m, r2m, p1m, p2m,reaction_rate])
            molecules.extend([r1,r2,p1,p2])
        except KeyError: # not all reaction is used in the input reaction_rate_data
            continue


    all_molecules = [maps(x) for x in list(set(molecules))]
    
    
    
    if "" in all_molecules:
        all_molecules.remove("")
        
    return all_molecules,reactions

def read_reaction_t(reaction_rate_data,input_reaction_file,species_file):

    item = "T"
    mapping = pd.read_excel(species_file)
    reaction_data = np.genfromtxt(input_reaction_file)
    
    def maps(item):
        return "" if item == 0 else mapping[mapping["No."]==item]["Species Name"].values[0]
    
    molecules = []
    reactions = []
    for i,line in enumerate(reaction_data):
        r1m, p1m, p2m= [maps(x) for x in line]
        r1, p1, p2 = line
        
        try:
            reaction_number = "%s%s"%(item,i+1)
            reaction_rate = reaction_rate_data[reaction_number]
            reactions.append([reaction_number,r1m, p1m, p2m,reaction_rate])
            molecules.extend([r1,p1,p2])
        except KeyError: # not all reaction is used in the input reaction_rate_data
            continue


    all_molecules = [maps(x) for x in list(set(molecules))]
    
    
    
    if "" in all_molecules:
        all_molecules.remove("")
        
    return all_molecules,reactions

def read_reaction_p(reaction_rate_data,input_reaction_file,species_file):

    item = "P"
    mapping = pd.read_excel(species_file)
    reaction_data = np.genfromtxt(input_reaction_file)
    
    def maps(item):
        try:
            return mapping[mapping["No."]==item]["Species Name"].values[0]
        except:
            return "" 

    
    molecules = []
    reactions = []
    for i,line in enumerate(reaction_data):
        r1m, p1m, p2m, p3m, *ct = [maps(x) for x in line]
        r1, p1, p2, p3, *ct = line
        
        try:
            reaction_number = "%s%s"%(item,i+1)
            reaction_rate = reaction_rate_data[reaction_number]
            reactions.append([reaction_number,r1m, p1m, p2m, p3m,reaction_rate])
            molecules.extend([r1,p1,p2,p3])
        except KeyError: # not all reaction is used in the input reaction_rate_data
            continue


    all_molecules = [maps(x) for x in list(set(molecules))]
    
    
    if "" in all_molecules:
        all_molecules.remove("")
        
    return all_molecules,reactions

def read_reactions(reaction_rate_data,input_reaction_files,species_file):
    
    
    # load reaction_r and return list of molecules and list of reaction entries
    molecules_r, reaction_r = read_reaction_r(reaction_rate_data,input_reaction_files["R"],species_file)
    molecules_m, reaction_m = read_reaction_m(reaction_rate_data,input_reaction_files["M"],species_file)
    molecules_t, reaction_t = read_reaction_t(reaction_rate_data,input_reaction_files["T"],species_file)
    molecules_p, reaction_p = read_reaction_p(reaction_rate_data,input_reaction_files["P"],species_file)
    
    
    molecules = list(set(molecules_r+molecules_m+molecules_t+molecules_p))
    reactions = reaction_r+reaction_m+reaction_t+reaction_p
    
    return molecules,reactions
    
def compute_molecule_info_rank(molecules,reactions):

    molecule_data = {}
    molecule_rank = []
    
    for blm,molecule in enumerate(molecules):
        if molecule == "":
            print("null molecule")
            sys.exit()
        
        loss_reaction = []
        """
        loss_max = -1
        loss_max_tag = ""
        loss_max_value = ""
        """
        total_loss = 0
        for reaction in reactions:
            
            if "R" in reaction[0]:
                reaction_number,r1, r2, p1, p2, p3,reaction_rate = reaction
            elif "M" in reaction[0]:
                reaction_number,r1, r2, p1, p2,reaction_rate = reaction
                p3 = ""
            elif "T" in reaction[0]:
                reaction_number,r1, p1, p2, reaction_rate = reaction
                r2 = ""
                p3 = ""
            elif "P"in reaction[0]:
                reaction_number,r1, p1, p2, p3,reaction_rate = reaction
                r2 = ""
        
            
            if molecule == r1 or molecule == r2:
                
                
                
                loss_reaction.append(reaction)
                reaction_rate = float(reaction[-1])
                
                if molecule == r1 and molecule == r2:
                    total_loss += 2*reaction_rate
                else:
                    total_loss += reaction_rate
                """
                if reaction_rate > loss_max:
                    loss_max = reaction_rate
                    loss_max_tag = reaction[0]
                    loss_max_value = reaction[-1]
                """
        
        production_reaction = []
        """
        production_max = -1
        production_max_value = ""
        production_max_tag = ""
        """
        total_production = 0
        for reaction in reactions:
            
            if "R" in reaction[0]:
                reaction_number,r1, r2, p1, p2, p3,reaction_rate = reaction
            elif "M" in reaction[0]:
                reaction_number,r1, r2, p1, p2,reaction_rate = reaction
                p3 = ""
            elif "T" in reaction[0]:
                reaction_number,r1, p1, p2, reaction_rate = reaction
                r2 = ""
                p3 = ""
            elif "P"in reaction[0]:
                reaction_number,r1, p1, p2, p3,reaction_rate = reaction
                r2 = ""
        
            
            if molecule == p1 or molecule == p2 or molecule == p3:
                production_reaction.append(reaction)
                reaction_rate = float(reaction[-1])
                
                if (molecule == p1 and molecule == p2 and molecule == p3):
                    sys.exit()
                elif (molecule == p1 and molecule == p2) or (molecule == p1 and molecule == p3) or (molecule == p2 and molecule == p3):
                    
                    total_production += 2*reaction_rate
                else:
                    total_production += reaction_rate
                
                
                
                
                
                
                """
                if reaction_rate > production_max:
                    production_max = reaction_rate
                    production_max_tag = reaction[0]
                    production_max_value = reaction[-1]
                """
        
        production_header = "%s\t\t\t\tPRODUCTION RXS\t\tINT RX RATE\tTP = %.2e"%(molecule,total_production)
        loss_header = "%s\t\t\t\tLOSS RXS\t\tINT RX RATE\tTL = %.2e"%(molecule,total_loss)
        
        
        molecule_rank.append([molecule,total_production,total_loss])  

        molecule_data[molecule] = [production_header,
                                   production_reaction,
                                   loss_header,
                                   loss_reaction]
    
    return molecule_rank,molecule_data

def print_output(molecule,reaction,type):

    if "R" in reaction[0]:
        reaction_number,r1, r2, p1, p2, p3,reaction_rate = reaction
    elif "M" in reaction[0]:
        reaction_number,r1, r2, p1, p2,reaction_rate = reaction
        p3 = ""
    elif "T" in reaction[0]:
        reaction_number,r1, p1, p2, reaction_rate = reaction
        r2 = ""
        p3 = ""
    elif "P"in reaction[0]:
        reaction_number,r1, p1, p2, p3,reaction_rate = reaction
        r2 = ""
        
    if r2 == "":
        reactant = "%s\t\t"%(r1)
        

    else:
        
        if len(r2) > 5:
            reactant = "%s\t+ %s"%(r1,r2)
        else: 
            reactant = "%s\t+ %s\t"%(r1,r2)
            

        
    if p3 == "":
        if p2 == "":
            product  = "%s\t\t"%(p1)
        else:
            product  = "%s\t+ %s\t"%(p1,p2)
    else:
        product  = "%s\t+ %s\t+ %s"%(p1,p2,p3)


    double = False
    
    if type == "loss":
        if r1 == r2 and r1 == molecule:
            double= True    
    elif type == "production":
        if (p1 == p2 and p1 == molecule) or (p2 == p3 and p2 == molecule) or (p1 == p3 and p3 == molecule):
            double = True
        
    if double:
        print("  %s\t%s\t= %s\t%s (x2)"%(reaction_number,reactant,product,reaction_rate))
    else:
        print("  %s\t%s\t= %s\t%s"%(reaction_number,reactant,product,reaction_rate))
        
        
        
        
        
        

def generate_int_rates_out(name):
    
    # toggle output to file or print
    sys.stdout = open(name+'int.rates.out3.dat', 'w')
    
    # input files
    species_file         = "hu-code-sr/Documentation/SpeciesName.xlsx"
    input_rate_file      = (name+"colintrxnrates.dat")
    input_reaction_files  = {"R":"hu-code-sr/Data/Reaction_R.txt",
                             "M":"hu-code-sr/Data/Reaction_M.txt",
                             "T":"hu-code-sr/Data/Reaction_T.txt",
                             "P":"hu-code-sr/Data/Reaction_P.txt"}
    
    # load input reaction rates from colintrxnrates.dat into a python dictionary
    reaction_rate_data   = read_rate_input(input_rate_file)
    
    # load reaction_r and return list of molecules and list of reaction entries
    molecules, reactions = read_reactions(reaction_rate_data,input_reaction_files,species_file)

    # compute the maximum reaction rate and sort molecule based on production rate
    molecule_rank,molecule_data = compute_molecule_info_rank(molecules, reactions)
    
    # switch x:x[1] to x:x[2] if you want to sort by loss rate
    # reverse = True is descending    
    sorted_molecule_rank = sorted(molecule_rank,key=lambda x:x[1],reverse=True) 
    sorted_molecule = [x[0] for x in sorted_molecule_rank]
    
    # output text to int.rates.out.dat
    print("Total_Molecules: %s"%len(molecules))
    print("Total_Reactions: %s"%len(reactions))
    print("Sorting by     : TP Descending")
    
    for molecule in sorted_molecule:
        production_header,production_reaction,loss_header,loss_reaction = molecule_data[molecule]
        
        print("*"*60)
        print()
        print(production_header)
        print()
        production_reaction.sort(key= lambda x:float(x[-1]),reverse=True)
        for reaction in production_reaction:
            print_output(molecule,reaction,"production")
            
        print()
        print(loss_header)
        print()
        loss_reaction.sort(key= lambda x:float(x[-1]),reverse=True)
        for reaction in loss_reaction:
            print_output(molecule,reaction,"loss")
        print()
        
def analyze_colint_rxn_rates(name):
    """
    Takes:
        --string from which it constructs addresses for a concentrationSTD.dat and ChemicalRate.dat file for a simulation run from the Hu code, including filepath.
    Returns: 
        --file giving column-integrated reaction rates for that simulation run (colintrxnrates.dat file)
        --File giving ordered production and loss rates for each reaction from that simulation run (int.rates.out3.dat file)
    """
    
    concSTD=name+'ConcentrationSTD.dat'
    ChemReac=name+'ChemicalRate.dat'
    
    ###Step 1: Generate column-integrated reaction rate file. 
    return_colint_reaction_rates(concSTD, ChemReac, name)
    
    ###Step 2: Using newly-generated column-integrated reaction rate file, generate prod/loss file. 
    generate_int_rates_out(name)
    


# return_colint_reaction_rates('hu-code-sr/scenario_library/CO2_CH4/fco2_1e-1/ConcentrationSTD.dat','hu-code-sr/scenario_library/CO2_CH4/fco2_1e-1/ChemicalRate.dat','outputs/')
generate_int_rates_out('hu-code-sr/scenario_library/CO2_CH4/fco2_1e-1/')
