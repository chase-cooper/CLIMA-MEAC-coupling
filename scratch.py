import matplotlib.pyplot as plt
import numpy as np

from generate_photochemistry_rate_per_molecule_v4 import generate_int_rates_out
from compare_model_outputs_v2 import *
from parameters import *


def printMEACztpFromCLAST():
    f = open('cloudy_clima/CLIMA/IO/clima_last.tab','r')
    data = f.readlines()
    f.close()

    alt,t,p = [],[],[]
    for line in data[1:]:
        line = line.split()
        alt.append(float(line[0]))
        p.append(np.log10(float(line[1])*101325))
        t.append(float(line[2]))
    alt = alt[::-1]
    p = p[::-1]
    t = t[::-1]

    zgap = alt[-1]-alt[-2]
    pgap = p[-1]-p[-2]

    z = alt[-1]
    while z < 100:
        z += zgap
        alt.append(z)
        p.append(p[-1]+pgap)
        t.append(t[-1])

    for i in range(len(p)):
        aa = np.round(alt[i],decimals=6)
        pp = np.round(p[i],decimals=6)
        tt = np.round(t[i],decimals=6)
        print(f"{aa}    {pp}      {tt}")

def compMEACandPhotochemCompositions(meacpath,pcpath):
    ################
    ###   MEAC   ###
    ################

    #Indexes, corrected for Hu 1-indexing vs Python 0-indexing
    ind_o=1-1 #O
    ind_h=3-1 #H
    ind_cho2=63-1
    ind_ch3=60-1
    ind_ch2=59-1
    ind_ch21=80-1
    ind_ch=58-1
    ind_cho=61-1
    ind_c2h4=29-1
    ind_so2=43-1
    ind_so=42-1
    ind_h2s=45-1
    ind_no=12-1
    ind_n2o=11-1
    ind_s8=79-1
    ind_s8a=111-1
    ind_ch4o=24-1
    ind_c2h2=27-1
    ind_ocs=49-1
    ind_o1d=56-1
    ind_c2h6=31-1
    ind_h2o2=6-1
    ind_h2so4=73-1
    ind_h2so4a=78-1
    ind_ch2o=22-1
    ind_ho2=5-1

    # MEAC species indices, not including 5-column offset due to other atmospheric parameters
    m_ind_ch4=21-1
    m_ind_h2=53-1
    m_ind_h2o=7-1
    m_ind_co2=52-1
    m_ind_n2=55-1
    m_ind_co=20-1
    m_ind_o2=54-1
    m_ind_o3=2-1
    m_ind_oh=4-1

    # Data file
    m_raw = np.loadtxt(meacpath+'ConcentrationSTD.dat',skiprows=2)
    m_species = m_raw[:,5:]
    m_p = m_raw[:,4]
    m_z = m_raw[:,0]
    m_t = m_raw[:,3]
    total_counts = np.sum(m_raw[:,5:],axis=1)

    # Number densities
    n_n2    = m_species[:,m_ind_n2]
    n_o2    = m_species[:,m_ind_o2]
    n_o3    = m_species[:,m_ind_o3]
    n_h2    = m_species[:,m_ind_h2]
    n_h2o   = m_species[:,m_ind_h2o]
    n_ch4   = m_species[:,m_ind_ch4]
    n_co    = m_species[:,m_ind_co]
    n_co2   = m_species[:,m_ind_co2]
    n_oh    = m_species[:,m_ind_oh]
    
    # Mixing ratios
    m_mr_n2     = np.divide(n_n2,total_counts)
    m_mr_o2     = np.divide(n_o2,total_counts)
    m_mr_o3     = np.divide(n_o3,total_counts)
    m_mr_h2     = np.divide(n_h2,total_counts)
    m_mr_h2o    = np.divide(n_h2o,total_counts)
    m_mr_ch4    = np.divide(n_ch4,total_counts)
    m_mr_co     = np.divide(n_co,total_counts)
    m_mr_co2    = np.divide(n_co2,total_counts)
    m_mr_oh     = np.divide(n_oh,total_counts)

    #################
    ### Photochem ###
    #################

    # Indexes
    p_ind_n2    = 47 -1
    p_ind_o2    = 22-1
    p_ind_o3    = 61-1
    p_ind_h2    = 18-1
    p_ind_h2o   = 19-1
    p_ind_ch4   = 31-1
    p_ind_co    = 23-1
    p_ind_co2   = 24-1
    p_ind_oh    = 20-1

    # Data file
    p_raw = np.loadtxt(pcpath+'atmosphere.txt',skiprows=1)
    p_species = p_raw[:,5:]
    p_p = p_raw[:,1]*1e5
    p_t = p_raw[:,3]
    p_z = p_raw[:,0]

    # Mixing ratios (photochem documents species concentrations in terms of mixing ratios)
    p_mr_n2     = p_species[:,p_ind_n2]
    p_mr_o2     = p_species[:,p_ind_o2]
    p_mr_o3     = p_species[:,p_ind_o3]
    p_mr_h2     = p_species[:,p_ind_h2]
    p_mr_h2o    = p_species[:,p_ind_h2o]
    p_mr_ch4    = p_species[:,p_ind_ch4]
    p_mr_co     = p_species[:,p_ind_co]
    p_mr_co2    = p_species[:,p_ind_co2]
    p_mr_oh     = p_species[:,p_ind_oh]

    fig,ax = plt.subplots(ncols=2,width_ratios=[1,3])
    m_ls,p_ls = ['-',':']

    # TP profiles
    ax[0].plot(m_t,m_p,ls=m_ls,c='black',label="MEAC")
    ax[0].plot(p_t,p_p,ls=p_ls,c='black',label="Photochem")
    ax[0].legend()
    ax[0].invert_yaxis()
    ax[0].set_yscale('log')
    
    # Meac
    ax[1].plot(m_mr_n2,m_p,ls=m_ls,color='forestgreen',label="$N_2$")
    ax[1].plot(m_mr_o2,m_p,ls=m_ls,color='blue',label="$O_2$")
    ax[1].plot(m_mr_o3,m_p,ls=m_ls,color='navy',label="$O_3$")
    ax[1].plot(m_mr_h2,m_p,ls=m_ls,color='purple',label="$H_2$")
    ax[1].plot(m_mr_h2o,m_p,ls=m_ls,color='cyan',label="$H_2O$")
    ax[1].plot(m_mr_ch4,m_p,ls=m_ls,color='orange',label="$CH_4$")
    ax[1].plot(m_mr_co,m_p,ls=m_ls,color='maroon',label="$CO$")
    ax[1].plot(m_mr_co2,m_p,ls=m_ls,color='red',label="$CO_2$")
    ax[1].plot(m_mr_oh,m_p,ls=m_ls,color='pink',label="$OH$")

    # Photochem
    ax[1].plot(p_mr_n2,p_p,ls=p_ls,color='forestgreen')
    ax[1].plot(p_mr_o2,p_p,ls=p_ls,color='blue')
    ax[1].plot(p_mr_o3,p_p,ls=p_ls,color='navy')
    ax[1].plot(p_mr_h2,p_p,ls=p_ls,color='purple')
    ax[1].plot(p_mr_h2o,p_p,ls=p_ls,color='cyan')
    ax[1].plot(p_mr_ch4,p_p,ls=p_ls,color='orange')
    ax[1].plot(p_mr_co,p_p,ls=p_ls,color='maroon')
    ax[1].plot(p_mr_co2,p_p,ls=p_ls,color='red')
    ax[1].plot(p_mr_oh,p_p,ls=p_ls,color='pink')
    
    ax[1].legend()
    ax[1].invert_yaxis()
    ax[1].set_xlim(1e-15,1)
    ax[1].set_ylim(max(m_p),min(m_p))
    ax[0].set_ylim(max(m_p),min(m_p))
    ax[1].set_xscale('log')
    ax[1].set_yscale('log')

    fig.set_figwidth(8)
    plt.savefig('comp',dpi=300)
    plt.show()

def writePhotochemTPasMEACTP():
    #   Photochem TP profile
    ztp_pc = open('../../pc/ex_evoatmosphere/atmosphere.txt','r')
    z,t,p = [],[],[]
    for line in ztp_pc.readlines()[1:]:
        line = line.split()
        z.append(float(line[0]))
        t.append(float(line[3]))
        p.append(np.log10(float(line[1])*1e5))

    #   Fill missing altitudes up to 100km by drawing isotherm
    zz = z[-1]
    pp = p[-1]
    dz = z[-1]-z[-2]
    dp = p[-1]-p[-2]
    while zz<100:
        zz += dz
        pp += dp
        z.append(zz)
        t.append(t[-1])
        p.append(pp)

    #   MEAC format TP profile
    f = open('TP.dat','w')
    for i in range(len(t)):
        f.write(np.format_float_positional(z[i],precision=6,min_digits=6)+'\t')
        f.write(np.format_float_positional(p[i],precision=6,min_digits=6)+'\t')
        f.write(np.format_float_positional(t[i],precision=6,min_digits=6)+'\n')


    f.close()

def binPhotochemSpectrum():
    # Bin the solar spectrum from Photochem into bins with the same centers and
    #   widths as the bins in the MEAC solar spectrum. Only 100nm-400nm.

    fig,ax = plt.subplots(nrows=2)

    #   Photochem
    #   Spectral bins are 1nm-wide up to 120nm, then 0.05nm-wide thru 400nm and beyond
    spec_pc = open('../../pc/ex_evoatmosphere/Sun_now.txt','r')
    wvln,flux = [],[]
    for line in spec_pc.readlines()[1:]:
        line = line.split()
        wvln.append(float(line[0]))
        flux.append(float(line[1])/1e3)
    ax[0].plot(wvln,flux,c='grey',label="Photochem")

    new_wvln,new_flux = [],[]
    for i in range(100,401):
        bin_total = 0
        bin_count = 0
        for j in range(len(wvln)):
            if np.floor(wvln[j])==i:
                bin_total += flux[j]
                bin_count += 1
        new_wvln.append(i+0.5)
        new_flux.append(bin_total/bin_count)
    # ax[0].plot(new_wvln,new_flux,c='red',label="Photochem, binned")

    #   MEAC
    #   The solar spectrum employed by MEAC uses 1nm-wide bins up til 630, then goes to 2nm-wide bins
    spec_meac = open('hu-code-sr/Data/solar0.txt','r')
    wvln,flux = [],[]
    for line in spec_meac.readlines()[:301]:
        line = line.split()
        wvln.append(float(line[0]))
        flux.append(float(line[1]))
    ax[0].plot(wvln,flux,c='black',label="CLIMA/MEAC")

    #   Relative difference of Photochem from CLIMA/MEAC
    ax[1].hlines(1,xmin=100,xmax=400,colors='black')
    new_flux = np.asarray(new_flux)
    flux = np.asarray(flux)
    ax[1].plot(wvln,new_flux/flux,c='red',label="Photochem, binned")
    ax[1].set_xlim(100,160)
    
    ave = np.mean(new_flux/flux)
    print(ave)

    ax[0].set_ylabel("Flux [W/m^2/nm]")
    ax[0].set_yscale('log')
    ax[0].legend()
    ax[0].set_xlim(100,160)
    ax[1].set_xlabel("Wavelength [nm]")
    ax[1].set_ylabel("Relative flux")
    ax[1].legend()

    fig.set_figwidth(12)
    plt.tight_layout()
    plt.savefig('fig',dpi=300)
    plt.show()

def plotCH4GainLossMechanisms():

    fig,ax = plt.subplots(ncols=2)

    ### Production
    # scenario gain variables
    gains = []
    rTotals = []
    m57s = []
    r71s = []

    for i in range(1,5):
        # Define some variables
        totalGain = 0
        rTotal = 0
        m57 = 0
        r71 = 0

        # Open the reaction rates file and get all lines
        f = open(f'hu-code-sr/scenario_library/{MSCENARIONAME}/int.rates.out3.dat','r')
        lines = f.readlines()
        f.close()

        # Discard all lines except the CH4 loss section
        start,end = 0,0
        # Find start tediously
        i = 0
        while i < 1e6:
            if "CH4\t\t\t\tPRODUCTION" in lines[i]: 
                lines = lines[i:]
                break
            i+=1
        # Find end tediously
        i = 0
        while i < 1e6:
            if "CH4\t\t\t\tLOSS" in lines[i]:
                lines = lines[:i]
                break
            i += 1

        # Split each line
        for i in range(len(lines)):
            lines[i] = lines[i].split()

        # Set variables
        totalGain = float(lines[0][-1])
        for line in lines[1:]:
            if line==[]: continue
            if "R" in line[0]: 
                rTotal += float(line[-1])
            # Specific reactions
            if line[0]=="R71": r71=float(line[-1])
            if line[0]=="M57": m57=float(line[-1])
        
        # Save variables

        print(rTotal)
        gains.append(totalGain)
        rTotals.append(rTotal)
        m57s.append(m57)
        r71s.append(r71)

    # General trends
    ax[0].plot(gains,'ko-',label="CH4 Production")
    ax[0].plot(rTotals,'bo-',label="R Rxs")
    # Specific reactions
    ax[0].plot(r71s,color='navy',ls=':',label="CH3 + CHO = CH4 + CO")
    ax[0].plot(m57s,color='red',ls=':',label="H + CH3 = CH4")

    ax[0].set_xlabel(r"$f_{\rm CO2}$",size='x-large')
    ax[0].set_ylabel("Production Rate",size='x-large')
    ax[0].set_xticks(ticks=[0,1,2,3],labels=[r"$10^{-1}$",r"$10^{-2}$",r"$10^{-3}$",r"$10^{-4}$"],minor=False)
    ax[0].set_yscale('log')
    ax[0].legend()

    ### Losses
    # lists of scenario loss variables
    losses = []
    rTotals = []
    pTotals = []
    r501s = []
    r157s = []
    p33s = []
    p32s = []

    for i in range(1,5):
        # Define some variables
        totalLoss = 0
        rTotal = 0
        pTotal = 0
        r501 = 0
        r157 = 0
        p33 = 0
        p32 = 0

        # Open the reaction rates file and get all lines
        f = open(f'hu-code-sr/scenario_library/Sun/N2_CO2_1e-{i}-Full/int.rates.out3.dat','r')
        lines = f.readlines()
        f.close()

        # Discard all lines except the CH4 loss section
        start,end = 0,0
        # Find start tediously
        i = 0
        while i < 1e6:
            if "CH4\t\t\t\tLOSS" in lines[i]: 
                lines = lines[i:]
                break
            i+=1
        # Find end tediously
        i = 0
        while i < 1e6:
            if "*****" in lines[i]:
                lines = lines[:i]
                break
            i += 1

        # Split each line
        for i in range(len(lines)):
            lines[i] = lines[i].split()

        # Set variables
        totalLoss = float(lines[0][-1])
        for line in lines[1:]:
            if line==[]: continue
            # Simple -- just R vs P
            if "P" in line[0]: pTotal += float(line[-1])
            elif "R" in line[0]: rTotal += float(line[-1])
            # Specific reactions
            if line[0]=="R501": r501=float(line[-1])
            if line[0]=="R157": r157=float(line[-1])
            if line[0]=="P33": p33=float(line[-1])
            if line[0]=="P32": p32=float(line[-1])
        
        # Save variables
        losses.append(totalLoss)
        pTotals.append(pTotal)
        rTotals.append(rTotal)
        r501s.append(r501)
        r157s.append(r157)
        p33s.append(p33)
        p32s.append(p32)
    
    # Make the data pretty
    
    # General trends
    ax[1].plot(losses,'ko-',label="CH4 Loss")
    ax[1].plot(pTotals,'go-',label="P Rxs")
    ax[1].plot(rTotals,'bo-',label="R Rxs")

    # Specific reactions
    ax[1].plot(r501s,color='navy',ls=':',label="OH + CH4 = CH3 + H2O")
    ax[1].plot(r157s,color='cornflowerblue',ls=":",label="H + CH4 = CH3 + H2")
    ax[1].plot(p33s,color="seagreen",ls=":",label="CH4 = CH21 + H2")
    ax[1].plot(p32s,color="lawngreen",ls=":",label="CH4 = CH3 + H")

    ax[1].set_xlabel(r"$f_{\rm CO2}$",size='x-large')
    ax[1].set_ylabel("Loss Rate",size='x-large')
    ax[1].set_xticks(ticks=[0,1,2,3],labels=[r"$10^{-1}$",r"$10^{-2}$",r"$10^{-3}$",r"$10^{-4}$"],minor=False)
    ax[1].set_yscale('log')
    ax[1].legend()
    # ax.set_ylim(1e8,1e11)
    fig.set_figwidth(12)
    plt.tight_layout()
    # plt.savefig("outputs/plots/ch4_gain_loss_mechanisms",dpi=250)
    plt.show()
    
def plotProductionAndLossRates(species:str='CH4'):

    # declare file names
    vert_rates      = "hu-code-sr/scenario_library/guzman-marmolejo/fco2_1e-1/ChemicalRate.dat"
    prod_loss_rates = 'hu-code-sr/scenario_library/guzman-marmolejo/fco2_1e-1/int.rates.out3.dat'
    colint_rates    = "hu-code-sr/scenario_library/guzman-marmolejo/fco2_1e-1/colintrxnrates.dat"

    # get vertical chemical rates for later
    vr = open(vert_rates,'r')
    lines = vr.readlines()[1:]
    vr.close()
    rate_dict = {}
    for l in lines:
        l = l.split()
        rate_dict[l[0]] = l[1:]
    
    # Get column-integrated reaction rates
    colr = open(colint_rates,'r')
    lines = colr.readlines()[1:]
    colr.close()
    colint_rate_dict = {}
    for l in lines:
        l = l.split()
        colint_rate_dict[l[0]] = float(l[1])

    # Open the reaction rates file and get all lines
    f = open(prod_loss_rates,'r')
    lines = f.readlines()
    f.close()

    # Discard all lines except the species loss section
    prodstart,prodend,lossend = 0,0,0
    start,mid,end = 0,0,0
    # Find start tediously
    i = 0
    while i < 1e6:
        if f"{species}\t\t\t\tPRODUCTION" in lines[i]: 
            prodstart = i
            break
        i+=1
    # Find mid tediously
    i = 0
    while i < 1e6:
        if f"{species}\t\t\t\tLOSS" in lines[i]:
            prodend = i
            break
        i += 1
    prodlines = lines[prodstart:prodend]
    # Find end tediously
    i = prodend+1
    while i < 1e6:
        if "LOSS" in lines[i]:
            lossend = i
            break
        i += 1
    losslines = lines[prodend:lossend]
    print(prodstart,prodend,lossend)
    # Split each line and get reactions that contribute to production and loss, and plot them
    fig,ax = plt.subplots(ncols=2)

    z = np.linspace(0,100,50)
    for l in prodlines[1:11]:
        l = l.split()
        if l: 
            reac = l[0]
            rates = np.array(rate_dict[reac],dtype=np.float32)
            ax[0].plot(rates,z,label=reac)
    ax[0].legend()
    ax[0].set_xlim(1e-15,1e5)
    ax[0].set_xscale('log')
    for l in losslines[1:11]:
        l = l.split()
        if l:
            reac = l[0]
            rates = np.array(rate_dict[reac],dtype=np.float32)
            ax[1].plot(rates,z,label=reac)
    ax[1].legend()
    ax[1].set_xlim(1e-15,1e5)
    ax[1].set_xscale('log')

    ax[0].set_title(f"{species} Production Rates",size='x-large')
    ax[1].set_title(f"{species} Loss Rates",size='x-large')
    plt.show()

def compareCH4ProductionAndLoss(meac:str,photochem:str):
    ### Plot reaction profiles
    fig,ax = plt.subplots(nrows=2,ncols=2)
    ax = ax.flatten()
    ##### Production

    ### Photochem

    # Production
    p = np.loadtxt(f'{photochem}/atmosphere.txt',skiprows=1)
    p = p[:,1]*1e5       # bar to Pa    
    rxn = np.zeros((5,100))
    f = open(f'{photochem}/production.dat','r')
    lines = f.readlines()[:10]
    f.close()

    labels = []
    for i in range(5):
        line = lines[i].strip('\n').split("\t")
        labels.append(' '.join(line[:-101]))
        for j in range(100):
            print(len(line))
            rxn[i,j] = float(line[-101+j])
    
    for i in range(5):
        ax[0].plot(rxn[i,:],p,ls=':',label=labels[i])

    ### MEAC
    int_rates       = f"{meac}/int.rates.out3.dat"
    vert_rates      = f"{meac}/ChemicalRate.dat"
    conc            = f"{meac}/ConcentrationSTD.dat"

    # Get names of reactions with greatest rates
    f = open(int_rates,'r')
    lines = f.readlines()
    f.close()

    prod,loss = 0,0
    for i in range(len(lines)):
        line = lines[i]
        if "CH4\t\t\t\t" in line:
            if prod==0: prod = i
            else:
                loss = i
                break
    
    rxns = []
    labels = []
    for i in range(prod+2,prod+7):
        line = lines[i].split()
        rxns.append(line[0])
        labels.append(' '.join(line[1:-1]))
                

    # get vertical chemical rates for later
    vr = open(vert_rates,'r')
    lines = vr.readlines()[1:]
    vr.close()
    rate_dict = {}
    for l in lines:
        l = l.split()
        rate_dict[l[0]] = l[1:]

    # Pressure
    f = np.loadtxt(conc,skiprows=2)
    pres = f[:,4]
    
    for i in range(5):
        vals = rate_dict[rxns[i]]
        vals = np.array(vals,dtype=np.float32)
        ax[2].plot(vals,pres,label=labels[i])
    
    for a in ax:
        a.set_xscale('log')
        a.set_xlim(1e-15,1e5)
        a.set_ylim(1e5,1e-1)
        a.set_yscale('log')
    ax[0].set_title("CH4 Production",size='x-large')
    ax[0].set_xlabel(r"Production rate [$cm^{-2} s^{-1}$ ?]",size='large')
    ax[0].set_ylabel("Pressure [Pa]",size='large')
    ax[0].legend()

    ##### Loss

    ### Photochem

    # Loss reactions
    #   CH4 + O => CH3 + OH     ->  R423/R450
    #   CH4 + OH => CH3 + H2O   ->  R501
    #   CH4 + hv => 1CH2 + H2   ->  P33
    #   CH4 + hv => CH3 + H     ->  P32
    #   CH4 + hv => CH + H2 + H ->  P34

    rxn = np.zeros((10,100))
    f = open('photochem_rates/loss.dat','r')
    lines = f.readlines()
    f.close()
    for i in range(10):
        line = lines[i].split("\t")
        for j in range(100):
            rxn[i,j] = float(line[1+j])

    ax[1].plot(rxn[0,:],p,c='red',ls=':',label="CH4 + O => CH3 + OH")
    ax[1].plot(rxn[1,:],p,c='orange',ls=':',label="CH4 + OH => CH3 + H2O")
    ax[1].plot(rxn[2,:],p,c='forestgreen',ls=':',label="CH4 + hv => 1CH2 + H2")
    ax[1].plot(rxn[3,:],p,c='cornflowerblue',ls=':',label="CH4 + hv => CH3 + H")
    ax[1].plot(rxn[4,:],p,c='pink',ls=':',label="CH4 + hv => CH + H2 + H")
    ax[1].plot(0,0,c='black',label='Coupled CLIMA/MEAC')
    ax[1].plot(0,0,c='black',ls=':',label="Photochem")

    rxns = ['R423','R501','P33','P32','P34']
    cols = ['red','orange','forestgreen','cornflowerblue','pink']
    for i in range(5):
        vals = rate_dict[rxns[i]]
        vals = np.array(vals,dtype=np.float32)
        ax[1].plot(vals,pres,c=cols[i])
    
    ax[1].set_xscale('log')
    ax[1].set_xlim(1e-15,1e5)
    ax[1].set_ylim(1e5,1e-1)
    ax[1].set_yscale('log')
    ax[1].set_title("CH4 Loss",size='x-large')
    ax[1].set_xlabel(r"Loss rate [$cm^{-2} s^{-1} $?]",size='large')
    ax[1].set_ylabel("Pressure [Pa]",size='large')
    ax[1].legend()

    fig.set_figwidth(12)
    fig.set_figheight(5)
    plt.tight_layout()
    plt.savefig(f"outputs/{NAME}/prod_loss_rates",dpi=250)
    plt.show()


# print('test')
compMEACandPhotochemCompositions('hu-code-sr/scenario_library/CO2_CH4/fco2_1e-1/','../../pc/fco2_1e-1/')
# compareCH4ProductionAndLoss('hu-code-sr/scenario_library/CO2_CH4/fco2_1e-1/','../../pc/fco2_1e-1/')
# plotProductionAndLossRates()
# printMEACztpFromCLAST()
# writePhotochemTPasMEACTP()
# compMEACandPhotochem()
# binPhotochemSpectrum()

