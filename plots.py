import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import re

from parameters import *

# plot parameters
figsize   = (8,6)
dpi       = 300
paper     = True # False -> Presentation
if paper:
  linewidth  = 2.0 #linewidth
  tlinewidth = 1.5 #tick linewidth
  tfontsize  = 16  #tick label fontsize
  lfontsize  = 16  #internal labeling fontsize
  afontsize  = 16  #axis label fontsize
  axmajor    = 7   #major tick length
  axminor    = 4   #minor tick length
  pad        = 8   #tick label padding
else:
  linewidth  = 2.0 #linewidth
  tlinewidth = 1.5 #linewidth
  tfontsize  = 18  #tick fontsize
  lfontsize  = 18  #labeling fontsize
  afontsize  = 18  #axis fontsize
  axmajor    = 7   #major tick length
  axminor    = 4   #minor tick length
  pad        = 8   #tick label padding
plt.rc('font',  family='sans-serif')
plt.rc('xtick', labelsize=tfontsize)
plt.rc('ytick', labelsize=tfontsize)
plt.rc('axes',  linewidth=linewidth)
plt.rcParams['xtick.major.size'] = axmajor
plt.rcParams['ytick.major.size'] = axmajor
plt.rcParams['xtick.minor.size'] = axminor
plt.rcParams['ytick.minor.size'] = axminor
plt.rcParams['xtick.major.width'] = tlinewidth
plt.rcParams['ytick.major.width'] = tlinewidth
plt.rcParams['xtick.minor.width'] = tlinewidth
plt.rcParams['ytick.minor.width'] = tlinewidth
cmap = plt. get_cmap('tab20b')

def plotAtmosphericComposition(conc_file:str,id:str,ref_file:str='',out_dir:str=''):
    """
    Plot relevant mixing ratios from MEAC concentration file
    """
    f = open(conc_file,'r')
    data = f.read().replace('#','').split()[121:]
    data = np.asarray(data,dtype=np.float32).reshape(len(data)//116,116)
    f.close()

    # TOA to surface
    alts = data[::-1,0]

    # Total number densities
    nd_all  =   np.sum(data[::-1],axis=1)
    
    # number densities                      mixing ratios
    ndC2H6  =   data[::-1,35];              mrC2H6  =   np.divide(ndC2H6,nd_all)    # ABSOLUTE
    ndCH4   =   data[::-1,25];              mrCH4   =   np.divide(ndCH4,nd_all)     # ABSOLUTE
    ndCO2   =   data[::-1,56];              mrCO2   =   np.divide(ndCO2,nd_all)     # ABSOLUTE
    ndH2    =   data[::-1,57];              mrH2    =   np.divide(ndH2,nd_all)      # ABSOLUTE
    ndH2O   =   data[::-1,11];              mrH2O   =   np.divide(ndH2O,nd_all)     # ABSOLUTE
    ndN2    =   data[::-1,59];              mrN2    =   np.divide(ndN2,nd_all)      # ABSOLUTE
    ndO2    =   data[::-1,58];              mrO2    =   np.divide(ndO2,nd_all)      # ABSOLUTE      
    ndO3    =   data[::-1,6];               mrO3    =   np.divide(ndO3,nd_all)      # ABSOLUTE
    
    # Writing mixing ratio profiles
    nd_profiles =   [ndC2H6,ndCH4,ndCO2,ndH2,ndH2O,ndN2,ndO2,ndO3]
    mr_profiles =   [mrC2H6,mrCH4,mrCO2,mrH2,mrH2O,mrN2,mrO2,mrO3]
    colors      =   ['gold','orange','red','purple','cyan','green','cornflowerblue','blue']
    labels      =   ['C2H6','CH4','CO2','H2','H2O','N2','O2','O3']

    # Plot outputs
    fig,ax = plt.subplots()
    for j in range(len(nd_profiles)):
        # ax[0].plot(nd_profiles[j],alts,c=colors[j],label=labels[j])
        ax.plot(mr_profiles[j],alts,c=colors[j],label=labels[j])
    
    if ref_file:
        f = open(ref_file,'r')
        data = f.read().replace('#','').split()[121:]
        data = np.asarray(data,dtype=np.float32).reshape(len(data)//116,116)
        f.close()

        # TOA to surface
        alts = data[::-1,0]

        # Total number densities
        nd_all  =   np.sum(data[::-1],axis=1)
        
        # number densities                      mixing ratios
        ndC2H6  =   data[::-1,35];              mrC2H6  =   np.divide(ndC2H6,nd_all)    # ABSOLUTE
        ndCH4   =   data[::-1,25];              mrCH4   =   np.divide(ndCH4,nd_all)     # ABSOLUTE
        ndCO2   =   data[::-1,56];              mrCO2   =   np.divide(ndCO2,nd_all)     # ABSOLUTE
        ndH2    =   data[::-1,57];              mrH2    =   np.divide(ndH2,nd_all)      # ABSOLUTE
        ndH2O   =   data[::-1,11];              mrH2O   =   np.divide(ndH2O,nd_all)     # ABSOLUTE
        ndN2    =   data[::-1,59];              mrN2    =   np.divide(ndN2,nd_all)      # ABSOLUTE
        ndO2    =   data[::-1,58];              mrO2    =   np.divide(ndO2,nd_all)      # ABSOLUTE      
        ndO3    =   data[::-1,6];               mrO3    =   np.divide(ndO3,nd_all)      # ABSOLUTE
        
        # Writing mixing ratio profiles
        nd_profiles =   [ndC2H6,ndCH4,ndCO2,ndH2,ndH2O,ndN2,ndO2,ndO3]
        mr_profiles =   [mrC2H6,mrCH4,mrCO2,mrH2,mrH2O,mrN2,mrO2,mrO3]
        colors      =   ['gold','orange','red','purple','cyan','green','cornflowerblue','blue']
        labels      =   ['C2H6','CH4','CO2','H2','H2O','N2','O2','O3']

        # Plot outputs
        for j in range(len(nd_profiles)):
            # ax[0].plot(nd_profiles[j],alts,ls=':',c=colors[j])
            ax.plot(mr_profiles[j],alts,ls=':',c=colors[j])

    # ax[0].set_xlim(left=1e8)
    # ax[0].set_ylim(bottom=0,top=max(alts))
    # ax[0].set_xscale('log')
    # ax[0].set_xlabel("Number density",size='x-large')
    # ax[0].set_ylabel("Altitude [km]",size='x-large')
    # ax[0].legend()
    ax.set_xlim(left=1e-25)
    ax.set_ylim(bottom=0,top=max(alts))
    ax.set_xscale('log')
    ax.set_xlabel("Mixing ratio",size='x-large')
    ax.set_ylabel("Altitude [km]",size='x-large')
    ax.legend()
    if ref_file: fig.suptitle('New -> solid, ref -> dotted')
    fig.set_figwidth(11)
    plt.tight_layout()
    plt.savefig(f"{out_dir}/mr/MEAC_mixing_ratios_{id}",dpi=300)
    # plt.show()
    plt.close()

def plotAtmosphericEvolution(scen_name:str='',out_dir:str=''):         # WIP
    """
    
    """
    fig,axes = plt.subplots(nrows=2,ncols=4,sharey=True,gridspec_kw={'wspace':0})
    axes = axes.flatten()
    
    files       =   os.listdir(f"outputs/{scen_name}/meac-out/")
    files.sort(key = lambda x: int(re.search(r"[0-9]+",x)[0]))
    colors      =   ['Greys','Oranges','Reds','Purples','Blues','Greens','GnBu','Blues']
    labels      =   ['C2H6','CH4','CO2','H2','H2O','N2','O2','O3']

    for i in range(len(files)):
        file = files[i]
        l = len(files)
        if file=='.DS_Store':continue
        
        f = open(f"outputs/{scen_name}/meac-out/{file}",'r')
        data = np.fromstring(''.join(f.readlines()[1:-1]),dtype=np.float32,sep=' ')
        # print(data)
        data = data.reshape(len(data)//19,19)
        alts = data[:,1]

        for j in range(len(labels)):
            cmap = matplotlib.colormaps[colors[j]]
            axes[j].plot(data[:,4+2*j],alts,c=cmap(i/l))
            axes[j].set_title(labels[j],size='x-large')
        
    for k in range(len(axes)):
        ax = axes[k]
        ax.set_xscale('log')
        ax.set_xticks([],minor=True)
        ax.tick_params(rotation=-90)
    axes[2].set_xticks(ticks=[2e-5,5e-5,1e-4],labels=[r"$2\times10^{-5}$",r"$5\times10^{-5}$",r"$10^{-4}$"],minor=False)         # CO2
    axes[3].set_xticks(ticks=[4e-4,5e-4],labels=[r"$4\times10^{-4}$",r"$5\times10^{-4}$"],minor=False)  # H2
    axes[5].set_xticks(ticks=[0.99,0.995,1.0],labels=["",r"$9.95\times10^{-1}$",r"1.0"],minor=False)    # N2

    fig.set_figwidth(15)
    fig.set_figheight(8)
    plt.tight_layout()
    # plt.show()
    plt.savefig(f"{out_dir}/MEAC_mr_evolution.png",dpi=250)
    plt.show()
    
def plotTPprofile(clast:str,id:str,meac_conv:str='',out_dir:str=''):
    """
    stuff here
    """
    fig,ax = plt.subplots()

    if meac_conv != '':
        f_first = open(meac_conv,'r')
        data = f_first.read().replace('#','').split()[121:]
        data = np.asarray(data,dtype=np.float32).reshape(len(data)//116,116)
        f_first.close()
        t_first = data[::-1,3]  # TOA to surf
        p_first = data[::-1,4]  # TOA to surf
        ax.plot(t_first,p_first,c='cornflowerblue',ls=':',label="Initial MEAC TP profile")

    f_last = open(clast,'r')
    data = f_last.read().split()[9:]
    data = np.array(data,dtype=np.float32).reshape(ND,9)
    f_last.close()
    t_last = data[:,2]      # TOA to surf
    p_last = data[:,1]*atm2Pa   # TOA to surf

    ax.plot(t_last,p_last,c='navy',label="Final CLIMA TP profile")

    # ax.set_ylim(bottom=max(max(p_last),max(p_first)),top=min(min(p_first),min(p_last)))
    ax.set_ylim((1e5,1))
    ax.set_yscale('log')
    ax.set_xlabel("Temperature [K]",size='x-large')
    ax.set_ylabel("Pressure [Pa]",size='x-large')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'{out_dir}/TPprofile',dpi=200)
    # plt.show()
    plt.close()

def plotTPevolution(scen_name:str,out_dir:str):

    files = os.listdir(f"outputs/{scen_name}/clima-out/")
    files.sort(key = lambda x: int(re.search(r"[0-9]+",x)[0]))

def plotChemTP(conc_file:str,clast:str,ref_meac_file:str,out_dir:str=''):
    fig,axes = plt.subplots(ncols=2)
    ### TP Profile
    f_first = open(ref_meac_file,'r')
    data = f_first.read().replace('#','').split()[121:]
    data = np.asarray(data,dtype=np.float32).reshape(len(data)//116,116)
    f_first.close()
    t_first = data[::-1,3]  # TOA to surf
    p_first = data[::-1,4]  # TOA to surf

    f_last = open(clast,'r')
    data = f_last.read().split()[9:]
    data = np.array(data,dtype=np.float32).reshape(ND,9)
    f_last.close()
    t_last = data[:,2]      # TOA to surf
    p_last = data[:,1]*atm2Pa   # TOA to surf

    # Set TP profile plot parameters
    axes[0].plot(t_first,p_first,c='cornflowerblue',ls=':',label="Initial MEAC TP profile")
    axes[0].plot(t_last,p_last,c='navy',label="Final CLIMA TP profile")

    axes[0].set_ylim(bottom=max(max(p_last),max(p_first)),top=min(min(p_first),min(p_last)))
    axes[0].set_yscale('log')
    axes[0].set_xlabel("Temperature [K]",size='x-large')
    axes[0].set_ylabel("Pressure [Pa]",size='x-large')

    ###
    ### Mixing ratios
    ###
    f = open(conc_file,'r')
    data = f.read().replace('#','').split()[121:]
    data = np.asarray(data,dtype=np.float32).reshape(len(data)//116,116)
    f.close()

    # TOA to surface
    alts = data[::-1,0]

    # Total number densities
    nd_all  =   np.sum(data[::-1],axis=1)
    
    # number densities                      mixing ratios
    ndC2H6  =   data[::-1,35];              mrC2H6  =   np.divide(ndC2H6,nd_all)    # ABSOLUTE
    ndCH4   =   data[::-1,25];              mrCH4   =   np.divide(ndCH4,nd_all)     # ABSOLUTE
    ndCO2   =   data[::-1,56];              mrCO2   =   np.divide(ndCO2,nd_all)     # ABSOLUTE
    ndH2    =   data[::-1,57];              mrH2    =   np.divide(ndH2,nd_all)      # ABSOLUTE
    ndH2O   =   data[::-1,11];              mrH2O   =   np.divide(ndH2O,nd_all)     # ABSOLUTE
    ndN2    =   data[::-1,59];              mrN2    =   np.divide(ndN2,nd_all)      # ABSOLUTE
    ndO2    =   data[::-1,58];              mrO2    =   np.divide(ndO2,nd_all)      # ABSOLUTE      
    ndO3    =   data[::-1,6];               mrO3    =   np.divide(ndO3,nd_all)      # ABSOLUTE
    
    # Writing mixing ratio profiles
    nd_profiles =   [ndC2H6,ndCH4,ndCO2,ndH2,ndH2O,ndN2,ndO2,ndO3]
    mr_profiles =   [mrC2H6,mrCH4,mrCO2,mrH2,mrH2O,mrN2,mrO2,mrO3]
    colors      =   ['gold','orange','red','purple','cyan','green','cornflowerblue','blue']
    labels      =   ['C2H6','CH4','CO2','H2','H2O','N2','O2','O3']

    # Plot outputs
    for j in range(len(nd_profiles)):
        # ax[0].plot(nd_profiles[j],alts,c=colors[j],label=labels[j])
        axes[1].plot(mr_profiles[j],alts,c=colors[j],label=labels[j])
    
    if ref_meac_file:
        f = open(ref_meac_file,'r')
        data = f.read().replace('#','').split()[121:]
        data = np.asarray(data,dtype=np.float32).reshape(len(data)//116,116)
        f.close()

        # TOA to surface
        alts = data[::-1,0]

        # Total number densities
        nd_all  =   np.sum(data[::-1],axis=1)
        
        # number densities                      mixing ratios
        ndC2H6  =   data[::-1,35];              mrC2H6  =   np.divide(ndC2H6,nd_all)    # ABSOLUTE
        ndCH4   =   data[::-1,25];              mrCH4   =   np.divide(ndCH4,nd_all)     # ABSOLUTE
        ndCO2   =   data[::-1,56];              mrCO2   =   np.divide(ndCO2,nd_all)     # ABSOLUTE
        ndH2    =   data[::-1,57];              mrH2    =   np.divide(ndH2,nd_all)      # ABSOLUTE
        ndH2O   =   data[::-1,11];              mrH2O   =   np.divide(ndH2O,nd_all)     # ABSOLUTE
        ndN2    =   data[::-1,59];              mrN2    =   np.divide(ndN2,nd_all)      # ABSOLUTE
        ndO2    =   data[::-1,58];              mrO2    =   np.divide(ndO2,nd_all)      # ABSOLUTE      
        ndO3    =   data[::-1,6];               mrO3    =   np.divide(ndO3,nd_all)      # ABSOLUTE
        
        # Writing mixing ratio profiles
        nd_profiles =   [ndC2H6,ndCH4,ndCO2,ndH2,ndH2O,ndN2,ndO2,ndO3]
        mr_profiles =   [mrC2H6,mrCH4,mrCO2,mrH2,mrH2O,mrN2,mrO2,mrO3]
        colors      =   ['gold','orange','red','purple','cyan','green','cornflowerblue','blue']
        labels      =   ['C2H6','CH4','CO2','H2','H2O','N2','O2','O3']

        # Plot outputs
        for j in range(len(nd_profiles)):
            # ax[0].plot(nd_profiles[j],alts,ls=':',c=colors[j])
            axes[1].plot(mr_profiles[j],alts,ls=':',c=colors[j])
        
    # Set mixing ratio plot params
    axes[1].set_xlim(left=1e-10)
    axes[1].set_ylim(bottom=0,top=max(alts))
    axes[1].set_xscale('log')
    axes[1].set_xlabel("Mixing ratio",size='x-large')
    axes[1].set_ylabel("Altitude [km]",size='x-large')
    axes[1].legend()

    fig.set_figwidth(9)
    plt.tight_layout()
    plt.savefig(f"{OUTPUT}/mr_tp_combined",dpi=250)
    # plt.show()

def plotSurfaceTemperature(temps_file:str,runBreaks:list[int]=[],out_dir:str=''):
    """
    """
    fig,ax = plt.subplots(gridspec_kw={'left':0.1,'right':0.85})

    f = open(temps_file,'r')
    t = ','.join(f.readlines())
    surfTemps = np.fromstring(t,sep=',')
    f.close()
    
    n = np.linspace(1,len(surfTemps),len(surfTemps))
    ax.plot(n,surfTemps,'r')
    ax.scatter(n[-1],surfTemps[-1],s=30,c='red')

    t_last = str(np.round(surfTemps[-1],2))
    ax.text(n[-1]*1.06,surfTemps[-1],r'$T_{\rm Surf}\approx$'+t_last,size='large',c='red')

    ax.set_xlabel("# of CLIMA steps",size='x-large')
    ax.set_ylabel("Surface Temperature [K]",size='x-large')
    fig.set_figwidth(9)
    fig.set_figheight(6)
    plt.savefig(f"{out_dir}/Surface Temperature",dpi=200)
    # plt.show()

def compareScenarios():
    """
    Compare final TP profiles and mixing ratios of the four N2CO2 scenarios
    """
    fig,axes = plt.subplots(ncols=3,nrows=3,gridspec_kw={'wspace':0.1,'hspace':0.45})
    axes = axes.flatten()


    colors = ['mediumspringgreen','darkturquoise','cornflowerblue','slateblue']
    labels = ['C2H6','CH4','CO2','H2','H2O','N2','O2','O3']
    for j in range(len(labels)):
        a = axes[j+1]
        a.set_title(r"$f_{\rm "+labels[j]+r"}$",size='x-large')

    maxt,maxp = 0,0
    mint,minp = 300,100

    for i in range(4):
        scen_name = f"N2_CO2_1e-{i+1}"
        color = colors[i]

        # Temperature-pressure profiles
        temps = os.listdir(f"outputs/{scen_name}/meac-in/")
        temps.sort(key = lambda x: int(re.search(r"[0-9]+",x)[0]))
        temp = open(f"outputs/{scen_name}/meac-in/{temps[-1]}",'r')
        t, p= [],[]
        data = temp.readlines()[2:-1]
        for line in data:
            line = line.split()
            t.append(float(line[3]))
            p.append(float(line[4]))
        axes[0].plot(t,p,c=color,label="pCO2 = $10^{-"+str(i+1)+"}$")
        
        maxt = max(max(t),maxt)
        mint = min(min(t),mint)
        maxp = max(max(p),maxp)
        minp = min(min(p),minp)
        # Concentrations
        
        files       =   os.listdir(f"outputs/{scen_name}/meac-out/")
        files.sort(key = lambda x: int(re.search(r"[0-9]+",x)[0]))

        latest = files[-1]
        f = open(f'outputs/{scen_name}/meac-out/{latest}','r')
        data = np.fromstring(''.join(f.readlines()[1:-1]),dtype=np.float32,sep=' ')
        data = data.reshape(len(data)//19,19)
        alts = data[:,2]

        for j in range(8):
            a = axes[j+1]
            a.plot(data[:,4+2*j],alts,c=color)
            a.set_ylim(maxp,minp)
            a.set_xscale('log')
            a.set_yscale('log')
            if j not in [2,5]: a.set_yticks(ticks=[],minor=False)
            else:
                a.set_yticks(ticks=[1e5,1e2,1e-1],minor=False)

    axes[0].invert_yaxis()
    axes[0].set_yscale('log')
    axes[0].legend()
    axes[0].set_yticks(ticks=[1e5,1e2,1e-1],minor=False)
    axes[0].set_title(r"Temperature [K]")
    for i in [0,3,6]:
        axes[i].set_ylabel(r"Pressure [Pa]",size='x-large')

    # for a in axes[1:]:
    #     a.invert_yaxis()
    axes[6].set_xticks([],minor=True)
    axes[6].set_xticks(ticks=[0.9,1],labels=['0.9','1.0'],minor=False)

    fig.set_figwidth(14)
    fig.set_figheight(8)
    plt.tight_layout()
    plt.savefig('test',dpi=200)
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
        f = open(f'hu-code-sr/scenario_library/Sun/N2_CO2_1e-{i}-Full/int.rates.out3.dat','r')
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
    
def plotH2OGainLossMechanisms():

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
        f = open(f'hu-code-sr/scenario_library/Sun/N2_CO2_1e-{i}-Full/int.rates.out3.dat','r')
        lines = f.readlines()
        f.close()

        # Discard all lines except the CH4 loss section
        start,end = 0,0
        # Find start tediously
        i = 0
        while i < 1e6:
            if "H2O\t\t\t\tPRODUCTION" in lines[i]: 
                lines = lines[i:]
                break
            i+=1
        # Find end tediously
        i = 0
        while i < 1e6:
            if "H2O\t\t\t\tLOSS" in lines[i]:
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
            # if line[0]=="R71": r71=float(line[-1])
            # if line[0]=="M57": m57=float(line[-1])
        
        # Save variables

        print(rTotal)
        gains.append(totalGain)
        # rTotals.append(rTotal)
        # m57s.append(m57)
        # r71s.append(r71)

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
            if "H2O\t\t\t\tLOSS" in lines[i]: 
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
            # if line[0]=="R501": r501=float(line[-1])
            # if line[0]=="R157": r157=float(line[-1])
            # if line[0]=="P33": p33=float(line[-1])
            # if line[0]=="P32": p32=float(line[-1])
        
        # Save variables
        losses.append(totalLoss)
        pTotals.append(pTotal)
        rTotals.append(rTotal)
        # r501s.append(r501)
        # r157s.append(r157)
        # p33s.append(p33)
        # p32s.append(p32)
    
    # Make the data pretty
    
    # General trends
    ax[1].plot(losses,'ko-',label="CH4 Loss")
    ax[1].plot(pTotals,'go-',label="P Rxs")
    ax[1].plot(rTotals,'bo-',label="R Rxs")

    # Specific reactions
    # ax[1].plot(r501s,color='navy',ls=':',label="OH + CH4 = CH3 + H2O")
    # ax[1].plot(r157s,color='cornflowerblue',ls=":",label="H + CH4 = CH3 + H2")
    # ax[1].plot(p33s,color="seagreen",ls=":",label="CH4 = CH21 + H2")
    # ax[1].plot(p32s,color="lawngreen",ls=":",label="CH4 = CH3 + H")

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
 

plotTPprofile(CLAST,id='init_tp',out_dir='gmtest')
