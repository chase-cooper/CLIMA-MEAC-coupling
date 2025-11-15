import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os

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

def plotAtmosphericComposition(conc_file:str,ref_file:str='',out_dir:str=''):
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
    colors      =   ['yellow','orange','red','purple','cyan','green','cornflowerblue','blue']
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
        colors      =   ['yellow','orange','red','purple','cyan','green','cornflowerblue','blue']
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
    ax.set_xlim(left=1e-10)
    ax.set_ylim(bottom=0,top=max(alts))
    ax.set_xscale('log')
    ax.set_xlabel("Mixing ratio",size='x-large')
    ax.set_ylabel("Altitude [km]",size='x-large')
    ax.legend()
    if ref_file: fig.suptitle('New -> solid, ref -> dotted')
    fig.set_figwidth(5)
    plt.tight_layout()
    # plt.show()
    plt.savefig(f"{out_dir}/MEAC_mixing_ratios",dpi=300)
    plt.close()

def plotAtmosphericEvolution(scen_name:str='',out_dir:str=''):         # WIP
    """
    
    """
    fig,axes = plt.subplots(nrows=2,ncols=4,sharey=True,gridspec_kw={'wspace':0})
    axes = axes.flatten()
    
    files       =   os.listdir(f"outputs/{scen_name}/meac-out/")
    files.sort()
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
        if k not in [2,5,6]: ax.set_xscale('log')
    fig.set_figwidth(15)
    fig.set_figheight(8)
    plt.tight_layout()
    # plt.show()
    plt.savefig(f"{out_dir}/MEAC_mr_evolution.png",dpi=250)
    plt.show()
    
def plotTPprofile(clast:str,meac_conv:str='',out_dir:str=''):
    """
    stuff here
    """
    fig,ax = plt.subplots()

    f_first = open(meac_conv,'r')
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

    ax.plot(t_first,p_first,c='cornflowerblue',ls=':',label="Initial MEAC TP profile")
    ax.plot(t_last,p_last,c='navy',label="Final CLIMA TP profile")

    ax.set_ylim(bottom=max(max(p_last),max(p_first)),top=min(min(p_first),min(p_last)))
    ax.set_yscale('log')
    ax.set_xlabel("Temperature [K]",size='x-large')
    ax.set_ylabel("Pressure [Pa]",size='x-large')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'{out_dir}/TPprofile',dpi=200)
    plt.close()

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
    colors      =   ['yellow','orange','red','purple','cyan','green','cornflowerblue','blue']
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
        colors      =   ['yellow','orange','red','purple','cyan','green','cornflowerblue','blue']
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
    plt.savefig(f"{out_dir}/surfTemps",dpi=200)
    # plt.show()

# plotSurfaceTemperature('outputs/test/surftemps.dat',out_dir=OUTPUT)
# plotAtmosphericComposition(MCONC,MCONV,OUTPUT)
# plotChemTP(MCONC,CLAST,MCONV,OUTPUT)
