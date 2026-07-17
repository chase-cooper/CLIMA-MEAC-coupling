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
    ndOH    =   data[::-1,8];               mrOH    =   np.divide(ndOH,nd_all)
    
    # Writing mixing ratio profiles
    nd_profiles =   [ndC2H6,ndCH4,ndCO2,ndH2,ndH2O,ndN2,ndO2,ndO3,ndOH]
    mr_profiles =   [mrC2H6,mrCH4,mrCO2,mrH2,mrH2O,mrN2,mrO2,mrO3,mrOH]
    colors      =   ['gold','orange','red','purple','cyan','green','cornflowerblue','blue','pink']
    labels      =   ['C2H6','CH4','CO2','H2','H2O','N2','O2','O3','OH']

    # Plot outputs
    fig,ax = plt.subplots()
    for j in range(len(nd_profiles)):
        if j == 0: continue
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
        ndOH    =   data[::-1,8];               mrOH    =   np.divide(ndOH,nd_all)
        
        # Writing mixing ratio profiles
        nd_profiles =   [ndC2H6,ndCH4,ndCO2,ndH2,ndH2O,ndN2,ndO2,ndO3,ndOH]
        mr_profiles =   [mrC2H6,mrCH4,mrCO2,mrH2,mrH2O,mrN2,mrO2,mrO3,mrOH]
        colors      =   ['gold','orange','red','purple','cyan','green','cornflowerblue','blue','pink']
        labels      =   ['C2H6','CH4','CO2','H2','H2O','N2','O2','O3','OH']

        ax.plot(0,0,color='black',ls=':',label="Old (Sat. reduction of 0.2)")

        # Plot outputs
        for j in range(len(nd_profiles)):
            if j == 0: continue
            # ax[0].plot(nd_profiles[j],alts,ls=':',c=colors[j])
            ax.plot(mr_profiles[j],alts,ls=':',c=colors[j])

    # ax[0].set_xlim(left=1e8)
    # ax[0].set_ylim(bottom=0,top=max(alts))
    # ax[0].set_xscale('log')
    # ax[0].set_xlabel("Number density",size='x-large')
    # ax[0].set_ylabel("Altitude [km]",size='x-large')
    # ax[0].legend()
    ax.set_xlim(left=1e-15)
    ax.set_ylim(bottom=0,top=max(alts))
    ax.set_xscale('log')
    ax.set_xlabel("Mixing ratio",size='x-large')
    ax.set_ylabel("Altitude [km]",size='x-large')
    ax.legend()
    if ref_file: fig.suptitle(r'$f_{\rm CO_2}=10^{-4}$')
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
    # axes[2].set_xticks(ticks=[2e-5,5e-5,1e-4],labels=[r"$2\times10^{-5}$",r"$5\times10^{-5}$",r"$10^{-4}$"],minor=False)         # CO2
    # axes[3].set_xticks(ticks=[1.5e-5,2e-5,2.5e-5],labels=[r"$1.5\times10^{-5}$",r"$2\times10^{-5}$",r"$2.5\times10^{-5}$"],minor=False)  # H2
    # axes[5].set_xticks(ticks=[0.99,0.995,1.0],labels=["",r"$9.95\times10^{-1}$",r"1.0"],minor=False)    # N2

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

#####################
### Data Products ###
#####################

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
    p_ind_n2    = 47
    p_ind_o2    = 22
    p_ind_o3    = 61
    p_ind_h2    = 18
    p_ind_h2o   = 19
    p_ind_ch4   = 31
    p_ind_co    = 23
    p_ind_co2   = 24
    p_ind_oh    = 20

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

    plt.show()

def getMEACintRates(meac:str,mol:str='CH4',mode:str='loss'):
    """
    Docstring for getMEACintRates
    
    :param meac: MEAC scenario folder path
    :type meac: str
    :param mode: specify production or loss rates
    :type mode: str
    """

    int_rates       = f"{meac}/int.rates.out3.dat"

    # Open rates file
    f = open(int_rates,'r')
    lines = f.readlines()
    f.close()

    # Get section start
    index = -1
    for i in range(len(lines)):
        line = lines[i]

        key = f"{mol.upper()}\t\t\t\t{mode.upper()}"
        if key in line:
            index = i
            break

    # Build dictionary of reaction rates
    rxns = {}
    j = i+2
    while '=' in lines[j]:
        line = lines[j].split()
        reaction = ' '.join(line[1:-1])
        rxns[reaction] = float(line[-1])
        j += 1
    
    print(f"{mol.upper()} {mode.lower()} rates:")
    for k in rxns.keys():
        print(k,' '*(30-len(k)),rxns[k])

    return rxns

def getMEACcolRates(meac:str,rxn:str="M57"):
    
    rates = meac+'ChemicalRate.dat'
    f = open(rates,'r')
    lines = f.readlines()
    f.close()

    for line in lines:
        line = line.split()
        if line[0]==rxn:
            vals = line[1:]
            break
    return np.array(vals,dtype=np.float32)

def getPCintRates(pc:str,mode:str='loss'):
    """
    Docstring for getPCintRates
    Only CH4 rxns
    
    :param pc: Description
    :type pc: str
    :param mode: Description
    :type mode: str
    """
    f = open(pc+"int_"+mode.lower()+'.dat','r')
    lines = f.readlines()
    f.close()

    rxns = {}
    for line in lines:
        line = line.split()
        rxn = ' '.join(line[:-1])
        rxns[rxn] = float(line[-1])
    
    return rxns

def getPCcolRates(pc:str,mol:str='CH4',mode:str='loss'):
    file = pc+mol.upper()+"/"+mode.lower()+'.dat'
    f = open(file,'r')
    lines = f.readlines()
    f.close()

    rxns = {}
    for line in lines:
        line = line.split()
        rxns[' '.join(line[:-100])] = np.array(line[-100:],dtype=np.float32)
    return rxns

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
        scen_name = f"Sun_fCO2_1e-{i+1}"
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

def compareCH4ProductionAndLoss(meac:str,photochem:str):
    ### Plot reaction profiles
    fig,ax = plt.subplots(nrows=2,ncols=2)
    ax = ax.flatten()

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
            rxn[i,j] = float(line[-101+j])
    
    for i in range(5):
        ax[0].plot(rxn[i,:],p,ls=':',label=labels[i])
    ax[0].set_title("CH4 Production",size='x-large')
    ax[0].set_xlabel(r"Production rate [$cm^{-2} s^{-1}$ ?]",size='large')
    ax[0].set_ylabel("Photochem",size='x-large')
    ax[0].legend()

    # Loss
    rxn = np.zeros((5,100))
    f = open(f'{photochem}/loss.dat','r')
    lines = f.readlines()[:10]
    f.close()

    labels = []
    for i in range(5):
        line = lines[i].strip('\n').split("\t")
        labels.append(' '.join(line[:-101]))
        for j in range(100):
            rxn[i,j] = float(line[-101+j])
    
    for i in range(5):
        ax[1].plot(rxn[i,:],p,ls=':',label=labels[i])
    ax[1].set_title("CH4 Loss",size='x-large')
    ax[1].set_xlabel(r"Production rate [$cm^{-2} s^{-1}$ ?]",size='large')
    ax[1].legend()


    ### MEAC

    int_rates       = f"{meac}/int.rates.out3.dat"
    vert_rates      = f"{meac}/ChemicalRate.dat"
    conc            = f"{meac}/ConcentrationSTD.dat"

    # Pressure
    f = np.loadtxt(conc,skiprows=2)
    pres = f[:,4]

    # get vertical chemical rates for later
    vr = open(vert_rates,'r')
    lines = vr.readlines()[1:]
    vr.close()
    rate_dict = {}
    for l in lines:
        l = l.split()
        rate_dict[l[0]] = l[1:]

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
    
    # Production
    rxns = []
    labels = []
    for i in range(prod+2,prod+7):
        line = lines[i].split()
        rxns.append(line[0])
        labels.append(' '.join(line[1:-1]))

    for i in range(5):
        vals = rate_dict[rxns[i]]
        vals = np.array(vals,dtype=np.float32)
        ax[2].plot(vals,pres,label=labels[i])
    ax[2].legend()
    ax[2].set_ylabel("MEAC",size='x-large')

    # Loss
    rxns = []
    labels = []
    for i in range(loss+2,loss+7):
        line = lines[i].split()
        rxns.append(line[0])
        labels.append(' '.join(line[1:-1]))
    
    for i in range(5):
        vals = rate_dict[rxns[i]]
        vals = np.array(vals,dtype=np.float32)
        ax[3].plot(vals,pres,label=labels[i])
    ax[3].legend()

    for a in ax:
        a.set_xscale('log')
        a.set_xlim(1e-15,1e5)
        a.set_ylim(1e5,1e-1)
        a.set_yscale('log')

    fig.set_figwidth(12)
    fig.set_figheight(5)
    plt.tight_layout()
    # plt.savefig(f"outputs/{NAME}/prod_loss_rates",dpi=250)
    plt.show()

def pCH4_vs_pCO2():

    fig,ax = plt.subplots()
    ### MEAC
    meac_co2,meac_ch4 = [],[]
    for i in range(1,5):
        data = np.loadtxt(f"hu-code-sr/scenario_library/CO2_CH4/fco2_1e-{i}/ConcentrationSTD.dat",skiprows=2)
        nd_total = np.sum(data[:,5:])
        nd_ch4 = np.sum(data[:,25])
        nd_co2 = np.sum(data[:,56])
        meac_co2.append(nd_co2/nd_total)
        meac_ch4.append(nd_ch4/nd_total)
    
    ### Photochem
    pc_co2,pc_ch4 = [],[]
    for i in range(1,5):
        data = np.loadtxt(f"../../pc/scenarios/fco2_1e-{i}/atmosphere.txt",skiprows=1)
        nd_total = data[:,2]
        mr_co2 = np.sum(data[:,29]*nd_total)/np.sum(nd_total)
        mr_ch4 = np.sum(data[:,36]*nd_total)/np.sum(nd_total)

        pc_co2.append(mr_co2)
        pc_ch4.append(mr_ch4)
    
    ax.plot(meac_co2,meac_ch4,'rs-',label='MEAC')
    ax.plot(pc_co2,pc_ch4,'bs-',label="Photochem")
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r"$f_{\rm CO_{\rm 2}}$",size='x-large')
    ax.set_ylabel(r"$f_{\rm CH_{\rm 4}}$",size='x-large')
    plt.tight_layout()
    plt.legend()
    plt.savefig("figs/fco2_vs_fch4",dpi=250)
    plt.show()

    pass

def H2Ophotolysis():
    fig,ax = plt.subplots()

    ### MEAC H2O photolysis rates
    data = np.loadtxt('hu-code-sr/scenario_library/CO2_CH4/fco2_1e-1/ConcentrationSTD.dat',skiprows=2)
    pres = data[:,4]
    x1 = getMEACcolRates('hu-code-sr/scenario_library/CO2_CH4/fco2_1e-1/',rxn='P6')
    y1 = pres
    ax.plot(x1,y1,label="MEAC")

    ### Photochem H2O photolysis rates
    data = np.loadtxt('../../pc/fco2_1e-1/atmosphere.txt',skiprows=1)
    pres = data[:,1]*1e5
    rates = getPCcolRates('../../pc/fco2_1e-1/',mol='H2O',mode='loss')
    x1 = rates['H2O + hv => OH + H']
    y1 = pres
    ax.plot(x1,y1,label="Photochem")

    ax.set_xlabel("H2O Photolysis Rxn Rate",size='x-large')
    ax.set_ylabel("Pressure [Pa]",size='x-large')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.invert_yaxis()
    plt.legend()
    plt.tight_layout()
    plt.show()

def compareMEACandPClossRates():
    ### MEAC
    totals = []
    R501,R157,R423,P32,P33,P34 = [],[],[],[],[],[]
    for i in range(1,5):
       rxns = getMEACintRates(f"hu-code-sr/scenario_library/CO2_CH4/fco2_1e-{i}/")
       totals.append(sum([items[1] for items in rxns.items()]))
       R501.append(rxns['R501'])
       R423.append(rxns['R423'])
       R157.append(rxns['R157'])
       P32.append(rxns['P32'])
       P33.append(rxns['P33'])
       P34.append(rxns['P34'])
    
    fig,ax = plt.subplots(ncols=2)
    fs = [1e-1,1e-2,1e-3,1e-4]
    ax[0].plot(fs,totals,marker='s',c='black',label="Total")
    ax[0].plot(fs,R501,marker='s',label=r'$OH+CH_4 -> CH_3+H_2O$')
    ax[0].plot(fs,R423,marker='s',label=r"$O+CH_4 -> CH_3 + OH$")
    ax[0].plot(fs,R157,marker='s',label=r"$H+CH_4 -> CH_3 + H_2$")
    ax[0].plot(fs,P32,marker='s',label=r"$CH_4 -> CH_3 + H$")
    ax[0].plot(fs,P33,marker='s',label=r"$CH_4 -> ^*CH_2 + H_2$")
    ax[0].plot(fs,P34,marker='s',label=r"$CH_4 -> CH + H_2 + H$")

    ax[0].set_title("CLIMA-MEAC",size='x-large')
    ax[0].set_xlabel(r"$fCO_2$",size='x-large')
    ax[0].set_ylabel("Column-integrated rxn rate",size='x-large')
    ax[0].set_ylim(1e5,2e10)
    ax[0].set_yscale('log')
    ax[0].set_xscale("log")
    ax[0].legend()

    ### Photochem
    totals = []
    R501,R157,R423,P32,P33,P34 = [],[],[],[],[],[]
    for i in range(1,5):
        rxns = getPCintRates(f'../../pc/scenarios/fco2_1e-{i}/')
        totals.append(rxns['Total'])
        R501.append(rxns["CH4 + OH => CH3 + H2O"])
        R423.append(rxns["CH4 + O => CH3 + OH"])
        R157.append(rxns["CH4 + H => CH3 + H2"])
        P32.append(rxns["CH4 + hv => CH3 + H"])
        P33.append(rxns["CH4 + hv => 1CH2 + H2"])
        P34.append(rxns["CH4 + hv => CH + H2 + H"])
    
    ax[1].plot(fs,totals,marker='s',c='black',label="Total")
    ax[1].plot(fs,R501,marker='s',label=r'$OH+CH_4 -> CH_3+H_2O$')
    ax[1].plot(fs,R423,marker='s',label=r"$O+CH_4 -> CH_3 + OH$")
    ax[1].plot(fs,R157,marker='s',label=r"$H+CH_4 -> CH_3 + H_2$")
    ax[1].plot(fs,P32,marker='s',label=r"$CH_4 -> CH_3 + H$")
    ax[1].plot(fs,P33,marker='s',label=r"$CH_4 -> ^*CH_2 + H_2$")
    ax[1].plot(fs,P34,marker='s',label=r"$CH_4 -> CH + H_2 + H$")

    ax[1].set_title("Photochem",size='x-large')
    ax[1].set_xlabel(r"$fCO_2$",size='x-large')
    ax[1].set_ylabel("Column-integrated rxn rate",size='x-large')
    ax[1].set_ylim(1e5,2e10)
    ax[1].set_yscale('log')
    ax[1].set_xscale("log")
    ax[1].legend()

    fig.set_figwidth(12)
    fig.set_figheight(6)
    plt.tight_layout()
    plt.savefig('figs/rate_comp',dpi=250)
    plt.show()

def compareRReactionRates():
    """
    Docstring for compareRReactionRates
    """

    fig,ax = plt.subplots(ncols=3)
    t = np.linspace(280,320,100)

    # Reaction 1: O + CH4 => CH3 + OH
    rate_m  = (2.26e-12) * (t/298)**2.2 * np.exp(-3820/t)
    rate_pc = (1.08e-21) * t**2.75 * np.exp(-1600/t)
    ax[0].plot(t,rate_m,c='cornflowerblue',label="CLIMA-MEAC")
    ax[0].plot(t,rate_pc,c='salmon',label="Photochem")
    ax[0].legend()
    # ax[0].set_xscale('log')
    ax[0].set_yscale('log')
    ax[0].set_title(r"$O + CH_4 -> CH_3 + OH$",size='x-large')
    ax[0].set_xlabel("Temperature [K]",size='x-large')
    ax[0].set_ylabel("Reaction rate",size='x-large')

    # Reaction 2: H + CH4 => CH3 + H2
    rate_m  = (5.83e-13) * (t/298)**3 * np.exp(-4040/t)
    rate_pc = (1.14e-20) * t**2.74 * np.exp(-4700/t)
    ax[1].plot(t,rate_m,c='cornflowerblue',label="CLIMA-MEAC")
    ax[1].plot(t,rate_pc,c='salmon',label="Photochem")
    ax[1].legend()
    # ax[1].set_xscale('log')
    ax[1].set_title(r"$H + CH_4 -> CH_3 + H_2$",size='x-large')
    ax[1].set_yscale('log')
    ax[1].set_xlabel("Temperature [K]",size='x-large')
    # ax[1].set_ylabel("Reaction rate",size='x-large')

    # Reaction 3: OH + CH4 => CH3 + H2O
    rate_m  = (2.45e-12) * np.exp(-1775/t)
    rate_pc = (2.58e-17) * t**1.83 * np.exp(-1396/t)
    ax[2].plot(t,rate_m,c='cornflowerblue',label="CLIMA-MEAC")
    ax[2].plot(t,rate_pc,c='salmon',label="Photochem")
    ax[2].legend()
    # ax[2].set_xscale('log')
    ax[2].set_title(r"$OH + CH_4 -> CH_3 + H_2O$",size='x-large')
    ax[2].set_yscale('log')
    ax[2].set_xlabel("Temperature [K]",size='x-large')
    # ax[2].set_ylabel("Reaction rate",size='x-large')

    fig.set_figwidth(12)
    plt.tight_layout()
    plt.savefig('figs/nonphot_rates')
    plt.show()

def importPhotochemWaterProfile(root:str = 'fco2_1e-1'):
    """
    This function takes the H2O mixing ratio profile from a specified Photochem
    run and formats it to to serve as a MEAC model fixed profile.
    """


    # Get PC altitudes and H2O mixing ratios
    PCPATH = '../../pc/'
    data = np.loadtxt(PCPATH+'/scenarios/'+root+'/atmosphere.txt',skiprows=1)
    alts_pc = data[:,0]
    ndens_pc= data[:,2]
    h2o_pc  = data[:,24]*ndens_pc
    co2_pc  = data[:,29]*ndens_pc   # not necessary but idk

    # Interpolate values
    alts_new    = np.linspace(1,99,50)
    h2o_new     = np.interp(alts_new,alts_pc,h2o_pc)        # Mixing ratios
    co2_new     = np.interp(alts_new,alts_pc,co2_pc)

    # Write to ConstantMixing.dat
    f = open('hu-code-sr/Data/ConstantMixing.dat','w')
    f.write('z\t\tH2O\t\tCO2\n')
    for i in range(len(alts_new)):
        f.write(np.format_float_positional(alts_new[i],precision=6,min_digits=6)+'\t')
        f.write(np.format_float_scientific(h2o_new[i],precision=6,min_digits=6)+'\t')
        f.write(np.format_float_scientific(co2_new[i],precision=6,min_digits=6)+'\n')
    f.close()


# conc_file = 'hu-code-sr/scenario_library/CO2_CH4/fco2_1e-4/ConcentrationSTD.dat'
# ref_file  = 'hu-code-sr/scenario_library/CO2_CH4/fco2_1e-4_satred/ConcentrationSTD.dat'
# plotAtmosphericComposition(conc_file=conc_file,id='test4',ref_file=ref_file,out_dir='outputs/test')

# plotAtmosphericEvolution(scen_name='co2_1e-6_ch4_1e8',out_dir='outputs/co2_1e-6_ch4_1e8')

