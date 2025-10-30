import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os

def plotAtmosphericComposition(conc_file:str,ref_file:str=None,out_dir:str=''):
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
    ndC2H6  =   data[::-1,35];              mrC2H6  =   np.divide(ndC2H6,nd_all)     # relative
    ndCH4   =   data[::-1,25];              mrCH4   =   np.divide(ndCH4,nd_all)      # relative
    ndCO2   =   data[::-1,56];              mrCO2   =   np.divide(ndCO2,nd_all)     # ABSOLUTE
    ndH2    =   data[::-1,57];              mrH2    =   np.divide(ndH2,nd_all)       # relative
    ndH2O   =   data[::-1,11];              mrH2O   =   np.divide(ndH2O,nd_all)     # ABSOLUTE
    ndN2    =   data[::-1,59];              mrN2    =   np.divide(ndN2,nd_all)       # relative
    ndO2    =   data[::-1,58];              mrO2    =   np.divide(ndO2,nd_all)       # relative      
    ndO3    =   data[::-1,6];               mrO3    =   np.divide(ndO3,nd_all)       # relative
    
    # Writing mixing ratio profiles
    nd_profiles =   [ndC2H6,ndCH4,ndCO2,ndH2,ndH2O,ndN2,ndO2,ndO3]
    mr_profiles =   [mrC2H6,mrCH4,mrCO2,mrH2,mrH2O,mrN2,mrO2,mrO3]
    colors      =   ['yellow','orange','red','purple','cyan','green','cornflowerblue','blue']
    labels      =   ['C2H6','CH4','CO2','H2','H2O','N2','O2','O3']

    # Plot outputs
    fig,ax = plt.subplots(ncols=2)
    for j in range(len(nd_profiles)):
        ax[0].plot(nd_profiles[j],alts,c=colors[j],label=labels[j])
        ax[1].plot(mr_profiles[j],alts,c=colors[j],label=labels[j])
    
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
        ndC2H6  =   data[::-1,35];              mrC2H6  =   np.divide(ndC2H6,nd_all)     # relative
        ndCH4   =   data[::-1,25];              mrCH4   =   np.divide(ndCH4,nd_all)      # relative
        ndCO2   =   data[::-1,56];              mrCO2   =   np.divide(ndCO2,nd_all)     # ABSOLUTE
        ndH2    =   data[::-1,57];              mrH2    =   np.divide(ndH2,nd_all)       # relative
        ndH2O   =   data[::-1,11];              mrH2O   =   np.divide(ndH2O,nd_all)     # ABSOLUTE
        ndN2    =   data[::-1,59];              mrN2    =   np.divide(ndN2,nd_all)       # relative
        ndO2    =   data[::-1,58];              mrO2    =   np.divide(ndO2,nd_all)       # relative      
        ndO3    =   data[::-1,6];               mrO3    =   np.divide(ndO3,nd_all)       # relative
        
        # Writing mixing ratio profiles
        nd_profiles =   [ndC2H6,ndCH4,ndCO2,ndH2,ndH2O,ndN2,ndO2,ndO3]
        mr_profiles =   [mrC2H6,mrCH4,mrCO2,mrH2,mrH2O,mrN2,mrO2,mrO3]
        colors      =   ['yellow','orange','red','purple','cyan','green','cornflowerblue','blue']
        labels      =   ['C2H6','CH4','CO2','H2','H2O','N2','O2','O3']

        # Plot outputs
        for j in range(len(nd_profiles)):
            ax[0].plot(nd_profiles[j],alts,ls=':',c=colors[j])
            ax[1].plot(mr_profiles[j],alts,ls=':',c=colors[j])

    ax[0].set_xlim(left=1e8)
    ax[0].set_ylim(bottom=0,top=max(alts))
    ax[0].set_xscale('log')
    ax[0].set_xlabel("Number density",size='x-large')
    ax[0].set_ylabel("Altitude [km]",size='x-large')
    ax[0].legend()
    ax[1].set_xlim(left=1e-10)
    ax[1].set_ylim(bottom=0,top=max(alts))
    ax[1].set_xscale('log')
    ax[1].set_xlabel("Mixing ratio",size='x-large')
    ax[1].set_ylabel("Altitude [km]",size='x-large')
    ax[1].legend()
    if ref_file: fig.suptitle('New -> solid, ref -> dotted')
    fig.set_figwidth(9)
    plt.tight_layout()
    plt.savefig(f"{out_dir}/mixing ratios",dpi=300)
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
        data = data.reshape(len(data)//11,11)
        alts = data[:,1]

        for j in range(len(labels)):
            cmap = matplotlib.colormaps[colors[j]]
            axes[j].plot(data[:,3+j],alts,c=cmap(i/l))
            axes[j].set_title(labels[j],size='x-large')
        
    for k in range(len(axes)):
        ax = axes[k]
        if k not in [2,5,6]: ax.set_xscale('log')
    fig.set_figwidth(15)
    fig.set_figheight(8)
    plt.tight_layout()
    plt.show()
    
def plotTPprofile(io_dir:str,out_dir:str=''):
    """
    stuff here
    """
    fig,ax = plt.subplots()

    f_in = open(f"{io_dir}/clima_first.tab",'r')
    f_out = open(f"{io_dir}/clima_last.tab",'r')
    data_in = np.fromstring(''.join(f_in.readlines()[1:]),sep=' ').reshape(101,9)
    data_out = np.fromstring(''.join(f_out.readlines()[1:]),sep=' ').reshape(101,9)
    
    ax.plot(data_in[:,2],data_in[:,1],ls=':',c='darkgrey',label="Before")
    ax.plot(data_out[:,2],data_out[:,1],c='black',label='After')

    ax.invert_yaxis()
    plt.yscale('log')
    plt.legend()
    plt.savefig(out_dir+'/TPprofiles',dpi=250)
    # plt.show()
    plt.close()



