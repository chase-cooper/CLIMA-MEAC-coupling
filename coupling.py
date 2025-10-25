import matplotlib.pyplot as plt
import numpy as np
import os
# from scipy.optimize import curve_fit
import subprocess

from compare_model_outputs_v2 import plot_comparison          # modified from original
from make_seed_concSTDfile import *

###########################################################################################################
###########################################################################################################

#######################
###    VARIABLES    ###
#######################

# Coupling parameters
ND              =   101                     # CLIMA layer variables DONT TOUCH
NLOOPS          =   10                       # Number of CLIMA-MEAC loops
NMINSTEPS       =   500                      # Minimum number of CLIMA steps per loop
NMAXSTEPS       =   500                      # Max Number of CLIMA steps per loops
TCONV           =   5E-5                    # CLIMA convergence criterion, delta temperature

# Planet + atmosphere parameters
MASS            =   5.9722e24               # Planet mass in kg
RAD             =   6371.0                  # Planet radius in km
G               =   980                     # grav*M/R**2; Surface gravity, cgs
A               =   1.0                     # Semimajor axis in AU
SURFALB         =   0.25                     # Planet surface albedo
INSTELL         =   1.0                     # Planet instelation, relative to Earth
P0              =   1e-5                    # Top-of-atmosphere pressure [atm]
PSURF           =   1e+0                    # Surface pressure [atm]
T0              =   200                     # Top-of-atmosphere temperature [K]
TSURF           =   280                     # Surface temperature [K]
TROPOPAUSE      =   22                      # CLIMA tropopause layer, default 22 (of 101)
AR              =   2e-2                    # Argon mixing ratio
RELHUM          =   0.7                     # Surface relative humidity
FIXH2O          =   1                       # Fixed H2O flag 
IO3             =   1                       # Ozone flag -- IO3=0 means ozone isn't read in
IME             =   1                       # Methane/ethane flag
ZMAX            =   -1                      # Altitude of highest layer

# Eddysed variables
DOEDDY          =   'true'                  # eddysed flag; if false, nothing below this matters
FCLOUD          =   0.0                     # fractional cloud coverage   
KZ_MIN          =   100000.0                # Minimum eddydiffusivity
CRAINF          =   0.01                    #
CSIG            =   2.0                     #
SUPERSAT        =   0.0                     #
COLDTRAPMINMIX  =   4e-10                   # Minimum coldtrap water mass mixing ratio
FCMINF          =   0.01                    # Minimum fraction for FC for upper atmosphere

# Misc
RESUMERUN       =   0                      # Continue from most recent run
ICONSERVE       =   1                      # CLIMA energy conservation flag, PLEASE KEEP ON
RAINOUT         =   1                      # controls rainout in 'main.c'. not used yet
ADIABATIC       =   0                      # tells CLIMA that atmosphere is dry adiabat  
atm2Pa          =   101_325

#######################
###  Paths & Files  ###
#######################

# Paths
PATH            =   os.getcwd()                                     # Path to this python file
CLIMAPATH       =   f'{PATH}/cloudy_clima'                          # Path to the cloudy-CLIMA folder
MEACPATH        =   f'{PATH}/hu-code-sr'                            # Path to the MEAC folder
MSCENARIONAME   =   "Earth"                                         # Name of MEAC scenario folder

# CLIMA files
CINOUT          =   f'{CLIMAPATH}/CLIMA/IO'
CINCLUDE        =   f'{CLIMAPATH}/CLIMA/INCLUDE/header.inc'         # file with step count, conv. criterion
CINPUT          =   f'{CINOUT}/input_clima.dat'                     # input file
CMIXING         =   f'{CINOUT}/mixing_ratios.dat'                   # mixing ratios file
CLAST           =   f'{CINOUT}/clima_last.tab'                      # output file
CALLOUT         =   f'{CINOUT}/clima_allout.tab'                    # allout file
CTEMPIN         =   f'{CINOUT}/TempIn.dat'                          # temperature and water profile input
CTEMPOUT        =   f'{CINOUT}/TempOut.dat'                         # temperature and water profile output
C_O2            =   f'{CINOUT}/Profiles/O2.dat'                     # O2   vertical profile
C_N2            =   f'{CINOUT}/Profiles/N2.dat'                     # N2   vertical profile
C_H2            =   f'{CINOUT}/Profiles/H2.dat'                     # H2   vertical profile
C_CO2           =   f'{CINOUT}/Profiles/CO2.dat'                    # CO2  vertical profile
C_CH4           =   f'{CINOUT}/Profiles/CH4.dat'                    # CH4  vertical profile
C_C2H6          =   f'{CINOUT}/Profiles/C2H6.dat'                   # C2H6 vertical profile
C_H2O           =   f'{CINOUT}/Profiles/H2O.dat'                    # H2O  vertical profile
C_O3            =   f'{CINOUT}/Profiles/O3.dat'                     # O3   vertical profile

# MEAC files
MSCENARIOPATH   =   f'scenario_library/{MSCENARIONAME}'             # Path to MEAC scenario folder
MCONV           =   f"{MEACPATH}/{MSCENARIOPATH}/ConcentrationSTD_base_Earth.dat"    # Conc. file from last convergence
MSCENARIO       =   f'{MSCENARIOPATH}/planet_Earth_Full_T1986.h'    # MEAC scenario file with planet parameters
MSPECIES        =   f'{MSCENARIOPATH}/species_Earth_Full.dat'       # MEAC atmosphere species file
MCONC           =   f'{MEACPATH}/{MSCENARIOPATH}/ConcentrationSTD.dat'         # MEAC concentrations file
MZTP            =   f'{MSCENARIOPATH}/TP.dat'                       # MEAC ztp profile

#######################
###   Bookkeeping   ###
#######################

CLIMAstepIntervals      = []                   # Step counts where MEAC runs occur
surfTemps     = []               # Surface temperatures after each CLIMA step (not each loop)
toaTemps         = []                  # TOA temperatures after each CLIMA step

##############################################################################################################
##############################################################################################################

# Auxiliary functions

def format_scenario_file():                     # Used in updateMEAC function
    """
    mass            --> {1}
    radius          --> {2}
    semimajor axis  --> {3}
    surface albedo  --> {4}
    MZTP            --> {5}
    ND              --> {6}
    MSPECIES        --> {7}
    MSCENARIONAME   --> {8}
    MSCENARIOPATH   --> {9}
    """
    f = open("templates/meac_scenario.txt",'r')
    template = f.read()
    f.close()

    # Do the writing to the scenario file
    template = template.replace('{1}',str(MASS))
    template = template.replace('{2}',str(RAD*1000))
    template = template.replace('{3}',str(A))
    template = template.replace('{4}',str(SURFALB))
    template = template.replace('{5}',MZTP)
    template = template.replace('{6}',MSPECIES)
    template = template.replace('{7}',MSCENARIOPATH)
    template = template.replace('{8}',str(ZMAX))

    f = open(f"{MEACPATH}/{MSCENARIO}",'w')
    f.write(template)
    f.close()

def format_clima_input():
    f = open('templates/clima_input.txt','r')
    template = f.read()
    f.close()

    # Do the writing
    template = template.replace('{1}',str(NMAXSTEPS))
    template = template.replace('{2}',str(RELHUM))
    template = template.replace('{3}',str(P0))
    template = template.replace('{4}',str(PSURF))
    template = template.replace('{5}',str(G))
    template = template.replace('{6}',str(IO3))
    template = template.replace('{7}',str(ICONSERVE))
    template = template.replace('{8}',str(SURFALB))
    template = template.replace('{9}',str(INSTELL))
    template = template.replace('{10}',"1.0")
    template = template.replace('{11}',str(IME))
    template = template.replace('{12}',str(DOEDDY).lower())
    template = template.replace('{13}',str(FCLOUD))
    template = template.replace('{14}',str(KZ_MIN))
    template = template.replace('{15}',str(CRAINF))
    template = template.replace('{16}',str(CSIG))
    template = template.replace('{17}',str(SUPERSAT))
    template = template.replace('{18}',str(COLDTRAPMINMIX))
    template = template.replace('{19}',str(FCMINF))

    f = open(CINPUT,'w')
    f.write(template)
    f.close()

def saveParameters():
    """
    Write all* initial parameters to a file.
    """
    f = open('outputs/parameters.txt','w')
    f.write("###    MODEL PARAMETERS\n")
    f.write(f"ND                    =   {ND}\n")
    f.write(f"NLOOPS                =   {NLOOPS}\n")
    f.write(f"NMINSTEPS             =   {NMINSTEPS}\n")
    f.write(f"NMAXSTEPS             =   {NMAXSTEPS}\n")
    f.write(f"TCONV                 =   {TCONV}\n")
    f.write(f"Resuming from previous?   {bool(RESUMERUN)}\n\n")

    f.write("### PLANET/ATMOSPHERE PARAMETERS\n")
    f.write(f"Radius [km]           =   {RAD}\n")
    f.write(f"G                     =   {G}\n")
    f.write(f"Separation (AU)       =   {A}\n")
    f.write(f"Surface Albedo        =   {SURFALB}\n")
    f.write(f"Instellation [Solar]  =   {INSTELL}\n")
    f.write(f"TOA Pressure [atm]    =   {P0}\n")
    f.write(f"Surface Pressure[atm] =   {PSURF}\n")
    f.write(f"TOA Temperature [K]   =   {T0}\n")
    f.write(f"Surface Temp [K]      =   {TSURF}\n")
    f.write(f"Tropopause layer      =   {TROPOPAUSE}\n")
    f.write(f"Argon Mixing Ratio    =   {AR}\n")
    f.write(f"Ozone Flag            =   {bool(IO3)}\n")
    f.write(f"Methane Flag          =   {bool(IME)}\n")
    f.write(f"Surface rel.humidity  =   {RELHUM}\n")
    f.write(f"Energy cons. flag     =   {bool(ICONSERVE)}\n")
    f.write(f"Rainout Flag          =   {bool(RAINOUT)}\n")
    f.write(f"Dry Adiabat Flag      =   {bool(ADIABATIC)}\n")
    f.write(f"Fractional cloudiness =   {FCLOUD}\n\n")

    f.write("###    EDDYSED PARAMETERS\n")
    f.write(f"eddysed flag          =   {DOEDDY}\n")
    f.write(f"Fracional cloudiness  =   {str(FCLOUD)}\n")
    f.write(f"Min. eddy diffusivity =   {str(KZ_MIN)}\n")
    f.write(f"CRAINF                =   {str(CRAINF)}\n")
    f.write(f"CSIG                  =   {str(CSIG)}\n")
    f.write(f"SUPERSAT              =   {str(SUPERSAT)}\n")
    f.write(f"Min cold trap fH2O    =   {str(COLDTRAPMINMIX)}\n")
    f.write(f"Min upper atm. FC frac=   {str(FCMINF)}\n\n")

    f.write("###    FILE PATHS\n")
    f.write(f"Python File Path      =   {PATH}\n")
    f.write(f"Cloudy CLIMA Path     =   {CLIMAPATH}\n")
    f.write(f"MEAC Path             =   {MEACPATH}\n\n")

    f.write("###    CLIMA FILE PATHS\n")
    f.write(f"CLIMA IO Folder       =   {CINOUT}\n")
    f.write(f"CLIMA header file     =   {CINCLUDE}\n")
    f.write(f"CLIMA Input File      =   {CINPUT}\n")
    f.write(f"CLIMA Mixing Ratios   =   {CMIXING}\n")
    f.write(f"CLIMA Last Step Output=   {CLAST}\n")
    f.write(f"CLIMA Allout File     =   {CALLOUT}\n")
    f.write(f"CLIMA TempIn File     =   {CTEMPIN}\n")
    f.write(f"CLIMA Tempout File    =   {CTEMPOUT}\n")
    f.write(f"CLIMA Profiles Folder =   {CINOUT}/Profiles\n\n")

    f.write("###    MEAC FILE PATHS\n")
    f.write(f"MEAC Scenario Name    =   {MSCENARIONAME}\n")
    f.write(f"MEAC Scenario Path    =   {MSCENARIOPATH}\n")
    f.write(f"MEAC Converged File   =   {MCONV}\n")
    f.write(f"MEAC Scenario File    =   {MSCENARIO}\n")
    f.write(f"MEAC Species File     =   {MSPECIES}\n")
    f.write(f"MEAC ConcentrationFile=   {MCONC}\n")
    f.write(f"MEAC T-P Profile      =   {MZTP}\n")

    f.close()

def CLIMA_ratio_comp():
    fig,ax = plt.subplots()

    colors = {
        'O2':'cyan',
        'N2':'green',
        'H2':'purple',
        'CO2':'red',
        'CH4':'orange',
        'C2H6':'gold',
        'H2O':'blue',
        'O3':'cornflowerblue'
    }

    # Get altitudes by layer
    f = open(CLAST)
    lines = f.read().split('\n')[1:-1]
    f.close()
    alts = [np.float32(line.split()[0]) for line in lines]
    ax2 = ax.twinx()

    # Pre-CLIMA mixing ratio profiles
    files = os.listdir(f"{CINOUT}/Profiles/")
    for fname in files:
        f = open(f"{CINOUT}/Profiles/{fname}")
        data = f.read().split('\n')[:-1]
        f.close()
        layers = np.linspace(1,len(data),len(data))[::-1]
        c = colors[fname[:-4]]
        ax.plot(np.float32(data),layers,c=c,ls=':')
    
    files = os.listdir(f"{CINOUT}/Profiles_out/")
    for fname in files:
        f = open(f"{CINOUT}/Profiles_out/{fname}")
        data = f.read().split('\n')[:-1]
        f.close()
        layers = np.linspace(1,len(data),len(data))[::-1]
        c = colors[fname[:-4]]
        ax.plot(np.float32(data),layers,c=c,label=fname[:-4])
    
    ax.set_ylim(1,ND)
    ax2.set_ylim(0,max(alts))
    ax2.set_ylabel("Altitude [km]",size='x-large')

    plt.xscale('log')
    ax.set_xlabel('Mixing ratio*',size='x-large')
    ax.set_ylabel('Layer',size='x-large')
    plt.xlim(1e-12,1)
    ax.legend()
    plt.tight_layout()
    plt.savefig('outputs/CLIMA ratios b4 vs after',dpi=200)
    plt.close()

def writeIncludeFile():
    f = open(CINCLUDE,'w')
    f.write(f"""c---------------------------------------------------
c Include file to contain common declarations
c
c JHM, 06-16-06
c--------------------------------------------------
        PARAMETER(ND={ND})
        PARAMETER(RAD={RAD})
        PARAMETER(TCONV={TCONV})
        PARAMETER(MINSTEPS={NMINSTEPS})
        PARAMETER(FIXH20={FIXH2O})
        implicit real*8(A-H,O-Z)""")
    f.close()

##############################################################################################################
##############################################################################################################

def updateCLIMA(firstloop:bool):
    """
    Trying this
    """
    debug = open('outputs/debug.txt','w')   # Debug log

    # Getting mixing ratios from MEAC    
    if firstloop:
        debug = open('outputs/debug_MCONV.txt','w')   # Debug log                   
        file = open(MCONV,'r')
    else:
        debug = open('outputs/debug_MCONC.txt','w')   # Debug log
        file = open(MCONC,'r')
    
    data = file.read().replace('#','').split()[121:]
    data = np.asarray(data,dtype=np.float32).reshape(len(data)//116,116)
    file.close()
    debug.write("Imported data array size: " + str(np.shape(data))+'\n\n')

    # TOA to surface
    alts = data[::-1,0]
    pressures = data[::-1,4]

    # Total number densities
    nd_all  =   np.sum(data[::-1],axis=1)
    nd_nc   =   nd_all - data[::-1,56]- data[::-1,11]   # Non-condensible number densities
    
    # number densities                      mixing ratios
    ndC2H6  =   data[::-1,35];              mrC2H6  =   np.divide(ndC2H6,nd_nc)     # relative
    ndCH4   =   data[::-1,25];              mrCH4   =   np.divide(ndCH4,nd_nc)      # relative
    ndCO2   =   data[::-1,56];              mrCO2   =   np.divide(ndCO2,nd_all)     # ABSOLUTE
    ndH2    =   data[::-1,57];              mrH2    =   np.divide(ndH2,nd_nc)       # relative
    ndH2O   =   data[::-1,11];              mrH2O   =   np.divide(ndH2O,nd_all)     # ABSOLUTE
    ndN2    =   data[::-1,59];              mrN2    =   np.divide(ndN2,nd_nc)       # relative
    ndO2    =   data[::-1,58];              mrO2    =   np.divide(ndO2,nd_nc)       # relative      
    ndO3    =   data[::-1,6];               mrO3    =   np.divide(ndO3,nd_nc)       # relative
    
    fH2O    =   []      # For later
    # Writing mixing ratio profiles
    mr_profiles =   [mrC2H6,mrCH4,mrCO2,mrH2,mrH2O,mrN2,mrO2,mrO3]
    mr_files    =   [C_C2H6,C_CH4,C_CO2,C_H2,C_H2O,C_N2,C_O2,C_O3]

    # Get pressure range
    # Don't think this is the issue
    if firstloop:
        pressure_range = np.logspace(np.log10(P0*atm2Pa),np.log10(PSURF*atm2Pa),ND)
        debug.write("Calculated pressure range manually\n")
    else:
        # Pull from clima_last otherwise
        f = open(CLAST)
        data = f.read().split()[9:]
        data = np.array(data,dtype=np.float32).reshape(ND,9)
        pressure_range = data[:,1]*atm2Pa
        f.close()
        debug.write("Copied pressure range from clima_last.tab\n")
    for i in range(len(pressure_range)):
        debug.write("Layer: "+str(i+1)+"\tPressure [Pa]: " + str(pressure_range[i])+'\n')

    # Write mixing ratio profiles
    for i in range(len(mr_files)):
        debug.write("\n\nWriting to "+mr_files[i]+':\n')
        file = open(mr_files[i],'w')
        profile = mr_profiles[i]
        mr = np.interp(pressure_range,pressures,profile)

        # !!! Capping per-layer ozone to 1e-5
        excessO3flag = 0
        for j in range(ND):
            if (i==7) and (mr[j] > 1e-5):
                excessO3flag += 1
                val = np.format_float_scientific(1e-5,precision=15,trim='k',unique=True,exp_digits=2,min_digits=15)
            else:
                val = np.format_float_scientific(mr[j],precision=15,trim='k',unique=True,exp_digits=2,min_digits=15)
            if (i==4): # H2O
                fH2O.append(val)
            file.write(val+'\n')
            debug.write(val+'\n')

        if excessO3flag and (i==7):
            debug.write(f"WARNING: excess O3 detected in {excessO3flag} layer(s)\n")
            
        debug.write("Mixing ratio profile complete!")
        file.close()
    
    # whole-atmosphere mixing ratios for clima input
    fC2H6   = (1-AR) * np.sum(ndC2H6)/np.sum(nd_nc)  # relative
    fCH4    = (1-AR) * np.sum(ndCH4)/np.sum(nd_nc)   # relative
    fCO2    = (1-AR) * np.sum(ndCO2)/np.sum(nd_all)  # ABSOLUTE
    fH2     = (1-AR) * np.sum(ndH2)/np.sum(nd_nc)    # relative
    # fH2O    = (1-AR) * np.sum(ndH2O)/np.sum(nd_all)  # ABSOLUTE, not used
    fN2     = (1-AR) * np.sum(ndN2)/np.sum(nd_nc)    # relative
    fO2     = (1-AR) * np.sum(ndO2)/np.sum(nd_nc)    # relative
    # fO3     = (1-AR) * np.sum(ndO3)/np.sum(nd_nc)    # relative, not used

    # Writing to 'mixing_ratios.dat'                       
    text = f"""\
{np.format_float_scientific(AR,precision=3,trim='k',unique=True,exp_digits=2,min_digits=3)}     ! Argon
{np.format_float_scientific(fCH4,precision=3,trim='k',unique=True,exp_digits=2,min_digits=3)}   ! Methane
{np.format_float_scientific(fC2H6,precision=3,trim='k',unique=True,exp_digits=2,min_digits=3)}  ! Ethane
{np.format_float_scientific(fCO2,precision=3,trim='k',unique=True,exp_digits=2,min_digits=3)}   ! Carbon Dioxide
{np.format_float_scientific(fN2,precision=3,trim='k',unique=True,exp_digits=2,min_digits=3)}    ! Nitrogen
{np.format_float_scientific(fO2,precision=3,trim='k',unique=True,exp_digits=2,min_digits=3)}    ! Oxygen
{np.format_float_scientific(fH2,precision=3,trim='k',unique=True,exp_digits=2,min_digits=3)}    ! Hydrogen
1.000e-72          !Nitrogen Dioxide
{TROPOPAUSE}                 !Tropopause layer
"""
    f = open(CMIXING,"w")
    f.write(text)
    f.close()

    # The rest of necessary CLIMA input
    format_clima_input()

    # Writing TempIn.dat
    # Doesn't seem bad either...
    temps = []
    if firstloop and not RESUMERUN:                         
        l = np.linspace(T0,TSURF,ND)
        for val in l:
            temps.append(np.format_float_positional(val,precision=12,trim='k',unique=True,min_digits=12))
    else:
        # Copy CLIMA TempOut to TempIn
        out = open(CTEMPOUT,'r')
        lines = out.read().split('\n')
        out.close()
        for i in range(ND):
            line = lines[i].split()
            temps.append(line[0])

    inn = open(CTEMPIN,'w')
    for i in range(ND):
        inn.write(' '*3)
        inn.write(temps[i])
        inn.write(' '*8)
        inn.write(fH2O[i]+'\n')
        print(fH2O[i])
    inn.close()

    # Plot outputs
    fig,ax = plt.subplots(ncols=2)
    ax[0].plot(ndC2H6,alts,c='yellow',label='C2H6')
    ax[0].plot(ndCH4,alts,c='orange',label='CH4')
    ax[0].plot(ndCO2,alts,c='red',label='CO2')
    ax[0].plot(ndH2,alts,c='purple',label='H2')
    ax[0].plot(ndH2O,alts,c='cyan',label='H2O')
    ax[0].plot(ndN2,alts,c='green',label='N2')
    ax[0].plot(ndO2,alts,c='cornflowerblue',label='O2')
    ax[0].plot(ndO3,alts,c='blue',label='O3')
    ax[0].set_xlim(left=1e10)
    ax[0].set_xscale('log')
    ax[0].set_xlabel("Number density",size='x-large')
    ax[0].set_ylabel("Altitude [km]",size='x-large')
    ax[0].legend()

    ax[1].plot(mrC2H6,alts,c='yellow',label='C2H6')
    ax[1].plot(mrCH4,alts,c='orange',label='CH4')
    ax[1].plot(mrCO2,alts,c='red',label='CO2')
    ax[1].plot(mrH2,alts,c='purple',label='H2')
    ax[1].plot(mrH2O,alts,c='cyan',label='H2O')
    ax[1].plot(mrN2,alts,c='green',label='N2')
    ax[1].plot(mrO2,alts,c='cornflowerblue',label='O2')
    ax[1].plot(mrO3,alts,c='blue',label='O3')
    ax[1].set_xlim(left=1e-10)
    ax[1].set_xscale('log')
    ax[1].set_xlabel("Mixing ratio",size='x-large')
    ax[1].set_ylabel("Altitude [km]",size='x-large')
    ax[1].legend()

    fig.set_figwidth(9)
    plt.tight_layout()
    plt.savefig("outputs/imported mixing ratios",dpi=300)

    debug.close()

def runCLIMA():
    atm2Pa = 101_325

    # Compile and run the cloudy-CLIMA model.
    os.chdir(CLIMAPATH)
    subprocess.run(["make","clean","-f","ClimaMake"])
    subprocess.run(["make","-f","ClimaMake"])
    subprocess.run("./clima.run")
    os.chdir(PATH)

    # Compare profiles before vs after CLIMA
    CLIMA_ratio_comp()

    ### SANITY CHECK -- TempIn vs TempOut ###
    fig,axes = plt.subplots(ncols=2)

    f = open(f'{CINOUT}/Profiles/H2O.dat','r')
    h2o_in = f.read().split('\n')[:-1]
    h2o_in = np.array(h2o_in,dtype=float)
    f.close()

    temps_out,p_out,h2o_out = [],[],[]
    f = open(CLAST)
    data = f.read().split()[9:]
    data = np.array(data,dtype=np.float32).reshape(ND,9)
    f.close()
    temps_out = data[:,2]
    p_out = data[:,1] * atm2Pa
    h2o_out = data[:,3]

    axes[0].plot(temps_out,p_out,color='black')
    axes[1].plot(h2o_in,p_out,c='blue',ls=':')
    axes[1].plot(h2o_out,p_out,color='blue')

    axes[0].set_title("Temp-Pressure Profile",size='x-large')
    axes[0].set_xlabel("Temperature [K]",size='large')
    axes[1].set_title("Water Vapor Profile",size='x-large')
    axes[1].set_xlabel("H2O Mixing Ratio",size='large')
    axes[1].set_xscale('log')
    for a in axes:
        a.set_ylabel("Pressure [Pa]",size='large')
        a.set_ylim(bottom=min(p_out),top=max(p_out))
        a.set(yscale='log')
        a.invert_yaxis()
    fig.set_figwidth(7)
    plt.tight_layout()
    plt.savefig('outputs/temp_h2o_profiles',dpi=250)
    # plt.show()
    plt.close()

    ### Plotting water quantities and the like
    fig,ax = plt.subplots()
    ax.plot(h2o_out,p_out,color='blue',label='H2O mixing ratio')
    
    # fsatur
    f = open(f'{CINOUT}/extras/fsatur.dat','r')
    data = f.read().split('\n')[:-1]
    f.close()
    data = np.array(data,dtype=np.float32)
    ax.plot(data,p_out,color='cornflowerblue',ls=':',label='fsatur(J)')

    # relhum (Relative Humidity)
    f = open(f'{CINOUT}/extras/relhum.dat','r')
    data = f.read().split('\n')[:-1]
    f.close()
    data = np.array(data,dtype=float)
    ax.plot(data,p_out,color='darkgrey',ls='--',label='RELHUM(J)')

    # relhum_vec (Relative Humidity but different)
    f = open(f'{CINOUT}/extras/relhumvec.dat','r')
    data = f.read().split('\n')[:-1]
    f.close()
    data = np.array(data,dtype=float)
    ax.plot(data,p_out,color='slategrey',label='relhum_vec(J)')
    plt.xscale('log')
    plt.yscale('log')
    ax.invert_yaxis()
    plt.legend()
    plt.savefig('outputs/water quantities',dpi=250)
    # plt.show()
    plt.close()

def updateMEAC(firstloop:bool):
    """
    Write MEAC input files, and generate new ConcentrationSTD.dat file (code thanks to Sukrit)
    """

    # Open and read zTP profile generated by CLIMA
    f = open(CLAST)
    data = f.read().split()[9:]
    data = np.array(data,dtype=np.float32).reshape(ND,9)
    f.close()
    
    # ZTP values (values in CLAST go from TOA to surface, so they get flipped)
    p = data[:,1][::-1]*atm2Pa             # atm to Pa                          
    z = data[:,0][::-1]
    t = data[:,2][::-1]
    zmax = int(z[-1])
    global ZMAX; ZMAX=zmax

    ### SANITY CHECK -- Surface Temperature over time ###
    fig,axes = plt.subplots(ncols=2)

    # updating temperature records
    f = open(f"{CINOUT}/extras/surftemp.dat",'r')   # new surface temp values
    newSurfTemps = f.read().split('\n')[:-1]
    newSurfTemps = np.array(newSurfTemps,dtype=np.float32)
    for te in newSurfTemps: surfTemps.append(te)
    f.close()
    ff = open("outputs/surftemps.dat",'w')
    for val in surfTemps:
        ff.write(str(val)+'\n')
    ff.close()

    g = open(f"{CINOUT}/extras/toatemp.dat",'r')    # new TOA temp values
    newTOATemps = g.read().split('\n')[:-1]
    newTOATemps = np.array(newTOATemps,dtype=np.float32)
    for te in newTOATemps: toaTemps.append(te)
    g.close()
    gg = open("outputs/toatemps.dat",'w')
    for val in toaTemps:
        gg.write(str(val)+'\n')
    gg.close()

    # plotting temperatures
    n = np.linspace(1,len(surfTemps),len(surfTemps))
    axes[0].plot(n,surfTemps,color='limegreen')
    axes[1].plot(n,toaTemps,color='palevioletred')

    # plotting temp differences
    ax2s = axes[0].twinx()
    ax2t = axes[1].twinx()
    n = np.linspace(2,len(surfTemps),len(surfTemps)-1)
    ysurf,ytoa = [],[]
    for i in range(1,len(n)+1):
        diff_surf = (surfTemps[i]-surfTemps[i-1])/surfTemps[i-1]
        diff_toa = (toaTemps[i]-toaTemps[i-1])/toaTemps[i-1]
        ysurf.append(diff_surf*100)
        ytoa.append(diff_toa*100)
    ax2s.plot(n,ysurf,c='darkgreen',ls='-.')
    ax2t.plot(n,ytoa,c='pink',ls='-.')

    # Plotting where coupling occurs
    CLIMAstepIntervals.append(len(surfTemps)+0.5)
    axes[0].vlines(CLIMAstepIntervals,ymin=np.min(surfTemps)-1,\
               ymax=np.max(surfTemps)+1,colors='grey',linestyles=':')
    axes[1].vlines(CLIMAstepIntervals,ymin=np.min(toaTemps)-1,\
               ymax=np.max(toaTemps)+1,colors='grey',linestyles=':')
    
    # Plot params
    axes[0].set_xlim(left=1)
    axes[0].set_ylim(np.min(surfTemps)-1,np.max(surfTemps)+1)
    axes[0].set_xlabel("CLIMA step")
    axes[0].set_ylabel("Temperature")
    ax2s.set_ylabel("% Temp Change")
    axes[1].set_xlim(left=1)
    axes[1].set_ylim(np.min(toaTemps)-1,np.max(toaTemps)+1)
    axes[1].set_xlabel("CLIMA step")
    axes[1].set_ylabel("Temperature")
    ax2t.set_ylabel("% Temp Change")
    fig.set_figwidth(9)
    plt.tight_layout()
    plt.savefig("outputs/surface temperatures",dpi=300)
    plt.close()
    # input()

    # Update global variables (to be written back into input_clima)
    global P0; P0 = p[-1]/atm2Pa           
    global PSURF; PSURF = p[0]/atm2Pa       
    global TSURF; TSURF = surfTemps[-1]           

    # Write temperature-pressure profile
    f = open(f"{MEACPATH}/{MZTP}","w")
    for i in np.linspace(0,zmax,ND):
        f.write(f"{i:.6f} {np.interp(i,z,np.log10(p)):.6f} {np.interp(i,z,t):.6f}\n")
    f.close()

    # Update scenario file using a format file
    format_scenario_file()

    # Update concentration file
    if firstloop:
        old_conc = MCONV
    else:
        old_conc = MCONC

    newzTP = f"{MEACPATH}/{MZTP}"
    zmin,zmax,zbin = 0,zmax,ND
    new_conc = MCONC
    
    generate_new_concentrationSTD(old_conc,newzTP,zmin,zmax,zbin,new_conc)

def runMEAC():
    #Compile and run the MEAC model.
    os.chdir(MEACPATH)                                      # Move to MEAC folder
    subprocess.run(["gcc","-o","main","main.c"])            # Compile main.c
    os.chmod('main',0b111101101)                            # Ensure main has permissions to be executed
    subprocess.run('./main')                                # Run main
    os.chdir(PATH)                                          # Return to PATH

    ### SANITY CHECK -- Plot old (converged) and new concentration profiles ###
    # plot_comparison(MCONV,MCONC,saveas='Sanity Checks/runMEAC converged vs new conc')

##############################

saveParameters()
writeIncludeFile()
for i in range(NLOOPS):
    first = (i==0)
    updateCLIMA(first)
    # input()
    runCLIMA()
    # input()
    updateMEAC(first)
    # input()
    runMEAC()
    # input()
    os.system('clear')
