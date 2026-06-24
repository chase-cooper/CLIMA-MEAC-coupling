import matplotlib.pyplot as plt
import numpy as np
import os
from shutil import rmtree
import subprocess
import sys
import time

from compare_model_outputs_v2 import *
from files import *
from make_seed_concSTDfile import *
from parameters import *
from plots import *

#######################
###   Bookkeeping   ###
#######################

CLIMAstepIntervals      = []     # Step counts where MEAC runs occur
surfTemps     = []               # Surface temperatures after each CLIMA step (not each loop)
zmax = ZMAX

########################
### Elements of main ###
########################

def updateCLIMA(stepnum:int):
    """
    Prep CLIMA input files using the outputs of the previous MEAC run (if applicable).
    """
    firstloop = (stepnum==0)

    ### Getting mixing ratios from MEAC    
    # Open concentration file
    conc_file = MBASE if firstloop else MCONC
    f = open(conc_file,'r')
    data = f.read().replace('#','').split()[121:]
    data = np.asarray(data,dtype=np.float32).reshape(len(data)//116,116)
    f.close()

    # Arrange number densities from TOA to surface
    pressures = data[::-1,4]    
    nd_all  =   np.sum(data[::-1],axis=1)               # Total number densities
    nd_nc   =   nd_all - data[::-1,56]- data[::-1,11]   # Non-condensible number densities
    
    # number densities                      mixing ratios 
    ndC2H6  =   data[::-1,35];              mrC2H6  =   np.divide(ndC2H6,nd_nc)     # relative to all non-condensibles
    ndCH4   =   data[::-1,25];              mrCH4   =   np.divide(ndCH4,nd_nc)      # relative to all non-condensibles
    ndCO2   =   data[::-1,56];              mrCO2   =   np.divide(ndCO2,nd_all)     # ABSOLUTE
    ndH2    =   data[::-1,57];              mrH2    =   np.divide(ndH2,nd_nc)       # relative to all non-condensibles
    ndH2O   =   data[::-1,11];              mrH2O   =   np.divide(ndH2O,nd_all)     # not needed, H2O read in thru TempIn.dat
    ndN2    =   data[::-1,59];              mrN2    =   np.divide(ndN2,nd_nc)       # relative to all non-condensibles
    ndO2    =   data[::-1,58];              mrO2    =   np.divide(ndO2,nd_nc)       # relative to all non-condensibles      
    ndO3    =   data[::-1,6];               mrO3    =   np.divide(ndO3,nd_nc)       # relative to all non-condensibles

    # Writing mixing ratio(*) profiles
    mr_profiles =   [mrC2H6,mrCH4,mrCO2,mrH2,mrH2O,mrN2,mrO2,mrO3]
    mr_files    =   [C_C2H6,C_CH4,C_CO2,C_H2,C_H2O,C_N2,C_O2,C_O3]

    ### Get pressure range
    if firstloop:
        # If first loop, interpolate over given pressure boundary values
        pressure_range = np.logspace(np.log10(P0*atm2Pa),np.log10(PSURF*atm2Pa),ND)
    else:
        # Pull from clima_last otherwise
        f = open(CLAST)
        data = f.read().split()[9:]
        data = np.array(data,dtype=np.float32).reshape(ND,9)
        pressure_range = data[:,1]*atm2Pa
        f.close()

    ### Write mixing ratio profiles
    fH2O = []
    for i in range(len(mr_files)):
        file = writeOrCreate(mr_files[i])
        profile = mr_profiles[i]
        mr = np.interp(pressure_range,pressures,profile)

        for j in range(ND):
            if (i==7) and (mr[j] > 1e-5):
                val = np.format_float_scientific(1e-5,precision=15,trim='k',unique=True,exp_digits=2,min_digits=15)
            else:
                val = np.format_float_scientific(mr[j],precision=15,trim='k',unique=True,exp_digits=2,min_digits=15)
            if (i==4): # H2O
                fH2O.append(val)
            file.write(val+'\n')
        file.close()
    
    ### whole-atmosphere mixing ratios for clima input
    fC2H6   = (1-AR) * np.sum(ndC2H6)/np.sum(nd_nc)  # relative
    fCH4    = (1-AR) * np.sum(ndCH4)/np.sum(nd_nc)   # relative
    fCO2    = (1-AR) * np.sum(ndCO2)/np.sum(nd_all)  # ABSOLUTE
    fH2     = (1-AR) * np.sum(ndH2)/np.sum(nd_nc)    # relative
    fN2     = (1-AR) * np.sum(ndN2)/np.sum(nd_nc)    # relative
    fO2     = (1-AR) * np.sum(ndO2)/np.sum(nd_nc)    # relative
    # Writing to 'mixing_ratios.dat'                       
    text = f"""\
{np.format_float_scientific(AR,precision=3,trim='k',unique=True,exp_digits=2,min_digits=3)}     ! Argon
{np.format_float_scientific(fCH4,precision=3,trim='k',unique=True,exp_digits=2,min_digits=3)}     ! Methane
{np.format_float_scientific(fC2H6,precision=3,trim='k',unique=True,exp_digits=2,min_digits=3)}     ! Ethane
{np.format_float_scientific(fCO2,precision=3,trim='k',unique=True,exp_digits=2,min_digits=3)}     ! Carbon Dioxide
{np.format_float_scientific(fN2,precision=3,trim='k',unique=True,exp_digits=2,min_digits=3)}     ! Nitrogen
{np.format_float_scientific(fO2,precision=3,trim='k',unique=True,exp_digits=2,min_digits=3)}     ! Oxygen
{np.format_float_scientific(fH2,precision=3,trim='k',unique=True,exp_digits=2,min_digits=3)}     ! Hydrogen
1.000e-72     !Nitrogen Dioxide
{TROPOPAUSE}     !Tropopause layer
"""
    f = writeOrCreate(CMIXING)
    f.write(text)
    f.close()
    # The rest of necessary CLIMA input
    writeCLIMAinput(firstloop)

    ### Writing TempIn.dat
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

    inn = writeOrCreate(CTEMPIN)
    for i in range(ND):
        inn.write(' '*3)
        inn.write(temps[i])
        inn.write(' '*8)
        inn.write(fH2O[i]+'\n')
    inn.close()

    # Copy input files to outputs folder for bookkeeping
    subprocess.run(['cp',CINPUT,f"{OUTPUT}/clima-in/input_clima_{stepnum}.dat"])    # input_clima.dat
    subprocess.run(['cp',CMIXING,f"{OUTPUT}/clima-in/mixing_ratios_{stepnum}.dat"]) # mixing_ratios.dat

def runCLIMA(stepnum:int):
    # Compile and run the cloudy-CLIMA model.
    os.chdir(CLIMAPATH)
    subprocess.run(["make","clean","-f","ClimaMake"])
    subprocess.run(["make","-f","ClimaMake"])
    subprocess.run("./clima.run")
    os.chdir(PATH)

    # Save formatted outputs
    writeCLIMAout(CLAST,OUTPUT,str(stepnum))

def updateMEAC(stepnum:int):
    """
    Write MEAC input files, and generate new ConcentrationSTD.dat file (code thanks to Sukrit)
    """

    firstloop = (stepnum==0)

    # Open and read zTP profile generated by CLIMA
    f = open(CLAST)
    data = f.read().split()[9:]
    data = np.array(data,dtype=np.float32).reshape(ND,9)
    f.close()
    z = data[:,0][::-1]             # values in CLAST go from TOA to surface, so they get flipped
    p = data[:,1][::-1]*atm2Pa      # atm to Pa                          
    t = data[:,2][::-1]
    zmax = z[-1]
    print("zmax = ",zmax,' km')

    # Add each CLIMA step's surface temperature and add it to the log of surf. temperature values.
    f = open(f"{CINOUT}/extras/surftemp.dat",'r')   
    newSurfTemps = f.read().split('\n')[:-1]
    newSurfTemps = np.array(newSurfTemps,dtype=np.float32)
    for te in newSurfTemps: surfTemps.append(te)
    f.close()
    ff = writeOrCreate(f"{OUTPUT}/surftemps.dat")
    for val in surfTemps:
        ff.write(str(val)+'\n')
    ff.close()

    # Update global variables (to be written back into input_clima)
    global P0; P0 = p[-1]/atm2Pa           
    global PSURF; PSURF = p[0]/atm2Pa       
    global TSURF; TSURF = surfTemps[-1]   
    
    # Write temperature-pressure profile
    f = writeOrCreate(f"{MEACPATH}/{MZTP}")
    # First, interpolate the CLIMA ztp over the MEAC altitude layers
    for i in np.linspace(0,zmax,ND):
        f.write(f"{i:.6f} {np.interp(i,z,np.log10(p)):.6f} {np.interp(i,z,t):.6f}\n")
    # For altitudes above the highest CLIMA altitude, assume an isotherm
    i = zmax
    pres = np.interp(i,z,np.log10(p))
    delta_p = np.interp(i,z,np.log10(p)) - np.interp(i-zmax/100,z,np.log10(p))
    print("delta_p: ",delta_p)
    while i < 100:
        i += zmax/100
        pres += delta_p
        f.write(f"{i:.6f} {pres:.6f} {np.interp(i,z,t):.6f}\n")
    f.close()

    # write new concentration file, and save in outputs folder
    if firstloop:
        generate_new_concentrationSTD(MBASE,f"{MEACPATH}/{MZTP}",0,100,NBIN,MCONC)
    else:
        generate_new_concentrationSTD(MCONC,f"{MEACPATH}/{MZTP}",0,100,NBIN,MCONC)
    subprocess.run(['cp',MCONC,f"{OUTPUT}/meac-in/conc_{stepnum}.dat"])

    # Update scenario file with CLIMA zTP
    writeScenarioFile()       

    # Update species scenario file
    writeMEACspecies(tsurf=t[0])      

    # Plot surface temperature evolution
    plotSurfaceTemperature(f"{OUTPUT}/surftemps.dat",runBreaks=CLIMAstepIntervals,out_dir=OUTPUT)

def runMEAC(stepnum:int):
    
    writeMEACmain()                 # Update MEAC main.c to include scenario file

    #Compile and run the MEAC model.
    os.chdir(MEACPATH)                                      # Move to MEAC folder
    subprocess.run(["gcc","-o","main","main.c"])            # Compile main.c
    os.chmod('main',0b111101101)                            # Ensure main has permissions to be executed
    subprocess.run('./main')                                # Run main
    os.chdir(PATH)                                          # Return to PATH

    writeMEACout(MCONC,NAME,str(stepnum))
    plotAtmosphericComposition(conc_file=MCONC,id=str(stepnum),out_dir=OUTPUT)

########################
###       main       ###
########################

def main(name=None):

    start = time.time()

    warning = input(f"Warning: if an output folder with the name '{NAME}' exists, it will be overwritten. Press enter to continue.\n")

    # Clear any output subdirectory with the name NAME, then populate it
    if (NAME in os.listdir('outputs')):
        rmtree(OUTPUT)
    os.mkdir(OUTPUT)
    os.chdir(OUTPUT)
    os.mkdir('clima-in')
    os.mkdir('clima-out')
    os.mkdir('meac-in')
    os.mkdir('meac-out')
    os.mkdir('mr')
    os.chdir(PATH)

    # Make scenario folder for this run, if it doesn't exist yet
    if not (NAME in os.listdir(f"{MEACPATH}/scenario_library/")):
        scen_path = f"{MEACPATH}/scenario_library/{NAME}"
        os.mkdir(scen_path)
        # os.system(f'cp {PATH}/templates/Concentration_STD_Earth.dat {PATH}/{MBASE}')    # Not wokring

        # Update water vapor lower boundary condition
        writeMEACspecies(tsurf=TSURF)    

        # Write a flat Kzz profile. You can change this file later
        f = open(f'{scen_path}/Eddy.dat','w')
        f.write('0.000000 100000\n100.000000 100000\n')
        f.close()

        # os.chdir(PATH)

    writeParameters()
    writeIncludeFile()
    for i in range(NLOOPS):
        updateCLIMA(i)
        runCLIMA(i)
        updateMEAC(i)
        runMEAC(i)
        os.system('clear')
    
    end = time.time()
    print(f"Start:      {start}")
    print(f"End:        {end}")
    print(f"Duration:   {(end-start)/60} minutes")

main(NAME)
