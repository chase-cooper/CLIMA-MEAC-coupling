import matplotlib.pyplot as plt
import numpy as np
import os
import subprocess

from compare_model_outputs_v2 import *
from files import *
from make_seed_concSTDfile import *
from parameters import *
from plots import *

#######################
###   Bookkeeping   ###
#######################

CLIMAstepIntervals      = []                   # Step counts where MEAC runs occur
surfTemps     = []               # Surface temperatures after each CLIMA step (not each loop)
zmax = ZMAX

##############################################################################################################
##############################################################################################################

def writeScenarioFile(zmax:str='100.0'):
    f = open('templates/meac_scenario.txt','r')
    text = f.read()
    f.close()

    text = text.replace('{1}',str(SURFALB))
    text = text.replace('{2}',MZTP)
    text = text.replace('{3}',zmax)
    text = text.replace('{4}',str(NMAXT))
    text = text.replace('{5}',MSPECIES)
    text = text.replace('{6}',MSCENARIOPATH)
    text = text.replace('{7}',str(NBIN))

    o = open(MEACPATH + '/' + MSCENARIO,'w')
    o.write(text)
    o.close()

    print("ENSURE # OF BINS IS PROPERLY SET")

def updateCLIMA(firstloop:bool):
    """
    
    """
    subprocess.run(['cp',CLAST,f"{CINOUT}/clima_first.tab"])
    ### Getting mixing ratios from MEAC    
    # Open concentration file
    conc_file = MCONV if firstloop else MCONC
    f = open(conc_file,'r')
    data = f.read().replace('#','').split()[121:]
    data = np.asarray(data,dtype=np.float32).reshape(len(data)//116,116)
    f.close()

    # TOA to surface
    pressures = data[::-1,4]    
    nd_all  =   np.sum(data[::-1],axis=1)               # Total number densities
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

    ### Get pressure range
    if firstloop:
        pressure_range = np.logspace(np.log10(P0*atm2Pa),np.log10(PSURF*atm2Pa),ND)
    else:
        # Pull from clima_last otherwise
        f = open(CLAST)
        data = f.read().split()[9:]
        data = np.array(data,dtype=np.float32).reshape(ND,9)
        pressure_range = data[:,1]*atm2Pa
        f.close()

    ### Write mixing ratio profiles
    for i in range(len(mr_files)):
        file = writeOrCreate(mr_files[i])
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
{np.format_float_scientific(fCH4,precision=3,trim='k',unique=True,exp_digits=2,min_digits=3)}   ! Methane
{np.format_float_scientific(fC2H6,precision=3,trim='k',unique=True,exp_digits=2,min_digits=3)}  ! Ethane
{np.format_float_scientific(fCO2,precision=3,trim='k',unique=True,exp_digits=2,min_digits=3)}   ! Carbon Dioxide
{np.format_float_scientific(fN2,precision=3,trim='k',unique=True,exp_digits=2,min_digits=3)}    ! Nitrogen
{np.format_float_scientific(fO2,precision=3,trim='k',unique=True,exp_digits=2,min_digits=3)}    ! Oxygen
{np.format_float_scientific(fH2,precision=3,trim='k',unique=True,exp_digits=2,min_digits=3)}    ! Hydrogen
1.000e-72          !Nitrogen Dioxide
{TROPOPAUSE}                 !Tropopause layer
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

def runCLIMA(stepnum:int):
    # Compile and run the cloudy-CLIMA model.
    os.chdir(CLIMAPATH)
    subprocess.run(["make","clean","-f","ClimaMake"])
    subprocess.run(["make","-f","ClimaMake"])
    subprocess.run("./clima.run")
    os.chdir(PATH)

    writeCLIMAout(CLAST,OUTPUT,str(stepnum))
    plotTPprofile(CLAST,str(stepnum),MCONV,OUTPUT)

def updateMEAC():
    """
    Write MEAC input files, and generate new ConcentrationSTD.dat file (code thanks to Sukrit)
    """
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

    # updating temperature records
    f = open(f"{CINOUT}/extras/surftemp.dat",'r')   # new surface temp values
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
    for i in np.linspace(0,zmax,ND):
        f.write(f"{i:.6f} {np.interp(i,z,np.log10(p)):.6f} {np.interp(i,z,t):.6f}\n")
    
    i = zmax
    pres = np.interp(i,z,np.log10(p))
    delta_p = np.interp(i,z,np.log10(p)) - np.interp(i-zmax/100,z,np.log10(p))
    print("delta_p: ",delta_p)
    while i < 100:
        i += zmax/100
        pres += delta_p
        f.write(f"{i:.6f} {pres:.6f} {np.interp(i,z,t):.6f}\n")
    f.close()

    # Update scenario file with CLIMA zTP
    print(zmax)
    val = np.format_float_positional(zmax,precision=1,trim='k',unique=True,min_digits=1)
    print(val)
    writeScenarioFile()             

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
    plotAtmosphericComposition(conc_file=MCONC,id=str(stepnum),ref_file=MCONV,out_dir=OUTPUT)

def resetMEAC():
    os.chdir(PATH)
    subprocess.run(["cp","-rf","scp2000/",MEACPATH])

##############################

# make output folder for run if it doesn't exist yet
if not (NAME in os.listdir('outputs')):
    os.mkdir(OUTPUT)
    os.chdir(OUTPUT)
    os.mkdir('clima-out')
    os.mkdir('meac-out')
    os.mkdir('mr')
    os.chdir(PATH)
else:
    clearDirectory = input("Do you want to clear this directory? [y/n]\n")
    if clearDirectory=='y':
        subprocess.call(['rm','-r',OUTPUT])
        os.mkdir(OUTPUT)
        os.chdir(OUTPUT)
        os.mkdir('clima-out')
        os.mkdir('meac-out')
        os.mkdir('mr')
        os.chdir(PATH)

writeParameters()
writeIncludeFile()
for i in range(NLOOPS):
    first = ((i==0) and not RESUMERUN)
    updateCLIMA(first)
    # input()
    runCLIMA(i)
    # # input()
    updateMEAC()
    # input()
    runMEAC(i)
    # input()
    os.system('clear')
