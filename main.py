import matplotlib.pyplot as plt
import numpy as np
import os
from shutil import rmtree
import subprocess
import sys
import time

from compare_model_outputs_v2 import *
# from files import *
from make_seed_concSTDfile import *
from parameters import *
from plots import *

#######################
###   Bookkeeping   ###
#######################

CLIMAstepIntervals      = []     # Step counts where MEAC runs occur
surfTemps     = []               # Surface temperatures after each CLIMA step (not each loop)
zmax = ZMAX
atm2Pa  = 101_325

########################
###       main       ###
########################

def main(name='test'):

    NAME = name
    # Set relevant paths
    OUTPUT          =   "outputs/"+NAME
    MSCENARIOPATH   =   f'scenario_library/{NAME}'             # Path to MEAC scenario folder
    MZTP            =   f'{MSCENARIOPATH}/TP.dat'                       # MEAC ztp profile
    # MBASE           =   f"{MEACPATH}/{MSCENARIOPATH}/ConcentrationSTD_base.dat"    # Starting conc file
    MBASE           =   f"{MEACPATH}/Earth/ConcentrationSTD_base_Earth.dat"
    MSCENARIO       =   f'{MSCENARIOPATH}/planet_scenario_N2r.h'    # MEAC scenario file with planet parameters
    MSPECIES        =   f'{MSCENARIOPATH}/species_scenario_N2r.dat'       # MEAC atmosphere species file
    MCONC           =   f'{MEACPATH}/{MSCENARIOPATH}/ConcentrationSTD.dat'         # MEAC concentrations file

    ########################
    ###  File functions  ###
    ########################

    #   I have moved certain constants that used to be in 'parameters.py' to 'main.py' so that
    #   they can be initialized using a user-provided scenario name. 
    #   
    #   Because some file functions rely on these constants, I've moved them to 'main.py' as well.
    #
    #   Additionally, since these constants are defined after user input, which is retrieved
    #   in the main function call, these file functions, and the main component functions, are
    #   now moved to inside the definition of the main function.
    #
    #   Brief testing incidicates the coupled model still runs, the correct MEAC scenario folder
    #   is being filled,and outputs are saved to the right place.

    def writeOrCreate(file:str):    # dumb
        # This needs to go at some point
        try:
            f = open(file,'x')
        except FileExistsError:
            f = open(file,'w')
        return f

    def writeParameters():
        """
        Write all* initial parameters to a file.
        """
        f = writeOrCreate(f'{OUTPUT}/parameters.txt')
        f.write("###    MODEL PARAMETERS\n")
        f.write(f"ND                    =   {ND}\n")
        f.write(f"NBIN                  =   {NBIN}\n")
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
        f.write(f"MEAC Scenario Path    =   {MSCENARIOPATH}\n")
        f.write(f"MEAC Converged File   =   {MBASE}\n")
        f.write(f"MEAC Scenario File    =   {MSCENARIO}\n")
        f.write(f"MEAC Species File     =   {MSPECIES}\n")
        f.write(f"MEAC ConcentrationFile=   {MCONC}\n")
        f.write(f"MEAC T-P Profile      =   {MZTP}\n")

        f.close()

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

    def writeIncludeFile():
        f = writeOrCreate(CINCLUDE)
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

    def writeCLIMAinput(first:bool):
        f = open('templates/clima_input.txt','r')
        template = f.read()
        f.close()

        # Do the writing
        template = template.replace('{1}',str(NMAXSTEPS))           # number of CLIMA steps
        template = template.replace('{2}',str(RELHUM))              # relative humidity
        template = template.replace('{21}',str(TSURF))              # Surface temperature
        template = template.replace('{22}',str(T0))                 # Temperature at model top
        template = template.replace('{3}',str(P0))                  # TOA pressure [bar]
        template = template.replace('{4}',str(PSURF))               # Surface pressure [bar]
        template = template.replace('{5}',str(G))                   # Surface gravity [cgs]
        template = template.replace('{6}',str(IO3))                 # Ozone flag
        template = template.replace('{7}',str(ICONSERVE))           # Energy conservation flag
        template = template.replace('{8}',str(SURFALB))             # Surface albedo
        template = template.replace('{9}',str(INSTELL))             # Instellation [S_Earth]
        template = template.replace('{10}',"1.0")                   # Max CO2 mixing ratio
        template = template.replace('{11}',str(IME))                # Methane/ethane flag
        template = template.replace('{12}',str(DOEDDY).lower())     # Eddy flag
        template = template.replace('{13}',str(FCLOUD))             # Fractional cloudiness
        template = template.replace('{14}',str(KZ_MIN))             # Eddy diffusivity
        template = template.replace('{15}',str(CRAINF))             # Rainout parameter
        template = template.replace('{16}',str(CSIG))               # ???
        template = template.replace('{17}',str(SUPERSAT))           # ???
        template = template.replace('{18}',str(COLDTRAPMINMIX))     # Minimum mixing ratio above coldtrap
        template = template.replace('{19}',str(FCMINF))             # See above
        template = template.replace('{20}',str(int(first)))         # Using previous run outputs flag

        f = writeOrCreate(CINPUT)
        f.write(template)
        f.close()

    def writeMEACmain():
        f = open('templates/meac_main.txt','r')
        text = f.read()
        f.close()

        text = text.replace('{1}',MSCENARIO)
        text = text.replace('{2}',str(NMAXT))

        o = open(MEACPATH+'/main.c','w')
        o.write(text)
        o.close()

    def writeMEACspecies(tsurf:float):
        # This function updates the lower boundary flux of water only -- other 
        #   mixing ratios need to be changed manually!

        # The boundary condition is calculated from surface pressure and surface
        #   temperature using the below function. Surface values are sourced from
        #   CLIMA out file.
        def waterPressure(temp:float):
            if temp < 273.16:               # Murphy & Koop (2005)
                res = np.exp(9.550426 - 5723.265/temp + 3.53068*np.log(temp) - 0.00728332*temp)
                res /= 1e5  # Convert to bar
            else:                           # Seinfield & Pandis (2006)
                a = 1 - 373.15/temp
                res = np.exp(13.3185*a - 1.97*a*a - 0.6445*a*a*a - 0.1229*a*a*a*a)
            
            res = np.round(res,6)
            print(f"When surface temperature = {temp} K, surface water vapor partial pressure is {res}")
            return res
        
        val = waterPressure(tsurf)

        f = open('templates/meac_species.txt','r')
        text = f.read()
        f.close()

        text = text.replace("{1}",str(val))

        f = open(MEACPATH+'/'+MSPECIES,'w')
        f.write(text)
        f.close()

    def writeMEACout(conc_file:str,out_dir:str='',id:str=''):
        # holy shit my code is ass
        file = open(conc_file,'r')
        data = file.read().replace('#','').split()[121:]
        data = np.asarray(data,dtype=np.float32).reshape(len(data)//116,116)
        file.close()

        # TOA to surface
        alts = data[::-1,0]
        pressures = data[::-1,4]

        nd_all  =   np.sum(data[::-1],axis=1)

        # number densities                      mixing ratios
        ndC2H6  =   data[::-1,35];              mrC2H6  =   np.divide(ndC2H6,nd_all)     
        ndCH4   =   data[::-1,25];              mrCH4   =   np.divide(ndCH4,nd_all)      
        ndCO2   =   data[::-1,56];              mrCO2   =   np.divide(ndCO2,nd_all)     
        ndH2    =   data[::-1,57];              mrH2    =   np.divide(ndH2,nd_all)       
        ndH2O   =   data[::-1,11];              mrH2O   =   np.divide(ndH2O,nd_all)     
        ndN2    =   data[::-1,59];              mrN2    =   np.divide(ndN2,nd_all)       
        ndO2    =   data[::-1,58];              mrO2    =   np.divide(ndO2,nd_all)       
        ndO3    =   data[::-1,6];               mrO3    =   np.divide(ndO3,nd_all)

        f = writeOrCreate(f"outputs/{out_dir}/meac-out/mr_{id}.dat")
        f.write('Layer\tAltitude [km]\tPressure [Pa]\tC2H6\t\tC2H6_mr\t\tCH4\t\t\tCH4_mr\t\tCO2\t\t\tCO2_mr\t\t')
        f.write('H2\t\t\tH2_mr\t\tH2O\t\t\tH2O_mr\t\tN2\t\t\tN2_mr\t\tO2\t\t\tO2_mr\t\tO3\t\t\tO3_mr\n')
        for j in range(len(nd_all)):
            f.write(f"{j}\t\t")
            f.write(np.format_float_scientific(alts[j],precision=2,trim='k',unique=True,exp_digits=2,min_digits=2)+'\t\t')
            f.write(np.format_float_scientific(pressures[j],precision=3,trim='k',unique=True,exp_digits=2,min_digits=3)+'\t\t')
            f.write(np.format_float_scientific(ndC2H6[j],precision=3,trim='k',unique=True,exp_digits=2,min_digits=3)+'\t')
            f.write(np.format_float_scientific(mrC2H6[j],precision=3,trim='k',unique=True,exp_digits=2,min_digits=3)+'\t')
            f.write(np.format_float_scientific(ndCH4[j],precision=3,trim='k',unique=True,exp_digits=2,min_digits=3)+'\t')
            f.write(np.format_float_scientific(mrCH4[j],precision=3,trim='k',unique=True,exp_digits=2,min_digits=3)+'\t')
            f.write(np.format_float_scientific(ndCO2[j],precision=3,trim='k',unique=True,exp_digits=2,min_digits=3)+'\t')
            f.write(np.format_float_scientific(mrCO2[j],precision=3,trim='k',unique=True,exp_digits=2,min_digits=3)+'\t')
            f.write(np.format_float_scientific(ndH2[j],precision=3,trim='k',unique=True,exp_digits=2,min_digits=3)+'\t')
            f.write(np.format_float_scientific(mrH2[j],precision=3,trim='k',unique=True,exp_digits=2,min_digits=3)+'\t')
            f.write(np.format_float_scientific(ndH2O[j],precision=3,trim='k',unique=True,exp_digits=2,min_digits=3)+'\t')
            f.write(np.format_float_scientific(mrH2O[j],precision=3,trim='k',unique=True,exp_digits=2,min_digits=3)+'\t')
            f.write(np.format_float_scientific(ndN2[j],precision=3,trim='k',unique=True,exp_digits=2,min_digits=3)+'\t')
            f.write(np.format_float_scientific(mrN2[j],precision=3,trim='k',unique=True,exp_digits=2,min_digits=3)+'\t')
            f.write(np.format_float_scientific(ndO2[j],precision=3,trim='k',unique=True,exp_digits=2,min_digits=3)+'\t')
            f.write(np.format_float_scientific(mrO2[j],precision=3,trim='k',unique=True,exp_digits=2,min_digits=3)+'\t')
            f.write(np.format_float_scientific(ndO3[j],precision=3,trim='k',unique=True,exp_digits=2,min_digits=3)+'\t')
            f.write(np.format_float_scientific(mrO3[j],precision=3,trim='k',unique=True,exp_digits=2,min_digits=3)+'\n')
        f.close()

    def writeCLIMAout(clima_last:str,out_dir:str='',id:str=''):
        # Open clima_last which contains the final temperature-pressure profile
        f = open(clima_last,'r')
        data = ''.join(f.readlines()[1:])
        data = np.fromstring(data,dtype=np.float32,sep=' ').reshape((ND,9))
        f.close()

        alts = data[:,0]
        temps = data[:,2]
        pres = data[:,1]*atm2Pa

        # Open species mixing ratio profiles
        c2h6    = np.genfromtxt(CINOUT+'/Profiles_out/C2H6.dat')
        ch4     = np.genfromtxt(CINOUT+'/Profiles_out/CH4.dat')
        co2     = np.genfromtxt(CINOUT+'/Profiles_out/CO2.dat')
        h2      = np.genfromtxt(CINOUT+'/Profiles_out/H2.dat')
        h2o     = np.genfromtxt(CINOUT+'/Profiles_out/H2O.dat')
        n2      = np.genfromtxt(CINOUT+'/Profiles_out/N2.dat')
        o2      = np.genfromtxt(CINOUT+'/Profiles_out/O2.dat')
        o3      = np.genfromtxt(CINOUT+'/Profiles_out/O3.dat')

        f = writeOrCreate(f'{out_dir}/clima-out/ztp_{id}.dat')
        f.write("Altitude [km]\tTemperature [K]\tPressure [Pa]\tC2H6\t\tCH4\t\t\tCO2\t\t\tH2\t\t\tH2O\t\t\tN2\t\t\tO2\t\t\tO3\n")
        for j in range(len(alts)):
            p = np.format_float_scientific(pres[j],precision=4,trim='k',unique=True,min_digits=4)
            f.write(f"{alts[j]:.4f}\t\t\t{temps[j]:.4f}\t\t{p}\t\t")
            for arr in [c2h6,ch4,co2,h2,h2o,n2,o2,o3]:
                f.write(np.format_float_scientific(arr[j],precision=4,min_digits=4)+'\t')
            f.write('\n')

        f.close()

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
        pressures   = data[::-1,4]    
        nd_all      =   np.sum(data[::-1],axis=1)               # Total number densities
        nd_nc       =   nd_all - data[::-1,56]- data[::-1,11]   # Non-condensible number densities
        
        # number densities                      mixing ratios 
        ndC2H6  =   data[::-1,35];              mrC2H6  =   np.divide(ndC2H6,nd_nc)     # relative to all non-condensibles
        ndCH4   =   data[::-1,25];              mrCH4   =   np.divide(ndCH4,nd_nc)      # relative to all non-condensibles
        ndCO2   =   data[::-1,56];              mrCO2   =   np.divide(ndCO2,nd_all)     # ABSOLUTE
        ndH2    =   data[::-1,57];              mrH2    =   np.divide(ndH2,nd_nc)       # relative to all non-condensibles
        ndH2O   =   data[::-1,11];              mrH2O   =   np.divide(ndH2O,nd_all)     # not needed, H2O read in thru TempIn.dat
        ndN2    =   data[::-1,59];              mrN2    =   np.divide(ndN2,nd_nc)       # relative to all non-condensibles
        ndO2    =   data[::-1,58];              mrO2    =   np.divide(ndO2,nd_nc)       # relative to all non-condensibles      
        ndO3    =   data[::-1,6];               mrO3    =   np.divide(ndO3,nd_nc)       # relative to all non-condensibles

        # get surface pressure by linearly extrpolating MEAC pressures to 0km
        p = pressures[::-1]
        PSURF = (p[0] + 0.5*(p[0]-p[1]))/atm2Pa     # atm


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


    start = time.time()

    # Warning about pre-existing folders with the given name
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

    # Save model parameters and inputs
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

if __name__ == '__main__':
    name = sys.argv[1]
    main(name)
