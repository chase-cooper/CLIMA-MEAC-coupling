import numpy as np
import re

from parameters import *

def writeOrCreate(file:str):
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
    f.write(f"MEAC Scenario Name    =   {MSCENARIONAME}\n")
    f.write(f"MEAC Scenario Path    =   {MSCENARIOPATH}\n")
    f.write(f"MEAC Converged File   =   {MCONV}\n")
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

    print("ENSURE # OF BINS IS PROPERLY SET")

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

def writeMEACspecies(tsurf:float,psurf:float=101325):
    # This function updates the lower boundary flux of water only -- other 
    #   mixing ratios need to be changed manually!

    # The boundary condition is calculated from surface pressure and surface
    #   temperature using the below function. Surface values are sourced from
    #   CLIMA out file.
    def waterPressure(temp:float,psurf:float=101325):
        if temp < 273.16:               # Murphy & Koop (2005)
            res = np.exp(9.550426 - 5723.265/temp + 3.53068*np.log(temp) - 0.00728332*temp)
            res /= psurf
        else:                           # Seinfield & Pandis (2006)
            a = 1 - 373.15/temp
            res = (101325/psurf)*np.exp(13.3185*a - 1.97*a*a - 0.6445*a*a*a - 0.1229*a*a*a*a)
        
        res = np.round(res,4)
        # print(f"When surface temperature = {temp} K, surface water vapor mixing ratio is {res}")
        return res
    
    val = waterPressure(tsurf,psurf)

    f = open('templates/meac_species.txt','r')
    text = f.read()
    f.close()

    text = text.replace("{1}",str(val))

    f = open(MEACPATH+'/'+MSPECIES,'w')
    f.write(text)
    f.close()

    print("ENSURE MIXING RATIOS ARE PROPERLY SET IN SPECIES TEMPLATE FILE")

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
    f = open(clima_last,'r')
    data = ''.join(f.readlines()[1:])
    data = np.fromstring(data,dtype=np.float32,sep=' ').reshape((ND,9))
    f.close()

    alts = data[:,0]
    temps = data[:,2]
    pres = data[:,1]*atm2Pa

    f = writeOrCreate(f'{out_dir}/clima-out/ztp_{id}.dat')
    f.write("Altitude [km]\tTemperature [K]\tPressure [Pa]\n")
    for j in range(len(alts)):
        p = np.format_float_positional(pres[j],precision=6,trim='k',unique=True,min_digits=6)
        f.write(f"{alts[j]:.4f}\t\t\t{temps[j]:.4f}\t\t\t{p}\n")
    f.close()




