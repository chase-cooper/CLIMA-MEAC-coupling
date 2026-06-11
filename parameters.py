import os

# Coupling parameters
ND              =   101                     # CLIMA number of layers -- PLEASE DONT TOUCH
NBIN            =   50                      # MEAC number of layers
NLOOPS          =   15                       # Number of CLIMA-MEAC loops
NMINSTEPS       =   500                      # Minimum number of CLIMA steps per loop
NMAXSTEPS       =   500                      # Max Number of CLIMA steps per loops
NMAXT           =   1e100                    # Max MEAC cumulative integration time
TCONV           =   1E-1                    # CLIMA convergence criterion, delta temperature (not too important)

# Planet + atmosphere parameters
MASS            =   5.9722e24               # Planet mass in kg
RAD             =   6371.0                  # Planet radius in km
G               =   980                     # grav*M/R**2; Surface gravity, cgs
A               =   1.0                     # Semimajor axis in AU
SURFALB         =   0.3                     # Planet surface albedo
INSTELL         =   1.0                     # Planet instelation, relative to Earth

# CLIMA parameters
P0              =   1e-5                    # Top-of-atmosphere pressure [atm] considered by CLIMA. If this value is too low, temperatures might blow up :(
PSURF           =   1e+0                    # Surface pressure [atm]
T0              =   180                     # Top-of-atmosphere temperature initial guess [K]   * not used
TSURF           =   330                     # Surface temperature initial guess [K]             * not used
TROPOPAUSE      =   20                      # CLIMA tropopause layer, default 22 (of 101)
AR              =   0                       # Argon mixing ratio
RELHUM          =   1.0                     # Surface relative humidity
FIXH2O          =   1                       # Fixed H2O flag. I really recommend leaving this on. Also, doesn't actually fix H2O?
IO3             =   1                       # Ozone flag -- IO3=0 means ozone isn't read in
IME             =   1                       # Methane/ethane flag

# Eddysed parameters (also for CLIMA)
DOEDDY          =   'true'                  # eddysed flag; if false, nothing below this matters
FCLOUD          =   0.0                     # fractional cloud coverage   
KZ_MIN          =   100000.0                # Minimum eddydiffusivity
CRAINF          =   0.01                    #
CSIG            =   2.0                     #
SUPERSAT        =   0.0                     #
COLDTRAPMINMIX  =   4e-10                   # Minimum coldtrap water mass mixing ratio
FCMINF          =   0.01                    # Minimum fraction for FC for upper atmosphere

# Misc
RESUMERUN       =   0                       # Continue from most recent run. Set to 0 for coupling
ICONSERVE       =   1                       # CLIMA energy conservation flag, PLEASE KEEP ON 1
RAINOUT         =   1                       # controls rainout in 'main.c'. not used yet
ADIABATIC       =   0                       # tells CLIMA that atmosphere is dry adiabat  
ZMAX            =   50                      # Uppermost altitude; overwritten in main.py, not used

atm2Pa          =   101_325

#######################
###  Paths & Files  ###
#######################

# Paths
NAME            =   'fco2_1e-1'
OUTPUT          =   "outputs/"+NAME
PATH            =   os.getcwd()                                     # Path to this python file
CLIMAPATH       =   f'{PATH}/cloudy_clima'                          # Path to the cloudy-CLIMA folder
MEACPATH        =   f'{PATH}/hu-code-sr'                            # Path to the MEAC folder

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
MSCENARIOPATH   =   f'scenario_library/{NAME}'             # Path to MEAC scenario folder
MZTP            =   f'{MSCENARIOPATH}/TP.dat'                       # MEAC ztp profile
MBASE           =   f"{MEACPATH}/{MSCENARIOPATH}/ConcentrationSTD_base.dat"    # Starting conc file
MSCENARIO       =   f'{MSCENARIOPATH}/planet_scenario_N2r.h'    # MEAC scenario file with planet parameters
MSPECIES        =   f'{MSCENARIOPATH}/species_scenario_N2r.dat'       # MEAC atmosphere species file
MCONC           =   f'{MEACPATH}/{MSCENARIOPATH}/ConcentrationSTD.dat'         # MEAC concentrations file