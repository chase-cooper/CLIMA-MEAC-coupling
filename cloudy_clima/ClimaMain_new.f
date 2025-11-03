C     PROGRAM SURFT(INPUT,OUTPUT,TAPE1,TAPE2,TAPE3)
c
C  This program is a modified version of the climate model SURFTEM made 
c  by James Kasting. The program has been modified by Michael Mischna (mm),
c  Alex Pavlov (AP), Kara Krelove (KK), Hilary Justh (HJ), Ravi Kopparapu (RK), 
c  Ramses Ramirez(RR), and Antigona Segura (AS). Some changes are identified 
c  with the initials of the author.

c  The code is mostly written in f77 but is compiled in f90 and it 
c  contains some f90 features.

c  This code is a 1-D, cloud-free (it has clouds now), radiative-convective climate model.R
c  The calculation of temperature profiles begins with an initial 
c  temperature-pressure profile and a solar constant. 

c  The net absorbed solar radiation is calculated by a delta two-stream
c  approximation (Toon, et al. JGR Vol. 94, 16287-16301, 1989). It uses
c  4-term correlated k coefficients to parameterize absorption by O3,
c  CO2, H2O, O2 and CH4 in 38 spectral intervals.

c  The IR is calculated by the RRTM routine developed by Mlawer et. al
c  (JGR, Vol.102 (D14), 16663-16682, 1997). It uses 16 term sums in 
c  each of its spectral bands in which the k-coefficients are concentrated 
c  in areas of most rapidly changing absorption. The version 3.0 of RRTM 
c  was implemeted on August/2003 (www.rtweb.aer.com).

c  When the mixing ratio of CO2 is greater than CO2MAX, the maximum 
c  level of CO2 that RRTM can manage, the former IR subroutine is used.
c  (Pavlov et al. J. Geophys. Res. 105: 11,981-11,990, 2000).

c  Units in cgs unless otherwise is stated.
   
c  Temperature in each layer is calculated from:
c              dT/dt = - (g/c_p) dF/dp
c  in this case the derivates are partial. T= temperature, t= time, 
c  g= gravitational constant, F=Flux, c_p= Heat capacity, p=pressure.

c  Two types of reach convergence have been set up. One uses a non-strict
c  time stepping mode which is faster and better for high O2-low CO2 runs,
c  like present Earth. The other one is slower but needed on high CO2 
c  atmospheres. 

c  This model can work alone or coupled to a photochemical model. 
c  Modifications for the coupled mode were made by Kara Krelove.  

c Input data files required by the program are:
C     Unit   File
C      3     H2O_tables.pdat
C      4     solar_data_38.pdat (Read by 2-stream code)
C      8     nearir_expsums.pdat
c      9     CO2_tables.pdat
c     20     ir_expsums.pdat
c        21         BIG_DATAFILE.dat

C
C   THE VERTICAL GRID IS STAGGERED BETWEEN TEMPERATURE AND FLUX
C   GRID POINTS.  THE FLUX GRID IS DEFINED FROM THE VERY TOP OF THE
C   ATMOSPHERE (J=1) TO THE GROUND (J=ND).  THE TEMPERATURE GRID POINTS
C   ARE HALFWAY BETWEEN THE FLUX POINTS, EXCEPT FOR T(ND) WHICH IS
C   LOCATED AT THE GROUND.
C
C   PARAMETERS:
C   ND = # OF ALTITUDE POINTS  (J)
C   NF = # OF FREQUENCIES  (N)
C   NA = # OF ANGLES  (M)
C   NS = # OF CHEMICAL SPECIES  (I)
C   NT = # OF TEMPERATURES IN THE STEAM TABLE
C   NSOL = # OF SOLAR FREQUENCIES
C
C   T = TEMPERATURE (K)
C   P = PRESSURE (bar)
C   Z = LOG PRESSURE + A CONSTANT (ZCON)
C   PF = PRESSURE AT FLUX GRID POINTS
C   ZF = LOG P AT FLUX POINTS
C   ALT = ALTITUDE (KM)
C   GAM = DTDZ
C   BVK = PLANCK FUNCTION
C   LAM = WAVELENGTHS (MICRONS)
C   AV = FREQUENCIES (1/S)
C   TAU = SLANT OPTICAL DEPTH TO OTHER PRESSURE LEVELS
C   F = INTEGRATED NET FLUX
C   FS = INTEGRATED SOLAR FLUX
C   FI = SPECIES MIXING RATIOS   1 = water, 2 = co2, 3 = ch4, 4 = o3, 5 = ethane
C   FH2O - H2O MIXING RATIO
C   T,TN - TEMPERATURES
C   FLAGCONVEC - Tags for the type of convection
c                1. = Water moist adiabat
c                2. = Water dry diabat
c                3. = CO2 adiabat 
c                0. = Non convective layer
 
C-KK        NLAYERS is a translation parameter between this climate model
C-KK    and Mlawer's RRTM code. 
C_KK    SurfTem indexes from 1 at the top to ND at the ground, while 
C_KK    RRTM indexes from 0 at the ground to NLAYERS at the top.
C-KK        NZ is the number of layers being carried in atm_chem. 
      INCLUDE 'CLIMA/INCLUDE/header.inc'
      INCLUDE 'globals.h'
      INCLUDE 'prog_params' !included from EddySed JDW

      PARAMETER(NF=55, NA=1, NLAYERS=ND-1, NZ=200)
      PARAMETER(NS=3, NS1=NS+2, NS4=NS+5) !gna: changed NS1 from NS+1 to NS+2 to add ethane
      PARAMETER(NT=76, MT=36)
      PARAMETER(NSOL=38, NGS=8, IK=8)  ! Added IK=8 parameter and NGS is 7 now, 3/26/2012
      parameter(nrow=11)
      integer eJCOLD,k,zz

      CHARACTER*5 :: STARR   !Changed to make STARR hold up to 5 characters
      CHARACTER*11 :: AA
      CHARACTER*200 :: BB
      CHARACTER :: DIRINOUT*8,DIRDATA*10

      DIMENSION TRAD(ND),DZ(ND),Z(ND),ZF(ND)
      DIMENSION temp_alt(NZ), water(NZ), O3(NZ), PRESS(NZ), !EWS - temp_t(NZ) removed because it wasn't used
     & CH4(NZ), CO2(NZ), ethane(NZ)
      DIMENSION T(ND),TOLD(ND),FTOTAL(ND),FTIR(ND),
     & FTSO(ND),PF1(ND),DELT(ND),DELTRAD(ND),TN(ND),
     & DIVF(ND),TCOOL(ND),THEAT(ND),FLAGCONVEC(ND),
     & told1(ND),told2(ND),told3(ND),pold(ND),pold1(ND),
     & pold2(ND),pold3(ND)
      DIMENSION FSATURATION(ND),FSATUR(ND),FSAVE(ND) !Removed T(ND)!EWS dt(ND) removed, not used 8/18/2015
      DIMENSION HEATNET(ND),BETA(ND),FCO2V(ND),FH2O(ND)
      DIMENSION AVOLD(NF) !EWS - ALAM not used
      DIMENSION PSATCO2(ND) !EWS - PML(ND) removed because it wasn't used
      DIMENSION FNC(ND)        ! Added FNC array c-rr 6/7/2012
      DIMENSION FNC_cloudy(ND)

      double precision alt_convec
      double precision eFI(NS1,ND)
      double precision sat(NS1,ND)

c     vectors for gaussian zenith angles
      dimension fdnsoltot(nd), fupsoltot(nd)
      dimension fdnsoltot_cloudy(nd),fupsoltot_cloudy(nd) !JDW
      dimension fdnsoltot_clear(nd),fupsoltot_clear(nd) !JDW
      dimension fupir_combined(nd),fdnir_combined(nd)
      dimension xi(nrow,20), wi(nrow,20), ngauss(nrow)
      dimension iconvec(nd)
      dimension FI_cloudy1(1,ND)

      REAL*8 newalt ! removed extraneous kappa, kappa_ir, and FLAGCONVE(ND) declarations
      double precision mw_atmos,RSURF
      double precision bwni, wnoi, dwni, wlni

      logical doCloud
      logical doEddy
      logical do_highres
      logical laststep
      double precision fcloud
      double precision teff
      double precision kz_min
      double precision Crainf
      double precision Csig
      double precision supsat
      double precision relhum_holder(ND)
      double precision fc_minf

      double precision cloud_hum1,clear_hum1
      double precision Mean_mmr
      double precision new_relhum
      double precision cFI(NS1,ND)
      double precision relhum_vec(ND)
      double precision upatm_mix
      double precision AV_highres(1000)
      double precision ALAM_highres(1000)
      double precision LAM_highres(1000)

      logical couplesmart
      character(len=200) :: fileout,fileout_highres
      double precision srmix(10,nd)
      double precision Fnet_rhr(nd),pressure_rhr(nd),altitude_rhr(nd),
     & c_p_rhr(nd),rhr_solar(nd),rhr_thermal(nd),Fnet_rhr_s(nd),Fnet_rhr_t(nd)
      COMMON/fh2o_ed/FH2O_e(ND)
      common/molec_weight/mw_atmos
      COMMON/SPECTI/ BWNI(NSPC1IT),WNOI(NSPECIT),DWNI(NSPECIT),
     & WLNI(NSPECIT)

      common/fc_min/fc_minf
      common/relhum_eddy/relhum_eddy_july_new,foundloc !JDW

      COMMON/DIR/DIRINOUT,DIRDATA
      COMMON/WAVE/AV(NF),LAM(NF),W(NF)
      COMMON/ABLOK/LTYPE(NF,3),XH2O(NF),YH2O(NF),XCO2(NF),YCO2(NF),
     & AXH(NF),AYH(NF),BXH(NF),BYH(NF),AXC(NF),AYC(NF),BXC(NF),
     & BYC(NF),PDOP(NF),CPR(NF),TPR(NF),PATH(NS4),PATHP(NS4),
     & PATHT(NS4),P1,P2,T1,T2,TAU2,TAUP2,ALPHA(4),BETH2O(4,5,NF),
     & BETCO2(4,5,NF),CA(19),CB(19),CC(19),CD(19),CE(19),CH(19),CK(19)
    
      COMMON/CBLOK/FO2,FN2,FCO2,FAR,FCH4,FC2H6,FNO2, FI(NS1,ND),FH22,
     & FI_cloudy(NS1,ND)
      COMMON/CBLOK/AO2(ND),AN2(ND),ACO2(ND),ACH4(ND),AC2H6(ND),AH2(ND),AH2O(ND),AO3(ND)
      DIMENSION FIXEDWATER(ND)      ! initial water mixing ratio

      COMMON/ALTBLOK/DALT(ND-1),RADIUS(ND-1),PARTICLES(ND),RAER(ND),
     & ALT(ND)
      COMMON/EBLOK/PG,TG,PG0,IMW,RSURF,OMEGA,POCEAN,IMOIST,
     & BETA1,BETA2,FVDRY,PDRY
      COMMON/FBLOK/TTAB(NT),PVAP(NT),DPVAP(NT),SVC(NT),DSV(NT),DSC(NT),
     & RHOV(NT),DRHOV(NT),BETAM(70,75),TCP(75),PCP(70),DPDTL(70,75),
     & DRDTL(70,75)
      COMMON/GBLOK/TCTAB(MT),PCVAP(MT),BETASC(MT),DPCVAP(MT),
     & DRCVAP(MT),SVSC(MT),DSCC(MT),TKTAB(MT),TCC(25),PCC(36),
     & BETAMC(25,36),CPC(25,36),DVDTC(25,36),DVDPC(25,36),DSVC(MT)
      COMMON/SBLOK/P0P,T0P,R,SUBL
      COMMON/PRESSURE/P(ND),PLOG(ND)
      COMMON/PRESS/BETIR1(4,5,NSOL),BETIR2(4,5,NSOL),
     & kappa_solh2o(NSOL,8,8,IK), kappa_solco2(NSOL,8,8,IK) ! Added new kappa matricies for each of CO2 and H2O coefficients. 8/26/2012

      COMMON/AOZONE/BETAO3(nsol), BETAO2(2),WGHTO2(NSOL,2)
      COMMON/RSOL/ALPHAZ(4,2),BETAZ(4,2),NPROB(2),
     & NG(2),SIGG(4,2,NSOL)
      common/colblok/eJCOLD,upatm_mix
      COMMON/STEPS/NST
      COMMON/SOLARBLK/AMU0,SRFALB,OMG0A(NSOL,ND-1),
     & ASYA(NSOL,ND-1),TAUAER(NSOL),SIGERT(NSOL),FMA(NSOL),PF(ND),
     & ALAMBDA(NSOL),CGAS(ND,NGS),FUPSOL(ND),FDNSOL(ND),
     & NGAS(2,NSOL),WGHT(4,2,NSOL),NPR(2,NSOL),SOLINT(NSOL),
     & TAULAM(ND-1),ASY(ND-1),OMG0(ND-1),FMT(ND-1),QEXT(NSOL,ND-1),
     & fdnsol_cloudy(nd),fupsol_cloudy(nd),
     & fdnsol_clear(nd),fupsol_clear(nd),ASY_cloudy(ND-1),
     & OMG0_cloudy(ND-1),TAULAM_cloudy(ND-1),
     & ASY_clear(ND-1),OMG0_clear(ND-1),TAULAM_clear(ND-1),
     & ASYEDDY(ND-1),OMG0EDDY(ND-1),TAUEDDY(ND-1),CGAS_cloudy(ND,NGS),
     & relhum_eddy_july(nd) !JDW  !JDW,CGAS_cloudy(ND,NGS) !JDW

      common/eddyblok/eddyopdIR(NF,MAXNZ), eddyw0IR(NF,MAXNZ),
     & eddyopdSOL(NSOL,MAXNZ),eddyw0SOL(NSOL,MAXNZ),
     & eddyg0SOL(NSOL,MAXNZ),eddyg0IR(NF,MAXNZ),
     & eddyqt(MAXNZ,MAXNGAS),eddyqc(MAXNZ,MAXNGAS),JCOLD

      common/cinputs/doEddy,doCloud,fcloud,kz_min,Crainf,Csig,
     & supsat,nsub_max,cld_hum1,clr_hum1,new_relhum(ND)
      common/tclima/tclim(ND),alt_convec(ND)
      COMMON/IRDATA/WEIGHTCH4(6),xkappa(3,12,55,8),
     & CIA(7,NF), CPRW(ND,NF)
      COMMON/VARIR/kappa_irh2o(NF,8,8,IK), kappa_irco2(NF,8,8,IK)! Added kappa matrix in IR for kpsectrum Co2 and H2O coefficients 8/26/2012
      COMMON/weightsIR/ weightco2_h2oIR(IK)
      COMMON/IRBLK/FUPIR(ND),FDNIR(ND),SRFALBIR,OMG0AIR(NF,ND-1),
     & ASYAIR(NF,ND-1),IO3,QEXTIR(NF,ND-1),FUPIR_cloudy(ND),
     & FDNIR_cloudy(ND),FUPIR_clear(ND),FDNIR_clear(ND)

      COMMON/HYDROCARB/Qextirst(73,55),w0irst(73,55),
     &  girst(73,55),Qextsolst(73,38),w0solst(73,38),gsolst(73,38),
     &  radstand(73)
      COMMON/CH4BLOCK/ALPHACH4T188(4,17),BETACH4T188(4,17),
     & ALPHACH4T295(4,17),BETACH4T295(4,17),ALPHACH4Kark(4,21),
     & BETACH4Kark(4,21),GAMMAEXP188(17),GAMMAEXP295(17),
     & ALPHACH4NEW(6),BETACH4NEW(17,3,5,6),ALCH4(6,38)
      COMMON/CO2BLOK/betac1,betac2,PC0,TC0,VAPCL0,SUBCL0,DLVCDT,
     & DLSCDT,CCL,CCS
      COMMON/NO2BLOK/SIGNO2(NSOL)

      COMMON/ MLAWERI/  layers, numspec, newalt(ND), TempT(0:NLAYERS),
     & Pres(0:NLAYERS), gasses(7, 0:NLAYERS), COLDEP(ND)
      COMMON/CONSS/C,BK,G,GNEW(ND),PI,SM,DM,DM2   ! Adding DM2 to common block entry 5/3/2011. DM and DM2 are AMN and AMN2 respectively in CONVEC

      COMMON/CPHEAT/CPO2(ND),CPCO2(ND), CPN2(ND), CPH2O(ND),
     & CPN(ND), CPNT(ND), CPH2(ND)  ! Added CPH2 5/31/2012 c-rr
      common/CPHEAT_cloudy/CPN_cloudy(ND) ! Added 04/14/2021 JDW

      COMMON/smart_optics/sngas,sgas_name,snwave,swave,snrad,sradius,sdr,   !added for the clima to smart outputs. JDW
     & sqscat,sqext,scos_qscat
      common/smart/QextsIR,QExtsSOL,QscasIR,QscasSOL,g0sIR,g0sSOL
     & dtauextsIR,dtauextsSOL,Qexts,Qscas,g0s,dtauexts !JDW+TDR

      common/inverse/INVERSE

      DATA BETA/ND*1./
      DATA BETH2O/1100*0./
      DATA BETCO2/1100*0./
      DATA SIGNO2/1.1E-20, 5.0E-20, 9.5E-20, 2.23E-19, 3.36E-19,
     & 5.1E-19, 5.36E-19, 2.58E-19, 1.07E-19, 8.0E-20, 4.75E-20,
     & 2.65E-20, 1.25E-20, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
     & 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0./

C   FREQUENCIES AT ENDS OF SPECTRAL INTERVALS (1/CM)
      DATA AV/40., 100., 160., 220., 280., 330., 380., 440., 495.,
     & 545., 617., 667., 720., 800., 875., 940., 1000., 1065.,
     & 1108., 1200., 1275., 1350., 1450., 1550., 1650., 1750., 1850.,
     & 1950., 2050., 2200., 2397., 2494., 2796., 3087., 3425., 3760.,
     & 4030., 4540., 4950., 5370., 5925., 6390., 6990., 7650., 8315.,
     & 8850., 9350., 9650., 10400., 11220., 11870., 12790., 13300.,
     & 14470., 15000./

      DATA C,HP,BK,SIGMA,SM,PI/3.E10, 6.63E-27, 1.38E-16, 5.67E-5,
     & 1.67E-24,3.14159274d0/

c Names of the subdirectories for the data, inputs and outputs
      DIRINOUT = 'CLIMA/IO'
      DIRDATA =  'CLIMA/DATA'
!================    FILE SECTION    ==================

      open(unit=901,file='cloud_optics/init_optics_bwni.txt')
      open(unit=902,file='cloud_optics/init_optics_wnoi.txt')
      open(unit=903,file='cloud_optics/init_optics_dwni.txt')
      open(unit=904,file='cloud_optics/init_optics_wave.txt')

C  INPUT FILES
      OPEN (unit=1,file= DIRINOUT//'/input_clima.dat')
      OPEN (unit=3,file= DIRDATA//'/H2O_tables.pdat',status='old')
      OPEN (unit=4,file= DIRDATA//'/solar_data_38.pdat',status='old')
      OPEN (unit=8,file= DIRDATA//'/nearIR_expsums.pdat',status='old')
      OPEN (unit=1111,file=DIRINOUT//'/AV_highres_python.txt')
      OPEN (unit=2222,file=DIRINOUT//'/ALAM_highres.txt')

!====================================================================
C New k-coefficients for H2O and CO2 were calculated by Eric Wolf
C using HELIOS-K (https://github.com/exoclime/HELIOS-K), an
C ultrafast GPU-driven correlated-k sorting program
C (Grimm et al. 2015, doi.org/10.1088/0004-637X/808/2/182).
C For H2O we use the HITRAN2016 line-list, assuming 25 cm-1
C line cut-offs using Lorentz profiles and with the plinth removed.
C For CO2 we also use the HITRAN2016 database, but we assume
C 500 cm-1 line cut-offs using the Perrin and Hartman
C (1989, doi.org/10.1016/0022-4073(89)90077-0) sub-Lorentzian
C line profiles.  These conventions represent the current standard
C practices for the treatment of H2O and CO2 lines within coarse
C spectral resolution climate model radiation schemes.
C It is assumed that the H2O self and foreign broadening components,
C and CO2-CO2 CIA, are included elsewhere in the code,
C both of which are independent of the line treatment.
C For further discussions contact eric.wolf@colorado.edu.

      OPEN (unit=15, file=DIRDATA//'/Wolf_HITRAN2016_solar_38_H2O.dat',
     & status='old')  ! H2O solar coefficients
      OPEN (unit=16, file=DIRDATA//'/Wolf_HITRAN2016_solar_38_CO2.dat',
     & status='old')  ! CO2 solar coefficients
      OPEN (unit=17, file=DIRDATA//'/Wolf_HITRAN2016_ir_55_H2O.dat',
     & status='old')     ! H2O ir coefficients
      OPEN (unit=18, file=DIRDATA//'/Wolf_HITRAN2016_ir_55_CO2.dat',
     & status='old')     ! CO2 ir coefficients

!=============    MISC FILE SECTION    ==================

      OPEN(unit=9,file= DIRDATA//'/CO2_tables.pdat',status='old')
      OPEN(unit=21,file= DIRDATA//'/BIG_DATAFILE.DAT',status='old')
      OPEN(unit=66, file = DIRINOUT//'/weight_factors.txt')
      OPEN(unit=10, 
     & file=DIRDATA//'/STELLAR_SPECTRA_update.pdat',status='old')
      OPEN(unit=30, file=DIRDATA//'/FinalCIAcoeffs2.dat', status='old')
      OPEN(unit=90, file=DIRINOUT//'/FTIR.dat')
      OPEN(unit=91, file=DIRINOUT//'/FTSO.dat')
      OPEN(unit=11,file= DIRINOUT//'/TempIn.dat')
      OPEN(unit=12,file= DIRINOUT//'/TempOut.dat')

      OPEN(unit=113,file= 'COUPLE/fromPhoto2Clima.dat')
      OPEN(unit=222,file= 'COUPLE/fromPhoto2Clima.dat')
      OPEN(unit=116,file= 'COUPLE/fromClima2Photo.dat')   

      OPEN(UNIT=98,FILE= DIRINOUT//'/clima_allout.tab')
      OPEN(UNIT=96,FILE= DIRINOUT//'/SolarHeating.tab')
      OPEN(UNIT=97,FILE= DIRINOUT//'/clima_last.tab')
      OPEN(UNIT=80,FILE= DIRINOUT//'/IR_wavelength_grid.tab')

      OPEN(UNIT=2021,FILE= DIRINOUT//'/Fluxes_all.tab')
      OPEN(UNIT=20212,FILE= DIRINOUT//'/radiative_heating_output.tab')

!================= VERTICAL MIXING RATIOS ======================    CC2025
      OPEN(unit=500,file=DIRINOUT//'/Profiles/O2.dat')
      OPEN(unit=501,file=DIRINOUT//'/Profiles/N2.dat')
      OPEN(unit=502,file=DIRINOUT//'/Profiles/H2.dat')
      OPEN(unit=503,file=DIRINOUT//'/Profiles/CO2.dat')
      OPEN(unit=504,file=DIRINOUT//'/Profiles/CH4.dat')
      OPEN(unit=505,file=DIRINOUT//'/Profiles/C2H6.dat')
      OPEN(unit=506,file=DIRINOUT//'/Profiles/H2O.dat')
      OPEN(unit=507,file=DIRINOUT//'/Profiles/O3.dat')

      OPEN(unit=600,file=DIRINOUT//'/Profiles_out/O2.dat')
      OPEN(unit=601,file=DIRINOUT//'/Profiles_out/N2.dat')
      OPEN(unit=602,file=DIRINOUT//'/Profiles_out/H2.dat')
      OPEN(unit=603,file=DIRINOUT//'/Profiles_out/CO2.dat')
      OPEN(unit=604,file=DIRINOUT//'/Profiles_out/CH4.dat')
      OPEN(unit=605,file=DIRINOUT//'/Profiles_out/C2H6.dat')
      OPEN(unit=606,file=DIRINOUT//'/Profiles_out/H2O.dat')
      OPEN(unit=607,file=DIRINOUT//'/Profiles_out/O3.dat')

!=========== Altitude array ================ CC2025
      OPEN(unit=700,file=DIRINOUT//'/alt.dat')

!   IDK what these do
      read(901,*) bwni
      read(902,*) wnoi
      read(903,*) dwni
      read(904,*) wlni

!   IDK what these do either
      Idry = 0
      read(1111,*) AV_highres
      read(2222,*) LAM_highres

!=================== CLIMA INPUTS ========================
!   Reading in from input_clima.dat. See that file for variable descriptions
      READ(1,51)            !   removes leading characters
51    FORMAT(4/)

      READ(1,*) AA,NSTEPS       !step number
      READ(1,*) AA,IMW
      READ(1,*) AA,RSURF               
      READ(1,*) AA,zy
      READ(1,*) AA,DTAU0
      READ(1,*) AA,ZCON
      READ(1,*) AA,P0           !Pressure at the top
      READ(1,*) AA,PG0
      READ(1,*) AA,G            !Gravity (Mars=373., Earth=980.) 
      READ(1,*) AA,FAC
      READ(1,*) AA,IO3                !Ozone?
      READ(1,*) AA,IUP                    
      READ(1,*) AA,TG0                !Surface temperature for IUP=1   
      READ(1,*) AA,TSTRAT       !Stratospheric temperature for  IUP=1
      READ(1,*) AA,STARR        !What star?
      READ(1,*) AA,ICONSERV     !Type of energy conservation
      READ(1,*) AA,ICOUPLE      !Coupled(1) or not(0)
      READ(1,*) AA,SRFALB       !fixed planetary albedo (0.2)
      READ(1,*) AA,SOLCON       !SOLCON=S/So
      READ(1,*) AA,dtmax        !maximum time step allowed (seconds)
      READ(1,*) AA,CO2MAX
      READ(1,*) AA, IMET        ! IMET (flag 0 or 1)
      READ(1,*) AA, IMETETH     ! IMETETH (flag 0 or 1)
      READ(1,*) AA, nga
      READ(1,*) AA, IHAZE       ! IHAZE (flag 0 or 1)
      READ(1,*) AA, ihztype
      READ(1,*) AA, icealbedo
      READ(1,*) AA, INVERSE
      READ(1,*) AA, FRAK        !can get a fractal haze without being coupled now
!   cloudy-CLIMA variables added by JDW
      read(1,*) AA,doEddy
      read(1,*) AA,fcloud
      read(1,*) AA,kz_min
      read(1,*) AA,Crainf 
      read(1,*) AA,Csig
      read(1,*) AA,supsat
      read(1,*) AA,upatm_mix
      read(1,*) AA,fc_minf
!   getting initial mixing ratios
      IF (ICOUPLE.eq.0) THEN
            OPEN (unit=114,file= DIRINOUT//'/mixing_ratios.dat')
      ELSE 
            OPEN (unit=114,file= 'COUPLE/mixing_ratios.dat')
      END IF

!   reading in photo parameters (not needed for MEAC coupling) (removed)

C **** Read the Gauss points and weights for the solar zenith angle int
      call data_grabber(xi,wi,ngauss)

! Reading the atmospheric composition from mixing_ratios.dat
      READ(114,*) FAR                  !Argon
      READ(114,*) FCH4                 !Methane
      READ(114,*) FC2H6                !Ethane
      READ(114,*) FCO2                 !Carbon dioxide
      READ(114,*) FN2                  !Nitrogen - added Nitrogen mixing ratio c-rr 6/5/2012
      READ(114,*) FO2                  !Oxygen
      READ(114,*) FH22                 ! c-rr 5/29/2012 added H2 mixing ratio
      READ(114,*) FNO2                 !Nitrogen dioxide
      READ(114,*) Jcold                !Tropopause layer

!   Reading mixing ratio profiles     CC2025
      DO K=1, ND
            READ(500,*) AO2(K)                  ! Oxygen
            READ(501,*) AN2(K)
            READ(502,*) AH2(K)
            READ(503,*) FI(2,K)                ! CO2
c            ACO2(K) = FI(2,K)
            READ(504,*) FI(3,K)                ! CH4
c            ACH4(K) = FI(3,K)
            READ(505,*) FI(5,K)               ! C2H6
c            AC2H6(K) = FI(5,K)
            READ(506,*) FI(1,K)                ! H2O
c            AH2O(K) = FI(1,K)
c            FIXEDWATER(K) = FI(1,K)            ! Writing water mixing ratio
            READ(507,*) FI(4,K)                 ! Ozone
c            AO3(K) = FI(4,K)
      ENDDO

!   Calculate new FCO2 (removed)

!   methane/ethane flags
      IF ((IMET.eq.0).and.(IMETETH.eq.0)) THEN
            FCH4=1.e-60
            FC2H6=1.e-60
      ENDIF

      IF ((IMET.eq.1).and.(IMETETH.eq.0)) THEN
            FC2H6=1.e-60
      ENDIF

!   oxygen flag
      IO2 = 0
      IF (AO2(1).ge.1e-40)IO2 = 1

c-rr Noncondensible molecular weight of the atmosphere when CO2 is condensing (for a colder planet)          5/3/2011
      DM2 = 28.*FN2 + 32.*FO2 + 40.*FAR + 16.*FCH4 
     & + 46.*FNO2+ 2.*FH22 ! c-rr 5/29/2012 added H2 mixing ratio

c jfk DM is the noncondensible molecular weight when CO2 is not condensing
      DM = 44.*FCO2 + (1.-FCO2)*DM2

      LAST = 0
      AMU0 = COS(ZY * PI/180.)

C   CONSTANT FACTORS (cgs)
      BCON = 2.*HP/C/C
      HK = HP/BK
      BKM = BK/(SM*G)
      ND1 = ND - 1

      R = 1.9872
      P0P = 6.103E-3
      T0P = 273.15
      SUBL = 677.

c  TRIPLE POINT PARAMETERS FOR CO2
      PC0 = 5.179
      TC0 = 216.56
      VAPCL0 = 83.2765
      SUBCL0 = 130.893
      DLVCDT = - 0.4817
      DLSCDT = - 0.1732
      CCL = 0.5
      CCS = 0.3

C Read Solar Data
      CALL READSOL

c Choosing a star
      CALL pickstar(STARR,SOLINT)
C try to accelerate ir.f, von Paris, 21/04/2006
      CALL IREXPSUMS
c Reading an initial temperature and water profile
998   FORMAT(3x,F16.12,7x,E22.15)

      IF(IUP.EQ.0) THEN
        DO J = 1,ND
          READ(11,998) T(J), FSAVE(J)
        END DO
        TG=T(ND)
      ENDIF


c Reading the ozone and water from the photochemical model (removed)
 
c  Interpolate the grid from the photochemical model to the grid of the climate model (removed)
       
C  Initialize pressure grid
      IF(IUP.EQ.1) TG = TG0
      CALL GRID(P0,FAC,ZCON,Z,ZF,DZ)

c  Reading the US Standard Atmosphere ozone profile       
c      if(IO3.eq.1.and.ICOUPLE.eq.0) CALL OZONE(FI,P)

!   removed block of commented-out code

C   CONVERT FREQUENCIES TO UNITS OF 1/SEC AND COMPUTE WEIGHTING FACTORS
      DO N=1,NF
        AV(N) = C*AV(N)
      END DO
      W(1) = AV(1)
      DO N=2,NF
        W(N) = AV(N) - AV(N-1)
      END DO
C
C   CENTER FREQUENCIES IN MIDDLE OF INTERVALS AND COMPUTE WAVELENGTHS
      SAV = 0.
      DO N=1,NF
        AVOLD(N) = AV(N)/C
        SAV2 = AV(N)
        AV(N) = 0.5*(AV(N) + SAV)
        LAM(N) = 3.E14/AV(N)
        print*,N,3.E14/AV(N),ALAMBDA(N)
        SAV = SAV2
      END DO

c Constructing temperature and water profiles in case they are not provided
       IF(IUP.EQ.1) THEN
          JCOLD = 1
          CALL PROFILE(TSTRAT,P,T,DZ,FSAVE,FCO2V,BETA,JCOLD,
     &    IDRY,FLAGCONVEC)
       ENDIF

c Building the water profile
      if(ICOUPLE.eq.0)then
        DO J = 1,ND
!            FI(1,J)=FSAVE(J)!eddyqt(J,1)-eddyqc(J,1)!FSAVE(J)
            FI(1,J)=FSAVE(J)    ! saving initial water profile  CC2025
        END DO
       
c jfk 6/25/08 Added four lines below
      IF(IMW.EQ.2) THEN
            DO J=1,JCOLD
                FI(1,J) = 4.0e-6 !JDW Changed in July 2021
            END DO
      END IF

      else
            DO J=1,ND
                CALL SATRAT(T(J),PSAT)
                FSATURATION(J) = (PSAT/P(J))*RELHUM(P(J)) !1st water vapor build, before loop JDW 10/2022
            END DO

       !first try to make JCOLD more sensible - giada
            IF (IUP.EQ.0.and.P(ND).GE.0.93) THEN !it's starting w/ fresh JCOLD otherwise ! EWS- and pressure is high enough
                IF (FSATURATION(1).GT.1) THEN !test if > 1 at start of grid
                    JCOLD_NEW = -1
                    DO j =1, ND
                        IF ((JCOLD_NEW .EQ. -1).and.(FSATURATION(J).LT.1)) THEN
                            JCOLD_NEW = J
                        end if
                    end do
                    !Update JCOLD if needed
                    If (JCOLD_NEW.NE.-1) THEN
                        JCOLD = JCOLD_NEW
                        eJCOLD=JCOLD
                    end if
                end if
                JCOLD = max(JCOLD,13) !EWS ensure JCOLD isn't too small.
                eJCOLD=JCOLD
            end if
       
!   removed code block here

c jkf 6/26/08 Change H2O initialization in the stratosphere
            do j=1,jcold
                if (imw.eq.2) FI(1,j) = 4.0e-6!FI(1,JCOLD)
                !fi(1,j) = FI(1,JCOLD)!4.0e-6 !JDW Changed in July 2021  !2nd water volume mixing build 10/2022 JDW
            end do
      endif
 

      tclim=T
      DO J=1,ND                   !   why the 2?
            PF1(J) = PF(J)*1.E6      !PF1 in dyn/cm^2 PF is pressure at flux grid points
            TOLD(J) = T(J)
            IF(IUP.EQ.1) FI(2,J)=FCO2V(J)
      ENDDO

      do j=1,nd
            fsave(j) = FI(1,j) !Re-save the FI water
      end do

      
C *** Initial time step
      dt0=5.e3
      IFLAGTIME = 0
      TIME = 0.
      tclim=T
    
c  Altitude calculation
      CALL ALTITUDE(NST,T,FI,DZ)
      DO J=1,ND
            WRITE(700,*) ALT(J)
      ENDDO

c Reading the ozone and water from the photochemical model
     
      IF(ICOUPLE.EQ.1) THEN
            DO JREAD=1,NZ  !number of layers in photochem code
                READ(222,*) temp_alt(JREAD),PRESS(JREAD),O3(JREAD),
     &                 water(JREAD),CH4(JREAD), CO2(JREAD), 
     &                 ethane(JREAD)
                temp_alt(JREAD)=temp_alt(JREAD)/1.0e5
            END DO
            FH2O=water(1)
            FCH4=CH4(1)
            FO3=O3(1)
            FCO2=CO2(1)
            FC2H6 = ethane(1)
            IF(FC2H6.LT.1.e-60) FC2H6 = 1.e-60 !!! Debug to prevent memory underflow issues - Eddie (8/3/2015)
            IF(FCH4.LT.1.e-60) FCH4 = 1.e-60   !!! Note that these are read whether or not IMETH or IEMETH flags are set
            print *, 'FC2H6 is ', FC2H6

c  Interpolate the grid from the photochemical model to the grid of the
c  climate model

            CALL INPUT_INTERP(temp_alt, water, O3, CH4, CO2, ethane, Jcold,
     &          T, FI)
      ENDIF
      close(222)

c  Reading the US Standard Atmosphere ozone profile
      if(IO3.eq.1.and.ICOUPLE.eq.0) then
            CALL OZONE(FI,P)
      endif

c Aerosol calculation (commented when not used)
      CALL AERABSDATA(FRAK, ihztype)
      CALL GRIDAER(ICOUPLE, IHAZE)
      CALL INTERPAR1(RAER)

C***********************************************************
C ****************** START ITERATIVE LOOP *******************
      DO 40 NST=1,NSTEPS
C************************************************************
      print *, 'TIME STEP = ', NST
      ITROP = 1
      TIME = Time + dt0

c Saving the former temperature-pressure profile

!   block removed

      tclim=T

      do i=1,ND
            alt_convec(i)=FLAGCONVEC(i)*alt(i)
      enddo

!   block removed

C     Calculate the atmospheric molecular weight at the surface with the noncondensibles
C     Use this value and pass it to Eddysed inside of Calleddy.f JDW 2021
       !if(NST .eq. 1)then

!   block removed

      mw_atmos=FI(2,nd)*44.01 + (1.0-FI(2,nd))*(ACH4(nd)*16.04 +
     & AO3(nd)*48 +  AO2(nd)*32.0 + FNO2*46.0055 + AC2H6(nd)*30.06904
     & + AH2(nd)*2.01588 + FAR*39.948 +FN2*28.0134) ! Do i need to instead caluclate the moist molecular weight? JDW 2023
c       above is only calculated for nd=101?
C ==========================================================================
C     Calling A&M 2001 Cloud Model
C ==========================================================================
      do J=1,ND
            FI_cloudy(1,J)=AH2O(J)
      enddo

      ! In Cloud-free clima, the Inverse climate model never recalculates water vapor, from the initial altitude setup. 
      ! Here we recalculate it with EddySed. 
      
      eJCOLD=JCOLD    !The problem might be that JCOLD is broken when coupled to EddySed JDW 2022
      laststep=.false.
      call calleddy(JCOLD,laststep,ALT)
      call relhum_vect((eddyqt(:,1)-eddyqc(:,1)),T,
     & P,mw_atmos,relhum_vec)
      
!   block removed

C ===============================================================================================
C                 Re-caluclating water vapor, based on the relhum_vec data from EddySed. 
C ===============================================================================================

!   huuuuuuuuge block removed

      DO J=1,ND
            call satrat(T(J),psat)
      enddo

      DO J=1,ND-1
            if(doEddy)then
                FI_cloudy(1,J)=AH2O(J)!eddyqt(J,1)-eddyqc(J,1) !+ FI(1,J)
            else
                FI_cloudy(1,J)=AH2O(J)
            endif

            FI_cloudy(1,nd)=FI_cloudy(1,nd-1)

            FI_cloudy(2,J)=FI(2,J)
            FI_cloudy(3,J)=ACH4(J)
            FI_cloudy(4,J)=AO3(J)
            FI_cloudy(5,J)=AO2(J)
            FI_cloudy(6,J)=FNO2

            FI_cloudy(2,nd)=FI(2,nd)
            FI_cloudy(3,nd)=ACH4(nd)
            FI_cloudy(4,nd)=AO3(nd)
            FI_cloudy(5,nd)=AO2(nd)
            FI_cloudy(6,nd)=FNO2
        

!            do I=2,NS1
!                  FI_cloudy(I,J)=FI(I,J)
!                  FI_cloudy(I,nd)=FI(I,nd)
!           enddo
      enddo

!     Initializing FNC c-rr 6/7/2012
      do J = 1,ND
            FNC(J) = 0.0  
            !FNC_cloudy(J) =0.0
      enddo
      do J = 1, ND
            FNC(J) = 1. - AH2O(J) - FI(2,J)   ! Added initial FNC array c-rr 6/7/2012
      enddo
      
      
c Initial non-condensible mixing ratio at surface (used for write statement in output file) 6/7/2011    
      FNCI = FNC(ND) ! c-rr 5/29/2012 added H2 mixing ratio
      do J=1,ND
            FNC_cloudy(J) =0.0
      enddo
C Be sure to check huge bug here ^^ KLUGE. Especially when looking for convergence with same water mixing ratios. JDW
      do J=1,ND
            FNC_cloudy(J)=1. - FI_cloudy(1,J) - FI_cloudy(2,J)!-FI_cloudy(1,J)-FI_cloudy(2,J) !FI goes to a NaN eventually? JDW
      enddo

!   block removed

c       GASCON is only called once per loop
      CALL GASCON(T,PF,FO2,FH22,FI,FI_cloudy,FNC, !Gascon may need to be looked into JDW 2021
     & FNC_cloudy,CGAS,CGAS_cloudy,NST)
      
c-rr gna  Created IRME.F (IR clone with methane and ethane loops turned on). When there is methane call IRM instead of IR. 5/2/2011 

      IF (IMET.eq.0) THEN
            CALL IR(T,PF,P,FNC,FNC_cloudy,CGAS,CGAS_cloudy)
      ENDIF
      IF ((IMET.eq.1).and.(IMETETH.eq.0)) THEN
c            CALL IRM(T,PF,P,FNC,CGAS) ! Passes FNC to IRM c-rr 6/7/2012
            CALL IRM(T,PF,P,FNC,FNC_cloudy,CGAS,CGAS_cloudy)
      ENDIF
      IF (IMETETH.eq.1) THEN
c            CALL IRME(T,PF,P,FNC,CGAS) ! Passes FNC to IRM c-rr 6/7/2012
            CALL IRME(T,PF,P,FNC,FNC_cloudy,CGAS,CGAS_cloudy)
      ENDIF

!   block removed

      IF (NST .EQ. NSTEPS) LAST = 1
C =================================================================
C  Solar code
c =================================================================
c  approximating the solar zenith angle with a gaussian summation

      do j = 1, nd
            fdnsoltot_clear(j) = 0.
            fupsoltot_clear(j) = 0.
            fdnsoltot_cloudy(j) = 0.
            fupsoltot_cloudy(j) = 0.
      enddo

C Find the right row in the matrix
      do i=1,11
            isave = i
            if (ngauss(i).eq.nga) exit
      enddo
C  isave holds the correct row number for the matrix

      do k=1,nga
            amu0 = xi(isave,k)
            zy = acos(amu0)*180./3.14159
c-rr    setting zenith angle to 60 degrees when nga = 1. WHY JDW
            if (nga.eq.1)then
                amu0 = 0.5
                zy=60.
            endif

            weightt = wi(isave,k)

C Heat capacity calculation
            DO J=1,ND-1
c-rr 3/30/11 The new CPCO2 and CPN2 curve fit equations
                CPCO2(J) = 5.89 + 6.06E-3*T(J) + 2.39E-5*T(J)*T(J)
     &              -3.44E-8*T(J)*T(J)*T(J)
c        if(j.eq.1)print *, 'CPCO2new=', CPCO2(J), T(J)
        
                CPN2(J) = 6.76 + 6.06E-4*T(J) + 1.3E-7*T(J)*T(J)
                CPO2(J) = 7.47 -4.84E-3*T(J) + 1.38E-5*T(J)*T(J)
     &              -8.73E-9*T(J)*T(J)*T(J) - 1.76E-9/T(J)/T(J)
                CPH2(J) = 7.17e-11*T(J)*T(J)*T(J)*T(J)
     &              -1.0e-07*T(J)*T(J)*T(J) + 4.77E-05*T(J)*T(J)
     &              -8.10E-03*T(J) + 7.17
                CPH2O(J) = 7.46 +4.52E-3*T(J)-1.38E-5*T(J)*T(J)
     &              + 1.74E-08*T(J)*T(J)*T(J)

!   old curve fit code block removed
      
                CPO2(J) = AMAX1(CPO2(J),CPN2(J))
          
!   block removed

                CpNC = AN2(K)*CPN2(J) + AO2(J)*CPO2(J) + FAR*4.97 +FCH4*8.3

c Total heat capacity
                CPN(J) = AH2O(J)*CPH2O(J)+FI(2,J)*CPCO2(J) + FNC(J)*CPNC

                CPN_cloudy(J)=FI_cloudy(1,J)*CPH2O(J) + FI_cloudy(2,J)
     &              *CPCO2(J)+ FNC_cloudy(J)*CPNC !CPN Cloudy added for cloud atmosphere column JDW
    
!   block removed

C since CPN is in calories/mol/K we should convert them to erg/g/K
                CPNT(J) = CPN(J)*4.18*1.E7/DM
            ENDDO

C   Surface heat capacity (assumes a 50 cm deep ocean mixed layer)
c   Units erg/K/cm^2 

            CPNT(ND) = 50.* 4.18*1.E7
        
!   block removed

C     =================================================================================
C     Below I need to update the shortwave subroutines for the cloud treatment. Right now, 
C     only SOLOROX.f is the subroutine that is operational. JDW
C     ==================================================================================
      !print*,'Calling Solar' !Each Solar code is called six times for the weightt* array
            IF ((IMET.eq.1).and.(IO2.eq.1))THEN
                CALL SOLARMOX(T,LAST,FNC,FNC_cloudy,NST)
                print *, 'Called SOLARMOX'
            ELSEIF((IMET.eq.1).and.(IO2.eq.0))THEN
                CALL SOLARM(T,LAST,FNC,NST)             ! Work on this one JDW 2023
                print *, 'Called SOLARM'
            ELSEIF((IMET.eq.0).and.(IO2.eq.1))THEN
                CALL SOLAROX(T,LAST,FNC,FNC_cloudy,NST) ! This one works !jdw
                print *, 'Called SOLAROX'
            ELSEIF((IMET.eq.0).and.(IO2.eq.0))THEN
                CALL SOLAR(T,LAST,FNC,NST)
                print *, 'Called SOLAR'
            ENDIF
            print *, fdnsol_cloudy(1)
        
C     =====================================================================
            do j=1,nd !JDW

            fdnsoltot_cloudy(j) = fdnsol_cloudy(j)*weightt +
     &          fdnsoltot_cloudy(j)

            fupsoltot_cloudy(j) = fupsol_cloudy(j)*weightt +
     &          fupsoltot_cloudy(j)! + fupsol_cloudy(j)*weightt
          
            fdnsoltot_clear(j) = fdnsoltot_clear(j)
     &          + fdnsol_clear(j)*weightt

            fupsoltot_clear(j) = fupsoltot_clear(j)
     &          + fupsol_clear(j)*weightt

            enddo
      enddo

c =================================================================
c IR and SOLAR fluxes (erg/cm^2/s)
      
      J=0
      DO 31 J=1,ND
            FDNSOL(J) = ((1-fcloud)*SOLCON*0.5 * FDNSOLTOT_CLEAR(J))  !JDW !Also uncomment when no longer tidal locking JDW 2023
     &          + (fcloud*SOLCON*0.5 * FDNSOLTOT_CLOUDY(J))

            FUPSOL(J) = ((1-fcloud)*SOLCON*0.5 * FUPSOLTOT_CLEAR(J)) !Uncomment when no longer tidal locking JDW 2023
     &          + (fcloud*SOLCON*0.5 * FUPSOLTOT_CLOUDY(J))

            FDNIR(J)=(1-fcloud)*FDNIR_clear(J)
     &          + fcloud*FDNIR_cloudy(J)

            FUPIR(J)=(1-fcloud)*FUPIR_clear(J)        !The fluxes here are breaking.
     &          + fcloud*FUPIR_cloudy(J)

            FTOTAL(J) = FDNSOL(J)-FUPSOL(J)
     &          + FDNIR(J)-FUPIR(J)
     
      FTIR(J) = FDNIR(J)-FUPIR(J)
      FTSO(J) = FDNSOL(J)-FUPSOL(J)

  31  CONTINUE

      ALBP = FUPSOL(1)/FDNSOL(1)
      SEFF = abs(FTIR(1)/FTSO(1))      !to print out Seff, c-rr 4/21/2011
      print *, 'FTIR= ', FTIR(1)
      print *, 'FUPIR',FUPIR(1)
      print *, 'FDNIR',FDNIR(1)
      print *, 'FTSO= ', FTSO(1)
      PRINT *, 'Seff=',SEFF
      print *, 'JCOLD',JCOLD
      PRINT 166,ALBP

 166  FORMAT(/1X,"PLANETARY ALBEDO:  ALBP = ",F6.4)

C BEGIN INVERSE SKIPS


      IF(INVERSE.EQ.0) THEN !only do if not wanting inverse calculations

C New temperature calculation for all layers from radiative equilibrum !TN is layer, PF is layer.
            DO 41 J=1,ND-1
                TN(J)=T(J)-(FTOTAL(J+1)-FTOTAL(J))*dt0*GNEW(J)/CPNT(J)
     &              /(PF1(J+1)-PF1(J))

                TCOOL(J)=-(FTIR(J+1)-FTIR(J))*GNEW(J)/CPNT(J)
     &              /(PF1(J+1)-PF1(J))*86400.
                THEAT(J)=-(FTSO(J+1)-FTSO(J))*GNEW(J)/CPNT(J)
     &              /(PF1(J+1)-PF1(J))*86400.
  41  CONTINUE

c New surface temperature from radiative equilibrum
      select case (ICONSERV)
            case(1)
            TN(ND)=T(ND)+FTOTAL(ND)*dt0/CPNT(ND)
            TCOOL(ND)= FTIR(ND)*86400./CPNT(ND)
            THEAT(ND)= FTSO(ND)*86400./CPNT(ND)
            case(0)

            TN(ND-1) = T(ND-1)+FTOTAL(1)/(PF1(ND)-PF1(ND-1))*
     &          GNEW(ND-1)/CPNT(ND-1)*dt0
!   something may be going wrong here
            CALL SATRAT(TN(ND-1),PSAT)
            FI(1,ND-1)=relhum_vec(ND-1)*PSAT/P(ND-1)!*((PF(J)/PG)-0.02)/0.98
      end select

* Total heating rate
      do j=1,ND
            HEATNET(j)=THEAT(j)+TCOOL(j)
      enddo
c      print *, 'Hello5'
c-as TRAD is defined for printing and diagnostic purposes
      DO J=1,ND
            TRAD(J)=TN(J)
      ENDDO
      
cTEMPORARY DEBUGGING STATEMENT**************
C      print *,'Calling output early'    
C     GOTO 571
c Calculating tropospheric temperatures

      select case(ICONSERV)
*** Non strict time-stepping model
            case(0)

            DIVF(1) = FTOTAL(1)/FUPIR(1)
            TN(ND) = T(ND) * (1. + 0.1*DIVF(1))

            IF (TN(ND) .LT. T(ND)) GO TO 1400 !This will skip the call to convec JDW 2021 Why is this a thing?

            JCONV=ND
            ITROP=1

            DO J1=ND, 2, -1
                FLAGCONVEC(J1) =0.!. JDW2021
                T1=TN(J1)
                DZP = DZ(J1)
                P1 = P(J1)
                P2 = P(J1-1)
                FC1 =FI(2,J1)
                FH1 =AH2O(J1)
      
                CALL CONVEC(T1,T2,P1,P2,FH1,FH2,FC1,FC2,DZP,ITROP,cflag,
     &              Idry, imco2)


                IF(IO3.EQ.1 .AND. ALT(J1) .GT. 40.) GOTO 1401   ! Skip convection if ozone beyond 40km

                IF (TN(J1-1) .LE. T2 + 0.001*T2) THEN
                    TN(J1-1) = T2
!                    FI(2,J1-1) = FCO2
                    FLAGCONVEC(J1) = cflag
                    JCONV = J1
                END IF
            END DO

            FLAGCONVEC(ND) = 1.
            ND1 = ND-1
            DO J=ND1,1,-1
!                FI(2,J) = AMIN1(FI(2,J),FI(2,J+1))
            END DO

            GO TO 1401
c
c   If surface temperature is decreasing, then adjust all temperatures
c   below the cold trap downward by the same amount. This ensures that 
c   the upward IR flux will decrease as surface temperature decreases.

 1400       CONTINUE
            DTSURF = T(ND) - TN(ND)
            DO J=JCOLD, ND-1
                TN(J) = T(J) - DTSURF
            ENDDO
            JCONV = JCOLD    ! 5/23/2011 So that it knows what JCONV is when it is not convecting !Crap
 1401       CONTINUE
** End of the non strict time-stepping model
      
*** Here the temperatures are calculated conserving energy on each ***************************************************************************************
*** layer
            case(1)
            DO ITER=1,20        !starting convective adjustment
                ITROP = 1
                imco2=0
                JCONV=ND
                HC1=CPNT(ND-1)*(PF1(ND)-PF1(ND-1))/g
                DZP = DZ(ND)
                T1 = TN(ND)
                P1 = P(ND)
                P2 = P(ND-1)
                FH1 = AH2O(ND)
                FC1 = FI(2,ND)
                CALL CONVEC(T1,TadND1,P1,P2,FH1,FH2,FC1,FC2,DZP,1,cflag,
     &              Idry, imco2)
!                FI(1,ND-1)=FH2*relhum_vec(ND-1)
                TnewND=(CPNT(ND)*TN(ND)-HC1*(TadND1-TN(ND)-TN(ND-1)))/
     &              (HC1+CPNT(ND))
                TnewND1=TadND1-TN(ND)+TnewND
                TN(ND-1)=TnewND1
                TN(ND)=TnewND
                if(TN(J).lt. 0.0) stop '"TN" is negative 0'
                FLAGCONVEC(ND)= 1.
                FLAGCONVEC(1) = 0.
                imco2=0
c-as This part has been modified to consider energy balance in each
c-as  layer, as Hilary Justh did it (oct-2003)  
***** CONVECTIVE ADJUSTMENT (considering energy balance and convection)

                DO J1=ND-1,2,-1
                    T1 = TN(J1)
                    DZP = DZ(J1)
                    P1 = P(J1)
                    P2 = P(J1-1)
                    FH = AH2O(J1) !Water mixing ratio, bottom two layers are changed above!
                    FC1 = FI(2,J1)
        !print*,'calling cflag3',cflag
                    CALL CONVEC(T1,T2,P1,P2,FH,FH2,FC1,FC2,DZP,ITROP,cflag,
     &                  Idry, imco2)    ! something going wrong here
                    IF (TN(J1-1).LE.T2) THEN

                        FLAGCONVEC(J1) = cflag
                        IF(cflag.eq.1.or.cflag.eq.3) JCONV=J1
                        DELPCP1=(PF1(J1+1)-PF1(J1))*CPNT(J1)
                        DELPCP2=(PF1(J1)-PF1(J1-1))*CPNT(J1-1)
                        T2P=TN(J1-1)*(DELPCP2/(DELPCP1+DELPCP2))+(T2)*
     &                      (DELPCP1/(DELPCP1+DELPCP2))
                        T1P=T1-T2+T2P
                        print *, J1,TN(J1-1),T2P,TN(J1),T1P
                        TN(J1-1)=T2P
                        TN(J1)=T1P
                        if(TN(J).lt. 0.0) stop '"TN" is negative1'
!                        FI(2,J1-1) = FCO2  ! jiggered again

                    ELSE
                        ITROP=0
                        FLAGCONVEC(J1)=0.
                    ENDIF
                ENDDO

1403        CONTINUE
            ENDDO          ! End of convection adjustment loop
            ND1 = ND-1
            DO J=ND1,1,-1
!                FI(2,J) = AMIN1(FI(2,J),FI(2,J+1))
            END DO
      end select
      
c Water recalculation

      DO J=1,ND
            CALL SATRAT(TN(J),PSAT)
            FSATUR(J) = (PSAT/P(J))
      ENDDO
      FCT=FSATUR(ND)
      DO J=ND-1,1, -1
            FCT=AMIN1(FCT,FSATUR(J+1))
      ENDDO
      JCOLD = 1
c
c jfk 7/15/08 Simplify the logic for finding the cold trap

      DO J = ND-1,2,-1
            JCOLD = J
            IF (FSATUR(J-1) .GT. FSATUR(J)) GO TO 3100
          !IF (T(J-1) .GT. T(J)) GO TO 3100
      END DO
 3100 CONTINUE

      eJCOLD=JCOLD !JDW
c Water from the cold trap to the ground
      if(doEddy)then
            DO J = JCOLD, ND
!                FI(1,J) = FSATUR(J)*relhum_vec(J) !relhum_vec(J) is broken!! JDW
                !FI(1,J)=FSATUR(J)
                if(imw.eq.2) FI(1,J) = amax1(FI(1,J),4.e-6)
!                if(FI(1,J).lt.upatm_mix) FI(1,J)=upatm_mix
            enddo
      else
            DO J = JCOLD, ND
!                FI(1,J) = FSATUR(J)*RELHUM(P(J))
                if(imw.eq.2) FI(1,J) = amax1(FI(1,J),4.e-6)
            enddo

c Water from the cold trap to the top (if it is used in the coupled
c mode these values are given by the photochemical code)
            if(ICOUPLE.eq.0)then
                DO J = JCOLD-1, 1, -1
!                    FI(1,J)= 4.e-6!FI(1,JCOLD)
                END DO
            endif
      endif

C-KK To smooth over the profile around JCOLD.
      sum = 2*AH2O((JCOLD-1)) + AH2O((JCOLD+1)) + 2*AH2O(JCOLD) !JDW previously had commented out 2021
      AH2O(JCOLD) = sum/5.

C Smoothing over the previous temperatures
      DO J=1,ND
        TOLD(J) = T(J)
      ENDDO
      
C Smoothing of temperature profile conserving energy
      if(ICONSERV.eq.1) then
c jfk Replace the DO logic below to make sure that the smoothing does
c  not occur when CO2 is condensing (FLAGCONVEC=3) 
            DO J=2,JCOLD
                Tj1 = 0.5*TN(J) + 0.25*(TN(J-1) + TN(J+1))
                CPP0=(PF1(J)-PF1(J-1))*CPNT(J-1)
                CPP1=(PF1(J+1)-PF1(J))*CPNT(J)
                CPP2=(PF1(J+2)-PF1(J+1))*CPNT(J+1)
                En1=CPP0*TN(J-1)+CPP1*TN(J)+CPP2*TN(J+1)
                DELT1=Tj1-TN(J)
                DELT2=-(CPP1/(CPP0+CPP2))*DELT1
                TN(J+1)=TN(J+1)+DELT2
                TN(J-1)=TN(J-1)+DELT2
                TN(J)=Tj1
                print *, J, TN(J-1), TN(J), TN(J+1)
                if(TN(J).lt. 0.0) stop '"TN" is negative2'
            END DO
      endif

C  Diagnostics parameters
      print *, 'Printing delta temperature array'
      DO J=1,ND
            DELT(J) = (TN(J)-TOLD(J))
            DELTRAD(J) =TRAD(J)-TOLD(J)
            T(J) = TN(J)
            DIVF(J) = FTOTAL(J)/FUPIR(J)
!            print *, J, T(J), TN(J)
            print *, J, DELT(J)
      ENDDO

c Smoothing the temperature profile in the non-strict time step case
c-jdh **check on this when comparing "conserving" vs "non-conserving"**
      if(ICONSERV.eq.0) then      
c  Replace the DO logic below to make sure that the smoothing does
c  not occur when CO2 is condensing (FLAGCONVEC=3)
            DO J=2,JCOLD
                T(J) = 0.5*TN(J) + 0.25*(TN(J-1) + TN(J+1))
                tclim(J) = 0.5*TN(J) + 0.25*(TN(J-1) + TN(J+1))
!                tclim(J)=TN(J)           ! Duplicate line? cc2025
            END DO
      endif
      print*,'Surface temperature=',T(ND)

c adjust albedo based on ice-albedo feedback
c parameterization added by Giada based on Charnay et al 2014 
      if (icealbedo.eq.1) then 
            IF (T(ND).LT.240.) SRFALB = 0.65
            IF (T(ND).GT.290.) SRFALB = 0.30
            IF (T(ND).GE.240. .AND. T(ND).LE.290.) then
                SRFALB=0.65+(0.3-0.65)*( (T(ND)-240)/(290-240) )**0.37
            end if
            print *, 'Surface albedo=', SRFALB
      end if
c Adjusting the time stepper
      DTS = dt0
      CHG = 0.
      DO J=2,ND-1
            REL = ABS(DELT(J)/TOLD(J))
            CHG = AMAX1(CHG,REL)
      END DO
      IF (CHG.LT.0.01) dt0 = DTS*1.5
      IF (CHG.LT.0.001) dt0 = DTS*5.
      IF (CHG.GT.0.02) dt0 = DTS/2.
      IF (dt0.GE.dtmax) dt0 = dtmax
            CALL ALTITUDE(NST,T,FI,DZ)
      END IF !end skipping for inverse model

C     Calculate the relative heating rates for the Longwave and shortwave
C     Shortwave first
C     Build appropriate input vector from fluxes JDW 2021
      do zz=1,nd
            Fnet_rhr_s(zz)=(FUPSOL(zz)-FDNSOL(zz))
            pressure_rhr(zz)=P(zz)
            altitude_rhr(zz)=ALT(zz)!*100000 JDW 2023
            c_p_rhr(zz)=CPNT(zz)
      enddo
      call radiative_heating_rate(G,c_p_rhr,Fnet_rhr_s,
     & pressure_rhr,altitude_rhr,rhr_solar)

      do zz=1,nd
            Fnet_rhr_t(zz)=(FUPIR(zz)-FDNIR(zz))
            pressure_rhr(zz)=P(zz)
            altitude_rhr(zz)=ALT(zz)!*100000 !JDW 2023
            c_p_rhr(zz)=CPNT(zz)
      enddo
      call radiative_heating_rate(G,c_p_rhr,Fnet_rhr_t,
     & pressure_rhr,altitude_rhr,rhr_thermal)

!   Resetting water mixing ratio
      DO K=1,ND
!            AH2O(K) = FIXEDWATER(K) ! just in case it gets modified somewhere
      ENDDO

!     convergence stuff
      IF (NST.LT.MINSTEPS) GOTO 6789
      DO K=1,ND
            IF (ABS(DELT(K)).GT.TCONV) GOTO 6789
      ENDDO
      GOTO 99999
6789  CONTINUE
C***********************************************************
c***  WRITING OUTPUT FILES
************************************************************
      IF(NST.EQ.1) THEN
            WRITE(98,*)
            WRITE(98,*) "   OUTPUT FILES FOR THE ",STARR
            WRITE(98,*)
            WRITE(98,555) SOLCON,FCH4*FNCI,FCO2,FO2*FNCI,FN2*FNCI,
     &          FH22*FNCI,FAR*FNCI,IO3,IUP
 555        format(1x,"Solar Constant= ",F5.3,3x,"F_CH4= ",1pe10.4,2x,
     &          "F_CO2= ",1pe10.4,2x,"F_O2= ",1pe10.4,2x,"F_N2= ",1pe10.4,2x,
     &          "F_H2= ", 1pe10.4,2x,"F_AR= ", 1pe10.4,2x,
     &          "IO3 = ",I2,3x,"IUP= ",I2)
            WRITE(98,556) ICONSERV,FAC,ND,SRFALB,G,IMW, INVERSE
 556        format(1x,'ICONSERV=',I2,2X,'FAC=',F4.1,2X,'ND=',I3,2X,'SRFALB='
     &          ,F6.3,2X,'G=',F6.1,2X,'IMW=',I2, 2X, 'INVERSE=', I2)
            WRITE(98,557) FNO2
 557        FORMAT(/1x,'FNO2 =',1pe10.3)
            WRITE(98,*)
      ENDIF
      nsteps2 = nsteps-2
      nsteps3 = nsteps-3
      IF(INVERSE.EQ.1)then
            JCONV=foundloc! Make sure this pulls the value of the top of the flagconvec array JDW 2021
            CHG=0.
      endif
      DIVFrms = 0.
      JC1 = JCONV-1 !JCONV is undefined
      DO J=1,JC1
            DIVFrms = DIVFrms + DIVF(J)*DIVF(J)
      END DO
      DIVFrms = SQRT(DIVFrms/JC1)

      WRITE(98,966) NST,JCONV,CHG,dt0,DIVF(1),DIVFrms,DELT(ND),T(ND)
 966  FORMAT(1x,"NST=",I6,1X,'JCONV=',I3,1x,'CHG=',1pe8.2,1x,"dt0=",
     & 1pe8.2,1X,"DIVF(1)=",1PE9.2,1X,"DIVFrms=",1PE9.2,1x,
     & "DT(ND)=",1PE9.2,1x,"T(ND)=",1PE10.4)
      if(nst .eq. NSTEPS) goto 1234
      if(nst.eq.1 .or. abs(SEFF -1.0 ).le. 0.0001 .and. nst > 40) then

 1234 continue
c       WRITE(98,965) NST,dt0,DIVF(1),FTOTAL(ND-1),FTIR(ND-1),
c    & FTSO(ND-1),DELT(ND),T(ND) 
c 965   FORMAT(1x,"NST=",I3,2X,"dt0=",1PE9.3,
c     & 2X,"DIVF(1)=",1PE12.5,2X,"Ftot(ND-1)=",1pe11.4,2x,"FtIR(ND-1)="
c     & ,1pe11.4,2x,"FtSol(ND-1)=",1pe11.4,/,1x,"DT(ND)=",1PE10.3,2x,
c     & "T(ND)=",1PE10.4)
c
      TIMEDAYS = TIME/24./3600.
      WRITE(98,567) TIME,TIMEDAYS
 567  FORMAT(/1X,'TIME=',1PE10.3,2X,'TIME IN DAYS =',E10.3)
      WRITE (98,166) ALBP
      WRITE(98,683)
        DO J=1,ND
c      WRITE(98,680) J,P(J),ALT(J),T(J),FLAGCONVEC(J),
c     & DELT(J),TOLD(J),FI(1,J),HEATNET(J),TCOOL(J),THEAT(J)
      WRITE(98,680) J,P(J),ALT(J),T(J),FLAGCONVEC(J),
     & DELT(J),TOLD(J),AH2O(J),FSAVE(J),AO3(J),TCOOL(J),THEAT(J)
       ENDDO
      WRITE(98,1683)
      Write (98,1112)
 1112 FORMAT(/1x,'FCO2')
      WRITE(98,1111) (FI(2,J),J=1,ND)
 1111 FORMAT(1x,1p10e9.2)
      WRITE(98,685)
         DO J=1,ND
      WRITE(98,684) J,PF(J),ALT(J),FTOTAL(J),FTIR(J),FDNIR(J),
     & FUPIR(J),FTSO(J),FDNSOL(J),FUPSOL(J),DIVF(J)
         ENDDO
      WRITE(98,*)   
      END IF        
  683  FORMAT(/2x,"J",5X,"P",9X,"ALT",9X,"T",8X,"CONVEC",
     & 7X,"DT",10X,"TOLD",8x,"FH20",
     &  7x,'FSAVE',8x,'FO3',8x,'TCOOL',7x,'THEAT') ! top of file
 1683  FORMAT(2x,"J",5X,"P",9X,"ALT",9X,"T",8X,"CONVEC",
     & 7X,"DT",10X,"TOLD",8x,"FH20",
     &  7x,'FSAVE',8x,'FO3',8x,'TCOOL',7x,'THEAT')  !bottom of file
 680  FORMAT(I3,3(1x,1PE10.4),1X,1PE9.2,2X,1PE11.4,7(1X,1PE11.4))
 685  FORMAT(/2x,"J",4X,"PF",9X,"ALT",7X,"FTOTAL",7X,"FTIR",7X,"FDNIR",
     & 7X,"FUPIR",7X,"FTSOL",7X,"FDNSOL",7X,"FUPSOL",7X,"DIVF")
 684  FORMAT(I3,2(1x,1PE10.4),8(1X,1PE11.4))
C
      DO J=1,ND
      TSAT = T(J)
      CALL SATCO2(TSAT,PSAT)
      PSATCO2(J) = PSAT
      END DO


      do j=1,nd   ! FH2O becomes old FSAVE for next time step 4/23/2012
      fsave(j) = AH2O(j)
      end do

      do j = 1,nd  ! redefines FNC for next time step c-rr 6/7/2012
      FNC(J) = 1. - AH2O(J) - FI(2,J)   ! Added initial FNC array c-rr 6/7/2012
      enddo
      if(nst.eq.NSTEPS)goto 12341
      if(abs(SEFF -1.0 ).le. 0.0001 .and. nst > 450) then  !0.00001
      !if(NST.eq.NSTEPS)then
12341 continue
      print*,'seff',SEFF

      goto 401
      endif
      write(2021,9661) NST
      write(20212,9661) NST
 9661   FORMAT(1x,"NST=",I6)
      write(2021,2023)
 2023 format(/2x,"J",2X,"FUPIR_clear",2X,"FUPSOLTOT_clear",2X,"FDNIR_clear",2X,"FDNSOLTOT_clear",
     & 2X,"FUPIR_cloudy",2X,"FUPSOLTOT_cloudy",2X,"FDNIR_cloudy",2X,"FDNSOLTOT_cloudy")

      write(20212,20232)
20232 format (/2x,"J",2X,"RHR Solar (K/s)",2X,"RHR Thermal (K/s)")
      do J=1,ND
            write(2021,2022) J,FUPIR_clear(J),SOLCON*0.5*FUPSOLTOT_clear(J),
     & FDNIR_clear(J),SOLCON*0.5*FDNSOLTOT_clear(J),FUPIR_cloudy(J),
     & SOLCON*0.5*FUPSOLTOT_cloudy(J),
     & FDNIR_cloudy(J),SOLCON*0.5*FDNSOLTOT_cloudy(J)

            write(20212,2022) J,rhr_solar(J),rhr_thermal(J)
      enddo

 2022 format(I3,*(1X,1P1E11.3))
***************************************************************
C   End of iterative loop
  40  CONTINUE
  401 continue
***************************************************************

      laststep = .true.
      call calleddy(eJCOLD,laststep,ALT)
      write(2021,2022)

      !couplesmart = .true.
      !Need to build rmix, instead of FI 
      !if(couplesmart) then 

      !print*,FNC,FC
C       COMMON/CBLOK/FO2,FN2,FCO2,FAR,FCH4,FC2H6,FNO2,FH2more 

C         READ(114,*) FAR                  !Argon
C         READ(114,*) FCH4                 !Methane
C         READ(114,*) FC2H6                !Ethane        
C         READ(114,*) FCO2                 !Carbon dioxide
C         READ(114,*) FN2                  !Nitrogen - added Nitrogen mixing ratio c-rr 6/5/2012        
C         READ(114,*) FO2                  !Oxygen        
C         READ(114,*) FH22                 ! c-rr 5/29/2012 added H2 mixing ratio
C         READ(114,*) FNO2                 !Nitrogen dioxide
c C   FI = SPECIES MIXING RATIOS   1 = water, 2 = co2, 3 = ch4, 4 = o3, 5 = ethane
      do zz=1,nd !Make sure the overall mixing ratios are not exceeding unity. Maybe scale everything to (1-FI(1,zz))*FI(j,*)


            srmix(1,zz)=AH2O(zz) !water 1 !Numbered according to HITRAN
            !print*,'FI(1,zz)',FI(1,zz)
            srmix(2,zz)=FI(2,zz) !CO2 2
            srmix(3,zz)=AO3(zz) !O3  3
            srmix(4,zz)=ACH4(zz) !CH4 4
            srmix(5,zz)=AO2(zz)  !O2      5
            srmix(6,zz)=FNO2!Nitrogen Dioxide 6
!             srmix(7,zz)=(1-FI(3,zz)-FI(4,zz)-FO2-FNO2-FH22-
!      & FI(5,zz)-FAR)*FN2 !N2
            srmix(7,zz)=FNC(zz) -ACH4(zz)-AO3(zz) !JDW testing 2021
     & - AO2(zz) - FNO2 -AC2H6(zz)-AH2(zz)-FAR  !Nitrogen 7
            if(srmix(7,zz).lt. 0.0) srmix(7,zz) = 0.0
            ! srmix(7,zz) = FN2
            srmix(8,zz)=AC2H6(zz) !C2H6 8
            srmix(9,zz)=AH2(zz) !H2      9
            srmix(10,zz)=FAR  !Argon    10
      enddo

c       Testing: making sure CO2 transfers properly CC2025
c      DO K=1,ND
c            PRINT *, FI(2,K),FI(2,K)
c      ENDDO

c       Write final mixing ratio profiles
      DO J=1,ND
            WRITE(600,*) AO2(J)
            WRITE(601,*) AN2(J)
            WRITE(602,*) AH2(J)
            WRITE(603,*) FI(2,J)
            WRITE(604,*) ACH4(J)
            WRITE(605,*) AC2H6(J)
            WRITE(606,*) AH2O(J)
            WRITE(607,*) AO3(J)
      ENDDO

        

**************************** Smart Output *****************************************
**************************** Change output filename below**************************
99999 CONTINUE
      if(ICOUPLE.eq.1) then
       print *, 'output photo'
       CALL OUTPUT_PHOTO(T, FI, water, ALT, nzp)
      endif

        WRITE(97,466)
 466  FORMAT(5X,'ALT',10X,'P',10X,'T',10X,'FH2O',11X,'O3',11X,
     2  'THEAT',8X,'TCOOL',8X,'PSATCO2',8x,'FCO2')
      DO J=1,ND
c jkl 6/27/08 Print out H2O from the initial profile
        WRITE(12,998) T(J),AH2O(J)
c        WRITE(97,467) ALT(J),P(J),T(J),FI(1,J),FI(4,J),THEAT(J),
c     &  TCOOL(J)
        WRITE(97,467) ALT(J),P(J),T(J),FSAVE(J),AO3(J),THEAT(J),
     &  TCOOL(J),PSATCO2(J),FI(2,J)
 467  FORMAT(1PE10.4,2X,1PE10.4,1X,1PE10.4,3X,1PE11.4,3X,1PE11.4,
     & 3X,1PE11.4,3X,1PE11.4,3x,1pe11.4,3x,1pe11.4) 
      END DO

      STOP

      END                 !end of the main program
      
*********************************************************************
      SUBROUTINE ALTITUDE(NST,T,FI,DZ)
      INCLUDE 'CLIMA/INCLUDE/header.inc'
      PARAMETER(NS1=5)       !gna: changed ns1 from 4 to 5
      COMMON/CONSS/C,BK,G,GNEW(ND),PI,SM,DM,DM2
      COMMON/ALTBLOK/DALT(ND-1),RADIUS(ND-1),PARTICLES(ND),RAER(ND),
     & ALT(ND)
      DIMENSION T(ND),FI(NS1,ND),DZ(ND)    

c-as  This subroutine calculates the altitude.
c-as  The water vapor was eliminated the first time this subroutine
c-as  is called (before the NSTEPS DO loop) in order to make easier 
c-as  the parameter translation to the photochemical model

      ALT(ND) = 0.
      DO J=ND-1,1,-1
       GNEW(J) = G*(RAD**2)/(RAD + ALT(J))**2  ! Turned Gravity into arrays 10/12/2012
       BKM = BK/(SM*GNEW(J))
       TA = 0.5*(T(J) + T(J+1))
       FH2O = 0.5 * (FI(1,J) + FI(1,J+1)) !smoothes over water mmr JDW
       FCO2J = 0.5* (FI(2,J) + FI(2,J+1))
       
       FNCA = 1. - FH2O - FCO2J
       AM = 18.*FH2O + 44.*FCO2J + DM2*FNCA  ! AM is the weight of entire parcel (noncondensible + condensible)
            

       IF(NST.lt.1) AM = DM2     
       BMG = BKM/AM
       ALT(J) = ALT(J+1) + BMG*TA*DZ(J+1)*1.E-5
       
       DALT(J) = ALT(J) - ALT(J+1)
      ENDDO
          
      RETURN
      END

C=============================================================================
        subroutine data_grabber(xi,wi,n)
        parameter(nrow=11)                        
        dimension xi(nrow,20), wi(nrow,20), n(nrow)        
100     format(2x, I2)                
200     format(F7.5,1x,F7.5)
c300     format (20(f7.5,1x))
c400     format (/)
c500     format (11i3)                  
        do i = 1,nrow
        read (66,100) n(i)
          do j=1,n(i)                                
          read(66,200) xi(i,j), wi(i,j)
          enddo
        enddo        
C        print*, 'n='
C        print 500,n
        !print 400        
        !print*,'the Abscissas(xi)'
        !print 300,((xi(i,j),j=1,n2), i=1,nrow)                
        !print 400
        !print*, 'the weight factors(wi)'
        !print 300,((wi(i,j),j=1,n2), i=1,nrow)


        do i=1,nrow
        sum = 0.
           do j=1,n(i)                                
        !print 20,j,xi(i,j),wi(i,j)
c20      format(5x,'j=',i2,1x,2f8.5)
          sum = sum + xi(i,j)*wi(i,j)
          enddo
C        print*,'n=',n(i),' sum =', sum
        enddo
        
        end
        
