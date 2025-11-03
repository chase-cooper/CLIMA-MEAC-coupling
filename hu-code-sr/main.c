#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "./scenario_library/Earth/planet_Earth_Full_T1986.h"
// #include "./scenario_library/Mars/planet_Mars_Full.h"

// #include "./scenario_library/Sun/CO2/planet_scenario_CO2r_emi.h"
// #include "./scenario_library/Sun/H2/planet_scenario_H2N2r2.h"
// #include "./scenario_library/Sun/N2/planet_scenario_N2r.h"

// #include "./scenario_library/Sun/CO2-Full/planet_scenario_CO2r_emi.h"
// #include "./scenario_library/Sun/H2-Full/planet_scenario_H2N2r2.h"
// #include "./scenario_library/Sun/N2-Full/planet_scenario_N2r.h"

// #include "./scenario_library/TRAPPIST-1/CO2-Full/planet_scenario_CO2r_emi.h"

#include "constant.h"
#include "routine.h"
#include "global.h"
#include "GetData.c"
#include "chemequil.c"
#include "Interpolation.c"
#include "RadTransfer.c"
#include "nrutil.h"
#include "nrutil.c"
#include "ReactionRate.c"
#include "ReactionRateM.c"
#include "ReactionRateT.c"
#include "Photorate.c"
#include "ChemReaction.c"
#include "lubksb.c"
#include "ludcmp.c"
#include "Convert.c"
#include "BTridiagonal.c"
#include "TPPara.c"
#include "TPScale.c"
#include "printout.c"
#include "printoutrate_noreverse.c"
#include "fallvelocity.c"
#include "saturationpressure.c"
#include "MeanOpacity.c"
#include "GreyTemp_Iso.c"
#include "RefIdx.c"
#include "planckmean.c"
#include "NewPressure_Iso.c"
#include "GlobalBalance.c"
/* #include "RadConv.c" */

/* external (global) variables */

double thickl;
double zl[zbin+1];
double pl[zbin+1];
double tl[zbin+1];
double MM[zbin+1];
double MMZ[zbin+1];
double wavelength[WaveBin+1];
double solar[WaveBin+1];
double crossr[WaveBin+1], crossa[3][WaveBin+1], sinab[3][WaveBin+1], asym[3][WaveBin+1];

double **opacCO2, **opacO2, **opacSO2, **opacH2O, **opacOH, **opacH2CO;
double **opacH2O2, **opacHO2, **opacH2S, **opacCO, **opacO3, **opacCH4; 
double **opacNH3;
double MeanCO2[zbin+1], MeanO2[zbin+1], MeanSO2[zbin+1], MeanH2O[zbin+1], MeanOH[zbin+1], MeanH2CO[zbin+1];
double MeanH2O2[zbin+1], MeanHO2[zbin+1], MeanH2S[zbin+1], MeanCO[zbin+1], MeanO3[zbin+1], MeanCH4[zbin+1];
double MeanNH3[zbin+1];	

double rainoutrate[zbin+1][NSP+1];
double Vesc[NSP+1], VFall[zbin+1];
double nsH2O[zbin+1], nsH2SO4[zbin+1], nsS8[zbin+1], tcondfH2O[zbin+1], tcondfH2SO4[zbin+1], tcondfS8[zbin+1];

double kk[zbin+1][NKin+1], kkM[zbin+1][NKinM+1], kkT[zbin+1][NKinT+1];
int    ReactionR[NKin+1][7], ReactionM[NKinM+1][5], ReactionP[NPho+1][9], ReactionT[NKinT+1][4];
double **DM, **dl, KE[zbin+1];

int    numr=0, numm=0, numt=0, nump=0, numx=0, numc=0, numf=0, numa=0, waternum=0, waterx=0;

double xx[zbin+1][NSP+1];

double mkv[zbin+1], Tnew[zbin+1], Pnew[zbin+1];

int main()
{
	  int s,i,ii,j,jj,jjj,k,nn,qytype,stdnum,iradmax;
	  int numr1=1,numm1=1,nump1=1,numt1=1;
	  int nums;
	  int numx1=1, numf1=1, numc1=1;
	  char *temp;
	  char dataline[10000];
	  double temp1, wavetemp, crosstemp, DD, GA, DenZ;
	  double z[zbin+1], T[zbin+1], PP[zbin+1], P[zbin+1];
      double **JJ, **cross, **qy, *wavep, *crossp, *qyp, *qyp1, *qyp2, *qyp3, *qyp4, *qyp5, *qyp6, *qyp7;
	  double **crosst, **qyt, *crosspt, *qypt, *qyp1t, *qyp2t, *qyp3t, *qyp4t, *qyp5t, *qyp6t, *qyp7t;
	  FILE *fspecies, *fzone, *fhenry, *fp, *fp1, *fp2, *fp3;
      FILE *fout, *fout1, *fout21, *fout22, *fout3, *fout4, *fcheck, *ftemp, *fout5, *foutp;

      GA=GRAVITY*MASS_PLANET/RADIUS_PLANET/RADIUS_PLANET; /* Planet Surface Gravity Acceleration, in SI */	
	
      /* Set the wavelength for computation of photolysis rate */
      double refidx0[WaveBin+1], DenS;
	  i=0;
	  iradmax = 0;
	  DenS=101325.0/KBOLTZMANN/273.0*1.0E-6; 
      while (i<=WaveBin)
      {
		  wavelength[i]=WaveMin+i*(WaveMax-WaveMin)/WaveBin;
		  if (wavelength[i]<=WaveMax1) { iradmax = i; }
		  if (RefIdxType == 0) { refidx0[i]=AirRefIdx(wavelength[i]);}
		  if (RefIdxType == 1) { refidx0[i]=CO2RefIdx(wavelength[i]);}
		  if (RefIdxType == 2) { refidx0[i]=HeRefIdx(wavelength[i]);}
		  if (RefIdxType == 3) { refidx0[i]=N2RefIdx(wavelength[i]);}
		  if (RefIdxType == 4) { refidx0[i]=NH3RefIdx(wavelength[i]);}
		  if (RefIdxType == 5) { refidx0[i]=CH4RefIdx(wavelength[i]);}
		  if (RefIdxType == 6) { refidx0[i]=H2RefIdx(wavelength[i]);}
		  if (RefIdxType == 7) { refidx0[i]=O2RefIdx(wavelength[i]);}
		  if (refidx0[i] < 1.0) { refidx0[i] = 1.0; }
		  /* crossr[i]=8.0*pow(PI,3)*pow(pow(refidx0[i],2)-1,2)/3.0/pow(wavelength[i]*1.0E-7,4)*1.061/DenS/DenS; */
		  crossr[i]=1.061*8.0*pow(PI,3)*pow(pow(refidx0[i],2)-1,2)/3.0/pow(wavelength[i]*1.0E-7,4)/DenS/DenS;
		  printf("%s\t%f\t%s\t%e\n", "The reyleigh scattering cross-section at wavelength", wavelength[i], "nm is", crossr[i]);
		  i = i+1;
      }
	printf("%s %d\n", "The maximum index for the UV-Visible Rad Transfer and Photolysis is", iradmax);
      
	/* Obtain the stellar radiation */
	fp2=fopen(STAR_SPEC,"r");
	fp3=fopen(STAR_SPEC,"r");
	s=LineNumber(fp2, 1000);
	double swave[s], sflux[s];
	GetData(fp3, 1000, s, swave, sflux);
	fclose(fp2);
	fclose(fp3);
	Interpolation(wavelength, WaveBin+1, solar, swave, sflux, s, 0);
	for (i=0; i<=WaveBin; i++) {
		solar[i] = solar[i]/ORBIT/ORBIT*FaintSun;  /* convert from flux at 1 AU */
		if (IFUVMULT == 1) {
			if (wavelength[i] < 200.0) {
				solar[i] = solar[i] * FUVMULT;
			} else if (wavelength[i] < 300.0) {
				solar[i] = solar[i] * MUVMULT;
			} else if (wavelength[i] < 400.0) {
				solar[i] = solar[i] * NUVMULT;
			}
		}
	}
	printf("%s\n", "The stellar radiation data are imported.");
	
	/* Set up the P-T-z for calculation */
	if (TPMODE==1) {
		DD=(zmax-zmin)/zbin;  /* Thickness of each layer, in unit of km */
		for (j=0; j<=zbin; j++) { z[j] = zmin+j*DD; } /* Set up the grid boundary */
		fp=fopen(TPLIST,"r");
		fp1=fopen(TPLIST,"r");
		s=LineNumber(fp, 1000);
		double Height[s];
		double Temp[s];
		double Pre[s];
		GetData3(fp1, 1000, s, Height, Pre, Temp);
		fclose(fp);
		fclose(fp1);
		Interpolation(z, zbin+1, T, Height, Temp, s, 0);
		Interpolation(z, zbin+1, PP, Height, Pre, s, 0);
		for (i=0; i<=zbin; i++) { P[i]=pow(10,PP[i]+PPOFFSET); } /*unit: Pa*/
    }
	if (TPMODE==0) {
		TPPara(P,T,TINV,zbin+1,PTOP,TTOP,TSTR,PSTR,PMIDDLE,TMIDDLE,PBOTTOM,TBOTTOM);
		DD=TPScale(P,T,zbin+1,z);
    }
	thickl = DD*1.0E+5; /* Thickness of each layer, in unit of cm*/
	for (j=1; j<=zbin; j++) {
		zl[j] = (z[j]+z[j-1])/2.0; /* Altitude at the center of layer */
		tl[j] = (T[j]+T[j-1])/2.0; /* Temperature at the center of layer */
		pl[j] = sqrt(P[j]*P[j-1]); /* Pressure at the center of layer */
	}
	FILE *TPPrint;
	TPPrint=fopen("Data/TPCheck.dat","w");
	for (j=0; j<=zbin; j++) {
		MMZ[j] = P[j]/KBOLTZMANN/T[j]*1.0E-6;
	}
	for (j=1; j<=zbin; j++) {
		MM[j]=pl[j]/KBOLTZMANN/tl[j]*1.0E-6; /*unit: Molecule cm-3*/
        printf("%lf %lf %lf %e\n", zl[j], pl[j], tl[j], MM[j]);
		fprintf(TPPrint, "%lf %lf %lf %e\n", zl[j], pl[j], tl[j], MM[j]);
	} 
	printf("%s\n", "The Z-T-P data are imported/calculated.");
	fclose(TPPrint);
	
	/* Determine the time-step limit due to the diffusion instability */
	double tslimit;
	tslimit = thickl*thickl/KET/4.0*TSPEED;
	printf("%s %e %s\n", "The diffusion-limit stepping time is ", tslimit, "s");
	
	/* Calculate the rainout rates throughout the atmosphere */
	double HenryH[NSP+1], HenryT[NSP+1], Heff;
	fhenry=fopen("Data/henry.dat","r");
	for (i=1; i<=NSP; i++) {
		fscanf(fhenry, "%lf\t", &HenryH[i]);
		fscanf(fhenry, "%lf\t", &HenryT[i]);
		printf("%e %e\n", HenryH[i], HenryT[i]);
	}
	fclose(fhenry);
	printf("%s\n","The Henry Law's constants are imported");
	for (i=1; i<=zbin; i++) {
		for (j=1; j<=NSP; j++) {
			Heff=HenryH[j]*exp(HenryT[j]*(1/tl[i]-1/298.15)); /* temperature-dependent Henry's law constant */ /*- sign in front of HenryT removed by JC, SR -- bugfix */
			rainoutrate[i][j]=2.0E-6/55.0/NAVOGADRO/(CloudDen*1.0E-9 + 1.0/(Heff*82.05746*tl[i]))*RainF;
			/* These rates need to be multiplied by the number density of H2O to yeild a true rainout rates */
			if (j==20 || j==21 || j==53 || j==9 || j==55 || j==27 || j==29 || j==31 || j==54) {
				rainoutrate[i][j] = 2.0E-6/55.0/NAVOGADRO/(CloudDen*1.0E-9 + 1.0/(Heff*82.05746*tl[i]))*RainF; /* turn off rainout for abiotic planets for H2, CO, CH4, NH3, N2, C2H2, C2H4, C2H6, and O2 */
			}
		}
	}
	
	/* Calculate the falling velocity for aerosol species */
	fallvelocity(GA);
	for (i=1; i<=zbin; i++) {
		printf("%s %f %s %f %s\n", "The Fall velocity at", zl[i],"km is", VFall[i], "cm s-1");		
	}
	
	/* Calculate the saturation density and condensation timescale of H2O, H2SO4 and S8 */
	double psH2O[zbin+1], psH2SO4[zbin+1], psS8[zbin+1];
	waterpressure(psH2O);
	sulfuridpressure(psH2SO4);
	sulfurpressure(psS8);
	for (i=1; i<=zbin; i++) {
		nsH2O[i] = psH2O[i]/KBOLTZMANN/tl[i]*1.0E-6*SATURATIONREDUCTION;
		printf("%s %f %s %e %s %e %s\n", "The saturation pressure and density of H2O at", zl[i], "km is", psH2O[i], "Pa and", nsH2O[i], "cm-3");
	}
	for (i=1; i<=zbin; i++) {
		nsH2SO4[i] = psH2SO4[i]/KBOLTZMANN/tl[i]*1.0E-6;
		printf("%s %f %s %e %s %e %s\n", "The saturation pressure and density of H2SO4 at", zl[i], "km is", psH2SO4[i], "Pa and", nsH2SO4[i], "cm-3");
	}
	for (i=1; i<=zbin; i++) {
		nsS8[i] = psS8[i]/KBOLTZMANN/tl[i]*1.0E-6;
		printf("%s %f %s %e %s %e %s\n", "The saturation pressure and density of S8 at", zl[i], "km is", psS8[i], "Pa and", nsS8[i], "cm-3");
	}
	for (i=1; i<=zbin; i++) {
		tcondfH2SO4[i]=49.0*AMU/AERDEN*pow(8.0*KBOLTZMANN*tl[i]/PI/98.0/AMU,0.5)/AERSIZE*1.0E+6;
		tcondfS8[i]=128.0*AMU/AERDEN*pow(8.0*KBOLTZMANN*tl[i]/PI/256.0/AMU,0.5)/AERSIZE*1.0E+6;
		tcondfH2O[i]=9.0*AMU/1.0E+3*pow(8.0*KBOLTZMANN*tl[i]/PI/18.0/AMU,0.5)/AERSIZE*1.0E+6;
		printf("%s %f %s %e %s\n", "The H2SO4 condensation timescale factor at ", zl[i], "km is", tcondfH2SO4[i], "cm3 s-1");
		printf("%s %f %s %e %s\n", "The S8 condensation timescale factor at ", zl[i], "km is", tcondfS8[i], "cm3 s-1");
		printf("%s %f %s %e %s\n", "The H2O condensation timescale factor at ", zl[i], "km is", tcondfH2O[i], "cm3 s-1");
	}
	
	/* Import the molecular opacity */
	opacCO2 = dmatrix(1,zbin,0,WaveBin);
	opacO2 = dmatrix(1,zbin,0,WaveBin);
	opacSO2 = dmatrix(1,zbin,0,WaveBin);
	opacH2O = dmatrix(1,zbin,0,WaveBin);
	opacOH = dmatrix(1,zbin,0,WaveBin);
	opacH2CO = dmatrix(1,zbin,0,WaveBin);
	opacH2O2 = dmatrix(1,zbin,0,WaveBin);
	opacHO2 = dmatrix(1,zbin,0,WaveBin);
	opacH2S = dmatrix(1,zbin,0,WaveBin);
	opacCO = dmatrix(1,zbin,0,WaveBin);
	opacO3 = dmatrix(1,zbin,0,WaveBin);
	opacCH4 = dmatrix(1,zbin,0,WaveBin);
	opacNH3 = dmatrix(1,zbin,0,WaveBin);
	
	/* Compute the planck mean of each molecular at each height */
	planckmean(MeanCO2, "Cross/opacCO2.dat", opacCO2);
	printf("CO2 mean opacity in the infrared calculated!\n");	
	planckmean(MeanO2, "Cross/opacO2.dat", opacO2);
	printf("O2 mean opacity in the infrared calculated!\n");	
	planckmean(MeanSO2, "Cross/opacSO2.dat", opacSO2);
	printf("SO2 mean opacity in the infrared calculated!\n");	
	planckmean(MeanH2O, "Cross/opacH2O.dat", opacH2O);
	printf("H2O mean opacity in the infrared calculated!\n");	
	planckmean(MeanOH, "Cross/opacOH.dat", opacOH);
	printf("OH mean opacity in the infrared calculated!\n");
	planckmean(MeanH2CO, "Cross/opacH2CO.dat", opacH2CO);
	printf("H2CO mean opacity in the infrared calculated!\n");
	planckmean(MeanH2O2, "Cross/opacH2O2.dat", opacH2O2);
	printf("H2O2 mean opacity in the infrared calculated!\n");
	planckmean(MeanHO2, "Cross/opacHO2.dat", opacHO2);
	printf("HO2 mean opacity in the infrared calculated!\n");
	planckmean(MeanH2S, "Cross/opacH2S.dat", opacH2S);
	printf("H2S mean opacity in the infrared calculated!\n");
	planckmean(MeanCO, "Cross/opacCO.dat", opacCO);
	printf("CO mean opacity in the infrared calculated!\n");
	planckmean(MeanO3, "Cross/opacO3.dat", opacO3);
	printf("O3 mean opacity in the infrared calculated!\n");
	planckmean(MeanCH4, "Cross/opacCH4.dat", opacCH4);
	printf("CH4 mean opacity in the infrared calculated!\n");
	planckmean(MeanNH3, "Cross/opacNH3.dat", opacNH3);
	printf("NH3 mean opacity in the infrared calculated!\n");
	
	/* Get the species list */
	fspecies=fopen(SPECIES_LIST, "r");
	fout21=fopen(OUT_FILE1,"w");
	fout22=fopen(OUT_FILE2,"w");
    fprintf(fout21, "%s\t\t\t", "z");
    fprintf(fout22, "%s\t\t\t", "z");
	s=LineNumber(fspecies, 10000);
	printf("Species list: \n");
	fclose(fspecies);
	fspecies=fopen(SPECIES_LIST, "r");
	struct Molecule species[s];
	temp=fgets(dataline, 10000, fspecies); /* Read in the header line */
	i=0;
	while (fgets(dataline, 10000, fspecies) != NULL )
	{
		sscanf(dataline, "%s %s %d %d %lf %lf %d %lf %lf", (species+i)->name, (species+i)->type, &((species+i)->num), &((species+i)->mass), &((species+i)->mix), &((species+i)->upper), &((species+i)->lowertype), &((species+i)->lower), &((species+i)->lower1));
		printf("%s %s %d %d %lf %lf %d %lf %lf\n",(species+i)->name, (species+i)->type, (species+i)->num, (species+i)->mass, (species+i)->mix, (species+i)->upper, (species+i)->lowertype, (species+i)->lower, (species+i)->lower1);
        if (strcmp("X",species[i].type)==0) {numx=numx+1; fprintf(fout21, "%s\t\t\t", (species+i)->name);}
		if (strcmp("F",species[i].type)==0) {numf=numf+1; fprintf(fout22, "%s\t\t\t", (species+i)->name);}
		if (strcmp("C",species[i].type)==0) {numc=numc+1;}
		if (strcmp("A",species[i].type)==0) {numx=numx+1; fprintf(fout21, "%s\t\t\t", (species+i)->name); numa=numa+1;}
		i=i+1;
    }
    fclose(fspecies);
	fprintf(fout21, "%s\n\n", "Air");
	fprintf(fout22, "%s\n\n", "Air");
    fclose(fout21);
    fclose(fout22);
	nums=numx+numf+numc;
	nn=zbin*numx;
	double Con[nn+1], fvec[nn+1], ConC[zbin*numc+1], Conf[zbin*numf+1]; 
	printf("%s\n", "The species list is imported.");
	printf("%s %d\n", "Number of species in model:", nums);
	printf("%s %d\n", "Number of species to be solved in full:", numx);
	printf("%s %d\n", "In which the number of aerosol species is:", numa);
	printf("%s %d\n", "Number of species to be solved in photochemical equil:", numf);
	printf("%s %d\n", "Number of species assumed to be constant:", numc);
	
	
	/* Variable initialization */
	int labelx[numx+1], labelc[numc+1], labelf[numf+1], MoleculeM[numx+1], listFix[numx+1], listAER[numa+1], AERCount=1; /* Standard number list of species */
	double Upflux[numx+1], Loflux[numx+1], Depo[numx+1], ConFix[numx+1], mixtemp;
	FILE *fimport;
	FILE *fimportcheck;
	for (i=0; i<s; i++) {
		if (strcmp("X",species[i].type)==0 || strcmp("A",species[i].type)==0) {
			if (IMODE==1) {
				for (j=1; j<=zbin; j++) { Con[(j-1)*numx+numx1]=MM[j]*species[i].mix;} /* Initialized the variables */
			}
			labelx[numx1]=species[i].num;
			if (species[i].num==7) {
				waternum=numx1;
				waterx=1;
			}
			MoleculeM[numx1]=species[i].mass;
			Upflux[numx1]=species[i].upper/thickl;
			Depo[numx1]=species[i].lower/thickl;
			Loflux[numx1]=species[i].lower1/thickl;
			if (species[i].lowertype==1) {
				ConFix[numx1]=species[i].lower1*MM[1];
				Con[numx1]=ConFix[numx1];
				listFix[numx1]=1;
			} else {
				listFix[numx1]=0;
			}
			if (strcmp("A",species[i].type)==0) {
				listAER[AERCount]=numx1;
				AERCount = AERCount+1;
				printf("%s %d\n", "The aerosol species is", numx1);
			}
			numx1=numx1+1;
		}
		if (strcmp("F",species[i].type)==0) {
			labelf[numf1]=species[i].num;
			if (IMODE==1) {
				for (j=1; j<=zbin; j++) { Conf[(j-1)*numf+numf1]=MM[j]*species[i].mix;} /* Initialized the variables */
			}
			numf1=numf1+1;
		}
		if (strcmp("C",species[i].type)==0) {
			labelc[numc1]=species[i].num;
			for (j=1; j<=zbin; j++) {
				ConC[(j-1)*numc+numc1]=MM[j]*species[i].mix; /* Initialize the constants */
			}
			/* import constant mixing ratio list for H2O */
			if (IFIMPORTH2O == 1 && species[i].num == 7) {
				fimport=fopen("Data/ConstantMixing.dat", "r");
				fimportcheck=fopen("Data/ConstantMixingH2O.dat", "w");
				temp=fgets(dataline, 10000, fimport); /* Read in the header line */
				for (j=1; j<=zbin; j++) {
					fscanf(fimport, "%lf\t", &temp1);
					fscanf(fimport, "%le\t", &mixtemp);
					fscanf(fimport, "%le\t", &temp1);
					ConC[(j-1)*numc+numc1]=mixtemp * MM[j];
					fprintf(fimportcheck, "%f\t%e\t%e\n", zl[j], mixtemp, ConC[(j-1)*numc+numc1]);
				}
				fclose(fimport);
				fclose(fimportcheck);
			}
			/* import constant mixing ratio list for H2O */
			if (IFIMPORTCO2 == 1 && species[i].num == 52) {
				fimport=fopen("Data/ConstantMixing.dat", "r");
				fimportcheck=fopen("Data/ConstantMixingCO2.dat", "w");
				temp=fgets(dataline, 10000, fimport); /* Read in the header line */
				for (j=1; j<=zbin; j++) {
					fscanf(fimport, "%lf\t", &temp1);
					fscanf(fimport, "%le\t", &temp1);
					fscanf(fimport, "%le\t", &mixtemp);
					ConC[(j-1)*numc+numc1]=mixtemp * MM[j];
					fprintf(fimportcheck, "%f\t%e\t%e\n", zl[j], mixtemp, ConC[(j-1)*numc+numc1]);
				}
				fclose(fimport);
				fclose(fimportcheck);
			}
			numc1=numc1+1;
		}
    }
	if (IMODE==0) { /* Calculate the initial mixing ratio from chemical equilibrium */
		double **mixequil;
		int labels[numx+numf+1];
		for (i=1; i<=numx; i++) {labels[i]=labelx[i];}
		for (i=1; i<=numf; i++) {labels[numx+i]=labelf[i];}
		mixequil=dmatrix(1,zbin,1,numx+numf);
		chemquil(pl, tl, zbin+1, labels, numx+numf, mixequil); 
		for (j=1; j<=zbin; j++) {
            for (i=1; i<=numx; i++) {Con[(j-1)*numx+i]=MM[j]*mixequil[j][i]; printf("X %d %d %e\n", j, i, Con[(j-1)*numx+i]);}
            for (i=1; i<=numf; i++) {Conf[(j-1)*numf+i]=MM[j]*mixequil[j][i+numx]; printf("F %d %d %e\n", j, i, Conf[(j-1)*numf+i]);}
        }
        free_dmatrix(mixequil,1,zbin,1,numx+numf);
    }
    if (IMODE==3) {
        double Contemp[4];
        for (j=1; j<=zbin; j++) {
			for (i=1; i<=numx; i++) {
				Con[(j-1)*numx+i]=0.0;
			}
			for (i=1; i<=numf; i++) {
				Conf[(j-1)*numf+i]=0.0;
			}
            SimChemCarbon(Contemp, tl[j]);
            for (i=1; i<=numx; i++) {
            if (labelx[i]==20) {Con[(j-1)*numx+i]=Contemp[1];} /* CO  */
            if (labelx[i]==21) {Con[(j-1)*numx+i]=Contemp[2];} /* CH4 */
            if (labelx[i]==7)  {Con[(j-1)*numx+i]=Contemp[3];} /* H2O */ 
            }
        }
    }
    if (IMODE==2) {
        fimport=fopen(IMPORTFILEX, "r");
		fimportcheck=fopen("Data/Fimportcheck.dat","w");
        temp=fgets(dataline, 10000, fimport); /* Read in the header line */
        temp=fgets(dataline, 10000, fimport); /* Read in the header line */
        for (j=1; j<=zbin; j++) {
            fscanf(fimport, "%lf\t\t", &temp1);
			fprintf(fimportcheck, "%lf\t", temp1);
            for (i=1; i<=numx; i++) {
                fscanf(fimport, "%le\t\t", &Con[(j-1)*numx+i]);
                fprintf(fimportcheck, "%e\t", Con[(j-1)*numx+i]);
            }
			fscanf(fimport, "%lf\t\t", &temp1); /* column of air */
			fprintf(fimportcheck,"\n");
        }
        fclose(fimport);
        fimport=fopen(IMPORTFILEF, "r");
        temp=fgets(dataline, 10000, fimport); /* Read in the header line */
        temp=fgets(dataline, 10000, fimport); /* Read in the header line */
        for (j=1; j<=zbin; j++) {
            fscanf(fimport, "%lf\t\t", &temp1);
			fprintf(fimportcheck, "%lf\t", temp1);
            for (i=1; i<=numf; i++) {
                fscanf(fimport, "%le\t\t", &Conf[(j-1)*numf+i]);
                fprintf(fimportcheck, "%e\t", Conf[(j-1)*numf+i]);
            }
            fscanf(fimport, "%lf\t\t", &temp1); /* column of air */
			fprintf(fimportcheck,"\n");
        }
        fclose(fimport);
		fclose(fimportcheck);
    }
	/* Fix type-1 boundary */
	for (i=1; i<=numx; i++) {
		if (listFix[i]==1) {
				Con[i]=ConFix[i];
		}
	}
	/*Generate General Variables*/
	Convert1(Con, ConC, Conf, labelx, labelc, labelf);
	if (IMODE == 4) {
		fimport=fopen(OUT_STD, "r");
		fimportcheck=fopen("Data/Fimportcheck.dat","w");
		temp=fgets(dataline, 10000, fimport); /* Read in the header line */
        temp=fgets(dataline, 10000, fimport); /* Read in the header line */
		for (j=1; j<=zbin; j++) {
            fscanf(fimport, "%lf\t", &temp1);
			fprintf(fimportcheck, "%lf\t", temp1);
			fscanf(fimport, "%lf\t", &temp1);
			fscanf(fimport, "%lf\t", &temp1);
			fscanf(fimport, "%lf\t", &temp1);
			fscanf(fimport, "%le\t", &temp1);
            for (i=1; i<=NSP; i++) {
                fscanf(fimport, "%le\t", &xx[j][i]);
                fprintf(fimportcheck, "%e\t", xx[j][i]);
            }
			fprintf(fimportcheck,"\n");
        }
		fclose(fimport);
		fclose(fimportcheck);
		Convert2(Con, ConC, Conf, labelx, labelc, labelf);
	}
	printf("%s\n", "Variable initialization successful");
	
	/* Calculate the empirical Eddy and Molecular Diffusion Coefficient, in unit of cm2 s-1*/
	DM=dmatrix(1, zbin, 1, numx);
	dl=dmatrix(1, zbin, 1, numx);
	if (EDDYPARA == 1) {
		for (i=1; i<zbin; i++) {
			KE[i]=KET;
			if (z[i]<ZT) {
				DenZ=MMZ[i];
			}
			if (z[i]>ZT) {
				KE[i] *= pow(DenZ/MMZ[i],0.5);
				if (KE[i]>KEH) {KE[i]=KEH;}
			}
			KE[i]=KE[i]/thickl/thickl*MMZ[i]; 
		}
	}
	if (EDDYPARA == 2) {
		fp=fopen(EDDYIMPORT,"r");
		fp1=fopen(EDDYIMPORT,"r");
		s=LineNumber(fp, 1000);
		double Ealt[s];
		double Etem[s];
		GetData(fp1, 1000, s, Ealt, Etem);
		fclose(fp);
		fclose(fp1);
		Interpolation(z, zbin+1, KE, Ealt, Etem, s, 0);
		for (i=1; i<zbin; i++) {
			KE[i]=KE[i]/thickl/thickl*MMZ[i];
		}
	}
	for (i=1; i<zbin; i++) {
		printf("%s %f %e\n", "Eddy diffusion at altitude", z[i], KE[i]/MMZ[i]);
	}
	for (i=1; i<=zbin; i++) {
		for (j=1; j<=numx; j++) {
			DM[i][j]=0; /*Calculate the molecular diffusion coefficient, unit s-1*/
			if (labelx[j]==3) { /* H */
				DM[i][j]=MDIFF_H_1*1.0E+17*pow(T[i],MDIFF_H_2)/thickl/thickl;
			} /* H */
			if (labelx[j]==53) { /* H2 */
				DM[i][j]=MDIFF_H2_1*1.0E+17*pow(T[i],MDIFF_H2_2)/thickl/thickl;
			}
			dl[i][j]=DM[i][j]/2.0*((AIRM-MoleculeM[j])*AMU*GA*thickl*1.0E-2/KBOLTZMANN/T[i]+0.38/T[i]*(tl[i+1]-tl[i]));
		}
	}
	
	/* Diffusion-limited escape velocity */
	for (j=1; j<=NSP; j++) {
		Vesc[j]=0;
	}
	for (j=1; j<=numx; j++) {
		Vesc[labelx[j]] = DM[zbin][j]/MMZ[zbin]*thickl*(AIRM-MoleculeM[j])*AMU*GA*1.0E-2/KBOLTZMANN/T[zbin];
	}
	Vesc[53] = Vesc[53]*MDIFF_H2_F; /* correct for H2 photo-dissociation above the considered atmosphere */
	printf("%s %e\n", "The escape frequency of H is", Vesc[3]);
	printf("%s %e\n", "The escape frequency of H2 is", Vesc[53]);
	
	
      
	/* Get the reaction list */
	fzone=fopen(REACTION_LIST, "r");
	s=LineNumber(fzone, 10000);
	fclose(fzone);
	fzone=fopen(REACTION_LIST, "r");
	struct Reaction React[s];
	temp=fgets(dataline, 10000, fzone); /* Read in the header line */
	i=0;
	while (fgets(dataline, 10000, fzone) != NULL )
	{
		sscanf(dataline, "%d %s %d", &((React+i)->dum), (React+i)->type, &((React+i)->num));
		printf("%d %s %d\n", (React+i)->dum, React[i].type, React[i].num);
		if (strcmp("R",React[i].type)==0) {numr=numr+1;}
		if (strcmp("M",React[i].type)==0) {numm=numm+1;}
		if (strcmp("P",React[i].type)==0) {nump=nump+1;}
		if (strcmp("T",React[i].type)==0) {numt=numt+1;}
		i=i+1;
    }
	fclose(fzone);
	int zone_r[numr+1], zone_m[numm+1], zone_p[nump+1], zone_t[numt+1];
	for (i=0; i<s; i++) {
		if (strcmp("R",React[i].type)==0) {
			zone_r[numr1]=React[i].num;
			numr1=numr1+1;
		}
		;if (strcmp("M",React[i].type)==0) {
			zone_m[numm1]=React[i].num;
			numm1=numm1+1;
		}
		if (strcmp("P",React[i].type)==0) {
			zone_p[nump1]=React[i].num;
			nump1=nump1+1;
		}
		if (strcmp("T",React[i].type)==0) {
			zone_t[numt1]=React[i].num;
			numt1=numt1+1;
		}
	}
	printf("%s\n", "The reaction lists are imported.");
	printf("%s %d\n", "Number of bi-molecular reactions:", numr);
	printf("%s %d\n", "Number of tri-molecular reactions:", numm);
	printf("%s %d\n", "Number of photolysis:", nump);
	printf("%s %d\n", "Number of thermo-dissociations:", numt);
	
	  /* Load the standard reaction list and information */
	  GetReaction();
	  printf("%s\n", "The standard reaction databases are imported.");
	  ReactionRate(); /* Calculate Reaction Rate*/
	  ReactionRateM(); /* Calculate Reaction Rate*/
	  ReactionRateT(); /* Calculate Reaction Rate*/
      printf("%s\n", "Reaction rates initialization successful.");
	/* Print-out the reaction rates used in the model */
	FILE *frates;
	frates=fopen("Data/ReactonRateCheck.dat", "w");
	fprintf(frates, "%s\t\t%s\t\t%s\t\t%s\n", "Type", "STD", "Rate at Bottom", "Rate at Top");
	for (i=1; i<=numr; i++) {
		fprintf(frates, "%s\t\t%d\t\t%e\t\t%e\n", "R", zone_r[i], kk[1][zone_r[i]], kk[zbin][zone_r[i]]);}
	for (i=1; i<=numm; i++) {
		fprintf(frates, "%s\t\t%d\t\t%e\t\t%e\n", "M", zone_m[i], kkM[1][zone_m[i]], kkM[zbin][zone_m[i]]);}
	for (i=1; i<=numt; i++) {
		fprintf(frates, "%s\t\t%d\t\t%e\t\t%e\n", "T", zone_t[i], kkT[1][zone_t[i]], kkT[zbin][zone_t[i]]);}
	fclose(frates);
	
		/* get the cross sections and quantum yields of molecules */   
        cross=dmatrix(1,nump,0,WaveBin);
		crosst=dmatrix(1,nump,0,WaveBin);
        qy=dmatrix(1,nump,0,WaveBin);
		qyt=dmatrix(1,nump,0,WaveBin);
	    int stdcross[nump+1];
	    double qysum[nump+1];
        fcheck=fopen("Data/CrossSectionCheck.dat","w"); 
	for (i=1; i<=nump; i++) {
		stdcross[i]=ReactionP[zone_p[i]][1];
		qytype=ReactionP[zone_p[i]][8];
		qysum[i]=ReactionP[zone_p[i]][7];
		j=0;
		while (species[j].num != stdcross[i]) {j=j+1;}
		/* printf("%s\n",species[j].name); */
		fp=fopen(species[j].name, "r");
		fp1=fopen(species[j].name, "r");
		s=LineNumber(fp, 1000);
		/* printf("%d\n",s); */
		wavep=dvector(0,s-1);
		crossp=dvector(0,s-1);
		qyp=dvector(0,s-1);
		qyp1=dvector(0,s-1);
		qyp2=dvector(0,s-1);
		qyp3=dvector(0,s-1);
		qyp4=dvector(0,s-1);
		qyp5=dvector(0,s-1);
		qyp6=dvector(0,s-1);
		qyp7=dvector(0,s-1);
		crosspt=dvector(0,s-1);
		qypt=dvector(0,s-1);
		qyp1t=dvector(0,s-1);
		qyp2t=dvector(0,s-1);
		qyp3t=dvector(0,s-1);
		qyp4t=dvector(0,s-1);
		qyp5t=dvector(0,s-1);
		qyp6t=dvector(0,s-1);
		qyp7t=dvector(0,s-1);
        k=0;
        if (qytype==1) {
        while (fgets(dataline, 1000, fp1) != NULL ) {
            sscanf(dataline, "%lf %le %le %lf %lf", wavep+k, crossp+k, crosspt+k, qyp+k, qypt+k);
            k=k+1; }
        }
        if (qytype==2) {
        while (fgets(dataline, 1000, fp1) != NULL ) {
            sscanf(dataline, "%lf %le %le %lf %lf %lf %lf", wavep+k, crossp+k, crosspt+k, qyp1+k, qyp1t+k, qyp+k, qypt+k);
            k=k+1; }
        }
        if (qytype==3) {
        while (fgets(dataline, 1000, fp1) != NULL ) {
            sscanf(dataline, "%lf %le %le %lf %lf %lf %lf %lf %lf", wavep+k, crossp+k, crosspt+k, qyp1+k, qyp1t+k, qyp2+k, qyp2t+k, qyp+k, qypt+k);
            k=k+1; }
        }
		if (qytype==4) {
		while (fgets(dataline, 1000, fp1) != NULL ) {
			sscanf(dataline, "%lf %le %le %lf %lf %lf %lf %lf %lf %lf %lf", wavep+k, crossp+k, crosspt+k, qyp1+k, qyp1t+k, qyp2+k, qyp2t+k, qyp3+k, qyp3t+k, qyp+k, qypt+k);
			k=k+1; }
        }
		if (qytype==5) {
		while (fgets(dataline, 1000, fp1) != NULL ) {
			sscanf(dataline, "%lf %le %le %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", wavep+k, crossp+k, crosspt+k, qyp1+k, qyp1t+k, qyp2+k, qyp2t+k, qyp3+k, qyp3t+k, qyp4+k, qyp4t+k, qyp+k, qypt+k);
			k=k+1; }
        }
		if (qytype==6) {
		while (fgets(dataline, 1000, fp1) != NULL ) {
			sscanf(dataline, "%lf %le %le %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", wavep+k, crossp+k, crosspt+k, qyp1+k, qyp1t+k, qyp2+k, qyp2t+k, qyp3+k, qyp3t+k, qyp4+k, qyp4t+k, qyp5+k, qyp5t+k, qyp+k, qypt+k);
			k=k+1; }
        }
		if (qytype==7) {
		while (fgets(dataline, 1000, fp1) != NULL ) {
			sscanf(dataline, "%lf %le %le %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", wavep+k, crossp+k, crosspt+k, qyp1+k, qyp1t+k, qyp2+k, qyp2t+k, qyp3+k, qyp3t+k, qyp4+k, qyp4t+k, qyp5+k, qyp5t+k, qyp6+k, qyp6t+k, qyp+k, qypt+k);
			k=k+1; }
        }
		if (qytype==8) {
		while (fgets(dataline, 1000, fp1) != NULL ) {
			sscanf(dataline, "%lf %le %le %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", wavep+k, crossp+k, crosspt+k, qyp1+k, qyp1t+k, qyp2+k, qyp2t+k, qyp3+k, qyp3t+k, qyp4+k, qyp4t+k, qyp5+k, qyp5t+k, qyp6+k, qyp6t+k, qyp7+k, qyp7t+k, qyp+k, qypt+k);
			k=k+1; }
        }
		fclose(fp);
		fclose(fp1);
		Interpolation(wavelength, WaveBin+1, *(cross+i), wavep, crossp, s, 0);
		Interpolation(wavelength, WaveBin+1, *(qy+i), wavep, qyp, s, 0);
		Interpolation(wavelength, WaveBin+1, *(crosst+i), wavep, crosspt, s, 0);
		Interpolation(wavelength, WaveBin+1, *(qyt+i), wavep, qypt, s, 0);
		free_dvector(wavep,0,s-1);
		free_dvector(crossp,0,s-1);
		free_dvector(qyp,0,s-1);
		free_dvector(qyp1,0,s-1);
        free_dvector(qyp2,0,s-1);
		free_dvector(qyp3,0,s-1);
		free_dvector(qyp4,0,s-1);
		free_dvector(qyp5,0,s-1);
		free_dvector(qyp6,0,s-1);
		free_dvector(qyp7,0,s-1);
		free_dvector(crosspt,0,s-1);
		free_dvector(qypt,0,s-1);
		free_dvector(qyp1t,0,s-1);
        free_dvector(qyp2t,0,s-1);
		free_dvector(qyp3t,0,s-1);
		free_dvector(qyp4t,0,s-1);
		free_dvector(qyp5t,0,s-1);
		free_dvector(qyp6t,0,s-1);
		free_dvector(qyp7t,0,s-1);
		printf("%s %s %s\n", "The", species[j].name, "Cross section and quantum yield data are imported.");
		fprintf(fcheck, "%s %s %s\n", "The", species[j].name, "Cross section and quantum yield data are imported.");
		for (j=0;j<=WaveBin;j++) {fprintf(fcheck, "%lf %le %le %lf %lf\n", wavelength[j], cross[i][j], crosst[i][j], qy[i][j], qyt[i][j]);}
	}
	
	/* cross section of aerosols */
	double *crossp1, *crossp2, *crossp3;
	double crossw1[WaveBin+1], crossw2[WaveBin+1], crossw3[WaveBin+1];
	fp=fopen(AERRADFILE1,"r");
	fp1=fopen(AERRADFILE1,"r");
	s=LineNumber(fp, 1000);
	wavep=dvector(0,s-1);
	crossp1=dvector(0,s-1);
	crossp2=dvector(0,s-1);
	crossp3=dvector(0,s-1);
	k=0;
	while (fgets(dataline, 1000, fp1) != NULL ) {
		sscanf(dataline, "%lf %lf %lf %lf", wavep+k, crossp1+k, crossp2+k, crossp3+k);
		k=k+1; 
	}
	fclose(fp);
	fclose(fp1);
	Interpolation(wavelength, WaveBin+1, crossw1, wavep, crossp1, s, 0);
	Interpolation(wavelength, WaveBin+1, crossw2, wavep, crossp2, s, 0);
	Interpolation(wavelength, WaveBin+1, crossw3, wavep, crossp3, s, 0);
	free_dvector(wavep,0,s-1);
	free_dvector(crossp1,0,s-1);
	free_dvector(crossp2,0,s-1);
	free_dvector(crossp3,0,s-1);
	for (i=0; i<=WaveBin; i++) {
		crossa[1][i] = crossw1[i];
		sinab[1][i]  = crossw2[i]/(crossw1[i]+1.0e-24);
		asym[1][i]   = crossw3[i];
	}
	fp=fopen(AERRADFILE2,"r");
	fp1=fopen(AERRADFILE2,"r");
	s=LineNumber(fp, 1000);
	wavep=dvector(0,s-1);
	crossp1=dvector(0,s-1);
	crossp2=dvector(0,s-1);
	crossp3=dvector(0,s-1);
	k=0;
	while (fgets(dataline, 1000, fp1) != NULL ) {
		sscanf(dataline, "%lf %lf %lf %lf", wavep+k, crossp1+k, crossp2+k, crossp3+k);
		k=k+1; 
	}
	fclose(fp);
	fclose(fp1);
	Interpolation(wavelength, WaveBin+1, crossw1, wavep, crossp1, s, 0);
	Interpolation(wavelength, WaveBin+1, crossw2, wavep, crossp2, s, 0);
	Interpolation(wavelength, WaveBin+1, crossw3, wavep, crossp3, s, 0);
	free_dvector(wavep,0,s-1);
	free_dvector(crossp1,0,s-1);
	free_dvector(crossp2,0,s-1);
	free_dvector(crossp3,0,s-1);
	for (i=0; i<=WaveBin; i++) {
		crossa[2][i] = crossw1[i];
		sinab[2][i]  = crossw2[i]/(crossw1[i]+1.0e-24);
		asym[2][i]   = crossw3[i];
	}
	printf("%s\n", "Cross sections of the aerosol are imported.");
	fprintf(fcheck, "%s\n", "Cross sections of the aerosol are imported.");
	for (j=0;j<NLAMBDA;j++) {fprintf(fcheck, "%lf %e %e %f %f %f %f\n", wavelength[j], crossa[1][j], crossa[2][j], sinab[1][j], sinab[2][j], asym[1][j], asym[2][j]);}
	fclose(fcheck);
	  
	/* Generate the first radiation field */
	double **rad, **opt;
	rad=dmatrix(0, iradmax, 0, zbin);
	opt=dmatrix(0, iradmax, 0, zbin);
	RadTransfer(rad, opt, stdcross, qysum, cross, crosst, iradmax+1);
	fout=fopen("Data/Radiation0.dat", "w");
	fout1=fopen("Data/OpticalDepth0.dat","w");
	for(ii=0; ii<=zbin; ii++)
	{
		fprintf(fout, "%s %f\n", "The initial radiation at z=", z[ii]);
		fprintf(fout1, "%s %f\n", "The optical depth at z=", z[ii]);
		for (jj=0; jj<=iradmax; jj++) {
		  fprintf(fout, "%f %e\n", wavelength[jj], rad[jj][ii]); /*write the initial radiation to a file*/
		  fprintf(fout1, "%f %f\n", wavelength[jj], opt[jj][ii]);} /*write the initial radiation to a file*/
	}
	fclose(fout);
	fclose(fout1);
	printf("%s\n", "The first radiation field is generated.");
	MeanOpacity(stdcross,qysum,cross,crosst);
	for (j=1; j<=zbin; j++) {
		printf("%s %f %s %e %s\n", "The mean opacity at altitude", zl[j], "km is", mkv[j], "cm-1");
	}
	JJ=dmatrix(1, zbin, 1, nump);
	Photorate(rad, cross, crosst, qy, qyt, zone_p, nump, JJ);
	printf("%s\n", "The photolysis rate is calculated for the first time.");
	ftemp=fopen(OUT_PHOTORATE, "w");
	for (ii=1; ii<=nump; ii++) {
        fprintf(ftemp, "%s\t%d\t%s\t%e\n", "The photolysis STD number", zone_p[ii], "rate at the top layer is", JJ[zbin][ii]);}
	fclose(ftemp);
	fclose(foutp);
	GreyTemp(rad, iradmax);
	foutp=fopen(OUT_NEWTEMP,"w");
	for (jj=0; jj<=zbin; jj++) {
		fprintf(foutp, "%s %f %s %e %s\n", "The new temperature at altitude", z[jj], "km is", Tnew[jj], "K");
	}
	fclose(foutp);
	
	/* Calculate the concentration of different molecules*/
	
	double tstep;
	tstep = fmin(TSINI, tslimit);
	double tstepold;  /* Initial time step */
	int check=1;
	double tt, control, controlt, controlt_old, ddd, test, emax;
	int controls, alindex;
	double controlz;
	double **Jaco, **Jaco1;
	Jaco=dmatrix(1, nn, 1, nn);
	Jaco1=dmatrix(1, nn, 1, nn);
	int erradjust, errn, errn1, errn2, errn3;
	erradjust=0; /* set to numx if Do not consider the error in the top and bottom layer */
	if (FINE1!=1) { erradjust=numx; }
	errn=numx+1; /* choose a species that you don't want to be considered in error. set to >numx if want to consider all */
	errn1=numx+1;
	errn2=numx+1;
	errn3=numx+1;
	
	fout4=fopen(OUT_CONVERGENCE,"w");
	fout5=fopen(OUT_HISTORY, "w");

	for (j=0; j<NMAX; j++) {
        
        Photorate(rad, cross, crosst, qy, qyt, zone_p, nump, JJ);
		
	    ChemEqu(Con, ConC, Conf, fvec, Jaco, labelx, labelc, labelf, listAER,
					 zone_r, zone_m, zone_t, JJ, zone_p, Upflux, Loflux, Depo, listFix);
		
		emax=1;
		tstep=2*tstep;
		
		while (emax>0.3) {
			
			tstep=0.5*tstep; /* Decrease the timestep if error larger than 0.3 */
		
		for (ii=1; ii<=nn; ii++) {
			for (jj=1; jj<=nn; jj++) {
				Jaco1[ii][jj] = -tstep*Jaco[ii][jj];
			}
			Jaco1[ii][ii] += 1;
		}
		
		BTridiagonal(Jaco1, fvec, zbin, numx);
		
		/* Calculate the time step control parameter */
		control=0;
		controlz=zmin;
		controls=1;
		controlt_old=controlt;
		controlt=0;
		for (ii=1+erradjust; ii<=nn-erradjust; ii++) {
			if ( (ii<=numx & listFix[ii]!=1) || ii>numx ) {
			test=fabs(fvec[ii]/fmax(MINNUM, Con[ii]));
			if (test>control & fmod(ii,numx)!=errn & fmod(ii,numx)!=errn1 & fmod(ii,numx)!=errn2 & fmod(ii,numx)!=errn3) {
				controlt=control;
				control=test;
				/* if (fvec[ii]<0) {
					controlt=test;
				} */
				alindex=trunc((ii-1)/numx)+1;
				controlz=zl[alindex]; /* Get the fastest decreasing specie */
				controls=fmod(ii, numx);
				if (controls==0) {controls=numx;}
			}
			if (test<control & test>controlt & fmod(ii,numx)!=errn & fmod(ii,numx)!=errn1 & fmod(ii,numx)!=errn2 & fmod(ii,numx)!=errn3) {
				controlt=test;
			}
			}
		}
		if (control/controlt>100 & FINE2!=1) {
			control=controlt;
			printf("%s\n", "disregard the fastest varying point");
		}
		emax=control*tstep;
			
		}
		
		fprintf(fout4, "%f %e %e %f %d\n", tt, control, controlt, controlz, labelx[controls]);
		
		if (fmod(j, NPRINT)==0) { printout_timescale(Con, fvec);}
		
		if (j>100 & control<Tol2) {
			check=0;
			Convert1(Con, ConC, Conf, labelx, labelc, labelf);	
			printf("%s\n", "converged!");
			/* record the result at each step */
			if (fmod(j, NPRINT)==0 && HISTORYPRINT == 1) {
				fprintf(fout5,"%f\t",tt);
				for (ii=1; ii<=zbin; ii++) {
					for (jj=1; jj<=NSP; jj++) {
						fprintf(fout5, "%e\t", xx[ii][jj]);
					}
				}
				for (ii=1; ii<=zbin; ii++) {
					fprintf(fout5, "%f\t", tl[ii]);
				}
				for (ii=0; ii<=zbin; ii++) {
					fprintf(fout5, "%f\t", z[ii]);
				}
				fprintf(fout5,"\n");;
			}
			break; /* Converged! */
		}
		
		if (tt > NMAXT) {
			Convert1(Con, ConC, Conf, labelx, labelc, labelf);	
			printf("%s\n", "completed");
			/* record the result at each step */
			if (fmod(j, NPRINT)==0 && HISTORYPRINT == 1) {
				fprintf(fout5,"%f\t",tt);
				for (ii=1; ii<=zbin; ii++) {
					for (jj=1; jj<=NSP; jj++) {
						fprintf(fout5, "%e\t", xx[ii][jj]);
					}
				}
				for (ii=1; ii<=zbin; ii++) {
					fprintf(fout5, "%f\t", tl[ii]);
				}
				for (ii=0; ii<=zbin; ii++) {
					fprintf(fout5, "%f\t", z[ii]);
				}
				fprintf(fout5,"\n");;
			}
			break; /* Completed! */
		}
		
        for (ii=1; ii<=nn; ii++) { Con[ii] = fmax(Con[ii]+tstep*fvec[ii], 0.0);} /* stepping! */
		
		/* Correct the lower boundary with fixed mixing ratio */
		for (ii=1; ii<=numx; ii++) {
			if (listFix[ii]==1) {
				Con[ii]=ConFix[ii];
			}
		}
		
		if (j==0) {
			tt=tstep;
		}else {
			tt += tstep;
		}
				
		tstepold=tstep;
		if (emax>0.15) tstep=tstepold*0.9;
		if (emax>0.20) tstep=tstepold*0.7;
		if (emax<0.10) tstep=tstepold*1.1;
		if (emax<0.05) tstep=tstepold*1.3;
		if (emax<0.03) tstep=tstepold*1.5;
        if (emax<0.01) tstep=tstepold*2.0;
        if (emax<0.003) tstep=tstepold*5.0;
        if (emax<0.001) tstep=tstepold*10.0;
		
		if (tstep>TMAX) { tstep=TMAX; }
		if (tstep<TMIN) { tstep=TMIN; }
		if (tstep>tslimit) { tstep=tslimit; }
		
		Convert1(Con, ConC, Conf, labelx, labelc, labelf);	
		
		/* record the result at each step */
		if (fmod(j, NPRINT)==0 && HISTORYPRINT == 1) {
			fprintf(fout5,"%f\t",tt);
			for (ii=1; ii<=zbin; ii++) {
				for (jj=1; jj<=NSP; jj++) {
					fprintf(fout5, "%e\t", xx[ii][jj]);
				}
			}
			for (ii=1; ii<=zbin; ii++) {
				fprintf(fout5, "%f\t", tl[ii]);
			}
			for (ii=0; ii<=zbin; ii++) {
				fprintf(fout5, "%f\t", z[ii]);
			}
			fprintf(fout5,"\n");;
		}
		
		printf("%s %d %s %e %s %e %s\n", "finish loop", j+1, "at time", tt, "s with the timestep", tstepold, "s");
	
		RadTransfer(rad, opt, stdcross, qysum, cross, crosst, iradmax+1);
		
		if (fmod(j, NPRINT)==0) {
			printout(labelx, labelf, rad, iradmax, z);
			printout_c();
			printout_std(z);
			printoutrate(zone_r, zone_m, zone_t, JJ, zone_p);
			printf("%s\n","print out");
			
			GlobalBalance(Con, labelx, listAER, zone_r, zone_m, zone_t, JJ, zone_p, Upflux, Loflux, Depo);
			
			MeanOpacity(stdcross,qysum,cross,crosst);
			foutp=fopen(OUT_MEANOPAC,"w");
			for (jj=1; jj<=zbin; jj++) {
				fprintf(foutp, "%s %f %s %e %s\n", "The mean opacity at altitude", zl[jj], "km is", mkv[jj], "cm-1");
			}
			fclose(foutp);
			
			GreyTemp(rad, iradmax);
			NewPressure(P[0], z, GA);
			/* foutp=fopen(OUT_NEWTEMP,"w");
			for (jj=0; jj<=zbin; jj++) {
				fprintf(foutp, "%s %f %s %e %s\n", "The new temperature at altitude", z[jj], "km is", Tnew[jj], "K");
			}
			fclose(foutp); */
		}
		
	}
	
	fclose(fout4);
	fclose(fout5);
	
	/* General printout */
	printout(labelx, labelf, rad, iradmax, z);
	printout_c();
	printout_std(z);
	/* printoutrate; */
	printoutrate(zone_r, zone_m, zone_t, JJ, zone_p);
	/* printout global balance */
	GlobalBalance(Con, labelx, listAER, zone_r, zone_m, zone_t, JJ, zone_p, Upflux, Loflux, Depo);
	printf("%s\n","print out");
	
	free_dmatrix(Jaco, 1, nn, 1, nn);
	free_dmatrix(Jaco1, 1, nn, 1, nn);
	free_dmatrix(cross,1,nump,0,WaveBin);
	free_dmatrix(qy,1,nump,0,WaveBin);
	free_dmatrix(crosst,1,nump,0,WaveBin);
	free_dmatrix(qyt,1,nump,0,WaveBin);
	free_dmatrix(JJ,1, zbin, 1,nump);
	free_dmatrix(opacCO2,1,zbin,0,WaveBin);
	free_dmatrix(opacO2,1,zbin,0,WaveBin);
	free_dmatrix(opacSO2,1,zbin,0,WaveBin);
	free_dmatrix(opacH2O,1,zbin,0,WaveBin);
	free_dmatrix(opacOH,1,zbin,0,WaveBin);
	free_dmatrix(opacH2CO,1,zbin,0,WaveBin);
	free_dmatrix(opacH2O2,1,zbin,0,WaveBin);
	free_dmatrix(opacHO2,1,zbin,0,WaveBin);
	free_dmatrix(opacH2S,1,zbin,0,WaveBin);
	free_dmatrix(opacCO,1,zbin,0,WaveBin);
	free_dmatrix(opacO3,1,zbin,0,WaveBin);
	free_dmatrix(opacCH4,1,zbin,0,WaveBin);
	free_dmatrix(opacNH3,1,zbin,0,WaveBin);
	free_dmatrix(DM, 1, zbin, 1, numx);
	free_dmatrix(dl, 1, zbin, 1, numx);
}
