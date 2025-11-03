/*Function to calculate temperature at a given height with grey atmosphere assumption */

#include <math.h>
#include "gasheat.c"

void GreyTemp(double **radiation, int NW);

void GreyTemp(double **radiation, int NW)
{

	int j;
	
	double Tirr, Tnew4;
	Tirr = STAR_TEMP/pow(ORBIT,0.5);
	
	double Fint, Tint;
	Fint = 0.0;
	for (j=0; j<NW; j++) {
		Fint += radiation[j][0]*(wavelength[j+1]-wavelength[j])*(1-PAB)*(1-PSURFAB);
	}
	Tint=pow(Fint/SIGMA, 0.25);
	Tint=TINTSET;
	printf("%s %f\n", "Tint is", Tint);
	
	double miu, gamma;
	miu=1/sqrt(3.0); /* ISOTROPIC IRRADIATION, GUILLOT 2010 EQ 29 */
	
	double tau[zbin+1], kth[zbin+1];
	
	for (j=1; j<=zbin; j++) {
		kth[j]=0.0;
		kth[j] += MeanCO2[j]*xx[j][52]/MM[j];  /* CO2  */
		kth[j] += MeanO2[j]*xx[j][54]/MM[j];   /* O2   */
		kth[j] += MeanSO2[j]*xx[j][43]/MM[j];  /* SO2  */
		kth[j] += MeanH2O[j]*xx[j][7]/MM[j];   /* H2O  */
		kth[j] += MeanOH[j]*xx[j][4]/MM[j];    /* OH   */
		kth[j] += MeanH2CO[j]*xx[j][22]/MM[j]; /* H2CO */
		kth[j] += MeanH2O2[j]*xx[j][6]/MM[j];  /* H2O2 */
		kth[j] += MeanHO2[j]*xx[j][5]/MM[j];   /* HO2  */
		kth[j] += MeanH2S[j]*xx[j][45]/MM[j];  /* H2S  */
		kth[j] += MeanCO[j]*xx[j][20]/MM[j];   /* CO   */
		kth[j] += MeanO3[j]*xx[j][2]/MM[j];    /* O3   */
		kth[j] += MeanCH4[j]*xx[j][21]/MM[j];   /* CH4   */
		kth[j] += MeanNH3[j]*xx[j][9]/MM[j];    /* NH3   */
	}
	
	for (j=1; j<=zbin; j++) {
		printf("%s %f %s %e %e %f\n", "The mean opacity at altitude", zl[j], "km is", mkv[j], kth[j], mkv[j]/kth[j]);
	}
	
	double kvtotal, kthtotal, tauint;
	kvtotal=0.0;
	kthtotal=0.0;
	tauint=0.0;
	for (j=zbin; j>=1; j--) {
		if (tauint<10.0) {
			kvtotal += mkv[j]*thickl;
			kthtotal+= kth[j]*thickl;
			tauint += kth[j]*thickl*MM[j];
		}
	}
	gamma = kvtotal/kthtotal;
	printf("%s %f\n","The gamma factor is", gamma);
	
	tau[zbin]=0;
	j=zbin-1;
	while (j>=0) {
		tau[j] = tau[j+1]+kth[j+1]*thickl*MM[j+1];
		j=j-1;
	}
	
	for (j=0; j<=zbin; j++) {
		
		Tnew4=3.0*pow(Tirr,4.0)/4.0*(2.0/3.0+miu/gamma+(gamma/3.0/miu-miu/gamma)*exp(-gamma*tau[j]/miu))*FADV*(1.0-PAB); /* ISOTROPIC IRRADIATION, GUILLOT 2010 EQ 29 */
		/* printf("%d %e %f\n", j, tau[j], pow(Tnew4,0.25)); */
		Tnew4 += 3.0*pow(Tint,4.0)/4.0*(2.0/3.0+tau[j]);
		Tnew[j]=pow(Tnew4,0.25);
		
	}
	
	/* check adiabats */
	double lapse[zbin+1], dtdz, gasheat, GA;
	GA=GRAVITY*MASS_PLANET/RADIUS_PLANET/RADIUS_PLANET; /* Planet Surface Gravity Acceleration, in SI */
	for (j=1; j<=zbin; j++) {
		if (RefIdxType == 0) { gasheat=AirHeat(tl[j]);}
		if (RefIdxType == 1) { gasheat=CO2Heat(tl[j]);}
		if (RefIdxType == 2) { gasheat=HeHeat(tl[j]);}
		if (RefIdxType == 3) { gasheat=N2Heat(tl[j]);}
		if (RefIdxType == 4) { gasheat=NH3Heat(tl[j]);}
		if (RefIdxType == 5) { gasheat=CH4Heat(tl[j]);}
		if (RefIdxType == 6) { gasheat=H2Heat(tl[j]);}
		if (RefIdxType == 7) { gasheat=O2Heat(tl[j]);}
		lapse[j] = GA*AIRM/1000.0/gasheat;
		/* printf("%d %e\n",j,lapse[j]); */
	}
	
	for (j=zbin; j>=1; j--) {
		dtdz = (Tnew[j-1] - Tnew[j])/(thickl*1.0E-2);
		if (dtdz > lapse[j]) {
			Tnew[j-1] = Tnew[j] + thickl*1.0E-2*lapse[j];
		}
	}
	
	/* for (j=0; j<zbin; j++) {
		dtdz = (Tnew[j] - Tnew[j+1])/(thickl*1.0E-2);
		if (dtdz > lapse[j]) {
			Tnew[j+1] = Tnew[j] - thickl*1.0E-2*lapse[j] ;
		}
	} */
	
	
}
