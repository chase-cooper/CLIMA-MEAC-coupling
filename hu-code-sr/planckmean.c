/*Function to calculate the infrared Planck mean opacity at each level */
/* input cross section in m^2 */
/* output mean cross section in cm^2 */
/* also output the cross section in cm^2 at each level in each wavelength */

#include <math.h>
#include "constant.h"

void planckmean(double Mean[], char Fname[], double **xsc);

void planckmean(double Mean[], char Fname[], double **xsc)
{		
	int i, j, k;
	FILE *fim;
	double ***opac;
	opac = f3tensor(0,NLAMBDA-1,0,NTEMP-1,0,NPRESSURE-1);
	double wave[NLAMBDA]; /* in m */
	double wave1[NLAMBDA]; /* in nm */
	double temp[NTEMP]; /* in k */
	double pres[NPRESSURE]; /* in Pa */
	double presdummy, cross[NLAMBDA];
	double Fplanck[zbin+1];
	double Cplanck[zbin+1];
	double pfitting, h1, h2, tfitting;
	
	printf("%s\n", Fname);
	fim = fopen(Fname,"r");	
	/* Header Lines */
	for (i=0; i<NTEMP-1; i++) {
		fscanf(fim, "%lf", temp+i);
	}
	i=NTEMP-1;
	fscanf(fim, "%lf\n", temp+i);
	for (i=0; i<NPRESSURE-1; i++) {
		fscanf(fim, "%le", pres+i);
	}
	i=NPRESSURE-1;
	fscanf(fim, "%le\n", pres+i);
	/* Read in data */
	for (i=0; i<NLAMBDA; i++) {
		fscanf(fim, "%le\n", wave+i);
		for (j=0; j<NPRESSURE; j++) {
			fscanf(fim, "%le", &presdummy);
			for (k=0; k<NTEMP-1; k++) {
				fscanf(fim, "%le", &opac[i][k][j]);
			}
			k=NTEMP-1;
			fscanf(fim, "%le\n", &opac[i][k][j]);
		}
	}
	fclose(fim);
	
	/* get wave1 from wave */
	for (i=0; i<NLAMBDA; i++) {
		wave1[i] = wave[i]*1.0E+9; /* convert to nm */
	}
	
	/* Compute the total planck flux in this range */
	for (i=1; i<=zbin; i++) {
		Fplanck[i] = 0.0;
		for (j=0; j<NLAMBDA-1; j++) {
			Fplanck[i] += 2*HPLANCK*pow(CLIGHT,2)/pow(wave[j],5)/(exp(HPLANCK*CLIGHT/wave[j]/KBOLTZMANN/tl[i])-1)*(wave[j+1]-wave[j]);
		}
	}
	
	/* Calculate the cross section */
	for (i=1; i<=zbin; i++) {
		Cplanck[i] = 0.0;
		for (j=0; j<NLAMBDA; j++) {
			if (pl[i]>pres[NPRESSURE-1]) {
				pfitting=pres[NPRESSURE-1];
			} else {
				pfitting=pl[i];
			}
			if (tl[i]>temp[NTEMP-1]) {
				tfitting=temp[NTEMP-1];
			} else if (tl[i]<temp[0]) {
				tfitting=temp[0];
			} else {
				tfitting=tl[i];
			}

			cross[j] = fmax(Interpolation2D(tfitting,pfitting,temp, NTEMP, pres, NPRESSURE, *(opac+j))*1.0E+4,1E-60); /* convert to cm^2 */
			/* printf("%e\n", cross[j]); */
			/* Cplanck[i] += cross[j]*2*HPLANCK*pow(CLIGHT,2)/pow(wave[j],5)/(exp(HPLANCK*CLIGHT/wave[j]/KBOLTZMANN/tl[i])-1)*(wave[j+1]-wave[j]); */
		}
		for (j=0; j<NLAMBDA-1; j++) {
			h1=2*HPLANCK*pow(CLIGHT,2)/pow(wave[j],5)/(exp(HPLANCK*CLIGHT/wave[j]/KBOLTZMANN/tl[i])-1);
			h2=2*HPLANCK*pow(CLIGHT,2)/pow(wave[j+1],5)/(exp(HPLANCK*CLIGHT/wave[j+1]/KBOLTZMANN/tl[i])-1);
			if (fabs(cross[j+1]-cross[j])<1E-20) {
				Cplanck[i] += (wave[j+1]-wave[j])/2.0*(cross[j]*h1+cross[j+1]*h2);
				/* printf("%d %d %e %e\n", i, j, Cplanck[i], cross[j]); */
			} else {
				Cplanck[i] += (wave[j+1]-wave[j])*(cross[j+1]*h2-cross[j]*h1)/(log(cross[j+1]*h2*1.0E+30)-log(cross[j]*h1*1.0E+30));
				/* printf("%d %d %e %e %e\n", i, j, Cplanck[i], cross[j]*h1, (log((cross[j+1]*h2 )*1.0E+60)-log((cross[j]*h1 )*1.0E+60)) ); */
			}
		}
		Interpolation(wavelength, WaveBin+1, *(xsc+i), wave1, cross, NLAMBDA, 0);
	}
	
	/* Calculate the planck mean */
	for (i=1; i<=zbin; i++) {
		Mean[i] = Cplanck[i]/Fplanck[i];
		printf("%s %e\n", Fname, Mean[i]);
	}
	
	/* free_dmatrix(opacj,1,NTEMP,1,NPRESSURE);*/ 
	free_f3tensor(opac,0,NLAMBDA-1,0,NTEMP-1,0,NPRESSURE-1);

}
