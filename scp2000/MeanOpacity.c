/*Function to provide mean opacity at a given height*/

#include <math.h>

void MeanOpacity(int stdn[], double qysum[], double **cross, double **crosst);

void MeanOpacity(int stdn[], double qysum[], double **cross, double **crosst)
{

	int i, j, k;
	
	double flux, product, mole2dust;
	mole2dust = PI*pow(AERSIZE,3)*AERDEN/6.0/AMU*2.0558;
	
	/* process temperature for cross sections */
	double temperature[zbin+1], crossl;
	for (j=1; j<=zbin; j++) {
		temperature[j] = tl[j];
		if (temperature[j] > TDEPMAX) {
			temperature[j] = TDEPMAX; 
		}
		if (temperature[j] < TDEPMIN) {
			temperature[j] = TDEPMIN; 
		}
	}
	
	double cx[WaveBin+1];
	
	for (j=1; j<=zbin; j++) {
		
		flux=0.0;
		product=0.0;
		
		for (i=0; i<=WaveBin; i++) {
			
			cx[i] = 1e-80;
			flux += solar[i]*(wavelength[i+1]-wavelength[i]);
			
			/* Molecular absorption associated with photolysis */
			/* for (k=3; k<=4; k++) {
				crossl = cross[k][i] + crosst[k][i]* ( temperature[j] - 295.0 ) ;
				cx[i] += solar[i]*crossl*xx[j][stdn[k]]/MM[j]/qysum[k];
			} */
			if (IFGREYAER == 1 ) {
				/* H2SO4AER */
				product += solar[i]*crossa[1][i]*xx[j][78]/MM[j]*98.0/mole2dust*(1.0-sinab[1][i]);
				/* S8AER */
				product += solar[i]*crossa[2][i]*xx[j][111]/MM[j]*256.0/mole2dust*(1.0-sinab[2][i]);
			}
			/* Molecular absorption not associated with photolysis */
			cx[i] += solar[i]*opacCO2[j][i]*xx[j][52]/MM[j];
			cx[i] += solar[i]*opacO2[j][i]*xx[j][54]/MM[j];
			cx[i] += solar[i]*opacSO2[j][i]*xx[j][43]/MM[j];
			cx[i] += solar[i]*opacH2O[j][i]*xx[j][7]/MM[j];
			cx[i] += solar[i]*opacOH[j][i]*xx[j][4]/MM[j];
			cx[i] += solar[i]*opacH2CO[j][i]*xx[j][22]/MM[j];
			cx[i] += solar[i]*opacH2O2[j][i]*xx[j][6]/MM[j];
			cx[i] += solar[i]*opacHO2[j][i]*xx[j][5]/MM[j];
			cx[i] += solar[i]*opacH2S[j][i]*xx[j][45]/MM[j];
			cx[i] += solar[i]*opacCO[j][i]*xx[j][20]/MM[j];
			cx[i] += solar[i]*opacO3[j][i]*xx[j][2]/MM[j];
			cx[i] += solar[i]*opacCH4[j][i]*xx[j][21]/MM[j];
			cx[i] += solar[i]*opacNH3[j][i]*xx[j][9]/MM[j];
		}
		
		for (i=0; i<WaveBin; i++) {
			if (fabs(cx[i+1]-cx[i])<1E-30) {
				product += (wavelength[i+1]-wavelength[i])*(cx[i+1]+cx[i])/2.0;
			} else {
				product += (wavelength[i+1]-wavelength[i])*(cx[i+1]-cx[i])/(log(cx[i+1])-log(cx[i]));
				/* printf("%d %d %e %e %e\n", j, i, product, cx[i], (log(cx[i+1])-log(cx[i]))); */
			}
		}
		
		mkv[j] = product/flux;
		printf("%d %e\n", j, mkv[j]);
	}
	
}
