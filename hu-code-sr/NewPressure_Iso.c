/*Function to calculate temperature at a given height with grey atmosphere assumption */

#include <math.h>
#include "constant.h"

void NewPressure(double P0, double z[], double GA);

void NewPressure(double P0, double z[], double GA)
{

	int j;
	FILE *fp;
	double logp;
	double h;
	
	fp=fopen(OUT_NEWTEMP,"w");
	Pnew[0] = P0;
	logp = log10(Pnew[0]);
	fprintf(fp, "%f\t%f\t%f\n", z[0], logp, Tnew[0]);
	for (j=1; j<=zbin; j++) {
		h = KBOLTZMANN * (Tnew[j-1]+Tnew[j]) /2.0 / AIRM / AMU / GA /1000.0 ;
		Pnew[j] = Pnew[j-1] * exp(-(z[j]-z[j-1])/h);
		logp = log10(Pnew[j]);
		fprintf(fp, "%f\t%f\t%f\n", z[j], logp, Tnew[j]);
	}
	fclose(fp);
	
}
