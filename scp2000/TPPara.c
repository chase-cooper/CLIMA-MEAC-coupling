#include <math.h>

/* Calculating T-P Profile for Hot Jupiter
 
 Assuming the bottom of atmosphere is at 1 bar
 
 Output: P: Pressure profile in SI
         T: Temperature profile
 Input:  NL: number of layer
         P0: Pressure at the top of atmosphere in bar
		 T0: Temperature at the top of atmosphere
         T1: Temperature at the bottom of atmosphere or at the top of stratosphere
         Tinv: set to 1 if there is a temperature inversion
         P1: Pressure at the top of stratosphere
         P2: Pressure at the bottom of stratosphere
         T2: Temperature at the bottom of stratosphere
		 P3: Pressure at the bottom of atmosphere
		 T3: Temperature at the bottom of atmosphere			*/

void TPPara(double P[], double T[], int Tinv, int NL, 
			double P0, double T0, double T1, double P1, double P2, double T2, double P3, double T3)
{
	double grid[NL];
	double bar, gridi, a1, a2;
	int i;
	
	bar=pow(10,5);
	grid[0]=log(P0);
	T[0]=T0;
	P[0]=P0*bar;
	gridi=(log(P3)-grid[0])/(NL-1);
	if (Tinv==1) {
		i=1;
		a1=(log(P1)-grid[0])/sqrt(T1-T0);
		while (P[i-1]<P1*bar) {
			grid[i]=grid[i-1]+gridi;
			P[i]=exp(grid[i])*bar;
			T[i]=pow((grid[i]-grid[0])/a1,2)+T0;
			i++;
		}
		a2=(log(P2)-log(P1))/sqrt(T1-T2);
		while (i<NL) {
			grid[i]=grid[i-1]+gridi;
			P[i]=exp(grid[i])*bar;
			T[i]=pow((grid[i]-log(P2))/a2,2)+T2;
			i++;
		}
	} else {
		i=1;
		a1=(log(P1)-grid[0])/sqrt(T1-T0);
		while (P[i-1]<P1*bar) {
			grid[i]=grid[i-1]+gridi;
			P[i]=exp(grid[i])*bar;
			T[i]=pow((grid[i]-grid[0])/a1,2)+T0;
			i++;
		}
		a2=(log(P3)-log(P1))/sqrt(T3-T1);
		while (i<NL) {
			grid[i]=grid[i-1]+gridi;
			P[i]=exp(grid[i])*bar;
			T[i]=pow((grid[i]-log(P1))/a2,2)+T1;
			i++;
		}
		
	}

}
