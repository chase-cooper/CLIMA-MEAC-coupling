void printout(int lx[], int lf[], double **rad, int iradmax, double z[]);

void printout(int lx[], int lf[], double **rad, int iradmax, double z[])
{
     int i=1,j,jj,s;
     FILE *fout1, *fout2, *fout3;
     char header1[10000], header2[10000], *temp;
	
	 double ntotal, airtotal;
	 airtotal = 0.0;
	 for (j=1; j<=zbin; j++) {
		 airtotal = airtotal + MM[j]*thickl;
	 }
     
     /* Check if it is the first time exporting data */
     fout1=fopen(OUT_FILE1, "r");
     s=LineNumber(fout1, 10000);
     fclose(fout1);
     
     if (s<5) {
              
     fout1=fopen(OUT_FILE1, "a");
     fout2=fopen(OUT_FILE2, "a");
     
    while (i<=zbin) {
        fprintf(fout1, "%f\t\t", zl[i]);
        fprintf(fout2, "%f\t\t", zl[i]);
		for (j=1; j<=numx; j++) { fprintf(fout1, "%e\t\t", xx[i][lx[j]]); }
		for (j=1; j<=numf; j++) { fprintf(fout2, "%e\t\t", xx[i][lf[j]]); }
		fprintf(fout1, "%e\n", MM[i]);
	    fprintf(fout2, "%e\n", MM[i]);
		i++;
    }
    
    } else {
     fout1=fopen(OUT_FILE1, "r");
     fout2=fopen(OUT_FILE2, "r");
     temp=fgets(header1, 10000, fout1); 
     temp=fgets(header2, 10000, fout2);
   	 fclose(fout1);
	 fclose(fout2);     
     fout1=fopen(OUT_FILE1, "w");
     fout2=fopen(OUT_FILE2, "w");
     fprintf(fout1, "%s\n", header1);
     fprintf(fout2, "%s\n", header2);
     while (i<=zbin) {
        fprintf(fout1, "%f\t\t", zl[i]);
        fprintf(fout2, "%f\t\t", zl[i]);
		for (j=1; j<=numx; j++) { fprintf(fout1, "%e\t\t", xx[i][lx[j]]); }
		for (j=1; j<=numf; j++) { fprintf(fout2, "%e\t\t", xx[i][lf[j]]); }
		fprintf(fout1, "%e\n", MM[i]);
	    fprintf(fout2, "%e\n", MM[i]);
		i++;
    }
    
    }
	
	fout3=fopen(OUT_RADIATION,"w");
	for (i=0; i<=zbin; i++) {
		fprintf(fout3, "%s %f %s\n", "z=", z[i], "km");
		for (jj=0; jj<=iradmax; jj++) {
			fprintf(fout3, "%f %e\n", wavelength[jj], rad[jj][i]);
		}
	}

    fclose(fout3);
	fclose(fout1);
	fclose(fout2);
	
	/* Integration over the column and print-out the column density */
	fout1 = fopen(OUT_COLUMN, "w");
	for (i=1; i<=numx; i++) {
		fprintf(fout1, "%s", "X");
		fprintf(fout1, "%d\t", lx[i]);
		ntotal = 0.0;
		for (j=1; j<=zbin; j++) {
			ntotal = ntotal + xx[j][lx[i]]*thickl;
		}
		fprintf(fout1, "%e\t%e\n", ntotal, ntotal/airtotal);
	}
	for (i=1; i<=numf; i++) {
		fprintf(fout1, "%s", "F");
		fprintf(fout1, "%d\t", lf[i]);
		ntotal = 0.0;
		for (j=1; j<=zbin; j++) {
			ntotal = ntotal + xx[j][lf[i]]*thickl;
		}
		fprintf(fout1, "%e\t%e\n", ntotal, ntotal/airtotal);
	}
	fclose(fout1);
	
}

void printout_c();

void printout_c()
{
	int i;
	FILE *fp;
	fp=fopen("Illustration/ConstantMixing.dat","w");
	fprintf(fp, "%s\t\t%s\t\t%s\n", "z", "H2O", "CO2");
	for (i=1; i<=zbin; i++) {
		fprintf(fp, "%lf\t%e\t%e\n", zl[i], xx[i][7], xx[i][52]);
	}
	fclose(fp);
}



void printout_timescale(double Con[], double fvec[]);

void printout_timescale(double Con[], double fvec[])
{
	int i,j,temp;
	FILE *fout;
	fout=fopen(OUT_TIMESCALE,"w");
	
	for (i=1; i<=zbin; i++) {
		fprintf(fout, "%f\t\t", zl[i]);
		for (j=1; j<=numx; j++) { 
			temp=(i-1)*numx+j;
			fprintf(fout, "%e\t\t", Con[temp]/fvec[temp]); 
		}
		fprintf(fout, "\n");
	}
	
	fclose(fout);
}

void printout_std(double z[]);

void printout_std(double z[])
{
	int i, j;
	FILE *fp;
	fp=fopen(OUT_STD,"w");
	fprintf(fp, "%s\t\t%s\t\t%s\t\t%s\t\t%s\t\t", "z", "z0", "z1", "T", "P");
	for (i=1; i<=NSP; i++) { fprintf(fp, "%d\t\t", i); }
	fprintf(fp, "\n");
	fprintf(fp, "%s\t\t%s\t\t%s\t\t%s\t\t%s\t\t", "km", "km", "km", "K", "Pa");
	fprintf(fp, "\n");
	for (j=1; j<=zbin; j++) {
		fprintf(fp, "%f\t", zl[j]);
		fprintf(fp, "%f\t%f\t", z[j-1], z[j]);
		fprintf(fp, "%f\t%e\t", tl[j], pl[j]);
		for (i=1; i<=NSP; i++) {
			fprintf(fp, "%e\t", xx[j][i]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}