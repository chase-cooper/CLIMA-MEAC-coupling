#include <math.h>

/* cm3 molecule-1 s-1 */
/* M in molecule cm-3 */

// function to compute reaction rate of the reverse reaction when the forward reaction is known
// product/reactant refer to the forward reaction
double ReverseRate_2(int mu[], int n_prod, int n_react, double a1[], double a2[], double a3[], double a4[], double a5[], double a6[],double a7[], double T, double k_f) {
     int delta_n, n_total;
     double k;
     double delta_a1 = 0.0;
     double delta_a2 = 0.0;
     double delta_a3 = 0.0;
     double delta_a4 = 0.0;
     double delta_a5 = 0.0;
     double delta_a6 = 0.0;
     double delta_a7 = 0.0;
     double K_c, k_r;

     k=1.380658e-16; //boltzmann constant, erg/K

     // mu -- stoiochiometic coefficients, positive for products, negative for reactants
     // a -- thermodynamics coefficienets for each species
     // n_total -- total number of species

     // Estimate change in mols between products and reactants
     n_total = n_prod + n_react; // total number of reactants and products
     delta_n = n_prod - n_react; // difference in number of reactants and products

    // Calaulate delta_ai
     for (int j = 0; j < n_total; j++) {
          delta_a1 += mu[j]*a1[j];
          delta_a2 += mu[j]*a2[j];
          delta_a3 += mu[j]*a3[j];
          delta_a4 += mu[j]*a4[j];
          delta_a5 += mu[j]*a5[j];
          delta_a6 += mu[j]*a6[j];
          delta_a7 += mu[j]*a7[j];
     }

    // Calculate equilibrium constant (Rimmer&Helling+2016 and refernce therein)
     K_c = pow(k*1.0E-6*T, -delta_n)*exp(delta_a1*(log(T)-1)+delta_a2*T/2+delta_a3*pow(T,2)/6+delta_a4*pow(T,3)/12+delta_a5*pow(T,4)/20-delta_a6/T+delta_a7);

    // Calculate reverse reaction rate
    k_r = k_f/K_c;

    return k_r;
}

void ReactionRateM()
{
  int i, j;
  double RH, LH, ind, K0, K2, K3, Kf, Kinf, kint, k_inf, k_0_M, F_c, F, M_M_c, N;
  double n_prod, n_reac, k_f, k_f_0_M, k_f_inf;
	
  for (i=1; i<=zbin; i++) {
    
	k_0_M=5.46E-31*pow(tl[i]/298.0,-1.6)*MM[i]*THREEBODY;
	k_inf=2.16E-11;
    kkM[i][1]=k_0_M/(1.0+k_0_M/k_inf);
	
	k_0_M=6.89E-32*MM[i]*THREEBODY;
	k_inf=2.06E-11*exp(-56.5/tl[i]);
    kkM[i][2]=k_0_M/(1.0+k_0_M/k_inf);
	
	k_0_M=1.68E-24*pow(tl[i]/298.0,-7.0)*exp(-1390.0/tl[i])*MM[i]*THREEBODY;
	k_inf=6.0E-11;
	F_c=0.38*exp(-tl[i]/73.0)+0.62*exp(-tl[i]/1180.0);
	M_M_c=k_0_M/k_inf;
	N=0.75-1.27*log10(F_c);
	F=pow(F_c, 1.0/(1.0+pow((log10(M_M_c)/N),2.0)));
    kkM[i][3]=k_0_M*k_inf/(k_0_M+k_inf)*F;
    
	k_0_M=4.1E-31*pow(tl[i]/298.0, -3.6)*MM[i]*THREEBODY;
    k_inf=1.2E-12*pow(tl[i]/298.0, 1.1);
    ind=1.0/(1.0+pow(log10(k_0_M/k_inf),2.0));
    kkM[i][4]=((k_0_M*k_inf)/(k_inf+k_0_M))*pow(0.6,ind);
    
	k_0_M=4.1E-30*pow(tl[i]/298.0,-2.1)*MM[i]*THREEBODY;
	k_inf=3.00E-11*pow(tl[i]/298.0, -0.90);
    kkM[i][5]=k_0_M/(1.0+k_0_M/k_inf);
	
	k_0_M=7.83E-29*pow(tl[i]/298.0,-7.56)*exp(-5490.7/tl[i])*MM[i]*THREEBODY;
	k_inf=1.9E-13*pow(tl[i]/298.0, 2.25)*exp(-3033.4/tl[i]);
    kkM[i][6]=k_0_M/(1.0+k_0_M/k_inf);
	
	k_0_M=1.750E-16*pow(tl[i], -4.664)*exp(-1902.0/tl[i])*MM[i]*THREEBODY;
	k_inf=2.84E-14*pow(tl[i],1.266)*exp(-1363/tl[i]);
    kkM[i][7]=k_0_M/(1.0+k_0_M/k_inf); //3.31E-30*exp(-740.0/tl[i])*MM[i];
	
	k_0_M=1.38E-30*MM[i]*THREEBODY;
	k_inf=3.68E-12*pow(tl[i]/298.0, 1.61)*exp(-1321.8/tl[i]);
    kkM[i][8]=k_0_M/(1.0+k_0_M/k_inf);
	
	k_0_M=6.80E-23*pow(tl[i]/300.0, -6.20)*exp(-2394.6/tl[i])*MM[i]*THREEBODY;
	k_inf=8.62E-12*pow(tl[i]/300.0, 1.87)*exp(-586.4/tl[i]);
    kkM[i][9]=k_0_M/(1.0+k_0_M/k_inf);
	
	k_0_M=3.21E-30*pow(tl[i]/298.0,-2.57)*exp(-215.0/tl[i])*MM[i]*THREEBODY;
	k_inf=7.77E-14*exp(2280/tl[i]);
    kkM[i][10]=k_0_M/(1.0+k_0_M/k_inf);
	
    k_0_M=8.63E-30*pow(tl[i]/298.0, -2.20)*exp(-567/tl[i])*MM[i]*THREEBODY;
    k_inf=1.73E-10*pow(tl[i]/298,-0.50);
	M_M_c=k_0_M/k_inf;
    F_c=0.95-(1.0E-4)*tl[i];
    F=pow(F_c, 1.0/(1.0+pow(log10(M_M_c),2.0)));
	kkM[i][11]=k_0_M/(1.0+k_0_M/k_inf)*F;
    
	k_0_M=5.29E-34*exp(-370.0/tl[i])*MM[i]*THREEBODY;
	k_inf=1.96E-13*exp(-1370.0/tl[i]);
	kkM[i][12]=k_0_M/(1.0+k_0_M/k_inf);
	
	k_0_M=6.04E-33*pow(tl[i]/298.0, -1.0)*MM[i]*THREEBODY;
	k_inf=1.81E-13*exp(-754.0/tl[i]);
    kkM[i][13]=k_0_M/(1.0+k_0_M/k_inf);
	
	k_0_M=4.40E-30*pow(tl[i]/298.0, -1.76)*MM[i]*THREEBODY;
	k_inf=2.6E-10;
    kkM[i][14]=k_0_M/(1.0+k_0_M/k_inf);
    
    RH=2.44E-10*pow(tl[i]/298.0, -0.41);
    LH=1.34E-31*pow(tl[i]/298.0, -1.32)*exp(-370.5/tl[i])*THREEBODY;
    ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2));
    kkM[i][15]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
    
	k_0_M=4.36E-32*pow(tl[i]/298.0,-1.0)*MM[i]*THREEBODY;
	k_inf=1.0E-11;
    kkM[i][16]=k_0_M/(1.0+k_0_M/k_inf);
	
	k_0_M=5.3E-32*pow(tl[i]/298.0, -1.8)*MM[i]*THREEBODY;
    k_inf=9.5E-11*pow(tl[i]/298.0, 0.4);
    ind=1.0/(1.0+pow(log10(k_0_M/k_inf),2.0));
    kkM[i][17]=(k_0_M/(1+k_0_M/k_inf))*pow(0.6, ind);
	
	k_0_M=6.87E-31*pow(tl[i]/298.0, -2.0)*MM[i]*THREEBODY; 
    k_inf=4.17E-11*pow(tl[i],0.234)*exp(57.5/tl[i]); 
    kkM[i][18]=k_0_M/(1.0+k_0_M/k_inf);
	
	
	k_0_M=2.1E-33*exp(920.0/tl[i])*MM[i]*THREEBODY;
    k_inf=7.0E-12*exp(620.0/tl[i]);
    kkM[i][19]=k_0_M/(1.0+k_0_M/k_inf);
    
    RH=4.0E-12*pow(tl[i]/298.0, -0.3);
    LH=1.9E-31*pow(tl[i]/298.0, -3.4)*THREEBODY;
    ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
    kkM[i][20]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
    
	k_0_M=9.18E-34*pow(tl[i]/298.0, -1.69)*MM[i]*THREEBODY;
    k_inf=2.0E-10*pow(tl[i]/298.0,0.31);
    kkM[i][21]=k_0_M/(1.0+k_0_M/k_inf);
	
	k_0_M=9.4E-33*MM[i]*THREEBODY;
    k_inf=1.0E-10;
    kkM[i][22]=k_0_M/(1.0+k_0_M/k_inf);
	
	k_0_M=5.0E-32*MM[i]*THREEBODY;
    k_inf=1.0E-10;	
    kkM[i][23]=k_0_M/(1.0+k_0_M/k_inf);
	
	k_0_M=8.3E-38*MM[i]*THREEBODY;
    k_inf=1.94E-20;
    kkM[i][24]=k_0_M/(1.0+k_0_M/k_inf);

	k_0_M=1.38E-33*exp(502.7/tl[i])*MM[i]*THREEBODY; 
    k_inf=5.0E-16;
    kkM[i][25]=k_0_M/(1.0+k_0_M/k_inf);

	k_0_M=5.46E-33*exp(155.2/tl[i])*MM[i]*THREEBODY;
    k_inf=1.0E-10;
    kkM[i][26]=k_0_M/(1.0+k_0_M/k_inf);
    
    RH=3.8E-11*pow(tl[i]/298.0, -0.6);
    LH=2.3E-29*pow(tl[i]/298.0, -2.8)*THREEBODY;
    ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
    kkM[i][27]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
    
    RH=1.9E-11*pow(tl[i]/298.0, -1.8);
    LH=5.5E-29*pow(tl[i]/298.0, -4.4)*THREEBODY;
    ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
    kkM[i][28]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
    
    RH=1.6E-12*pow(tl[i]/298.0, 0.1);
    LH=2.4E-30*pow(tl[i]/298.0, -3.0)*THREEBODY;
    ind=1/(1.0+pow(log10(LH*MM[i]/RH),2.0));
    kkM[i][29]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
    
	k_0_M=2.0E-34*MM[i]*THREEBODY;
	k_inf=4.82E-15*pow(tl[i]/300.0,-1.0);
    kkM[i][30]=k_0_M/(1.0+k_0_M/k_inf);
	
	k_0_M=1.7E-33*exp(-1509.0/tl[i])*MM[i]*THREEBODY;
	k_inf=1.0E-14*exp(-1630.0/tl[i]);
    kkM[i][31]=k_0_M/(1.0+k_0_M/k_inf);
    
    RH=3.0E-11;
    LH=9.1E-32*pow(tl[i]/298.0, -1.5)*THREEBODY;
    ind=1.0/(1+pow(log10(LH*MM[i]/RH),2));
    kkM[i][32]=LH*MM[i]/(1+LH*MM[i]/RH)*pow(0.6,ind);
    
//    RH=2.3E-11*pow(298.0/tl[i], 0.2);
//    LH=3.4E-31*pow(298.0/tl[i], 1.6)*THREEBODY;
//    ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
	
	K0=3.4E-31*pow((298.0/tl[i]), 1.6)*THREEBODY;
    Kinf=2.3E-11*pow((298.0/tl[i]), 0.2);
    ind=1.0/(1.0+pow(log10(K0*MM[i]/Kinf),2));
    kkM[i][33]=(K0*MM[i]/(1.0+K0*MM[i]/Kinf))*pow(0.6, ind);
    
	k_0_M=5.21E-35*exp(900.0/tl[i])*MM[i]*THREEBODY;
	k_inf=1.21E-11*pow(tl[i]/300.0,-2.0);
    kkM[i][34]=k_0_M/(1.0+k_0_M/k_inf);

	k_0_M=6.0E-34*pow(tl[i]/298.0, -2.4)*MM[i]*THREEBODY;
	k_inf=2.8E-11;
    kkM[i][35]=k_0_M/(1.0+k_0_M/k_inf);
	
	k_0_M=2.8E-36*pow(tl[i]/298.0, -0.9)*MM[i]*THREEBODY;
	k_inf=3.4E-16;
    kkM[i][36]=k_0_M/(1.0+k_0_M/k_inf);
    
    // RH=8.3E-13*pow(tl[i]/298.0, 2.0);
    // LH=5.5E-30*THREEBODY;
    // ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
    // kkM[i][37]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
	kkM[i][37] = 0;
    
    RH=8.5E-12*pow(tl[i]/298.0, -1.75);
    LH=1.1E-28*pow(tl[i]/298.0, -3.5)*THREEBODY;
    ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2));
    kkM[i][38]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
    
//     RH=1.1E-12*pow(tl[i]/300.0, 1.3);
//     //LH=5.9E-33*pow(tl[i]/300.0, -1.4); //Update R39: OH + CO + M --> HCO2 + M from the Sander+2011 recommendation to the Burkholder+2015 recommendation
//     LH=5.9E-33*pow(tl[i]/300.0, -1.0); 
//     ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
//     kkM[i][39]=LH*MM[i]/(1+LH*MM[i]/RH)*pow(0.6,ind);
	
	K0=6.9E-33*pow((298.0/tl[i]), 2.1)*THREEBODY; //Update R39 OH + CO + M --> HCO2 + M to the Burkholder+2020 recommendation
	Kinf=1.1E-12*pow((298.0/tl[i]), -1.3); //JPL 2020
	ind=1.0/(1.0+pow(log10(K0*MM[i]/Kinf),2.0)); //JPL 2020
	kkM[i][39]=(K0*MM[i]/(1.0+K0*MM[i]/Kinf))*pow(0.6, ind); //JPL2020
     
	//kkM[i][40]=1.8E-27*pow(tl[i]/298.0,-3.85)*MM[i]*THREEBODY; //Origin of this rate constant unknown
    k_0_M=2.0E-23*pow(tl[i]/1.0, -1.3)*exp(-362.0/tl[i])*MM[i]*THREEBODY;
    k_inf=1.50E-10;
    kkM[i][40]=k_0_M/(1.0+k_0_M/k_inf);
	  
    K0=2.4E-14*exp(460.0/tl[i]);
    K2=2.7E-17*exp(2199.0/tl[i]);
    K3=6.5E-34*exp(1355.0/tl[i])*THREEBODY;
    kkM[i][41]=K0+K3*MM[i]/(1.0+K3*MM[i]/K2);
    
    RH=3.6E-11*pow(tl[i]/298.0, -0.1);
    LH=7.1E-31*pow(tl[i]/298.0, -2.6)*THREEBODY;
    ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
    kkM[i][42]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
    
    RH=2.8E-11;
    LH=1.8E-30*pow(tl[i]/298.0, -3.0)*THREEBODY;
    ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
    kkM[i][43]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
    
    RH=2.6E-11;
    LH=6.9E-31*pow(tl[i]/298.0, -1.0)*THREEBODY;
    ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
    kkM[i][44]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
    
    //kkM[i][45]=5.0E-32*exp(900.0/tl[i])*MM[i]*THREEBODY;
	k_0_M=1.98E-33*exp(206.0/tl[i])*MM[i]*THREEBODY;
	k_inf=2.26E-14*exp(415.0/tl[i]);
    kkM[i][45]=k_0_M/(1.0+k_0_M/k_inf);	
	
	k_0_M=4.82E-31*pow(tl[i]/298.0,-2.17)*MM[i]*THREEBODY;
    k_inf=5.30E-11;
    kkM[i][46]=k_0_M/(1.0+k_0_M/k_inf);

	//Forward reaction HOSO + M --> OH + SO + M
    k_f_0_M=1.92E22*pow(tl[i], -9.02)*exp(-26647/tl[i])*MM[i]*THREEBODY;
    k_f_inf=9.94E21*pow(tl[i], -2.54)*exp(-38190/tl[i]);
    k_f=k_f_0_M/(1.0+k_f_0_M/k_f_inf);

    n_prod=2;
    n_reac=1;
    int mu[3] = {-1, -1, 1}; //HOSO + M --> OH + SO + M is the FORWARD reaction we are reversing
    //thermo coefficents for [HOSO, OH, SO]. For 200-1000K. Via https://respecth.elte.hu/
    double a1[3] = {3.73732792, 3.99198424, 3.61859514};
    double a2[3] = {0.00930978241, -0.00240106655, -0.00232173768};
    double a3[3] = {-0.00000385404197, 0.00000461664033, 0.0000116462669};
    double a4[3] = {-0.00000000407782859, -0.00000000387916306, -0.000000014209251};
    double a5[3] = {0.0000000000030963234, 0.00000000000136319502, 0.0000000000056076537};
    double a6[3] = {-29826.0746, 3368.89836, -480.621641};
    double a7[3] = {10.1080861, -0.103998477, 6.36504115};

    kkM[i][47]=ReverseRate_2(mu, n_prod, n_reac, a1, a2, a3, a4, a5, a6, a7, tl[i], k_f);
	
	
    k_0_M=5.61E-30*pow(tl[i]/298.0,-5.19)*exp(-2271.0/tl[i])*MM[i]*THREEBODY;
	k_inf=7.58E-12*pow(tl[i]/298.0,1.59)*exp(-1243.6/tl[i]);
    kkM[i][48]=k_0_M/(1.0+k_0_M/k_inf);
	
    RH=4.1E-14*pow(tl[i]/298.0, 1.8);
    LH=1.8E-33*pow(tl[i]/298.0, 2.0)*THREEBODY;
    ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
    kkM[i][49]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
    
    RH=1.7E-12*pow(tl[i]/298.0,0.2);
    LH=2.9E-31*pow(tl[i]/298.0, -4.1)*THREEBODY;
    ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
    kkM[i][50]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
	  
	RH=7.5E-11;
	LH=5.7E-32*pow(tl[i]/300.0, -1.6)*THREEBODY;
	ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
	kkM[i][51]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
	  
	k_0_M=1.1E-30*pow(tl[i]/300.0, -2.0)*MM[i]*THREEBODY;
	k_inf=8.0E-11;
	kkM[i][52]=k_0_M/(1.0+k_0_M/k_inf);
	  
	k_0_M=1.1E-30*pow(tl[i]/300.0, -2.0)*MM[i]*THREEBODY;
	k_inf=8.0E-11;
	kkM[i][53]=k_0_M/(1.0+k_0_M/k_inf);

	k_0_M=4.0E-31*exp(900.0/tl[i])*MM[i]*THREEBODY;
	k_inf=3.0E-11;
	kkM[i][54]=k_0_M/(1.0+k_0_M/k_inf);
	
	k_0_M=4.0E-31*exp(900.0/tl[i])*MM[i]*THREEBODY;
	k_inf=1.0E-11;
	kkM[i][55]=k_0_M/(1.0+k_0_M/k_inf);
	  
	RH=2.37E-12*exp(523.0/tl[i]);
	LH=5.8E-30*exp(355.0/tl[i])*THREEBODY;
	ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
	kkM[i][56]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
	
	k_0_M=6.2E-29*pow(tl[i]/298.0,-1.8)*MM[i]*THREEBODY;
	k_inf=3.50E-10;
	kkM[i][57]=k_0_M/(1.0+k_0_M/k_inf);
	  
	  RH=1.5E-10;
	  LH=5.5E-23*pow(tl[i],-2.0)*exp(-1040/tl[i])*THREEBODY;
	  ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
	  kkM[i][58]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
	  
	  kkM[i][59]=4.42E-10*pow(tl[i]/298.0,0.22)*exp(43.0/tl[i]); // /THREEBODY;
	  kkM[i][60]=1.9E-10*pow(tl[i]/298.0,0.1)*exp(16.0/tl[i]); // /THREEBODY;
	  kkM[i][61]=1.1E-10*pow(tl[i]/298.0,0.21)*exp(86.6/tl[i]); // /THREEBODY;
	  kkM[i][62]=9.49E-12*pow(tl[i]/298.0,1.74)*exp(3872.8/tl[i]); // /THREEBODY;
	  kkM[i][63]=9.61E-12*exp(-1560.0/tl[i]); // /THREEBODY;
	  kkM[i][64]=9.49E-12*pow(tl[i]/298.0,1.74)*exp(3872.8/tl[i]); // /THREEBODY;
	  kkM[i][65]=1.01E-11*pow(tl[i]/298.0,0.69)*exp(1509.0/tl[i]); // /THREEBODY;
	  kkM[i][66]=2.64E-10*pow(tl[i]/298.0,0.18)*exp(62.5/tl[i]); // /THREEBODY;
	  kkM[i][67]=7.59E-12*pow(tl[i]/298.0,0.51)*exp(-1318.2/tl[i]); // /THREEBODY;
	  kkM[i][68]=3.32E-11; // /THREEBODY;
	  kkM[i][69]=4.23E-12*pow(tl[i]/298.0,2.54)*exp(-3400/tl[i]); // /THREEBODY;
	  
	  RH=1.39E-10*exp(-1184/tl[i]);
	  LH=1.0E-28*THREEBODY;
	  ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
	  kkM[i][70]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
	  
	  kkM[i][71]=6.68E-11*pow(tl[i]/298.0,0.31)*exp(93.8/tl[i])/THREEBODY;
	  
	  RH=1.39E-10*exp(-1184/tl[i]);
	  LH=1.0E-28*THREEBODY;
	  ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
	  kkM[i][72]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
	  
// 	  RH=2.37E-12*exp(523.0/tl[i]);
// 	  LH=5.8E-30*exp(355.0/tl[i])*THREEBODY;
// 	  ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
// 	  kkM[i][73]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
	k_0_M=9.0E-31*exp(550.0/tl[i])*MM[i]*THREEBODY;
    k_inf=8.55E-11*pow(tl[i],0.15);
	kkM[i][73]=k_0_M/(1.0+k_0_M/k_inf);
	  
	  kkM[i][74]=1.2E-10; // /THREEBODY;
	  kkM[i][75]=5.6E-11; // /THREEBODY;
	  kkM[i][76]=8.3E-12; // /THREEBODY;
	  kkM[i][77]=8.3E-12; // /THREEBODY;
	  kkM[i][78]=2.72E-11*pow(tl[i]/298.0,-0.32)*exp(66.1/tl[i]); // /THREEBODY;
	  kkM[i][79]=1.27E-14*pow(tl[i]/298.0,2.67)*exp(-3447.0/tl[i]); // /THREEBODY;
	  kkM[i][80]=5.16E-11*pow(tl[i]/298.0,-0.32); // /THREEBODY;
	  kkM[i][81]=2.09E-14*pow(tl[i]/298.0,1.9)*exp(-1060/tl[i]); // /THREEBODY;
	  
	  RH=8.2E-11;
	  LH=8.76E-6*pow(tl[i],-7.03)*exp(-1390/tl[i])*THREEBODY;
	  ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
	  kkM[i][82]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
	  
	  kkM[i][83]=6.5E-11;
	  kkM[i][84]=1.9E-11;
// 	  kkM[i][85]=1.96E-29*pow(tl[i]/298.0,-3.9)*MM[i]*THREEBODY;
	k_0_M=4.48E-14*pow(tl[i],-5.49)*exp(-1000.0/tl[i])*MM[i]*THREEBODY;
    k_inf=9.33E-10*pow(tl[i],-0.414)*exp(-33/tl[i]);
	kkM[i][85]=k_0_M/(1.0+k_0_M/k_inf);  
	  
	k_0_M=1.8E-27*pow(tl[i]/298.0,-3.85)*MM[i]*THREEBODY;
    k_inf=1.30E-10*pow(tl[i]/298.0,0.42);
	kkM[i][86]=k_0_M/(1.0+k_0_M/k_inf);
	  
// 	kkM[i][87]=2.2E-33*exp(-1780.0/tl[i])*MM[i]*THREEBODY; //THIS IS THE S+CO-->OCS RXN FROM ATMOS
	k_0_M=3.6E-34*pow(tl[i]/298.0,-0.57)*MM[i]*THREEBODY;
    k_inf=3.0E-14;
	kkM[i][87]=k_0_M/(1.0+k_0_M/k_inf);
	  
	  LH=1.22E-23*pow(tl[i],-3.0)*exp(-2900.0/tl[i])*THREEBODY;
	  RH=6.56E3*pow(tl[i],-5.0)*exp(-4000.0/tl[i]);
	  kkM[i][88]=MM[i]*LH/(1.0+(MM[i]*LH/RH)); //THIS IS BASED ON A FIT FROM WOGAN+2024 TO XU+2015 (in personal communication to SR 12/2023)
	  
//       LH=1.5E-13;
//       RH=2.1E9*pow(tl[i]/300.0, 6.1);
//       ind=1.0/(1.0+pow(log10(LH/(RH/MM[i])),2.0));
// 	  kkM[i][89]=LH/(1.0+LH/(RH/MM[i]))*pow(0.6,ind); //JPL2015
	  
	  K0=6.9E-33*pow((298.0/tl[i]), 2.1)*THREEBODY;
	  Kinf=1.1E-12*pow((298.0/tl[i]), -1.3);
      ind=1.0/(1.0+pow(log10(K0*MM[i]/Kinf),2.0));
      kint=1.85E-13*exp(-65.0/tl[i]);
      Kf=(K0*MM[i]/(1.0+K0*MM[i]/Kinf))*pow(0.6, ind);
	  kkM[i][89]=kint*(1.0-Kf/Kinf); //JPL2020 update

	  K0=3.4E-31*pow((298.0/tl[i]), 1.6)*THREEBODY;
	  Kinf=2.3E-11*pow((298.0/tl[i]), 0.2);
      ind=1.0/(1.0+pow(log10(K0*MM[i]/Kinf),2.0));
      kint=5.3E-12*exp(200.0/tl[i]);
      Kf=(K0*MM[i]/(1.0+K0*MM[i]/Kinf))*pow(0.6, ind);
	  kkM[i][90]=kint*(1.0-Kf/Kinf); //JPL2020 update

// 	  for (j=1; j<=NKinM; j++) {
// 		  kkM[i][j] = kkM[i][j]*THREEBODY;
// 	  }
  }
}











































