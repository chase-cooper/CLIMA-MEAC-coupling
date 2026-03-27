#include <math.h>

/* s-1 */
/* M in molecule cm-3 */

// function to compute reaction rate of the reverse reaction when the forward reaction is known
// product/reactant refer to the forward reaction
double ReverseRate(int mu[], int n_prod, int n_react, double a1[], double a2[], double a3[], double a4[], double a5[], double a6[],double a7[], double T, double k_f) {
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

void ReactionRateT()
{
     
     int i, n_prod, n_react;
	double k_inf, k_0_M, k_f, k_f_0_M, k_f_inf, F, F_c, M_M_c,N;
     double k;
     k=1.380658e-16; //boltzmann constant, erg/K
     for(i=1;i<=zbin;i++)
     {
       k_0_M=7.16E-10*exp(-11200.0/tl[i])*MM[i]*THREEBODY;
	   k_inf=7.60E12*exp(-12268.0/tl[i]);
       kkT[i][1]=k_0_M/(1.0+k_0_M/k_inf);
	   
	   k_0_M=2.41E-8*pow(tl[i]/298.0, -1.18)*exp(-24415.0/tl[i])*MM[i]*THREEBODY;
	   k_inf=5.82E11*pow(tl[i]/300.0, -2.18)*exp(-24400.0/tl[i]);
       kkT[i][2]=k_0_M/(1.0+k_0_M/k_inf);
	   
	   k_0_M=8.43E-6*pow(tl[i]/298.0, -2.30)*exp(-24536.9/tl[i])*MM[i]*THREEBODY;
	   k_inf=3.37E14*pow(tl[i]/298.0, 0.90)*exp(-24536.9/tl[i]);	   
	   kkT[i][3]=k_0_M/(1.0+k_0_M/k_inf);
	   
	   k_0_M=1.55E-9*exp(-30190.0/tl[i])*MM[i]*THREEBODY;
	   k_inf=7.94E10*exp(-30911.7/tl[i]);	 	   
	   kkT[i][4]=k_0_M/(1.0+k_0_M/k_inf);
	   
	   k_0_M=1.88E-4*pow(tl[i]/298.0, -3.37)*exp(-37645.0/tl[i])*MM[i]*THREEBODY;
	   k_inf=5.48E15*pow(tl[i]/298.0,-1.27)*exp(-36925.6/tl[i]);
       kkT[i][5]=k_0_M/(1.0+k_0_M/k_inf);

	   k_0_M=2.51E-14*exp(-1230.0/tl[i])*MM[i]*THREEBODY;
	   k_inf=2.5E6*exp(-6100.0/tl[i]);
       kkT[i][6]=k_0_M/(1.0+k_0_M/k_inf);
	   
	   k_0_M=1.3E-3*pow(tl[i]/300.0, -3.5)*exp(-11000.0/tl[i])*MM[i]*THREEBODY;
	   k_inf=9.7E14*pow(tl[i]/300.0, 0.1)*exp(-11080/tl[i]);
       kkT[i][7]=k_0_M/(1.0+k_0_M/k_inf);
	   
	   k_0_M=0.01*pow(tl[i], -1.61)*exp(-25585.0/tl[i])*MM[i]*THREEBODY;
	   k_inf=1.2E16*pow(tl[i], -0.43)*exp(-24922.0/tl[i]);
       kkT[i][8]=k_0_M/(1.0+k_0_M/k_inf);
	
	   k_0_M=5.0E6*pow(tl[i], -3.8)*exp(-25340.0/tl[i])*MM[i]*THREEBODY;
	   k_inf=1.2E19*pow(tl[i], -1.23)*exp(-25010.0/tl[i]);
       kkT[i][9]=k_0_M/(1.0+k_0_M/k_inf);

	   k_0_M=1.15E-6*exp(-23092.0/tl[i])*MM[i]*THREEBODY;
	   k_inf=9.33E15*exp(-24657.0/tl[i]);
       kkT[i][10]=k_0_M/(1.0+k_0_M/k_inf);
	   
       kkT[i][11]=0.0; //8.0E-2*pow(tl[i]/298.0, -6.55)*exp(-26099.0/tl[i])*MM[i];
	   
       kkT[i][12]=0.0; //6.31E17*exp(-13110.0/tl[i]);
       
       k_0_M=4.1E-5*exp(-10650.0/tl[i])*MM[i]*THREEBODY;
	   k_inf=4.8E15*exp(-11170.0/tl[i]);
       kkT[i][13]=k_0_M/(1.0+k_0_M/k_inf);

       k_0_M=7.3E14*pow(tl[i],-6.1)*exp(-47300.0/tl[i])*MM[i]*THREEBODY;
	   k_inf=3.7E13*exp(-36220.0/tl[i]);	   
       kkT[i][14]=k_0_M/(1.0+k_0_M/k_inf);
	   
       k_0_M=2.81E-9*exp(-25720.0/tl[i])*MM[i]*THREEBODY;
	   k_inf=4.46E13*exp(-34340.0/tl[i]);
       kkT[i][15]=k_0_M/(1.0+k_0_M/k_inf);

	   k_0_M=6.73E-9*exp(-26662.1/tl[i])*MM[i]*THREEBODY;
	   k_inf=7.5E14*exp(-34578.1/tl[i]);
       kkT[i][16]=k_0_M/(1.0+k_0_M/k_inf);
	   
	   k_0_M=8.7E7*pow(tl[i], -6.76)*exp(-16462.0/tl[i])*MM[i]*THREEBODY;
	   k_inf=4.31E19*pow(tl[i], -3.86)*exp(-18254.0/tl[i]);
       kkT[i][17]=k_0_M/(1.0+k_0_M/k_inf);
	   
       kkT[i][18]=0.0;
	   
	   k_0_M=4.04E-07*pow(tl[i]/300.0, 3.39)*exp(-50200.0/tl[i])*MM[i]*THREEBODY;
	   k_inf=1.19E7*pow(tl[i], 2.39)*exp(-50130.9/tl[i]);  
       kkT[i][19]=k_0_M/(1.0+k_0_M/k_inf);

	   k_0_M=1.1E-7*exp(-33075.0/tl[i])*MM[i]*THREEBODY;
	   k_inf=1.70E16*exp(-45706.0/tl[i]);  
       kkT[i][20]=k_0_M/(1.0+k_0_M/k_inf);

       kkT[i][21]=2.03E10*pow(tl[i]/298.0, 1.22)*exp(-43539.0/tl[i]);
       k_0_M=(8.03e4*pow(tl[i]/300,-10.2)*exp(-52454/tl[i]) + 6.245*pow(tl[i]/300,-6.577)*exp(-48007/tl[i]) )*MM[i]*THREEBODY;
       k_inf=9.443e15*pow(tl[i]/300,-1.0117)*exp(-46156/tl[i]);
       kkT[i][22]=k_0_M/(1.0+k_0_M/k_inf);
       k_0_M = (1.773*pow((tl[i]/298.0),-7.502)*exp(-23531/tl[i]))*MM[i]*THREEBODY;
       k_inf = (5.684e16*pow((tl[i]/298.0),-1.153)*exp(-22270/tl[i]));
       kkT[i][23] = k_0_M / (1.0 + k_0_M / k_inf);
       kkT[i][24]=4.15E13*exp(-22731.0/tl[i]);
       k_0_M = 4.3e3*pow(tl[i],-3.4)*exp(-18020/tl[i])*MM[i]*THREEBODY;
       k_inf = 3.9e8*pow(tl[i],-1.62)*exp(-18650/tl[i]);
       kkT[i][25]=k_0_M / (1.0 + k_0_M / k_inf);
       k_0_M = 3.4e-7*exp(-39390/tl[i])*MM[i]*THREEBODY;
       k_inf = pow(10,(12.9))*pow((tl[i]),0.44)*exp(-44700/tl[i]);
       kkT[i][26]=k_0_M / (1.0 + k_0_M / k_inf);
       k_0_M = 6.63e9*pow(tl[i],-4.99)*exp(-20130/tl[i])*MM[i]*THREEBODY;
       k_inf = 1.11e10*pow(tl[i],1.037)*exp(-18504/tl[i]);
       kkT[i][27]=k_0_M / (1.0 + k_0_M / k_inf);
       k_0_M = pow(10,42.838)*pow(tl[i],-6.431)*exp(-53938/tl[i])*MM[i]*THREEBODY;
       k_inf = pow(10,20.947)*pow(tl[i],-1.228)*exp(-51439/tl[i]);
       kkT[i][28]=k_0_M / (1.0 + k_0_M / k_inf);
       const double R = 8.314462618;
       const double alpha = 0.5757;
       const double T3 = 237.0;
       const double T1 = 1652.0;
       const double T2 = 5069.0;
       double Ti = tl[i];
       k_inf = 5.61818e19 * pow(Ti, -1.28) * exp(-309677.92/(R*Ti));
       k_0_M   = 3.02258e28 * pow(Ti, -5.02) * exp(-317778.14/(R*Ti)) * 1e6/6.022e23;
       double Pr = (k_0_M * MM[i]) / k_inf;
       double Fcent = (1 - alpha)*exp(-Ti/T3) + alpha*exp(-Ti/T1) + exp(-T2/Ti);
       double c = -0.4 - 0.67*log10(Fcent);
       double n = 0.75 - 1.27*log10(Fcent);
       double d = 0.14;
       double logPr_c = log10(Pr) + c;
       double denominator = n - d*logPr_c;
       F = pow(10, log10(Fcent) / (1 + pow(logPr_c/denominator, 2)));
       kkT[i][29] = k_inf * (Pr/(1 + Pr)) * F;
       k_0_M = 5.98e-9*exp(-29828.0/tl[i])*MM[i]*THREEBODY;
       k_inf = 3e14*exp(-35700.0/tl[i]);
       kkT[i][30]=k_0_M / (1.0 + k_0_M / k_inf);
       double a1[] = {0.808679682, 3.43126659, 2.41723661};
       double a2[] = {0.0233615762, 0.000631146866, 0.017671704};
       double a3[] = {-0.0000355172234, -0.00000192914359, -0.00000904883576};
       double a4[] = {0.0000000280152958, 0.00000000240618712, -0.00000000103230911};
       double a5[] = {-0.00000000000850075165, -0.000000000000866679361, 0.00000000000199106024};
       double a6[] = {26428.9808, -18508.5918, -8112.20974};
       double a7[] = {13.9396761, 1.07990541, 12.5095416};
       k_inf = 9e-13 * pow(tl[i]/300.0, 2.0);
       k_0_M = 5e-30 * pow(tl[i]/300.0, -1.5) * MM[i];
       k_f = k_0_M / (1.0 + k_0_M / k_inf);
       kkT[i][31] = ReverseRate((int[]){-1, -1, 1}, 1, 2, a1, a2, a3, a4, a5, a6, a7, tl[i], k_f);
       k_0_M = 6.82e-3*pow(tl[i]/298.0,-8.62)*exp(-11300.0/tl[i])*MM[i]*THREEBODY;
       k_inf = 2e10*exp(-7550.0/tl[i]);
       kkT[i][32]=k_0_M / (1.0 + k_0_M / k_inf);
       k_0_M = 9.62e9*pow((tl[i]/298.0),-11.92)*exp(-376000/8.314/tl[i])*MM[i]*THREEBODY;
       k_inf = 1.45e17*exp(-353000/8.314/tl[i]);
       kkT[i][33]=k_0_M / (1.0 + k_0_M / k_inf);
       k_inf = 3e14*exp(-351000/(8.314*tl[i]));
       k_0_M = 1.24e-5*pow((tl[i]/300),(-1))*exp(-42200/tl[i])*MM[i]*THREEBODY;
       kkT[i][34]=k_0_M / (1.0 + k_0_M / k_inf);
       k_inf = 1e15*exp(-42800/(tl[i]));
       k_0_M = 4.14e-5*pow((tl[i]/300),(-1))*exp(-42800/tl[i])*MM[i]*THREEBODY;
       kkT[i][35]=k_0_M / (1.0 + k_0_M / k_inf);
       k_inf = 1.33e15*pow((tl[i]/298.0),-2.02)*exp(-12749.0/tl[i]);
       k_0_M = 0.35*pow((tl[i]/298.0),-10.23)*exp(-13600/tl[i])*MM[i]*THREEBODY;
       kkT[i][36]=k_0_M / (1.0 + k_0_M / k_inf);
       k_inf = 9.94e16*pow((tl[i]/298.0),-6.54)*exp(-13567/tl[i]);
       k_0_M = 4.32E06*pow((tl[i]/298.0),1.51)*exp(-7639/tl[i])*k*MM[i]*THREEBODY*tl[i]/1333;
       kkT[i][37]=k_0_M / (1.0 + k_0_M / k_inf);
       k_inf = 1.25E14*pow((tl[i]/298),(-2.15))*exp(-92030/(8.314*tl[i]));
       k_0_M = 1.28E13*pow((tl[i]/298),(-7.58))*exp(-10754/(tl[i]))*k*MM[i]*THREEBODY*tl[i]/1333;
       kkT[i][38]=k_0_M / (1.0 + k_0_M / k_inf);
	   k_0_M=3.6E4*pow(tl[i], -3.1)*exp(-51280.0/tl[i])*MM[i]*THREEBODY;
	   k_inf=6.0E13*exp(-50228.0/tl[i]);	   
       kkT[i][39]=k_0_M/(1.0+k_0_M/k_inf);

	   k_0_M=3.95E-6*pow(tl[i]/300.0, -1.90)*exp(-30100.0/tl[i])*MM[i]*THREEBODY;
	   k_inf=9.54E13*pow(tl[i]/300.0, -2.90)*exp(-30100.0/tl[i]);
       kkT[i][40]=k_0_M/(1.0+k_0_M/k_inf);

       k_0_M=3.3E-9*pow(tl[i]/298, 0.65)*exp(-38607.7/tl[i])*MM[i]*THREEBODY;
       k_inf=2.85e19*pow(tl[i],-1.52)*exp(-43137.0/tl[i]);
       kkT[i][41]=k_0_M/(1.0+k_0_M/k_inf);

       k_0_M = 1.40E-8*exp(-29500.0/tl[i])*MM[i]*THREEBODY;          
       k_inf = 3.13E14*exp(-29587.0/tl[i]);      
       kkT[i][42]=k_0_M/(1.0+k_0_M/k_inf);

     //   k_0_M = 2.58e7*pow(tl[i], -4.53)*exp(-24750.0/tl[i])*MM[i]*THREEBODY;
     //   k_inf = 1.70e10*pow(tl[i],0.8)*exp(-23620/tl[i]);
     //   kkT[i][43]=k_0_M/(1.0+k_0_M/k_inf);

     // The reaction we are reversing: H + SO2 + M --> HOSO + M
          n_prod=1;
          n_react=2;
          int mu_6[3]={-1, -1, 1};

          //Thermo coefficients for [H, SO2, HOSO]
          double a1_6[3] = {0.25000000E+01, 3.67480752, 3.73732792};
          double a2_6[3] = {0.00000000E+00, 0.00228302107, 0.00930978241};
          double a3_6[3] = {0.00000000E+00, 0.0000084689304, -0.00000385404197};
          double a4_6[3] = {0.00000000E+00, -0.0000000136562039, -0.00000000407782859};
          double a5_6[3] = {0.00000000E+00, 0.00000000000576271873, 0.00000000000309632341};
          double a6_6[3] = {0.25473660E+05, -36945.5073, -29826.0746};
          double a7_6[3] = {-0.44668285E+00, 7.9686643, 10.1080861};
     
          // forward reaction rate from Hughes+2002: H + SO2 + M --> HOSO + M
          k_f_0_M=7.34e-10*pow(tl[i],-6.43)*exp(-5577./tl[i])*THREEBODY*MM[i];
          k_f_inf=5.18e-16*pow(tl[i], 1.61)*exp(-3606./tl[i]);
          k_f=k_0_M/(1.0+k_0_M/k_inf);
     
          kkT[i][43]=ReverseRate(mu_6, n_prod, n_react, a1_6, a2_6, a3_6, a4_6, a5_6, a6_6, a7_6, tl[i], k_f);

       k_0_M = 1.92e22*pow(tl[i], -9.02)*exp(-26647.0/tl[i])*MM[i]*THREEBODY;
       k_inf = 9.94e21*pow(tl[i], -2.54)*exp(-38190.0/tl[i]);
       kkT[i][44]=k_0_M/(1.0+k_0_M/k_inf);
       
       kkT[i][45]=0.00;//ZEROED BC Wrong reaction. This refers to thermal decay of HSOO, not HSO2. Zeroed to prevent accidental use; should be removed entirely eventually. //2.73E-8*pow(tl[i]/298.0, -2.82)*exp(3750.0/tl[i])*MM[i];
       
       // Compute k46 CHO + M --> CO + H + M using forward reaction rate (CO + H + M --> CHO + M) and equilibrium constant
       // Number of products, reactants, and total species
       n_prod = 1;
       n_react = 2;
       
       // Stoichiometric coefficients: negative for reactants, positive for products
       // Order: [CO, H, CHO]
       int mu[3] = {-1, -1, 1};
       
       // NASA polynomial thermodynamic coefficients (200–1000K) from Burcat&Ruscic+2005
       // Each array corresponds to a coefficient for [CO, H, CHO]
       double a1_5[3] = {3.5795335e+00, 2.5000000e+00, 4.2375461e+00};
       double a2_5[3] = {-6.1035369e-04, 0.0000000e+00, -3.32075257e-03};
       double a3_5[3] = {1.0168143e-06, 0.0000000e+00, 1.40030264e-05};
       double a4_5[3] = {9.0700586e-10, 0.0000000e+00, -1.34239995e-08};
       double a5_5[3] = {-9.0442449e-13, 0.0000000e+00, 4.37416208e-12};
       double a6_5[3] = {-1.4344086e+05, 2.5473660e+05, 3.87241185e+03};
       double a7_5[3] = {3.5084093e+00, -4.4668285e-01, 3.30834869e+00};
       
       // forward rate from Wagner and Bowman (1987) for low Pressure limit and Arai+1981 for high Pressure limit
     //   k_f_0_M = 1.40e-34*exp(-100/tl[i])*MM[i]*THREEBODY;
     //   k_f_inf = 1.96e-13*exp(-1370/tl[i]);
     //   k_f = k_f_0_M/(1.0+k_f_0_M/k_f_inf);

     //   kkT[i][46]=ReverseRate(mu,n_prod,n_react,a1,a2,a3,a4,a5,a6,a7,tl[i],k_f);

       k_0_M = 1.23E-2*pow(tl[i],-2.36)*exp(-9755/tl[i])*MM[i]*THREEBODY;
       k_inf = 4.93e+16*pow(tl[i],-0.93)*exp(-9927/tl[i]);
       F_c=0.897*exp(-tl[i]/139)+0.103*exp(-tl[i]/1.09E4)+exp(-4.55E3/tl[i]);
       M_M_c=k_0_M/k_inf;
       N=0.75-1.27*log10(F_c);
       F=pow(F_c, 1.0/(1.0+pow((log10(M_M_c)/N),2.0)));
       kkT[i][46]=k_0_M*k_inf/(k_0_M+k_inf)*F;
       
       k_0_M = 9.995E-8*pow(tl[i], -0.55)*exp(-9063.94/tl[i])*MM[i]*THREEBODY;
       k_inf = 1.13E10*pow(tl[i], 1.21)*exp(-12111.8944/tl[i]);
       F_c=0.659*exp(-tl[i]/28)+0.341*exp(-tl[i]/1000)+exp(-2339/tl[i]);
       M_M_c=k_0_M/k_inf;
       N=0.75-1.27*log10(F_c);
       F=pow(F_c, 1.0/(1.0+pow((log10(M_M_c)/N),2.0)));
       kkT[i][47] = k_0_M*k_inf/(k_0_M+k_inf)*F;
       
     //   k_0_M = 4.92e-29*pow(tl[i]/298.0,-2.4)*exp(-18862.0/tl[i])*MM[i]*THREEBODY;
     //   k_inf = 2.985e1*pow(tl[i]/298.0,0.13)*exp(-18349.0/tl[i]);
     //   kkT[i][48]= k_0_M / (1.0 + k_0_M / k_inf);

     k_0_M = pow(10,25.137)/6.022e23*pow(tl[i],-2.4)*exp(-18862.0/tl[i])*MM[i]*THREEBODY;
     k_inf = pow(10,14.074)*pow(tl[i], 0.132)*exp(-18349.0/tl[i]);
     M_M_c=k_0_M/k_inf;
     F_c = 0.729*exp(-513/tl[i])+exp(-tl[i]/540);
     N=0.75-1.27*log10(F_c);
     F=pow(F_c, 1.0/(1.0+pow((log10(M_M_c)/N),2.0)));
     kkT[i][48]=k_0_M*k_inf/(k_0_M+k_inf)*F;
       
     //   k_0_M = 7.21e-31*pow(tl[i]/298.0,-3.15)*exp(-18629.0/tl[i])*MM[i]*THREEBODY;
     //   k_inf = 1.25e2*pow(tl[i]/298.0,0.41)*exp(-17783.0/tl[i]);
     //   kkT[i][49] = k_0_M / (1.0 + k_0_M / k_inf);

     k_0_M = pow(10,26.775)/6.022e23*pow(tl[i],-3.15)*exp(-18629.0/tl[i])*MM[i]*THREEBODY;
     k_inf = pow(10, 11.915)*pow(tl[i], 0.413)*exp(-17783.0/tl[i]);
     M_M_c=k_0_M/k_inf;
     F_c = 1.049*exp(-2407/tl[i])+exp(-tl[i]/823);
     N=0.75-1.27*log10(F_c);
     F=pow(F_c, 1.0/(1.0+pow((log10(M_M_c)/N),2.0)));
     kkT[i][49]=k_0_M*k_inf/(k_0_M+k_inf)*F;

       k_0_M = 1.32E-6*exp(-17199.0/tl[i])*MM[i]*THREEBODY;
       k_inf = 1.06e15*exp(-18041.0/tl[i]);
       kkT[i][50]=k_0_M / (1.0 + k_0_M / k_inf);

       k_0_M = 1.18e-11*MM[i]*THREEBODY;
       k_inf = 3.98e13*exp(-19364.0/tl[i]);
       kkT[i][51]=k_0_M / (1.0 + k_0_M / k_inf);

       k_0_M = 4.14e-7 * (tl[i]/300) * exp(-16800 / tl[i])*MM[i]*THREEBODY;
       k_inf = 1.0e13 * exp(-16800 / tl[i]);
       kkT[i][52]=k_0_M / (1.0 + k_0_M / k_inf);

       k_0_M=3.16E-4*pow(tl[i]/298.0, -4.72)*exp(-13591.0/tl[i])*MM[i]*THREEBODY;
       k_inf = 4.34e14*exp(-12999/tl[i]);
       kkT[i][53]=k_0_M / (1.0 + k_0_M / k_inf);

          k_0_M = 1.0E-12 * MM[i] *THREEBODY + 1.5E3;
          k_inf = 1.0E10;
          kkT[i][54] = k_0_M / (1.0 + k_0_M / k_inf);

          k_0_M = 1.0E-11 * MM[i] *THREEBODY + 2.2E4;
          k_inf = 1.0E10;
          kkT[i][55] = k_0_M / (1.0 + k_0_M / k_inf);

          k_0_M = 1.5E-13 * MM[i] *THREEBODY + 1.13E3;
          k_inf = 1.0E10;
          kkT[i][56] = k_0_M / (1.0 + k_0_M / k_inf);
	
		 k_0_M=1.93E-4*pow(tl[i]/298.0,-2.44)*exp(-62782.1/tl[i])*MM[i]*THREEBODY;
           k_inf = 4.15E15*pow(tl[i]/298.0,-0.93)*exp(-62294/tl[i]);
           F_c = 0.875 - 0.5e-4*tl[i];
           kkT[i][57] = k_0_M / ((1.0 + k_0_M / k_inf)*F_c);

          // Compute k58 H2S + M --> H2 + S + M using forward reaction rate (H2 + S + M --> H2S + M) and equilibrium constant
          // Number of products, reactants, and total species
          n_prod = 1;
          n_react = 2;
          
          // Stoichiometric coefficients: negative for reactants, positive for products
          // Order: [H2, S, H2S]
          int mu_2[3] = {-1, -1, 1};
          
          // NASA polynomial thermodynamic coefficients (200–1000K) from Burcat&Ruscic+2005
          // Each array corresponds to a coefficient for [H2, S, H2S]
          double a1_2[3] = {3.50207268, 2.31725616, 4.12024455};
          double a2_2[3] = {0.0000865475654, 0.00478018342, -0.00187907426};
          double a3_2[3] = {-0.000000263683344, -0.0000142082674, 0.0000082142665};
          double a4_2[3] = {0.000000000337306621, 0.0000000156569538, -0.0000000070642573};
          double a5_2[3] = {-0.0000000000000292359706, -0.00000000000596588299, 0.0000000000021423486};
          double a6_2[3] = {-1046.31279, 32506.8976, -3682.15173};
          double a7_2[3] = {-4.25875759, 6.06242434, 1.53174068};
          
          // forward rate from Zahnle+2016
          k_f_0_M = 1.4e-31 * pow((tl[i]/298),-1.9) * exp(-8140 / tl[i]) * MM[i]*THREEBODY;
          k_f_inf = 1.0e-11;
          k_f = k_f_0_M / (1.0 + k_f_0_M / k_f_inf);

          kkT[i][58]=ReverseRate(mu_2,n_prod,n_react,a1_2,a2_2,a3_2,a4_2,a5_2,a6_2,a7_2,tl[i],k_f);

		// Compute k59 H2S + M --> HS + H + M using forward reaction rate (HS + H + M --> H2S + M) and equilibrium constant
          // Number of products, reactants, and total species
          n_prod = 1;
          n_react = 2;
          
          // Stoichiometric coefficients: negative for reactants, positive for products
          // Order: [HS, H, H2S]
          int mu_3[3] = {-1, -1, 1};
          
          // NASA polynomial thermodynamic coefficients (200–1000K) from Burcat&Ruscic+2005
          // Each array corresponds to a coefficient for [HS, H, H2S]
          double a1_3[3] = {3.68466877,  2.5, 4.12024455};
          double a2_3[3] = {0.00324608824, 0, -0.00187907426};
          double a3_3[3] = {-0.0000128635079, 0, 0.0000082142665};
          double a4_3[3] = {0.0000000169512196, 0, -0.0000000070642573};
          double a5_3[3] = {-0.00000000000707595387, 0, 0.0000000000021423486};
          double a6_3[3] = {15903.6477, 25473.66, -3682.15173};
          double a7_3[3] = {2.01781634, -0.44668285, 1.53174068};
          
          // forward rate from Zahnle+2016
          k_f_0_M = 1.4e-31 * pow((tl[i]/298),-2.5) * exp(500 / tl[i]) * MM[i]*THREEBODY;
          k_f_inf = 1.0e-10;
          k_f = k_f_0_M / (1.0 + k_f_0_M / k_f_inf);

          kkT[i][59]=ReverseRate(mu_3,n_prod,n_react,a1_3,a2_3,a3_3,a4_3,a5_3,a6_3,a7_3,tl[i],k_f);

		 kkT[i][60]=4.07E-10*exp(-30910.0/tl[i])*MM[i];
		 kkT[i][61]=3.98E+12*exp(-18402.0/tl[i]);
		 kkT[i][62]=1.40E+13*exp(-30188.3/tl[i]);
		 kkT[i][63]=1.26E+13*exp(-16838.0/tl[i]);
		 kkT[i][64]=1.18E+18*pow(tl[i]/298.0, -1.2)*exp(-49191.0/tl[i]);
		 kkT[i][65]=2.50E+15*exp(-43659.0/tl[i]);
		 kkT[i][66]=4.00E+13*exp(-40291.0/tl[i]);
		 kkT[i][67]=1.09E+13*pow(tl[i]/298.0, 0.17)*exp(-17921.6/tl[i]);
		 kkT[i][68]=1.31E+13*pow(tl[i]/298.0, 0.87)*exp(-15274.6/tl[i]);
		 kkT[i][69]=1.58E+16*exp(-49071.0/tl[i]);
		 kkT[i][70]=1.07E+12*pow(tl[i]/298.0, -15.74)*exp(-49672.0/tl[i])*MM[i];
		 kkT[i][71]=3.38E+10;
		 kkT[i][72]=1.66E-8*exp(-29707.0/tl[i])*MM[i];
		 kkT[i][73]=6.3E+13*exp(-43779.0/tl[i]);
		 kkT[i][74]=1.15E+15*exp(-43539.0/tl[i]);
		 kkT[i][75]=1.6E+14*exp(-20807.0/tl[i]);
		 kkT[i][76]=3.57E-8*pow(tl[i]/298.0, 0.7)*exp(-21288.0/tl[i])*MM[i];
		 kkT[i][77]=7.7E+14*exp(-44260.0/tl[i]);
		 kkT[i][78]=5.0E+15*exp(-38126.0/tl[i]);
		 kkT[i][79]=2.0E+15*exp(-43178.0/tl[i]);
		 kkT[i][80]=2.0E+15*exp(-37645.0/tl[i]);
		 kkT[i][81]=1.58E+16*exp(-55445.0/tl[i]);
		 kkT[i][82]=1.8E+13*exp(-42817.0/tl[i]);
		 kkT[i][83]=1.0E+14*exp(-37645.0/tl[i]);
		 kkT[i][84]=1.0E+16*exp(-36683.0/tl[i]);
		 kkT[i][85]=1.1E+13*pow(tl[i]/298.0, 0.25)*exp(-17921.6/tl[i]);
		 kkT[i][86]=7.7E+13*pow(tl[i]/298.0, 0.77)*exp(-15394.8/tl[i]);
		 kkT[i][87]=1.58E+16*exp(-49312.0/tl[i]);
		 kkT[i][88]=1.0E+17*exp(-42576.0/tl[i]);
		 kkT[i][89]=7.94E+16*exp(-40411.0/tl[i]);
		 kkT[i][90]=6.31E+14*exp(-31151.0/tl[i]);
		 
		 
		 k_0_M = 6.61E-9*exp(-20634.0/tl[i])*MM[i]*THREEBODY;          
		 k_inf = 7.94E13*exp(-27680.0/tl[i]);      
		 kkT[i][91]=k_0_M/(1.0+k_0_M/k_inf);
		 
		 k_0_M = 2.66E-11*exp(-24295.0/tl[i])*MM[i]*THREEBODY;          
		 k_inf = 6.43E8*pow(tl[i]/300.0, -1.00)*exp(-24300.0/tl[i]);     
		 kkT[i][92]=k_0_M/(1.0+k_0_M/k_inf);
		 
		// Compute T93 CH5N + M --> CH3 + NH2 + M using forward reaction rate (M86 CH3 + NH2 + M --> CH5N + M) and equilibrium constant
          // Number of products, reactants, and total species
          n_prod = 1;
          n_react = 2;
          
          // Stoichiometric coefficients: negative for reactants, positive for products
          // Order: [CH3, NH2, CH5N]
          int mu_4[3] = {-1, -1, 1};
          
          // NASA polynomial thermodynamic coefficients (200–1000K) from Burcat&Ruscic+2005
          // Each array corresponds to a coefficient for [CH3, NH2, CH5N]
          double a1_4[3] = {0.36571797E+01,  4.19198016E+00, 4.93595327E+00};
          double a2_4[3] = {0.21265979E-02, -2.04602827E-03, -1.06687240E-02};
          double a3_4[3] = {0.54583883E-05, 6.67756134E-06, 6.66595644E-05};
          double a4_4[3] = {-0.66181003E-08, -5.24907235E-09, -7.68165338E-08};
          double a5_4[3] = {0.24657074E-11, 1.55589948E-12, 2.88891949E-11};
          double a6_4[3] = {0.16422716E+05, 2.11864310E+04, -3.96311166E+03};
          double a7_4[3] = {0.16735354E+01, -9.04785244E-02, 1.01955189E+00};
          
          // forward rate from M86
          k_f_0_M = 1.8E-27*pow(tl[i]/298.0,-3.85)*MM[i]*THREEBODY;
          k_f_inf = 1.30E-10*pow(tl[i]/298.0,0.42);
          k_f = k_f_0_M/(1.0+k_f_0_M/k_f_inf);
		 kkT[i][93]=ReverseRate(mu_4,n_prod,n_react,a1_4,a2_4,a3_4,a4_4,a5_4,a6_4,a7_4,tl[i],k_f);
     }
} 
