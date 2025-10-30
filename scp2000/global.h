/*----------------------- global.h --------------------------------

Author: Renyu Hu (hury@mit.edu)
Last modified: July 20, 2011

--------------------------------------------------------------------- */

#ifndef __GLOBAL_H__
#define __GLOBAL_H__

/*---- External Variables ------------------------------------------- */

extern double thickl;
extern double zl[];
extern double pl[];
extern double tl[];
extern double MM[], MMZ[];
extern double wavelength[];
extern double solar[];
extern double crossr[], crossa[3][WaveBin+1], sinab[3][WaveBin+1], asym[3][WaveBin+1];
extern double **opacCO2, **opacO2, **opacSO2, **opacH2O, **opacOH, **opacH2CO; 
extern double **opacH2O2, **opacHO2, **opacH2S, **opacCO, **opacO3, **opacCH4; 
extern double **opacNH3;
extern double MeanCO2[], MeanO2[], MeanSO2[], MeanH2O[], MeanOH[], MeanH2CO[];
extern double MeanH2O2[], MeanHO2[], MeanH2S[], MeanCO[], MeanO3[], MeanCH4[];
extern double MeanNH3[];
extern double rainoutrate[zbin+1][NSP+1];
extern double Vesc[], VFall[];
extern double nsH2O[], nsH2SO4[], nsS8[], tcondfH2O[], tcondfH2SO4[], tcondfS8[];
extern double kk[zbin+1][NKin+1], kkM[zbin+1][NKinM+1], kkT[zbin+1][NKinT+1];
extern double Rkk[zbin+1][NKin+1], RkkM[zbin+1][NKinM+1], RkkT[zbin+1][NKinT+1];
extern int    ReactionR[NKin+1][7], ReactionM[NKinM+1][5], ReactionP[NPho+1][9], ReactionT[NKinT+1][4];
extern int    numr, numm, numt, nump, numx, numc, numf, numa, waternum, waterx;
extern double **DM, **dl, KE[];
extern double xx[zbin+1][NSP+1];
extern double mkv[], Tnew[], Pnew[];
extern double GibbsForm[NSP+1][zbin+1];

#endif /* !__GLOBAL_H__ */

/*---- end ------------------------ global.h ---------------------- */
