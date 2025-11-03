% read cross sections
% all in SI unless otherwised noted

clear all
close all

% name
NAME='H2O';

% T, P
T=500;
P=1;
P1=1E+5;


% fopen
f=fopen(['opac',NAME,'.dat'],'r');

% temp
nt=18;
for i=1:nt
    t(i)=fscanf(f,'%f',1);
end

% pres
np=13;
for i=1:np
    p(i)=fscanf(f,'%e',1);
end

% cross
for k=1:16000
    w(k)=fscanf(f,'%e',1);
    for j=1:np
        dum=fscanf(f,'%e',1);
        for i=1:nt
            c(j,i)=fscanf(f,'%e',1);
        end
    end
    s(k)=interp2(t,p,c,T,P)*1.0E+4;
    s1(k)=interp2(t,p,c,T,P1)*1.0E+4;
end

% fclose
fclose(f);

% Planck Mean

% Physical constant SI
HPLANCK = 6.626068E-34;
CLIGHT = 299792458;
KB = 1.3806503E-23;
SIGMA = 5.670373E-8;

h = 2*HPLANCK*CLIGHT^2./w.^5./(exp(HPLANCK*CLIGHT./w/KB/T)-1);
sm = trapz(w,h.*s)/trapz(w,h);
sm1 = trapz(w,h.*s1)/trapz(w,h);

% ultrafine wavelength resolution
l = linspace(w(1),w(end),2E+7);
h = 2*HPLANCK*CLIGHT^2./l.^5./(exp(HPLANCK*CLIGHT./l/KB/T)-1);
sl= exp(interp1(w,log(s+1e-50),l));
sl1= exp(interp1(w,log(s1+1e-50),l));
slm = trapz(l,h.*sl)/trapz(l,h);
slm1 = trapz(l,h.*sl1)/trapz(l,h);

% Interpolation and Integration
h = 2*HPLANCK*CLIGHT^2./w.^5./(exp(HPLANCK*CLIGHT./w/KB/T)-1);
dw=diff(w);
ds=diff(s.*h);
ds1=diff(s1.*h);
dslog=diff(log(s.*h+1e-80.*w));
dslog1=diff(log(s1.*h+1e-80.*w));
se=sum(dw.*ds./dslog)/trapz(w,h);
se1=sum(dw.*ds1./dslog1)/trapz(w,h);