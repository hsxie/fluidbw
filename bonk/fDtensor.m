% Hua-sheng XIE, huashengxie@gmail.com, ENN, 2021-08-31 17:03
% Hot magnetized plasma dispersion relation, only Maxwellian.
% Ref: Xie, H.; Ma, H. & Bai, Y., Plasma Waves Accessibility Diagrams: A 
%   Tutorial to Include the Fluid and Kinetic Thermal Effects, arXiv, 2021.
% Note the sqrt(2) of notations have been updated.
% 21-11-12 19:37 fixed a bug on Dxz & Dzx
% 22-06-06 14:19 update for YANG Hua (ASIPP) ICRF full wave
% after solve the kx, we can use this function to obtain all 3*3 tensor

% 2021-08-28 17:00 modified from previous (2013-06-17) pdrko.m
function [DD,Kxx,Kxy,Kxz,Kyx,Kyy,Kyz,Kzx,Kzy,Kzz]=fDtensor(...
    B0,qs,ms,ns0,Ts0,f,kz,kx) % DR, only Maxwellian
% all are SI unit, note Ts0 should also transform from eV to Kelvin
% note the tensor relation: K_ij=I-sigma_ij/(i*w*epsion0)
% note also lmx is the sum range of the Bessel function, increase it to
% check the convergence

w=2*pi*f;
S=length(qs);

% 1. default SI unit
c2=(2.99792458E8)^2; % speed of ligth c^2
epsilon0=8.854187817E-12;
% mu0=1/(c2*epsilon0);
kB=1.38064852e-23;
% qe=1.60217662e-19; % electron charge, coulombs
% mp=1.6726219e-27; % proton mass, kg
% me=9.1094e-31; % electron mass, kg

% Ts0K=Ts0*1.6022e-19/kB; % Ts, eV -> K
vts=sqrt(2*kB*Ts0./ms); % thermal velocity
% c=sqrt(c2);

wcs=qs*B0./ms; % the gyro frequency
wps=sqrt(ns0.*qs.^2./ms/epsilon0); % plasma frequency
wps2=wps.^2;
rhocs=sqrt(kB*Ts0./ms)./wcs; % cyclotron radius

nx=sqrt(c2)*kx/w; nz=sqrt(c2)*kz/w;
as=kx*rhocs;
bs=as.^2;

lmx=15; % change this to check convergence
%     lmx=200;
Kxx=1; Kxy=0; Kxz=0;
Kyx=0; Kyy=1; Kyz=0;
Kzx=0; Kzy=0; Kzz=1;
for j=1:S
    z0=w/(kz*vts(j));
    b=bs(j);
    a=as(j);
    for l=-lmx:lmx
        zn=(w-l*wcs(j))/(kz*vts(j));
        Xsn_xx=l^2*Gamman(l,b)/b;
        Xsn_xy=1i*l*Gammapn(l,b);
        Xsn_xz=sqrt(2)*zn*l*Gamman(l,b)/a;
        Xsn_yx=-Xsn_xy;
        Xsn_yy=l^2/b*Gamman(l,b)-2*b*Gammapn(l,b);
        Xsn_yz=-1i*sqrt(2)*zn*a*Gammapn(l,b);
        Xsn_zx=Xsn_xz;
        Xsn_zy=-Xsn_yz;
        Xsn_zz=2*zn^2*Gamman(l,b);
        
        Zn=zfun(zn);
        term1=wps2(j)/w^2*(z0*Zn);
        Kxx=Kxx+term1*Xsn_xx;
        Kxy=Kxy+term1*Xsn_xy;
        Kxz=Kxz+term1*Xsn_xz;
        Kyx=Kyx+term1*Xsn_yx;
        Kyy=Kyy+term1*Xsn_yy;
        Kyz=Kyz+term1*Xsn_yz;
        Kzx=Kzx+term1*Xsn_zx;
        Kzy=Kzy+term1*Xsn_zy;
        Kzz=Kzz+term1*Xsn_zz;
    end
    % Kzz=Kzz+2*z0^2; % 21-09-02 21:26 wrong?
    Kzz=Kzz+wps2(j)/w^2*2*z0^2;
end

Dxx=Kxx-nz^2; Dyx=Kyx; Dxy=Kxy;
% Dxz=Kxz-nx*nz; Dyy=Kyy-(nx^2+nz^2); Dzx=Kzx-nx*nz; % wrong
Dxz=Kxz+nx*nz; Dyy=Kyy-(nx^2+nz^2); Dzx=Kzx+nx*nz; % 21-11-12 19:37
Dzy=Kzy; Dyz=Kyz; Dzz=Kzz-nx^2;

DD=Dxx*Dyy*Dzz+Dyx*Dzy*Dxz+Dzx*Dyz*Dxy-...
    Dxz*Dyy*Dzx-Dyz*Dzy*Dxx-Dzz*Dyx*Dxy;
%     f=DD;


end
