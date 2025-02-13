% Hua-sheng XIE, huashengxie@gmail.com, ENN, 2022-06-08 08:02
% update 25-02-13 08:12
% function fD and fDtensor are traditional sum_n method, slow and bad
% convergent for large kperp
% fDR is fast method from Ronnmark82 WHAMP code
close all; 
clear; clc;

global S c2 wcs wps2 rhocs w kz vts;
% global S c2 wcs wps2 w kz vTzs;

% 1. default SI unit
c2=(2.99792458E8)^2; % speed of ligth c^2
epsilon0=8.854187817E-12;
mu0=1/(c2*epsilon0);
kB=1.38064852e-23;
qe=1.60217662e-19; % electron charge, coulombs
mp=1.6726219e-27; % proton mass, kg
me=9.1094e-31; % electron mass, kg

% 2. set parameters, modify here for your own case
jcase=1;
if(jcase==1)
    ms=[me,2*mp]; % speceis mass
    qs=[-1,1]*qe; % speceis charge
    ns_unit=[1,1]; % *den0, species density, m^-3
    Ts_unit=[1,1]; % eV
    S=length(ms);
    
    B0=6e0; % background B field, Tesla
    n00=2e20; % m^-3
    % f=28e9; % Hz
    
    T00=2*2e3; % eV
    % w=f*2*pi; % omega
    
    kz=2*5.0;
elseif(jcase==2)
    ms=me; % speceis mass
    qs=-1*qe; % speceis charge
    ns_unit=1; % *den0, species density, m^-3
    Ts_unit=1; % eV
    S=length(ms);

    B0=1e0; % background B field, Tesla
    n00=6.075e19; % m^-3
    % f=28e9; % Hz

    T00=1e3; % eV
    kz=0.0;

elseif(jcase==3)
    ms=[me,2*mp,3*mp]; % speceis mass
    qs=[-1,1,1]*qe; % speceis charge
    ns_unit=[1,0.1,0.9]; % *den0, species density, m^-3
    Ts_unit=[1,1,1]; % eV
    S=length(ms);
    
    B0=6e0; % background B field, Tesla
    n00=2e20; % m^-3
    % f=28e9; % Hz
    
    T00=2e3; % eV
    kz=1.0;
end

ns0=ns_unit*n00;
Ts0=Ts_unit*T00*1.6022e-19/kB; % Ts, eV -> K
vts=sqrt(2*kB*Ts0./ms); % thermal velocity
c=sqrt(c2);

wcs=qs*B0./ms; % the gyro frequency
wps=sqrt(ns0.*qs.^2./ms/epsilon0); % plasma frequency
wcs2=wcs.^2;
wps2=wps.^2;
wp2=sum(wps2);
rhocs=sqrt(kB*Ts0./ms)./wcs; % cyclotron radius

wcn=min(abs(wcs));
rhocn=max(abs(rhocs));
kn=1/rhocn;
% kn=1/abs(rhocs(1));
wn=wcn;

% kz=0.1*kn; % 25-02-02 01:54

%% 3. Test fD(kx) to solve kx with fixed kz and w
xx=(0.11:0.05:3.95).*wcn+0*wcn;
nx=length(xx);
fDr=zeros(nx,1);
runtime=cputime;
for jx=1:nx
    w=xx(jx);
    if(jx==1)
        % x01=10.0*kn;
        % x01=1e3+1i*1e3;
        x01=(0.6+0.0i)*kn;
        % x01=2.8e3;
    % else
    %     x01=kx1;
    end
    
    % x01=kxx(jx,5);
%     warning off;
    options = optimset('Display','off','TolFun',1e-16,'Tolx',1e-16);
%     options = optimset('Display','iter','TolFun',1e-100);
    % kx1=fsolve(@fD,x01,options)
    kx1=fsolve(@fDR,x01,options)
    nf(jx,1)=kx1;
end
runtime=cputime-runtime;
%% plot

if(jcase==1)

close all;
h=figure('unit','normalized','Position',[0.01 0.05 0.5 0.7],...
    'DefaultAxesLineWidth',2,...
    'DefaultAxesFontSize',14);
% subplot(2,2,[1,3]);

ax1=axes('position',[0.1,0.1,0.4,0.8]);
semilogy(xx/wcn,real(nf(:,1)),'.',...
    'linewidth',3,'Markersize',12);
hold on;
semilogy([1,1],[1e-2,1e4],'k:',[2,2],[1e-2,1e4],'k:',...
    [3,3],[1e-2,1e4],'k:',[4,4],[1e-2,1e4],'k:','linewidth',1);
ylim([1e0,1e4]);

xlabel('\omega/\omega_{ci}');
ylabel('Re(k_\perp) [m^{-1}]');
title(['S=',num2str(S),', B_0=',num2str(B0,3),'T, k_z=',num2str(kz,3),...
    'm^{-1}, n_{00}=',num2str(n00,3),'m^{-3}, T_{00}=',num2str(T00,3),'eV,',10, ...
    'n_s=[',num2str(ns_unit,3),'], m_s/m_e=[',num2str(ms/me,4),...
    '], q_s/q_e=[',num2str(qs/qe,3),...
    '], T_s/T_{00}=[',num2str(Ts_unit,3),']']);

% subplot(224);
ax3=axes('position',[0.62,0.1,0.35,0.3]);
semilogy(xx/wcn,abs(imag(nf(:,1))),'b.',...
    'linewidth',3,'Markersize',12);hold on;
% semilogy([1,1],[1e-2,1e4],'k:',[2,2],[1e-2,1e4],'k:',...
%     [3,3],[1e-2,1e4],'k:',[4,4],[1e-2,1e4],'k:','linewidth',1);
% ylim([1,1e3]);
xlabel('\omega/\omega_{cn}');ylabel('Im(k_\perp) [m]');
% 
% ax3=axes('position',[0.62,0.1,0.35,0.3]);
% semilogy(xx/wci,-1./imag(nf(:,2)),'r.',...
%     'linewidth',3,'Markersize',12);hold on;

title(['nx=',num2str(nx),', runtime=',num2str(runtime,3),'s']);
% semilogy([1,1],[1e-2,1e4],'k:',[2,2],[1e-2,1e4],'k:',...
%     [3,3],[1e-2,1e4],'k:',[4,4],[1e-2,1e4],'k:','linewidth',1);
% ylim([-1e1,1e1]);
% ylim([0.1,1e2]);
% xlabel('\omega/\omega_{ci}');ylabel('-Im(k_\perp)^{-1} [m]');

% save figure
% set(gcf,'Units','inches');
% screenposition = get(gcf,'Position');
% set(gcf,'PaperPosition',[0 0 screenposition(3:4)],...
%   'PaperSize',[screenposition(3:4)]);

print(gcf,'-dpng',['bonkR_S=',num2str(S),',B0=',num2str(B0,3),...
    ',n00=',num2str(n00,3),',T00=',num2str(T00,3),',ns=',num2str(ns_unit,3),...
    ',ms=',num2str(ms/me,4),',qs=',num2str(qs/qe,3),...
    ',Ts=',num2str(Ts_unit,3),',k_z=',num2str(kz,3),'.png']);

else

subplot(121);
plot(xx/wcn,abs(real(nf(:,1))),'r:x');
set(gca,'yscale','linear');
subplot(122);
plot(xx/wcn,abs(imag(nf(:,1))),'r:x');
set(gca,'yscale','log');
% ylim([0,5]);xlim([0,5]);
% xlabel('k_\perp');ylabel('\omega_r');
end
%% calculate fDtensor
% jx=10;
% f=xx(jx)/(2*pi); kx=nf(jx,1);
% [DD,Kxx,Kxy,Kxz,Kyx,Kyy,Kyz,Kzx,Kzy,Kzz]=fDtensor(...
%     B0,qs,ms,ns0,Ts0,f,kz,kx)

