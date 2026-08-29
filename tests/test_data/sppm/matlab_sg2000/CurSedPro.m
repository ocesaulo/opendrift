%   program CurSedPro.m
clear;

%         Last udate:  4-August-2017
%         This software at its present stage of development is not intended
%         for commercialization.  This software and any copies or derivatives 
%         is intended to be used for research and evaluation purposes only
%         and is provided as is WITHOUT ANY WARRANTY.
%         WARRANTIES OF MERCHANTABILITY AND OF FITNESS FOR A PARTICULAR
%         PURPOSE ARE EXPRESSLY DISCLAIMED.
%         The authors shall not be liable for any loss or damages arising
%         from any use, defect, omission, failure or the like of said 
%         software, nor shall they have any obligation to make available any
%         corrections, improvements, or other modifications or to provide 
%         any assistance or service of any kind.


% computes current and concentration profiles for data output from bblm02.m

% format for BBLMPRMS: Ro mu epsilon z1ozn z2ozn zroz1 zroz2 fofx kbs kbr ub ab ur znot
% format for SKNPRMS: Psi fwskn arg_ole

% gamma    - ratio of neutral eddy viscosity to neutral sediment diffusivity
% gammanot - resuspension coefficient
% kappa    - von Karman's constant
% cbed     - bed sediment concentration (NOTE default is to set it equal to 0.65
%            even in the presence of multiple size classes)
% z1000    - maximum height to compute profiles in the constant stress layer


% This script calls the function: neut_sed_fun.m

% OUTPUT
% SED    - non-structured array containing the concentration profiles for all size classes
% CUR    - non-structured array containing current profiles
% STRESS - non-structured array containing shear stress and other parameters
% SEDPRO - non-structured array containing the sediment transport for all size classes
% ZZ     - non-structured array containing the vertical coordinate for the above profiles

% Reference:  Styles, R., S. Glenn, and M. Brown 2017, "An optimized combined
% wave and current bottom boundary layer model for arbitrary bed roughness",
% U. S. Army Corps of Engineers, ERDC/CHL TR-17-11, 36 pp.


load model_output_file2017.mat;
kappa=0.4;
gamma=0.74;

% note bed concentration must be specified.
cbed=(d-d)+0.65;
gammanot=1.3e-3;
z1000=1000;

[m, n]=size(BBLMPRMS);
PSI=SKNPRMS(:,1);
for i=1:m
  Ro=BBLMPRMS(i,1);
  mu=BBLMPRMS(i,2);
  eps=BBLMPRMS(i,3);
  z1ozn=BBLMPRMS(i,4);
  z2ozn=BBLMPRMS(i,5);
  ub=BBLMPRMS(i,12);
  ab=BBLMPRMS(i,13);
  znot=BBLMPRMS(i,11);
  omega=ub./ab;
  ustarcw=Ro.*znot.*omega;
  ustarc=eps*ustarcw;
  ustarwm=mu*ustarcw;
  z1=z1ozn*znot;
  z2=z2ozn*znot;
% compute reference concentration for each size class
  [mm, nn]=size(d);
  a=0;
  b=2*pi./omega;
  dt=360;
  t=linspace(a,b,dt);
  z7d=7*d_median;
  for j=1:nn
    R=-gamma*wf(j)/kappa/ustarcw;
    psiopsicr=PSI(i)/psicr;
    F1=psiopsicr*sin(t*omega)-1;
    F2=F1;  Igt=find(t > mean(t));
    F2(Igt)=-F2(Igt)-2;
    Ilt0=find(F2 < 0);
    F2(Ilt0)=0;
    Sp(i)=trapz(t,F2)/b;
    C7d(i,j)=cbed(j)*gammanot*Sp(i);
    Cref(i,j)=C7d(i,j)*(znot/z7d)^R;
  end
% compute concentration profiles.
  cref=Cref(i,:);
  [SEDPRO,CURPRO,SEDTRANS,DTRNS,Z]=neutsed_fun(ustarcw,ustarc,znot,z1,z1000,d,cref,wf,kappa,gamma);
  SED(i).nsa=SEDPRO; CUR(i).nsa=CURPRO; ZZ(i).nsa=Z'; SCTRANS(i).nsa=SEDTRANS;
  DTRANS(i).nsa=DTRNS; STRESS(i).ustarcw=ustarcw; STRESS(i).ustarc=ustarc; STRESS(i).ustarwm=ustarwm;
  STRESS(i).znot=znot; STRESS(i).z1=z1; STRESS(i).z2=z2;
end

