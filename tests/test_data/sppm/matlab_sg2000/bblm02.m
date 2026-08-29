%  program bblm02
clear;
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

% Last update:  05-Dec-2002
% Last update:  13-Dec-2013 Note: matlab bessel functions have changed.  In this
% version 'bessel' is replaced with 'besselj'.  No more changes

% This code is designed for time series input.
% it relies on kb computed as described by Styles and Glenn
% In addition, if used in a muddy environment, (no wave generated sand ripples)
% then the roughness parameter should be modified.



% INPUT:
% Ub       - bottom wave orbital velocity  (cm/s)
% Ab       - bottom wave excursion amplitude (cm)
% Ur       - mean current at a known height above the bed zr  (cm/s)
% zr       - height above the bed the mean current Ur is measured  (cm)
% DEG      - angle between the wave and current (degrees)
% d        - sediment grain size (can be a scalar or vector) (cm)
% d_median - median (or mean) grain diameter used to compute
%            psicr and Psi (cm)
% s        - relative sediment density
% nu       - kinematic viscosity of the water (default is for seawater @ 15C) (cm^2/s)
% g        - acceleration due to gravity (cm/s^2)
% eta_def  - default value for ripple height (cm)
% lab_def  - default value for ripple wavelength (cm)
% alpha    - closure constant
% beta     - closure constant


% COMPUTED MODEL PARAMETERS
% znot     - hydraulic roughness  (cm)
% kb       - bottom roughness (=30*znot) (cm)
% ETA      - ripple height (cm)
% LAM      - ripple wavelength (cm)
% Psi      - Shields parameter based on skin friction
% psicr    - critical shear stress for initiation of sediment motion.
% mu       - ratio of the magnitude of the maximum wave shear velocity (ustarwm)
%            to the magnitude of the combined shear velocity (ustarcw).
% epsilon  - ratio of the magnitude of the time averaged shear velocity (ustarc) 
%            to the magnitude of the combined shear velocity (ustarcw).
% Ro       - internal friction Rossby number =utarcw/(omega*znot)
% sigma    - ub/ustarcw


% OUTPUT PRODUCTS
% SKNPRMS - skin friction parameters.
% BBLMPRMS - selected model parameters needed to compute shear stresses &
%            velocity profiles.
% These are saved in a mat file "model_output_file1.mat"


%addpath c:\mlab

global alpha kappa z1p mp
kappa=0.4; 
nu=0.0119; 
g=981;
mm=50;
beta=0.7;
d=[0.009 0.01 0.02 0.03 0.04];
d_median=0.04;
s=2.65;
tol=1.0e-4;
eta_def=1;
lam_def=15;
kbr_def=3;
star=d_median/(4*nu).*sqrt(g*d_median*(s-1));
psicrnorm=(s-1)*g*d_median;
den=0.011607+0.0744*d;
wf=(-3*nu+sqrt(9*nu^2+g*d.^2*(s-1).*...
   (0.003869+0.0248*d)))./den;
psinorm=(s-1)*g*d_median;
% calculate critcal shear stress for initiation of sediment motion

psicr=shldc(star);

% LOAD the data file with time series of Ub, Ab, Ur, Zr, Deg
load test2.mat;
INPUT=DATA;
Time=INPUT(:,1);
Ub=INPUT(:,2);
Ab=INPUT(:,3);
Ur=INPUT(:,4);
Zr=INPUT(:,5);
Deg=INPUT(:,6);
[m, n]=size(INPUT);

Omega=Ub./Ab;

% compute skin friction shear stress based on Ole Madsen's formula.
arg_ole=Ab/d_median;
fwcskn=exp(5.61*arg_ole.^(-0.109)-7.30);
ustarwmsknole2=0.5*fwcskn.*(1.42*Ub).^2;
Psi=ustarwmsknole2./psinorm;
SKNPRMS=[Psi fwcskn arg_ole];

% Alpha values based on work presented in Styles and Glenn, JGR, (2002)
% a value of Alpha = 0.3 and Con = 6.4 is presently recommended in
% the presence of waves over a sandy bed.
AlphaMatrix=[0.15 0.3 0.5];
AlphaCon=[9.8 6.4 4.3]; AlphaIdx=2;
Alpha=AlphaMatrix(AlphaIdx);
Con=AlphaCon(AlphaIdx);


ETA=zeros(m,1)+eta_def;
LAM=zeros(m,1)+lam_def;

CHI=4*nu*Ub.^2/(d_median*((s-1)*g*d_median)^(1.5));
Iclt=find(CHI < 2);
Icgt=find(CHI >= 2);
if isempty(Icgt) == 1
  ETA=Ab*0.30.*CHI.^(-0.39);
  LAM=Ab*1.96.*CHI.^(-0.28);
elseif isempty(Iclt) == 1
  ETA=Ab*0.45.*CHI.^(-0.99);
  LAM=Ab*2.71.*CHI.^(-0.75);
else
  ETA(Iclt)=Ab(Iclt)*0.30.*CHI(Iclt).^(-0.39);
  ETA(Icgt)=Ab(Icgt)*0.45.*CHI(Icgt).^(-0.99);
  LAM(Iclt)=Ab(Iclt)*1.96.*CHI(Iclt).^(-0.28);
  LAM(Icgt)=Ab(Icgt)*2.71.*CHI(Icgt).^(-0.75);
end

for i=1:m
  ub=Ub(i);
  ab=Ab(i);
  deg=Deg(i);
  zr=Zr(i);
  theta=deg*pi/180;
  ur=Ur(i);
  omega=ub/ab;
  if Psi(i)-psicr <= 0
    kbr=kbr_def;
    ETA(i)=eta_def;
    LAM(i)=lam_def;
  else
    kbr=Con*ETA(i);
  end
  ubour=ub/ur;
  ubokur=ubour/kappa;
  kbs=ab*0.0655*(ub^2/((s-1)*g*ab)).^(1.4);
  kb=d_median+kbr+kbs;
  znot=kb/30;
  zrozn=zr/znot;
  abozn=ab/znot;
  alpha=Alpha*(1+beta*kb/ab);
  z1p=alpha;
  delta=1/sqrt(2*z1p);
  mp=delta+sqrt(-1)*delta;
  t1=-alpha*kappa*z1p;
  a1=1.0e-6;
   
  INP=[abozn zrozn ubokur theta a1];
  OUT1=bstress2(INP);
  fofa=OUT1(end);
% compute pure wave limit for upper bound
  if abozn < 6.25
      ubouwmgs=1/abs(t1*mp);
  elseif abozn >= 6.25 && abozn < 10
      ubouwmgs=exp(1.488)*abozn^(-0.653)*...
             abozn^(0.185*log(abozn));
  elseif abozn >= 10 && abozn < 100
      ubouwmgs=exp(0.4599)*abozn^(0.1977)*...
               abozn^(0.0085*log(abozn));
  elseif abozn >= 100
      ubouwmgs=exp(0.13996)*abozn^(0.3539)*...
               abozn^(-0.0106*log(abozn));
  end
  ubouwm=pwave(abozn,ubouwmgs);
  b1=ubouwm;
  fofb=-fofa;
  c1=0.5*(a1+b1);
  INP(end)=c1;
  OUT1=bstress2(INP);
  fofc=OUT1(end);
   cnt1=0;
  for iii=1:mm
    cnt1=cnt1+1;
    sgn=fofb*fofc;
    if sgn < 0
       a1=c1;
    else
       b1=c1;
    end
    c1=0.5*(a1+b1);
    INP(end)=c1;
    OUT1=bstress2(INP);
    fofc=OUT1(end);
    if b1-c1 < tol; break; end
  end
  uboucw=c1;
% format for OUT1: Ro mu epsilon z1ozn z2ozn zroz1 zroz2 fofx
  BBLMPRMS(i,:)=[OUT1 kbs kbr znot ub ab ur];
end

save model_output_file2017 BBLMPRMS SKNPRMS psicr wf d d_median;



