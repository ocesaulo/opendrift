%     subprogram neutsed_fun.m

% Last update:  4-Aug-2017
% This is for neutral conditions.  It is a matlab function to solve
% the current, suspended sediment concentration and suspended sediment
% transport profiles.

% INPUT
% ustarcw - combined wave and current shear velocity
% ustarc  - time-average shear velocity in the presence of a current
% z1      - inner wave boundary layer height scale
% n       - total number of sediment classes
% znot    - hydraulic roughness
% g       - acceleration due to gravity
% nu      - kinematic viscosity of water at 10 degrees C
% gamma   - ratio of neutral eddy viscosity to neutral eddy diffusivity
% d       - particle diameter

% OUTPUT

% SEDPRO   - matrix containing concentration profiles for all sediment size classes
% CURPRO   - array containing the current profile
% SEDTRANS - matrix containing the sediment transport profiles for all sediment size classes


function [SEDPRO, CURPRO, SEDTRANS, DTRNS, Z]=neutsed_fun(ustarcw,ustarc,znot,z1,hr,d,cref,wf,kappa,gamma)

% nn - number of z-levels used to compute profiles.
nn=500;
z2=z1*ustarcw/ustarc;
uc2oucw=ustarc^2/kappa/ustarcw;
n=length(d);
% Make coordinate transformation
y1=log(z1/znot);
y2=log(z2/znot);
y3=log(hr/znot);
% calculate reference concentration at z1 and z2
arg1=-wf*gamma/kappa/ustarcw;
arg2=-wf*gamma/kappa/ustarc;
cref1=cref.*exp(arg1*y1);
cref2=cref1.*exp(arg1*(exp(y2-y1)-1));
% calculate concentration profiles
n1=y1/y3*nn;
n2=(y2-y1)/y3*nn;
n3=(y3-y2)/y3*nn;
dy1=y1/n1;
dy2=(y2-y1)/n2;
dy3=(y3-y2)/n3;
yy1=[0:dy1:y1];
yy2=[y1:dy2:y2];
yy3=[y2:dy3:y3];
% compute current profiles also
Upro1=uc2oucw*yy1';
Upro2=uc2oucw.*(exp(yy2'-y1)-1)+uc2oucw*y1;
U2=uc2oucw.*(exp(y2-y1)-1)+uc2oucw*y1;
Upro3=ustarc/kappa*(yy3'-y2)+U2;
for i=1:n
  spro1(:,i)=cref(1,i)*exp(yy1'*arg1(1,i));
  spro2(:,i)=cref1(1,i)*exp(arg1(1,i)*...
             (exp(yy2'-y1)-1));
  spro3(:,i)=cref2(1,i)*exp((yy3'-y2)*arg2(1,i));
end

z=[znot*exp(yy1');znot*exp(yy2');znot*exp(yy3')];
SEDPRO=[spro1;spro2;spro3;]; CURPRO=[Upro1;Upro2;Upro3];
for i=1:n
  SEDTRANS(:,i)=SEDPRO(:,i).*CURPRO;
  DTRNS(i)=trapz(z,SEDTRANS(:,i));
end
Z=z;
