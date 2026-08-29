%     program bstress2

% Last update:  19-Oct-99
% This program calculates mu for the 3-layer model
% in matlab, but neglects the region > z2 when calculating
% the wave shear.  It is the same as the fortran program
% fric2.f.


function OUT1=bstress2(INP)
global alpha kappa z1p mp

Abozn=INP(1);
zrozn=INP(2);
ubokur=INP(3);
theta=INP(4);
uboucw=INP(5);


Ro=Abozn/uboucw;
znotp=1/(kappa*Ro);
if z1p/znotp > 1
   phi=phi2_1(znotp);
else
   gammai=-kappa*z1p*mp;
   phi=abs(gammai);
end

mu=sqrt(uboucw*phi);
eps2=-mu^2*abs(cos(theta))+sqrt(1-mu^4*abs(sin(theta)^2));
epsilon=sqrt(eps2);
% if mu > 1; epsilon=1.0e-1; end;

z2p=z1p/epsilon;
Ror=Ro/zrozn;
zroz1=1./(alpha*kappa*Ror);
zroz2=epsilon*zroz1;
z1ozn=alpha*kappa*Ro;
z2ozn=z1ozn/epsilon;


if zroz2 > 1 & z1ozn > 1;
  fofx=ubokur*epsilon*(log(zroz2)+1-epsilon+...
       epsilon*log(z1ozn))-uboucw;
end

if zroz2 <= 1 & zroz1 > 1 & z1ozn > 1;
  fofx=ubokur*epsilon^2*(zroz1-1+log(z1ozn))-uboucw;
end

if zroz1 <= 1 & z1ozn > 1;
  fofx=ubokur*epsilon^2*log(zrozn)-uboucw;
end

if zroz2 > 1 & z1ozn <= 1 & z2ozn > 1;
  fofx=ubokur*epsilon*(log(zroz2)+1....
       -1/z2ozn)-uboucw;
end

if zroz2 <= 1 & zroz1 > 1 &...
  z1ozn <= 1 & z2ozn > 1;
  fofx=ubokur*epsilon^2*(zroz1-1./z1ozn)-uboucw;
end

if zroz2 > 1 & z2ozn <= 1;
  fofx=ubokur*epsilon*log(zrozn)-uboucw;
end

OUT1=[Ro mu epsilon z1ozn z2ozn zroz1 zroz2 fofx];


