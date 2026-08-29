% function pwave.m

% calculates nondimensional wave shear, phi
% input arguments: znotp,z1p,z2p
% oputput arguments: phi
% this program is for the 2-layer model when z1p/znotp > 1;
function ubouwm = pwave(Abozn,ubouwmgs);
global kappa z1p mp

tol=1.0e-4;
ubouwmn=0.5234511947*ubouwmgs;
ubouwm=ubouwmgs;
ubouwmo=ubouwmgs*.293847;
phio=.9;
cnt=0;
while abs((ubouwmn-ubouwm)/ubouwmn) > tol;
  if cnt > 0;
    ubouwmo=ubouwm;
    ubouwm=ubouwmn;
    phio=phi;
  end;
  if ubouwmn < 0; ubouwm=1.0e-8; end;
  cnt=cnt+1;
  if cnt > 40;
    cnt
      return;
  end;
%  ubouwm=ubouwmn;
  Ro=Abozn/ubouwm;
  znotp=1/(kappa*Ro);

  if z1p/znotp > 1;
    x=[2.*sqrt(znotp),2.*sqrt(z1p)];
    y=x*exp(3*pi*i/4);
    ber=real(besselj(0,y));
    bei=imag(besselj(0,y));
    ker=real(0.5*pi*i*besselh(0,y));
    kei=imag(0.5*pi*i*besselh(0,y));
    ber1=real(besselj(1,y));
    bei1=imag(besselj(1,y));
    berp=(ber1+bei1)/sqrt(2);
    beip=(-ber1+bei1)/sqrt(2);
    ker1=real(0.5*pi*i*besselh(1,y));
    kei1=imag(0.5*pi*i*besselh(1,y));
    kerp=(ker1+kei1)/sqrt(2);
    keip=(-ker1+kei1)/sqrt(2);
% calculate the coefficients for phi
    bnot=ber(1)+i*bei(1);
    knot=ker(1)+i*kei(1);
    bnotp=(berp(1)+i*beip(1))/sqrt(znotp);
    knotp=(kerp(1)+i*keip(1))/sqrt(znotp);

    b1=ber(2)+i*bei(2);
    k1=ker(2)+i*kei(2);
    b1p=(berp(2)+i*beip(2))/sqrt(z1p);
    k1p=(kerp(2)+i*keip(2))/sqrt(z1p);

    ll=mp*b1+b1p;
    nn=mp*k1+k1p;
    argi=bnotp*nn/(bnot*nn-knot*ll)+...
    knotp*ll/(knot*ll-bnot*nn);
    gammai=-kappa*znotp*argi;
    phi=abs(gammai);
  else
   gammai=-kappa*z1p*mp;
   phi=abs(gammai);
  end
  fofsigma=ubouwm-1/phi;
  fofsigmao=ubouwmo-1/phio;
  ubouwmn=ubouwm-fofsigma*(ubouwm-ubouwmo)/(fofsigma-fofsigmao);
end;

ubouwm=ubouwmn;
