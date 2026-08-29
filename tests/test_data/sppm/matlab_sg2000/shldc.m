function Psicr=shldc(Star);
% CALCULATE CRITICAL SHIELDS PARAMETER
% FOR INITIATION OF SEDIMENT MOTION
% FROM SHIELDS DIAGRAM.
% SCF IS CORRECTION FACTOR.  FORMULAT CAN HANDLE MULTIPLE SIZE CLASSES

scf=1;
for i=1:length(Star);
    star=Star(i);
  if star < 1.5;
    psicr=scf*0.0932*star^(-.707);
  elseif star < 4.0;
    psicr=scf*0.0848*star^(-.473);
  elseif star < 10.0;
    psicr=scf*0.0680*star^(-.314);
  elseif star < 34.0;
    psicr=scf*0.033;
  elseif star < 270;
    psicr=scf*0.0134*star^.255;
  elseif star >= 270;
    psicr=scf*0.056;
  end;
  Psicr(i)=psicr;
end;


