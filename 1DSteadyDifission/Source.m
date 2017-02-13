function [Su Sp]=Source(sStat)
%Return Source Term S=Sp*FIp*Su
global XC NIM Su0 Sp0 SU SP Sfunc DX S0
Su=zeros(NIM,1);  % Source Term
Sp=zeros(1,NIM);  % Source Term

switch sStat
    case 1
        %Do nothing
    case 2 %Constant Sourc Term;
        i=2:NIM;
        Su(i)=Su0.*DX(i).*Sfunc(S0,XC(i));
        Sp(i)=Sp0.*DX(i).*Sfunc(S0,XC(i));
    case 3  % Su or Sp is function of X
        i=2:NIM;
        Su(i)=SU(Su0,XC(i)).*DX(i).*Sfunc(S0,XC(i));
        Sp(i)=SP(Sp0,XC(i)).*DX(i).*Sfunc(S0,XC(i));
        
end
