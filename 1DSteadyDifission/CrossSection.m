function [Se Sw]=CrossSection
%Return cross Section S for each CV
% Se and Sw is the cross section for east & west Boundary, respectively

global X  NIM NI Sfunc S0
Sw=zeros(1,NI);
Se=zeros(1,NI);
for i=2:NIM
    Sw(i)=Sfunc(S0,X(i-1));
    Se(i)=Sfunc(S0,X(i));
end


end

