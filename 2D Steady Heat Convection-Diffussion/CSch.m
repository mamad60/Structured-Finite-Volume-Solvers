function [ A] = CSch(P)
% Returns the Coe. for discritisation of the Convective Method according to to the Selectet Scheme
%   P is the the point Peclect Number(ie. Pe,Ps), A is the output Coe.
global IC
%Return Coe. for Compact Impelemtation of Schemes
switch IC  % Discritize According to the choosen Scheme
    case 0 %CDS
        A=1-0.5*abs(P);     
    case 1 %First Order Upwind(Donor Cell) Scheme
        A=1;
    case 2 %Hybrid Scheme
        A=max(1-0.5*abs(P),0);   
    case 3 %Pwer Law Scheme
        A=max((1-0.1*abs(P))^5,0);   
end

end

