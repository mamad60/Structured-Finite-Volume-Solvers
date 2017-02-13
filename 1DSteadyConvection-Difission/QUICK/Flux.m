function [Fe,Fw,De,Dw]=Flux(i)
%Computes Convective Flux & Conductance
% Rertuns F(Convective Flux) and D(Conductance)

global rho u  Xe Xw Gammae Gammaw 
%==========================================================================
%Convective Fluxes
Fe=rho*u(i); %Dnsity- Velocity Product, Flux of u,East Boundary
Fw=rho*u(i-1); %Dnsity- Velocity Product, Flux of u,West Boundary
%==========================================================================
%Diffuse Fluxes
De=Gammae(i)/Xe(i); %Diffussion Coe., East Bounday
Dw=Gammaw(i)/Xw(i);  %Diffussion Coe., West oBunday

end



