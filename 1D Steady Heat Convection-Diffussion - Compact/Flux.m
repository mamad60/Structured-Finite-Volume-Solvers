function [Fe,Fw,De,Dw]=Flux(i)
%Computes Convective Flux & Conductance
% Rertuns F(Convective Flux) and D(Conductance)

global rho u  Xe Xw Gamma 
%==========================================================================
%Convective Fluxes
Fe=rho*u(i); %Dnsity- Velocity Product, Flux of u,East Boundary
Fw=rho*u(i-1); %Dnsity- Velocity Product, Flux of u,West Boundary
%==========================================================================
%Diffuse Fluxes
De=Gamma/Xe(i); %Diffussion Coe., East Bounday
Dw=Gamma/Xw(i);  %Diffussion Coe., West oBunday

end



