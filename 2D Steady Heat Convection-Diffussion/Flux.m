function [Fe,Fw,Fs,Fn,De,Dw,Dn,Ds]=Flux(i,j)
%Computes Convective Flux & Conductance
% Rertuns F(Convective Flux) and D(Conductance)
global rho u v Xe Xw Yn Ys Gamma 
%==========================================================================
%Convective Fluxes
Fe=rho*u(i,j); %Density- Velocity Product, Flux of u,East Boundary
Fw=rho*u(i-1,j); %Density- Velocity Product, Flux of u,West Boundary
Fs=rho*v(i,j-1); %Density- Velocity Product, Flux of u,South Boundary
Fn=rho*v(i,j); %Density- Velocity Product, Flux of u,North Boundary
%==========================================================================
%Diffuse Fluxes
De=Gamma/Xe(i); %Diffussion Coe., East Bounday
Dw=Gamma/Xw(i);  %Diffussion Coe., West oBunday
Dn=Gamma/Yn(j);  %Diffussion Coe., North ounday
Ds=Gamma/Ys(j); %Diffussion Coe., south Bounday

end



