function [Fe,Fw,Fs,Fn,De,Dw,Dn,Ds]=Flux(i,j)
%Computes Convective Flux & Conductance
% Rertuns F(Convective Flux) and D(Conductance)
global rho u v Xe Xw Yn Ys Gamma DX DY
%==========================================================================
%Convective Fluxes
Fe=rho*u(i,j)*DY(i); %Dnsity- Velocity Product, Flux of u,East Boundary
Fw=rho*u(i-1,j)*DY(i); %Dnsity- Velocity Product, Flux of u,West Boundary
Fs=rho*v(i,j-1)*DX(i); %Dnsity- Velocity Product, Flux of u,South Boundary
Fn=rho*v(i,j)*DX(i); %Dnsity- Velocity Product, Flux of u,North Boundary
%==========================================================================
%Diffuse Fluxes
De=Gamma/Xe(i)*DY(i); %Diffussion Coe., East Bounday
Dw=Gamma/Xw(i)*DY(i);  %Diffussion Coe., West oBunday
Dn=Gamma/Yn(j)*DX(i);  %Diffussion Coe., North ounday
Ds=Gamma/Ys(j)*DX(i); %Diffussion Coe., south Bounday

end



