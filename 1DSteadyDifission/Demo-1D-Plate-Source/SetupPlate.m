function SetupPlate
%Setups problem for a specific FV problem
% Customize According the problem solved
global Xmin Xmax EXX Gamma0 gStat S0 Sfunc sStat Su0 Sp0 SU SP 
global Bw Be FIw FIe qw qe
%Problem: A Plate with L=2(cm) and Heat Gerneration ans constant Temp. Bcs
%qdot=1000kW/m^3 Ta=100(C), Tb=200(C)
%==========================================================================
%Enter Domain Limits Here
Xmin=0; %Minimum X
Xmax=2e-2; %Maximum X
EXX=1; %Expansion Factor ---X Direction

%=========================================================================
%Enter Cross Section Data Here
S0=1;  %Base Cross-Section---A=1 in y-z Plane
Sfunc=@(S0,x) (S0); %Funtion for Calculation Cross-Sections

%==========================================================================
%Enter Material Properties Here
Gamma0=0.5;%(k=0.5 W/Mk) %Diffussion Coefficient
gStat=1; %1: Gamma is constant 2:Gamama varies with x

%==========================================================================
%Enter Source Terms Here
sStat=2; %1: No source Term 2: Constant Source Tem 3:Source Term is a function of  x
Su0=1000e3; %Base Su----qdot*DeltaV(ADX)---
Sp0=0; %Base Sp
SU=@(Su0,x) (0);  % Enter the Exprission of Su(x) Here
SP=@(Sp0,x) (0);  % Enter the Exprission of Sp(x) Here

%**************************************************************************
%Modify Boundary Conditions below
Bw=0; % 0:Fixed FI 1:Flux Bc
Be=0; % 0:Fixed FI 1:Flux Bc
FIw=100;   % FI at The Left Boundary, Applies only if NL=0
FIe=200;   % FI at The Right Boundary, Applies only if NR=0
qw=0; % Flux at The Left Boundary, Applies only if NL=1
qe=0; % Flux at The Right Boundary, Applies only if NR=1
%**************************************************************************

end

